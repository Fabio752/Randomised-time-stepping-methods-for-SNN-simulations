from brian2 import *
import tqdm
import os
import glob
import subprocess
import time
import numpy as np
import numpy as np_
import random as rd
from numba import jit

from src.model import *
import src.utils as utils
import src.measures as measures

def fixed_simulation(dt, rt, I_syn_bar_magnitude, method = 'euler', c = 0.5,
                     profile = True):

    I_syn_bar = I_syn_bar_magnitude * uA / cm**2 # Define coupling strength.

    # Initialise fixed dt and runtime.
    defaultclock.dt = dt * ms
    runtime = rt * ms

    # Initialise neuron group.
    G = NeuronGroup(neurons, eqs, threshold='V > theta', reset='V = Vl',
                    method = method, name = 'G')
    G.V = ''' Vl + I0 / gl * (1 - exp(-c * (i * T) / (N * tau))) '''
    G.s = 1

    # Initialise synapses and connect them all-to-all.     
    S = Synapses(G, G, on_pre= "s += 1", method = method, name = 'S')

    S.connect(condition = 'i != j')

    
    # Initialise spike and state monitor.
    sp_m = SpikeMonitor(G, name = 'sp_m') 
    st_m = StateMonitor(G,["V"], record = True, dt = 1 * ms, name = 'st_m')

    # Build and run network.
    net = Network(G)
    net.add(S)
    net.add(sp_m)
    net.add(st_m)
    
    net.run(runtime, profile = profile)
    return net



def srts_simulation(dt_fixed, dt_rad = 0.5, up_rad = None, low_rad = None, 
                    prj_folder = "SRTS", t_bin = "auto", v_bin = "auto"):

    if up_rad is None:
        up_rad = dt_rad
    if low_rad is None:
        low_rad = dt_rad

    dt = dt_fixed * 10**-3 
    upp_lim = dt * (1 + up_rad)  
    low_lim = dt * (1 - low_rad) 
    rand_dt = "\t\t\tdt = fRand("+ str(low_lim) + ", " + str(upp_lim) + ");\n"
    
    main = './SRTS/main.cpp'
    utils.write_file(main, 88, rand_dt)

    # use same I_syn_bars of fixeddt simulation for consistency.
    I_syn_bars = np.fromfile("bin/fixeddt_I_syn_bars", dtype=np.float64) * 10**-2
    coherences = []
    exec_time = 0

    # Find useful files with glob.
    for file in glob.glob('./SRTS/code_objects/G_stateupdater_*.cpp'):
        state_updater = file
    if t_bin == "auto":
        for file in glob.glob('./SRTS/results/_dynamic_array_st_m_t_*'):
            t_bin = file
    if v_bin == "auto":
        for file in glob.glob('./SRTS/results/_dynamic_array_st_m_V_*'):
            v_bin = file
    

    for I_syn_bar in tqdm.tqdm(I_syn_bars):

        # Plug I_syn_bar in _lio_2 variable in G_stateupdater_codeobject.cpp file
        new_I_syn_bar = "\tconst double _lio_2 = 1.0f*(" + str(I_syn_bar) + " * 0.01)/128; \n"
        utils.write_file(state_updater, 98, new_I_syn_bar)
       
        # Compile (make all) and run (./main) simulation
        cwd = os.path.abspath(prj_folder)
        process = subprocess.Popen(["make", "all"], stdout=subprocess.PIPE, cwd = cwd)
        process.wait()

        start_time = time.time()
        process = subprocess.Popen(["./main"], stdout=subprocess.PIPE, cwd = cwd)
        process.wait()
        exec_time += time.time() - start_time

        Ts = np.fromfile(t_bin, dtype=np.float64)
        Vs = np.reshape(np.fromfile(v_bin,  dtype=np.float64), 
                        (neurons, -1), order = 'F')

        # Find indexes closest to sampling points of 1ms multiples
        idxs = []
        instant_to_record = 5.0
        while instant_to_record < (duration - 0.001):
            exp_loc = int(Ts.shape[0] * instant_to_record / duration)
            idx = utils.find_idx(Ts, exp_loc, instant_to_record)
            idxs.append(idx)
            instant_to_record += 0.001
        
        # filter voltages based on sampling indexes, compute and append coherence.
        Vs = Vs[:, idxs]
        sigma = measures.sync_measure(Vs)
        coherences.append(sigma)

    return coherences, exec_time / len(I_syn_bars)


def nsrts_syngroup_simulation(rt, I_syn_bar_magnitude, method = 'euler', 
                              c = 0.5, dt = 0.25, rad = 0.5, g = 2):

    # SUPPRESS WARNINGS
    BrianLogger.suppress_hierarchy(
        'brian2.synapses.synapses.synapses_dt_mismatch'
    ) 

    start_scope()
    device.reinit()
    device.activate()    
    
    dt = dt * ms
    I_syn_bar = I_syn_bar_magnitude * uA / cm**2

    net = b2.Network()
    group_size = neurons / g

    for p in range(g):   
        rand_dt = rd.uniform(dt * (1 - rad), dt * (1 + rad)) 
        G = NeuronGroup(group_size, subgroup_eqs, 
                        threshold='V > theta', reset='V = Vl', 
                        method = method, name = 'G' + str(p+1), dt = rand_dt)
        
        G.V = '''Vl + I0 / gl * (1 - exp(-c * (p * group_size + i) * T / (g * N * tau)))'''
        G.s = 1
        net.add(G)
    
    idx = 1   
    for p in range(g):
        for q in range(g):
            ng1 = net.__getitem__('G' + str(p + 1))
            ng2 = net.__getitem__('G' + str(q + 1))
            S = Synapses(ng1, ng2, on_pre= "s += 1", 
                         method = method, name = 'S' + str(idx))
            if p == q:
                S.connect(condition='i!=j')
            else :
                S.connect()
            
            net.add(S)
            idx += 1
    
    for p in range(g):
        G = net.__getitem__('G' + str(p + 1))
        sp_m = b2.SpikeMonitor(G, name='sp_m' + str(p + 1)) # SPIKE 
        st_m = b2.StateMonitor(G, ["V", "I_syn", "s"], record = True,
                               dt = 1 * ms, name = 'st_m' + str(p + 1)) # STATE
        
        net.add(sp_m)
        net.add(st_m)
    
    net.run(rt * ms, profile = True)
    return net



@jit(nopython=True)
def nsrts_python_simulation(I_syn_bar_, dt_fixed = 0.25, c = 0.5, dt_rad = 0.5):
    c = 0.5
    t = np.zeros(neurons)
    I_syn_bar_ *= 10**-6 / 10**-4
    i = 0
    
    t_rec = np.zeros(neurons)
    v = np.zeros(neurons) 
    f = np.zeros(neurons)
    s = np.zeros(neurons)
    Vs = np.zeros((neurons, 5000)) 

    # Initialisations
    for n in range(neurons):
        v[n] =  Vl_ + I0_ / gl_ * (1 - np_.exp(-c * (n * T_) / (neurons * tau_))) 
        s[n] = 1                                                 

    while t[i] < rt:
        # Generate random dt
        dt = rd.uniform(dt_fixed * (1 - dt_rad), dt_fixed * (1 + dt_rad))
        if  v[i] > theta_:
            v[i] = Vl_ 
            for j in range(neurons):
                if i != j:
                    s[j] += 1
       
        curr_neuron_ts = int(t_rec[i])
        if t[i] >= 5000 + curr_neuron_ts:
            Vs[i, curr_neuron_ts] = v[i]
            t_rec[i] += 1

        t[i] += dt

        I_syn = (I_syn_bar_ / neurons * f[i] * tau_)
        v[i] += (dt * 10**-3) * (-gl_ * (v[i] - Vl_) + I_syn + I0_) / C_
        f[i] += (dt * 10**-3) * (k_ * s[i] - f[i]) / tau1_
        s[i] += (dt * 10**-3) * (- s[i] / tau2_)                     
        
        i = argmin(t)
    return Vs