import tqdm
import src.simulations as sim
import numpy as np
import src.utils as cutils
from src.model import *
import src.measures as measures
import time


gs = [2, 4, 8] # Number of NeuronGroups.
rads = [0.5] # Radius of random dt.

for g in gs:
    for rad in rads: 
        dts = np.fromfile("bin/fixeddt_dts", dtype=np.float64)
        I_syn_bars = np.fromfile("bin/fixeddt_I_syn_bars", dtype=np.float64)

        nsrts_neugroup_coherences = np.array([]) # Log for synch measures.
        nsrts_neugroup_execution_time = np.array([]) # Log for exec times.

        for dt in dts:
            coherences = np.array([])
            exec_time = 0
            for i_bar in tqdm.tqdm(I_syn_bars):
                
                # Run the simulation.
                net = sim.nsrts_syngroup_simulation(rt, i_bar, dt = dt, 
                                                    rad = rad, g = g)
                
                # Compute execution time from profiling info.
                exec_times = np.array(net.get_profiling_info())      
                
                for t in exec_times[:, 1]:
                    exec_time += float(t)
                        
                V_arr = cutils.concat_voltages(net, g)
                # Compute synchronisation measure and append it.
                sigma = measures.sync_measure(V_arr)
                coherences = np.append(coherences, sigma)
            
            # Log coherence and exec time for the full set of simulations.
            nsrts_neugroup_coherences = np.append(nsrts_neugroup_coherences, 
                                                  coherences)
            nsrts_neugroup_execution_time = np.append(nsrts_neugroup_execution_time, 
                                                exec_time / len(I_syn_bars))

        nsrts_neugroup_coherences.tofile("bin/nsrts_neugroup_coherences_g" \
                                            + str(g) + "_rad" + str(rad))
        nsrts_neugroup_execution_time.tofile("bin/nsrts_neugroup_execution_time_g" \
                                                + str(g) + "_rad" + str(rad))

# radanalysis on subgroup.
g = 2 # Number of NeuronGroups.
I_syn_bars = np.fromfile("bin/fixeddt_I_syn_bars", dtype=np.float64)

rads = np.fromfile("bin/srts_dt_rads", dtype=np.float64)

dt = 0.01

nsrts_neugroup_coherences = np.array([]) # Log for synch measures.
nsrts_neugroup_execution_time = np.array([]) # Log for exec times.

for rad in rads:
    coherences = np.array([])
    exec_time = 0
    for i_bar in tqdm.tqdm(I_syn_bars):
        
        # Run the simulation.
        net = sim.nsrts_syngroup_simulation(rt, i_bar, dt = dt, rad = rad, g = g)
        
        # Compute execution time from profiling info.
        exec_times = np.array(net.get_profiling_info())      
        
        for t in exec_times[:, 1]:
            exec_time += float(t)
                
        V_arr = cutils.concat_voltages(net, g)
        # Compute synchronisation measure and append it.
        sigma = measures.sync_measure(V_arr)
        coherences = np.append(coherences, sigma)
    
    # Log coherence and exec time for the full set of simulations.
    nsrts_neugroup_coherences = np.append(nsrts_neugroup_coherences, coherences)
    nsrts_neugroup_execution_time = np.append(nsrts_neugroup_execution_time, 
                                        exec_time / len(I_syn_bars))

nsrts_neugroup_coherences.tofile("bin/nsrts_neugroup_coherences_radanalysis_g" \
                                    + str(g) + "_dt" + str(dt))
nsrts_neugroup_execution_time.tofile("bin/nsrts_neugroup_execution_time_radanalysis_g" \
                                        + str(g) + "_dt" + str(dt))


# Each neuron has different random dt.
dt = 0.01

nsrts_neugroup_coherences = np.array([]) # Log for synch measures.
nsrts_neugroup_execution_time = np.array([]) # Log for exec times.

for rad in rads:
    coherences = np.array([])
    exec_time = 0
    for i_bar in tqdm.tqdm(I_syn_bars):
        
        # Run the simulation.
        net = sim.nsrts_syngroup_simulation(rt, i_bar, dt = dt, rad = rad, g = g)
        
        # Compute execution time from profiling info.
        exec_times = np.array(net.get_profiling_info())      
        
        for t in exec_times[:, 1]:
            exec_time += float(t)
                
        V_arr = cutils.concat_voltages(net, g)
        # Compute synchronisation measure and append it.
        sigma = measures.sync_measure(V_arr)
        coherences = np.append(coherences, sigma)
    
    # Log coherence and exec time for the full set of simulations.
    nsrts_neugroup_coherences = np.append(nsrts_neugroup_coherences, coherences)
    nsrts_neugroup_execution_time = np.append(nsrts_neugroup_execution_time, 
                                        exec_time / len(I_syn_bars))

nsrts_neugroup_coherences.tofile("bin/nsrts_neugroup_coherences_radanalysis_g" \
                                    + str(g) + "_dt" + str(dt))
nsrts_neugroup_execution_time.tofile("bin/nsrts_neugroup_execution_time_radanalysis_g" \
                                        + str(g) + "_dt" + str(dt))
dt_rad = 0.5
dts = np.fromfile("bin/fixeddt_dts", dtype=np.float64)
I_syn_bars = np.fromfile("bin/fixeddt_I_syn_bars", dtype=np.float64)

nsrts_pyvardt_coherences = np.array([])
nsrts_pyvardt_execution_time = np.array([])

for dt in dts: 
    coherences = []
    exec_time = 0
    for I_syn_bar in tqdm.tqdm(I_syn_bars):

        start_time = time.time()
        Vs = sim.nsrts_python_simulation(I_syn_bar, dt_fixed = dt, dt_rad = dt_rad) 
        exec_time += time.time() - start_time
        
        sigma = measures.sync_measure(Vs)
        coherences.append(sigma)
        
    # Log coherence and exec time for the full set of simulations.
    nsrts_pyvardt_coherences = np.append(nsrts_pyvardt_coherences, coherences)
    nsrts_pyvardt_execution_time = np.append(nsrts_pyvardt_execution_time, 
                                        exec_time / len(I_syn_bars))


nsrts_pyvardt_coherences.tofile("bin/nsrts_pyvardt_coherences_rad" + str(dt_rad))
nsrts_pyvardt_execution_time.tofile("bin/nsrts_pyvardt_execution_time_rad" + str(dt_rad))



# Radanalysis on individual.
dt = 0.01
# dt_rads = np.fromfile("bin/srts_dt_rads", dtype=np.float64)

I_syn_bars = np.fromfile("bin/fixeddt_I_syn_bars", dtype=np.float64)
nsrts_pyvardt_coherences = np.array([])
nsrts_pyvardt_execution_time = np.array([])
for dt_rad in rads: 
    coherences = []
    exec_time = 0
    for I_syn_bar in tqdm.tqdm(I_syn_bars):

        start_time = time.time()
        Vs = sim.nsrts_python_simulation(I_syn_bar, dt_fixed = dt, dt_rad = dt_rad) 
        exec_time += time.time() - start_time
        
        sigma = measures.sync_measure(Vs)
        coherences.append(sigma)

    # Log coherence and exec time for the full set of simulations.
    nsrts_pyvardt_coherences = np.append(nsrts_pyvardt_coherences, coherences)
    nsrts_pyvardt_execution_time = np.append(nsrts_pyvardt_execution_time, 
                                        exec_time / len(I_syn_bars))

nsrts_pyvardt_coherences.tofile("bin/nsrts_pyvardt_coherences_radanalysis_dt" + str(dt))
nsrts_pyvardt_execution_time.tofile("bin/nsrts_pyvardt_execution_time_radanalysis_dt" + str(dt))

