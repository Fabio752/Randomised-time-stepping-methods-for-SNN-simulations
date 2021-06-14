import numpy as np
import tqdm 
import matplotlib.pyplot as plt
import src.simulations as sim
import src.measures as measures
import src.utils as utils

img_folder = 'img/fixeddt/'

density = 50 # The number of different I_syn_bar simulated.
I_syn_bars = np.linspace(0.05, 1, num = density)
rt = 10000 # Simulation time = 10 seconds.

dts = [1, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01]

fixeddt_coherences = np.array([])
fixeddt_execution_time = np.array([])

for dt in dts:
    coherences = np.array([])
    exec_time = 0

    for i_bar in tqdm.tqdm(I_syn_bars):
        
        net = sim.fixed_simulation(dt, rt, i_bar, method='euler')

        # Compute execution time from profiling info.
        exec_times = np.array(net.get_profiling_info())      
        for t in exec_times[:, 1]:
            exec_time += float(t)
        
        st_m = net.__getitem__('st_m')

        # FROM PAPER: coherences analysed between 5 and 10s in a 10s simulation 
        # and sampled every millisecond.
        fifth_sec = int(st_m.t.shape[0] / 2)
        V_arr = np.array(st_m.V[:, fifth_sec:])
        sigma = measures.sync_measure(V_arr)
        coherences = np.append(coherences, sigma)
        
    fixeddt_coherences = np.append(fixeddt_coherences, coherences)
    fixeddt_execution_time = np.append(fixeddt_execution_time, 
                                    exec_time / len(I_syn_bars))
    
dts.tofile("bin/fixeddt_dts")
I_syn_bars.tofile("bin/fixeddt_I_syn_bars")
fixeddt_coherences.tofile("bin/fixeddt_coherences")
fixeddt_execution_time.tofile("bin/fixeddt_execution_time")

# random coupling 
fixeddt_randcoup_coherences = np.array([])
fixeddt_randcoup_execution_time = np.array([])

coupling = 0.5
new_synapses = \
    "\t\t\t\tconst char _cond = (i != _k) && (((double) std::rand() / (RAND_MAX)) < """ \
    +  str(coupling) + ");\n"

utils.write_file("SRTS/code_objects/S_synapses_create_generator_codeobject.cpp",
                171, new_synapses)
for dt in dts:
    # Define random ranges in main.cpp
    dt_string = "\t\t\tdt = " + str(dt * 10 **-3) + ";\n"
    utils.write_file("SRTS/main.cpp", 90, dt_string)

    coherences, exec_time = sim.srts_simulation(dt, dt_rad = 0)
    fixeddt_randcoup_coherences = np.append(fixeddt_randcoup_coherences, 
                                            coherences)
    fixeddt_randcoup_execution_time = np.append(fixeddt_randcoup_execution_time,
                                                exec_time)
                                            
fixeddt_randcoup_coherences.tofile("bin/fixeddt_randcoup_coherences")
fixeddt_randcoup_execution_time.tofile("bin/fixeddt_randcoup_execution_time")

dt = 0.0001
# Define random ranges in main.cpp
dt_string = "\t\t\tdt = " + str(dt * 10 **-3) + ";\n"
utils.write_file("SRTS/main.cpp", 90, dt_string)

coherences, exec_time = sim.srts_simulation(dt, dt_rad = 0)
coherences = np.array(coherences)
coherences.tofile("bin/randcoup_ground_truth")

# reset to normal
new_synapses = "\t\t\t\tconst char _cond = (i != _k);\n" 

utils.write_file("SRTS/code_objects/S_synapses_create_generator_codeobject.cpp",
                171, new_synapses)