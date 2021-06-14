import numpy as np
import tqdm 
import matplotlib.pyplot as plt
import src.simulations as sim
import src.measures as measures

img_folder = 'img/fixeddt/'

density = 50 # The number of different I_syn_bar simulated.
I_syn_bars = np.linspace(0.05, 1, num = density)
rt = 10000 # Simulation time = 10 seconds.

dts = [0.25, 0.1, 0.01, 0.001]
styles = ['k--', 'k:', 'k-.', 'k-']

fixed_coherences_log = np.array([])
fixed_execution_time = np.array([])

for dt, s in zip(dts, styles):
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
        
    fixed_coherences_log = np.append(fixed_coherences_log, coherences)
    fixed_execution_time = np.append(fixed_execution_time, 
                                     exec_time / len(I_syn_bars))
    
    # Plot results.
    plt.plot(I_syn_bars, coherences, s, label = 'dt = ' + str(dt))

plt.xlabel(r"$\bar{I_{syn}}$")
plt.ylabel(r"$\Sigma$")
plt.ylabel('y', rotation = 0)
plt.legend()
plt.grid()
plt.savefig(img_folder + 'replica_figure_1.pdf')

fixed_coherences_log = fixed_coherences_log.reshape((len(dts), len(I_syn_bars)))
ground_truth_coherences = fixed_coherences_log[3, :] 
ground_truth_coherences.tofile("bin/ground_truth")
