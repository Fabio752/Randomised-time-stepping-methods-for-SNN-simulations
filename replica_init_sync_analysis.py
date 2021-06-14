import numpy as np
import tqdm 
import matplotlib.pyplot as plt
import src.simulations as sim
import src.measures as measures

img_folder = 'img/fixeddt/'

density = 50 # The number of different I_syn_bar simulated.
I_syn_bars = np.linspace(0.05, 1, num = density)
rt = 10000 # Simulation time = 10 seconds.

dt = 0.25
styles = ['k--', 'k-', 'k:']
cs = [0.1, 0.5, 0.9]

for c, s in zip(cs, styles):
    coherences = np.array([])
    exec_time = 0

    for i_bar in tqdm.tqdm(I_syn_bars):
        
        net = sim.fixed_simulation(dt, rt, i_bar, method='euler', c = c)

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
        
    # Plot results.
    plt.plot(I_syn_bars, coherences, s, label = 'c = ' + str(c))

plt.xlabel(r"$\bar{I_{syn}}$")
plt.ylabel(r"$\Sigma$")
plt.ylabel('y', rotation = 0)
plt.legend()
plt.grid()
plt.savefig(img_folder + 'replica_figure_2.pdf')

