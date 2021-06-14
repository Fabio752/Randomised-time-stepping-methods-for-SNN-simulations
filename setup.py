from brian2 import *
import os
import src.simulations as sim
import src.utils as utils

if os.path.exists('SRTS/results'):
    filelist = [f for f in os.listdir('SRTS/results')]
    for f in filelist:
        os.remove(os.path.join('SRTS/results', f))

# Create bin folder.
if not os.path.exists('bin'):
    os.makedirs('bin')

# Create img folders.
if not os.path.exists('img'):
    os.makedirs('img')

if not os.path.exists('img/fixeddt'):
    os.makedirs('img/fixeddt')

if not os.path.exists('img/srts'):
    os.makedirs('img/srts')

if not os.path.exists('img/nsrts'):
    os.makedirs('img/nsrts')

if not os.path.exists('img/comparisons'):
    os.makedirs('img/comparisons')

# Density is the number of different I_syn_bar simulated.
density = 50
# lower bound of 0.05 is estimated visually from original paper.

# Simulation time = 10 seconds.
dt = 1

set_device('cpp_standalone', directory='SRTS')
sim.fixed_simulation(dt, dt, 1, method='euler', profile=False)

# Change files.
main = './SRTS/main.cpp'

rand_func = '''
#include <random>
double fRand(double fMin, double fMax){
    double f = (double)std::rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
'''

sim_loop = '''
        double duration = 10.0;
        double t = 0.0;
        double dt;

        srand (time(NULL)); 
        while (t < duration)
        {
            dt = 1e-05;
            _array_defaultclock_dt[0] = dt;
            _array_defaultclock_t[0] = t;
            network.run(dt, NULL, 10.0);
            t += dt;
        } 
'''
utils.write_file(main, 25, rand_func)
utils.write_file(main, 80, sim_loop)

clock = './SRTS/brianlib/clocks.h'
new_clock = "\t\tt[0] += dt[0];\n"
utils.write_file(clock, 25, new_clock)

new_synapses = """
                srand (1);
                const char _cond = (i != _k);
"""

utils.write_file("SRTS/code_objects/S_synapses_create_generator_codeobject.cpp",
                 169, new_synapses)