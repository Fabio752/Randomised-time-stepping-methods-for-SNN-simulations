from brian2 import *
# Parameters
neurons   =     128                           # number of neurons
gl        =     0.1   *    mS / cm**2         # conductance of the voltage-independent leak current
Vl        =     -60   *       mV              # reversal potential of the voltage-independent leak current
C         =       1   *  ufarad / cm**2       # membrane capacitance
tau       =      10   *       ms          
theta     =     -40   *       mV              # membrane potential threshold
tau1      =       3   *       ms              
tau2      =       1   *       ms    
I0        =     2.3   *     uA / cm**2        # I0 / C -> get 2.3 kHz, time that by tau to get 23 mvolts. 

# FOR INITIAL CONDITION 
T = -tau * log(1 - gl/I0 * (theta - Vl))

k = 1 / tau2

eqs = '''
dV/dt = (-gl * (V - Vl) + I_syn + I0) / C : volt
I_syn = (I_syn_bar / N) * f * tau         : amp / metre**2
df/dt = (k * s - f) / tau1                : 1 / second
ds/dt = -s / tau2                         : 1
'''

rt = 10000
duration = rt / 1000

gl_ = 0.1 * 10 **-3 / 10 ** -4
Vl_ =-60 * 10 **-3
C_ = 1 * 10 **-6 / 10 **-4
tau_ = 10 * 10 **-3
theta_ = -40 * 10 **-3
tau1_ = 3 * 10 **-3
tau2_ = 1 * 10 **-3
I0_ = 2.3 * 10 **-6 / 10 **-4
T_ = -tau_ * log(1 - gl_/I0_ * (theta_ - Vl_))
k_ = 1 / tau2_

subgroup_eqs = '''
dV/dt = (-gl * (V - Vl) + I_syn + I0) / C : volt
I_syn = (I_syn_bar / (g * N)) * f * tau   : amp / metre**2
df/dt = (k * s - f) / tau1                : 1 / second
ds/dt = -s / tau2                         : 1
'''