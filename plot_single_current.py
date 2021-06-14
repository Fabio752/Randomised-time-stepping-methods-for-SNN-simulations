import brian2 as b2
from brian2 import *
from src.model import *

img_folder = 'img/'
dt = 0.01
rt = 20
c = 0.5
method = 'euler'
I_syn_bar = 1 * uA / cm**2

# Initialise fixed dt and runtime.
defaultclock.dt = dt * ms
runtime = rt * ms

eqs = '''
dV/dt = (-gl * (V - Vl) + I_syn + I0) / C : volt
I_syn = I_syn_bar * f * tau               : amp / metre**2
df/dt = (k * s - f) / tau1                : 1 / second
ds/dt = -s / tau2                         : 1
'''
# Initialise neuron group.
G = NeuronGroup(1, eqs, threshold='V > theta', reset='V = Vl',
                method = method, name = 'G')
G.V = ''' Vl + I0 / gl * (1 - exp(-c * (i * T) / (1 * tau))) '''
G.s = 1

# Initialise synapses and connect them all-to-all.     
S = Synapses(G, G, on_pre= "s += 1", method = method, name = 'S')

S.connect(condition = 'i != j')

# Initialise state monitor.
st_m = b2.StateMonitor(G, ["I_syn"], record = True, name = 'st_m')

# Build and run network.
net = b2.Network(G)
net.add(S)
net.add(st_m)
    
net.run(runtime, profile = True)

# Plot results.
plt.plot(st_m.t / ms, st_m.I_syn[0] * cm**2 / uA, 'k-')

plt.xlabel(r"time ($ms$)")
plt.ylabel(r"$I_{syn}$ ($\frac{\mu A}{cm^2}$")
plt.grid()
plt.savefig(img_folder + '/single_current.pdf')