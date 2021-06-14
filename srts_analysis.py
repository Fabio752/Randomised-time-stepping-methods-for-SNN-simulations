import src.simulations as sim
import numpy as np
import src.utils as utils
# dt analysis.
dts = np.fromfile("bin/fixeddt_dts", dtype=np.float64)

# srts_coherences = np.array([])
# srts_execution_time = np.array([])

# dt_rad = 0.5
# for dt in dts:
#     print(dt)
#     coherences, exec_time = sim.srts_simulation(dt, dt_rad = dt_rad)
#     srts_coherences = np.append(srts_coherences, coherences)
#     srts_execution_time = np.append(srts_execution_time, exec_time)

# srts_coherences.tofile("bin/srts_coherences_rad" + str(dt_rad))
# srts_execution_time.tofile("bin/srts_execution_time_rad" + str(dt_rad))


# radius analysis.
dt_rads = np.array([0.00001, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 0.5])
dt = 0.01
srts_coherences = np.array([])
srts_execution_time = np.array([])

for dt_rad in dt_rads:
    coherences, exec_time = sim.srts_simulation(dt, dt_rad = dt_rad)
    srts_coherences = np.append(srts_coherences, coherences)
    srts_execution_time = np.append(srts_execution_time, exec_time)

dt_rads.tofile("bin/srts_dt_rads")
srts_coherences.tofile("bin/srts_coherences_radanalysis_dt" + str(dt))
srts_execution_time.tofile("bin/srts_execution_time_radanalysis_dt" + str(dt))


# up_lim = dt_fixed analysis.
# dt_rad = 0.5
# srts_coherences = np.array([])
# srts_execution_time = np.array([])

# for dt in dts:
#     coherences, exec_time = sim.srts_simulation(dt, dt_rad = dt_rad, up_rad=0)
#     srts_coherences = np.append(srts_coherences, coherences)
#     srts_execution_time = np.append(srts_execution_time, exec_time)

# srts_coherences.tofile("bin/srts_coherences_lowrad" + str(dt_rad))
# srts_execution_time.tofile("bin/srts_execution_time_lowrad" + str(dt_rad))

# # up_lim = dt_fixed radanalysis.
# dt = 0.01
# srts_coherences = np.array([])
# srts_execution_time = np.array([])

# for dt_rad in dt_rads:
#     coherences, exec_time = sim.srts_simulation(dt, dt_rad = dt_rad, up_rad=0)
#     srts_coherences = np.append(srts_coherences, coherences)
#     srts_execution_time = np.append(srts_execution_time, exec_time)

# dt_rads.tofile("bin/srts_dt_lowrads")
# srts_coherences.tofile("bin/srts_coherences_lowradanalysis_dt" + str(dt))
# srts_execution_time.tofile("bin/srts_execution_time_lowradanalysis_dt" + str(dt))

# # random coupling 
# srts_randcoup_coherences = np.array([])
# srts_randcoup_execution_time = np.array([])

# coupling = 0.5
# new_synapses = \
#     "\t\t\t\tconst char _cond = (i != _k) && (((double) std::rand() / (RAND_MAX)) < """ \
#     +  str(coupling) + ");\n"

# utils.write_file("SRTS/code_objects/S_synapses_create_generator_codeobject.cpp",
#                 171, new_synapses)
# for dt in dts:
#     # Define random ranges in main.cpp
#     dt_string = "\t\t\tdt = " + str(dt * 10 **-3) + ";\n"
#     utils.write_file("SRTS/main.cpp", 90, dt_string)

#     coherences, exec_time = sim.srts_simulation(dt, dt_rad = 0.5)
#     srts_randcoup_coherences = np.append(srts_randcoup_coherences, 
#                                         coherences)
#     srts_randcoup_execution_time = np.append(srts_randcoup_execution_time,
#                                             exec_time)
                                            
# srts_randcoup_coherences.tofile("bin/srts_randcoup_coherences")
# srts_randcoup_execution_time.tofile("bin/srts_randcoup_execution_time")

# # reset to normal
# new_synapses = "\t\t\t\tconst char _cond = (i != _k);\n" 
# utils.write_file("SRTS/code_objects/S_synapses_create_generator_codeobject.cpp",
#                     171, new_synapses)