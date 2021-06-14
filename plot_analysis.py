import numpy as np
from pyparsing import line
import src.graphics as customplt

## Fixeddt. ##
img_folder = 'img/fixeddt/'
plot_dts = [0.25, 0.1, 0.025, 0.01] # Change this to specify which dts to plot.

customplt.coherence_plotter(plot_dts, "fixeddt_dts", "fixeddt_coherences", 
                            "fixeddt_analysis", img_folder)

customplt.mse_plotter("fixeddt_dts", "fixeddt_coherences", 
                      "fixeddt_execution_time", "fixeddt_mse", img_folder,
                      lab = 'fixeddt method')
##############

##  srts. ##
img_folder = 'img/srts/'

# srts analysis
dt_rad = 0.5
customplt.coherence_plotter(plot_dts, "fixeddt_dts", "srts_coherences_rad" + str(dt_rad), 
                       "srts_analysis_rad" + str(dt_rad), img_folder, 
                       low_lim = 0.5, up_lim = 0.5)

customplt.mse_plotter("fixeddt_dts", "srts_coherences_rad" + str(dt_rad),
                      "srts_execution_time_rad" + str(dt_rad), 
                      "srts_mse_rad" + str(dt_rad), img_folder, 
                      low_lim = 0.5, up_lim = 0.5, lab = 'srts method')

# srts radanalysis.
dt = 0.01
plot_rads = [0.01, 0.1, 0.5, 0.9]

customplt.radanalysis_plotter(dt, plot_rads, "srts_dt_rads", 
                              "srts_coherences_radanalysis_dt" + str(dt), 
                              "srts_radanalysis_dt" + str(dt), img_folder)

customplt.mse_radanalysis_plotter(dt, "srts_dt_rads", 
                                  "srts_coherences_radanalysis_dt" + str(dt), 
                                  "srts_execution_time_radanalysis_dt" + str(dt), 
                                  "srts_radanalysis_mse_dt" + str(dt), 
                                  img_folder, lab = 'SRTS', 
                                  yoff=0.00015)


# up_lim = dt_fixed analysis
dt_rad = 0.5
customplt.coherence_plotter(plot_dts, "fixeddt_dts", 
                            "srts_coherences_lowrad" + str(dt_rad), 
                            "srts_analysis_lowrad" + str(dt_rad), img_folder, 
                            low_lim = 0.5, up_lim = 0)

customplt.mse_plotter("fixeddt_dts", "srts_coherences_lowrad" + str(dt_rad),
                      "srts_execution_time_rad" + str(dt_rad), 
                      "srts_mse_lowrad" + str(dt_rad), img_folder, 
                      low_lim = 0.5, up_lim = 0, lab = 'srts method')

                
# up_lim = dt_fixed radanalysis

customplt.radanalysis_plotter(dt, plot_rads, "srts_dt_rads", 
                            "srts_coherences_lowradanalysis_dt" + str(dt), 
                            "srts_lowradanalysis_dt" + str(dt), img_folder,
                            up_lim=0)

customplt.mse_radanalysis_plotter(dt, "srts_dt_rads", 
                                  "srts_coherences_lowradanalysis_dt" + str(dt), 
                                  "srts_execution_time_lowradanalysis_dt" + str(dt), 
                                  "srts_lowradanalysis_mse_dt" + str(dt),
                                  img_folder, up_lim = 0, lab = 'SRTS')


# Random Coupling #
# Fixeddt.
img_folder = 'img/fixeddt/'

customplt.coherence_plotter(plot_dts, "fixeddt_dts", 
                            "fixeddt_randcoup_coherences",
                            "fixeddt_randcoup_analysis", img_folder, 
                            gt_file="randcoup_ground_truth", gt_dt = 0.0001)

customplt.mse_plotter("fixeddt_dts", "fixeddt_randcoup_coherences", 
                      "fixeddt_randcoup_execution_time", 
                      "fixeddt_randcoup_mse", img_folder,
                      gt_file="randcoup_ground_truth",
                      lab = 'FTS with 0.5 coupling')

# srts.
img_folder = 'img/srts/'

customplt.coherence_plotter(plot_dts, "fixeddt_dts", 
                            "srts_randcoup_coherences",
                            "srts_randcoup_analysis", img_folder, 
                            gt_file="randcoup_ground_truth", gt_dt = 0.0001)

customplt.mse_plotter("fixeddt_dts", "srts_randcoup_coherences", 
                      "srts_randcoup_execution_time", 
                      "srts_randcoup_mse", img_folder,
                      gt_file="randcoup_ground_truth",
                      lab = 'SRTS(0.5) with 0.5 coupling')
##########


## nsrts. ##
img_folder = 'img/nsrts/'
rad = 0.5
for g in [2, 4, 8]:
    # neugroup shared
    customplt.coherence_plotter(plot_dts, "fixeddt_dts", 
                                "nsrts_neugroup_coherences_g" + str(g) +"_rad" + str(rad),
                                "nsrts_neugroup_analysis_g" + str(g) +"_rad" + str(rad), 
                                img_folder, up_lim=rad, low_lim=rad)

    customplt.mse_plotter("fixeddt_dts", 
                        "nsrts_neugroup_coherences_g" + str(g) +"_rad" + str(rad),
                        "nsrts_neugroup_execution_time_g" + str(g) +"_rad" + str(rad),
                        "nsrts_neugroup_mse_g" + str(g) +"_rad" + str(rad), img_folder, 
                        lab = 'SRTS (' + str(g) + ', ' + str(rad) + ')')


# neugroup radanalysis
dt = 0.01
plot_rads = [0.1, 0.5, 0.9]
g = 2

customplt.radanalysis_plotter(dt, plot_rads, "srts_dt_rads", 
                            "nsrts_neugroup_coherences_radanalysis_g" + str(g) \
                            + '_dt' + str(dt),
                            "nsrts_neugroup_radanalysis_g" + str(g) + "_dt" \
                            + str(dt), img_folder)

customplt.mse_radanalysis_plotter(dt, "srts_dt_rads", 
                                  "nsrts_neugroup_coherences_radanalysis_g" \
                                  + str(g) + '_dt' + str(dt),
                                  "nsrts_neugroup_execution_time_radanalysis_g" \
                                  + str(g) + '_dt' + str(dt), 
                                  "nsrts_neugroup_radanalysis_mse_g" + str(g) \
                                  + "_dt" + str(dt),  
                                  img_folder, lab = 'SSRTS(2, x)')



# individual
customplt.coherence_plotter(plot_dts, "fixeddt_dts", 
                            "nsrts_pyvardt_coherences_rad" + str(rad),
                            "nsrts_pyvardt_analysis_rad" + str(rad), img_folder)

customplt.mse_plotter("fixeddt_dts", "nsrts_pyvardt_coherences_rad" + str(rad), 
                      "nsrts_pyvardt_execution_time_rad" + str(rad), 
                      "nsrts_pyvardt_mse_rad" + str(rad), img_folder,
                       lab = 'indiv. nsrts (rad = ' + str(rad) + ')')


# individual radanalysis
dt = 0.01
plot_rads = [0.1, 0.5, 0.9]

customplt.radanalysis_plotter(dt, plot_rads, "srts_dt_rads", 
                            "nsrts_pyvardt_coherences_radanalysis", 
                            "nsrts_pyvardt_radanalysis", img_folder)

customplt.mse_radanalysis_plotter(dt, "srts_dt_rads", 
                                  "nsrts_pyvardt_coherences_radanalysis", 
                                  "nsrts_pyvardt_execution_time_radanalysis", 
                                  "nsrts_pyvardt_radanalysis_mse", 
                                  img_folder, lab = 'NSRTS')
dt = 0.1

customplt.radanalysis_plotter(dt, plot_rads, "srts_dt_rads", 
                            "nsrts_pyvardt_coherences_radanalysis_dt" + str(dt), 
                            "nsrts_pyvardt_radanalysis_dt" + str(dt), img_folder)

customplt.mse_radanalysis_plotter(dt, "srts_dt_rads", 
                                  "nsrts_pyvardt_coherences_radanalysis_dt" + str(dt), 
                                  "nsrts_pyvardt_execution_time_radanalysis_dt" + str(dt), 
                                  "nsrts_pyvardt_radanalysis_mse_dt" + str(dt), 
                                  img_folder, lab = 'NSRTS')



# comparison studies:
# FTS only
img_folder = 'img/fixeddt/'
mse_files = [
    "fixeddt_mse",
]
exec_files = [
    "fixeddt_execution_time",
]
labels = [
    "FTS",
]

customplt.multiple_mses_plotter("fixeddt_dts", mse_files, exec_files, labels,
                                "mse_fixeddt_labelled", img_folder, yoff = 0.75,
                                line = True, cols=['k'])

# FTS only
img_folder = 'img/fixeddt/'
mse_files = [
    "fixeddt_randcoup_mse",
]
exec_files = [
    "fixeddt_randcoup_execution_time",
]
labels = [
    "FTS with 0.5 random coupling",
]

customplt.multiple_mses_plotter("fixeddt_dts", mse_files, exec_files, labels,
                                "mse_fixeddt_randcoup_labelled", img_folder, 
                                xoff = 0.0006, yoff = -0.07, line = True, cols=['k'])


# SRTS only 
img_folder = 'img/srts/'
mse_files = [
    "srts_randcoup_mse",
]
exec_files = [
    "srts_randcoup_execution_time",
]
labels = [
    "SRTS with 0.5 random coupling",
]

customplt.multiple_mses_plotter("fixeddt_dts", mse_files, exec_files, labels,
                                "mse_srts_randcoup_labelled", img_folder, 
                                yoff = 0.15, line = True, cols=['k'])




# FTS vs SRTS
img_folder = 'img/comparisons/'
dt_rad = 0.5
skip = 0
mse_files = [
    "fixeddt_mse",
    "srts_mse_rad" + str(dt_rad), 
]
exec_files = [
    "fixeddt_execution_time",
    "srts_execution_time_rad" + str(dt_rad), 
]
labels = [
    "FTS",
    "SRTS(" + str(dt_rad) + ')'
]

customplt.multiple_mses_plotter("fixeddt_dts", mse_files, exec_files, labels,
                                "mse_fixeddt_vs_SRTS", img_folder, skip=skip)


# FTS vs SSRTS
img_folder = 'img/comparisons/'
dt_rad = 0.5
skip = 2
for g in [2, 4, 8] :
    mse_files = [
        "fixeddt_mse",
        "nsrts_neugroup_mse_g" + str(g) + "_rad" + str(dt_rad), 
    ]
    exec_files = [
        "fixeddt_execution_time",
        "nsrts_neugroup_execution_time_g" + str(g) + "_rad" + str(dt_rad),     ]
    labels = [
        "FTS",
        "SSRTS(" + str(g) + ", " + str(dt_rad) + ")",
    ]

    customplt.multiple_mses_plotter("fixeddt_dts", mse_files, exec_files, labels,
                                    "mse_fixeddt_vs_SSRTS_g" + str(g), 
                                    img_folder, skip=skip, yoff= g * 2, on_m=1)


# FTS vs NSRTS
img_folder = 'img/comparisons/'
dt_rad = 0.5
skip = 2

mse_files = [
    "fixeddt_mse",
    "nsrts_pyvardt_mse_rad" + str(dt_rad), 
]
exec_files = [
    "fixeddt_execution_time",
    "nsrts_pyvardt_execution_time_rad" + str(dt_rad),     ]
labels = [
    "FTS",
    "NSRTS(" + str(dt_rad) + ")",
]

customplt.multiple_mses_plotter("fixeddt_dts", mse_files, exec_files, labels,
                                "mse_fixeddt_vs_NSRTS", img_folder, skip=skip)


# gstudy

img_folder = 'img/nsrts/'
dt = 0.01
dt_rad = 0.5
plot_gs = [2, 4, 8]
coh_files = []
labels = []
for g in plot_gs :
    coh_files.append("nsrts_neugroup_coherences_g" + str(g) + "_rad" + str(dt_rad))
    labels.append("SSRTS(" + str(g) + ", " + str(dt_rad) + ")")

customplt.gstudy_plotter(dt, dt_rad, plot_gs, coh_files, "gstudy_coherences",
                         img_folder, labels)


# mse gstudy
img_folder = 'img/nsrts/'
mse_files = []
exec_files = []

for g in plot_gs :
    mse_files.append("nsrts_neugroup_mse_g" + str(g) + "_rad" + str(dt_rad))
    exec_files.append("nsrts_neugroup_execution_time_g" + str(g) + "_rad" + str(dt_rad))
    labels.append("SSRTS(" + str(g) + ", " + str(dt_rad) + ")")


customplt.multiple_mses_plotter("fixeddt_dts", mse_files, exec_files, labels,
                                "mse_gstudy" + str(g), img_folder, 
                                cols = ['k', 'b', 'skyblue'], connect = False, 
                                yoff = 12, size_mult = 250, skip = 2)

