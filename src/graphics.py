import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from numpy.core.arrayprint import dtype_short_repr

import src.measures as measures


def coherence_plotter(plot_dts, dts_file, coh_file, img_name, img_folder,
                      gt_file = "ground_truth", gt_dt = 0.001, up_lim = None, 
                      low_lim = None, styles = ['k--', 'k:', 'k-.', 'k-']):
    plt.clf()
    ground_truth = np.fromfile("bin/" + gt_file,  dtype=np.float64)
    I_syn_bars = np.fromfile("bin/fixeddt_I_syn_bars", dtype=np.float64)
    
    dts = np.fromfile("bin/" + dts_file, dtype=np.float64)
    coherences = np.fromfile("bin/" + coh_file, dtype=np.float64)
    coherences = coherences.reshape((len(dts), len(I_syn_bars)))

    plt.plot(I_syn_bars, ground_truth, 'b-', label = 'gt: dt = ' + str(gt_dt))

    styles_idx = 0
    for idx, dt in enumerate(dts):
        if dt in plot_dts:
            curr_coherences = coherences[idx, :]
            l = 'dt = ' + str(dt)
            if (up_lim is not None) and (low_lim is not None):
                l = 'dt = [' + str(round(dt * (1 - low_lim), 5)) + \
                    ', ' + str(round(dt * (1 + up_lim), 5)) + ']'
               
            plt.plot(I_syn_bars, curr_coherences, styles[styles_idx], label = l)
            styles_idx += 1

    plt.xlabel(r"$\bar{I_{syn}}$")
    plt.ylabel(r"$\Sigma$", rotation = 0)
    plt.legend()
    plt.grid()
    plt.savefig(img_folder + img_name + '.pdf')


def radanalysis_plotter(dt, plot_rads, rads_file, coh_file, img_name, img_folder,
                        gt_file = "ground_truth", gt_dt = 0.001, up_lim = None, 
                        low_lim = None, styles = ['k--', 'k:', 'k-.', 'k-']):
    plt.clf()

    ground_truth = np.fromfile("bin/" + gt_file,  dtype=np.float64)
    I_syn_bars = np.fromfile("bin/fixeddt_I_syn_bars", dtype=np.float64)
    
    dt_rads = np.fromfile("bin/" + rads_file, dtype=np.float64)
    coherences = np.fromfile("bin/" + coh_file, dtype=np.float64)
    coherences = coherences.reshape((-1, len(I_syn_bars)))

    plt.plot(I_syn_bars, ground_truth, 'b-', 
             label = 'gt: dt = ' + str(gt_dt))

    styles_idx = 0
    for idx, dt_rad in enumerate(dt_rads):
        u = str(round(dt * (1 + dt_rad), 5))
        if up_lim is not None:
            u =  str(round(dt * (1 + up_lim), 5))
        b = str(round(dt * (1 - dt_rad), 5))
        if low_lim is not None:
            b = str(round(dt * (1 - low_lim), 5))
        
        l = 'dt = [' + b +  ', ' + u + ']'
        
        if dt_rad in plot_rads:
            curr_coherences = coherences[idx, :]
            plt.plot(I_syn_bars, curr_coherences, styles[styles_idx], 
                    label = l)
            styles_idx += 1

    plt.xlabel(r"$\bar{I_{syn}}$")
    plt.ylabel(r"$\Sigma$", rotation = 0)
    plt.legend()
    plt.grid()
    plt.savefig(img_folder + img_name + '.pdf')


def mse_plotter(dts_file, coh_file, exec_file, img_name, img_folder, lab = None, 
                gt_file = "ground_truth", col = 'k', size_mult = 250,
                log_scale = False, line = True, up_lim = None, low_lim = None):

    plt.clf()
    execution_time = np.fromfile("bin/" + exec_file, dtype=np.float64)
    ground_truth = np.fromfile("bin/" + gt_file,  dtype=np.float64)
    I_syn_bars = np.fromfile("bin/fixeddt_I_syn_bars", dtype=np.float64)    
    dts = np.fromfile("bin/" + dts_file, dtype=np.float64)
    coherences = np.fromfile("bin/" + coh_file, dtype=np.float64)
    coherences = coherences.reshape((len(dts), len(I_syn_bars)))

    mses = np.array([]) 

    ylab = "Execution time (sec)"
    if log_scale:
        execution_time = np.log10(execution_time)
        ylab = "Log execution time (log10(sec))"

    for idx, dt in enumerate(dts):
        curr_coherences = coherences[idx, :]
        mse_score = measures.mse(curr_coherences, ground_truth)
        mses = np.append(mses, mse_score)
        exec_time = execution_time[idx]
        l = 'dt = ' + str(dt)
        if (up_lim is not None) and (low_lim is not None):
            l = 'dt = [' + str(round(dt * (1 - low_lim), 5)) + \
                ', ' + str(round(dt * (1 + up_lim), 5)) + ']'
        plt.scatter(mse_score, exec_time, color = col, 
                    s = 10 + dt * size_mult)    
    
    if line :
        plt.plot(mses, execution_time, col + '--', label = lab)
    
    plt.legend()
    plt.grid()
    plt.xlabel('MSE')
    plt.ylabel(ylab, rotation = 90)
    plt.savefig(img_folder + img_name + '.pdf')

    mses.tofile('bin/' + img_name)


def mse_radanalysis_plotter(dt, rads_file, coh_file, exec_file, img_name, 
                            img_folder, gt_file = "ground_truth", col = 'k', 
                            size_mult = 150, log_scale = False, line = True,
                            up_lim = None, low_lim = None, lab = None, 
                            range = None, xoff = 0, yoff = 0.75):

    plt.clf()

    dts = np.fromfile("bin/fixeddt_dts", dtype=np.float64)
    idx_dt = np.where(dts == dt)
    fixed_mses = np.fromfile('bin/fixeddt_mse', dtype=np.float64)
    fixed_exec_times = np.fromfile('bin/fixeddt_execution_time', dtype=np.float64)


    execution_time = np.fromfile("bin/" + exec_file, dtype=np.float64)
    ground_truth = np.fromfile("bin/" + gt_file,  dtype=np.float64)
    I_syn_bars = np.fromfile("bin/fixeddt_I_syn_bars", dtype=np.float64)    
    dt_rads = np.fromfile("bin/" + rads_file, dtype=np.float64)
    coherences = np.fromfile("bin/" + coh_file, dtype=np.float64)
    coherences = coherences.reshape(-1, len(I_syn_bars))

    if range is not None:
        execution_time = execution_time[range]
        dt_rads = dt_rads[range]
        coherences = coherences[range, :]


    mses = np.array([]) 
    ylab = "Execution time (sec)"
    if log_scale:
        execution_time = np.log10(execution_time)
        ylab = "Log execution time (log10(sec))"

    for idx, dt_rad in enumerate(dt_rads):
        curr_coherences = coherences[idx, :]
        mse_score = measures.mse(curr_coherences, ground_truth)
        mses = np.append(mses, mse_score)
        exec_time = execution_time[idx]
        
        u = str(round(dt * (1 + dt_rad), 5))
        if up_lim is not None:
            u =  str(round(dt * (1 + up_lim), 5))
        b = str(round(dt * (1 - dt_rad), 5))
        if low_lim is not None:
            b = str(round(dt * (1 - low_lim), 5))
        
        l = 'dt = [' + b +  ', ' + u + ']'
       
        plt.scatter(mse_score, exec_time, color = col, 
                    s = 30 + dt_rad * size_mult) 
  
    
    # plt.scatter(fixed_mses[idx_dt], fixed_exec_times[idx_dt],
    #             color = 'b', s = 100, label = "FTS")
    
    if line :
        plt.plot(mses, execution_time, col + '--', label = lab)        

    for mse_score, exec_time, dt_rad in zip(mses, execution_time, dt_rads):
        plt.text(xoff + mse_score, yoff + exec_time, 
                 str(dt_rad), fontsize = 10, ha="center")
    plt.legend()
    plt.grid()
    plt.xlabel('MSE')
    plt.ylabel(ylab, rotation = 90)
    plt.savefig(img_folder + img_name + '.pdf')

    mses.tofile('bin/' + img_name)


def multiple_mses_plotter(dts_file, mse_files, exec_files, labels, img_name, 
                          img_folder, line = False, cols = ['b', 'k'],
                          size_mult = 150, log_scale = False, skip = 0,
                          xoff = 0, yoff = 0.5, on_m = 0, connect = True):

    plt.clf()
    patches = []
    all_mses = np.array([])
    all_exec_times = np.array([])
    
    dts = np.fromfile("bin/" + dts_file, dtype=np.float64)
    for mse_file, exec_file, lab, c in zip(mse_files, exec_files, labels, cols):
        mse_scores = np.fromfile("bin/" + mse_file, dtype=np.float64)
        exec_times = np.fromfile("bin/" + exec_file, dtype=np.float64)
        
        ylab = "Execution time (sec)"
        if log_scale:
            exec_times = np.log10(exec_times)
            ylab = "Log execution time (log10(sec))"

        for dt, exec_time, mse_score in \
            zip(dts[skip:], exec_times[skip:], mse_scores[skip:]):

            plt.scatter(mse_score, exec_time, color = c, s = 30 + dt * size_mult)    
        
        if line :
            plt.plot(mse_scores[skip:], exec_times[skip:], c + '--', label = lab)
        else :
            p = mpatches.Patch(color=c, label=lab)
            patches.append(p)

        all_mses = np.append(all_mses, mse_scores[skip:])
        all_exec_times = np.append(all_exec_times, exec_times[skip:])

    if line:   
        plt.legend()
    else:
        plt.legend(handles = patches)
    
    all_mses = np.reshape(all_mses, (len(mse_files), -1))
    all_exec_times = np.reshape(all_exec_times, (len(mse_files), -1))
    for i in range(all_mses.shape[1]):
        if connect:
            plt.plot(all_mses[:, i], all_exec_times[:, i], 'k--')
        plt.text(xoff + all_mses[on_m, i], yoff + all_exec_times[on_m, i], 
                    str(dts[skip + i]), fontsize = 10, ha="center")

    plt.grid()
    plt.xlabel('MSE')
    plt.ylabel(ylab, rotation = 90)
    plt.savefig(img_folder + img_name + '.pdf')


def gstudy_plotter(dt, dt_rad, plot_gs, coh_files, img_name, img_folder,
                   labels, gt_file = "ground_truth", gt_dt = 0.001,
                   styles = ['k--', 'k:', 'k-.', 'k-']):
    plt.clf()
    dts = np.fromfile("bin/fixeddt_dts", dtype=np.float64)
    idx_dt = np.where(dts == dt)
    ground_truth = np.fromfile("bin/" + gt_file,  dtype=np.float64)
    I_syn_bars = np.fromfile("bin/fixeddt_I_syn_bars", dtype=np.float64)
    
    plt.plot(I_syn_bars, ground_truth, 'b-', 
        label = 'gt: dt = ' + str(gt_dt))
    for g, coh_file, l in zip(plot_gs, coh_files, labels):
        coherences = np.fromfile("bin/" + coh_file, dtype=np.float64)
        coherences = coherences.reshape((-1, len(I_syn_bars)))

        styles_idx = 0
        curr_coherences = coherences[idx_dt, :].squeeze()

        plt.plot(I_syn_bars, curr_coherences, styles[styles_idx], 
                label = l)
        styles_idx += 1

    plt.xlabel(r"$\bar{I_{syn}}$")
    plt.ylabel(r"$\Sigma$", rotation = 0)
    plt.legend()
    plt.grid()
    plt.savefig(img_folder + img_name + '.pdf')

