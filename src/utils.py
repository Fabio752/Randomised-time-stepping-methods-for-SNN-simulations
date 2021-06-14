from brian2 import *
import numpy as np

def write_file(path, line_num, line_content):
    f = open(path, "r")
    list_of_lines = f.readlines()

    list_of_lines[line_num] = line_content

    f = open(path, "w")
    f.writelines(list_of_lines)
    f.close()



def find_idx(Ts, idx, cut_off): 
    if idx >= len(Ts):
        return idx - 1
    elif Ts[idx] > cut_off :
        while(Ts[idx] > cut_off):
            idx -= 1
    else: 
        while(Ts[idx] <= cut_off):
            idx += 1  
            if idx == len(Ts):
                return idx - 1       
        # go one back
        idx -= 1  
    return idx

def concat_voltages(net, g):
    # Sample voltages in range [5s, 10s] with sampling freq of 1ms.
    stm_1 = net.__getitem__('st_m1')
    fifth_sec1 = int(stm_1.t.shape[0] / 2)  
    V_arr = np.array(stm_1.V[:, fifth_sec1:])

    for m in range(1, g):
        stm = net.__getitem__('st_m' + str(m + 1))
        fifth_sec = int(stm.t.shape[0] / 2)
        curr_V_arr = np.array(stm.V[:, fifth_sec:])
        
        # Concatenate in V_arr voltages sampled for each NeuronGroup.
        V_arr = np.concatenate([V_arr, curr_V_arr], axis=0)
    
    return V_arr