import math
import numpy as np

# Materials are passed or updated by the optimizer script
tau_max_f = 480e6 
tau_max_l = 300e6
sigma_y_lug = 300e6
ultimate_bending_lug = 450e6

# Metric bolt list [mm]
metric_bolt_d_i = [2.459, 2.850, 3.242, 3.688, 4.134, 4.917, 5.917, 6.647] # [1.567, 1.713, 2.013, 2.459, 2.850, 3.242, 3.688, 4.134, 4.917, 5.917, 6.647]
metric_bolt_name = ['M3', 'M3.5', 'M4', 'M4.5', 'M5', 'M6', 'M7', 'M8'] # ['M2', 'M2.2', 'M2.5', 'M3', 'M3.5', 'M4', 'M4.5', 'M5', 'M6', 'M7', 'M8']
metric_bolt_d_o = [3, 3.5, 4, 4.5, 5, 6, 7, 8] #[2, 2.2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8]

def fastener_diameter(F_x, F_z, n_l, n_f, tau_max_f_in):
    F_s = F_z / (n_l * n_f)   # shear
    F_n = F_x / (n_l * n_f)   # axial
    D = math.sqrt((4 / (math.pi * tau_max_f_in)) * math.sqrt(F_n**2 + 3*F_s**2))
    
    for i, d in enumerate(metric_bolt_d_i):
        if d >= D*1000:
            D_metric_i = metric_bolt_d_i[i]*0.001       # meters
            D_metric_o = metric_bolt_d_o[i]*0.001
            M_size = metric_bolt_name[i]
            break
    else:
        raise ValueError("No metric bolt large enough")
    
    A = math.pi*D_metric_i**2/4
    sigma_vm = math.sqrt((F_n/A)**2 + 3*(F_s/A)**2)
    return (tau_max_f_in / sigma_vm - 1), D_metric_o, M_size

def Ixx(w_1, t_1):
    return max((t_1**3 * w_1) / 12.0, 1e-15)

def bending(I_xxh, w_2, F_z, t_2, n_l):
    M = (F_z / n_l) * (w_2 / 2)
    bending_applied = (M * (t_2 / 2)) / I_xxh
    bending_margin = (ultimate_bending_lug / bending_applied) - 1
    return t_2, w_2, bending_margin, M

def length_lug(F_z, t_2, tau_max_l_in, sigma_y_lug_in, n_l, w_1, M, w_2):
    I_xxv = max((t_2**3 * w_1) / 12.0, 1e-15)
    bending_stress = (M * (t_2 / 2)) / I_xxv
    shear_stress = (F_z / n_l) / (w_1 * t_2)
    combined_stress = math.sqrt(bending_stress**2 + 3 * shear_stress**2)
    
    combined_margin = (sigma_y_lug_in / combined_stress) - 1.0
    shear_capacity = (8.0 * t_2 * tau_max_l_in) / (3.0 * w_1)
    shear_margin = shear_capacity / (F_z / n_l) - 1.0
    return combined_margin, shear_margin

def CreateFastenerList_v(n_f, l, w_1, D_fo):
    array = []
    n_side = n_f // 2
    if n_side <= 1:
        y_positions = [0.0]
    else:
        d = (l - 3*D_fo) / (n_side - 1)
        y_positions = [-(2*D_fo + i*d) + l/2 for i in range(n_side)]
    for y in y_positions:
        array.append([+w_1/2 - 2*D_fo, y, 0, 0, 0])
        array.append([-w_1/2 + 2*D_fo, y, 0, 0, 0])
    return array

def ForceAtFastener_v(fasteners, w_2, n_f, n_l, F_z, F_y):
    ri_sq_tot = sum(f[0]**2 + f[1]**2 for f in fasteners)
    for i in range(len(fasteners)):
        fasteners[i][2] = F_y / n_l / n_f
        M_x = F_z * w_2 / 2
        fasteners[i][2] += (M_x * fasteners[i][0] / ri_sq_tot) if ri_sq_tot > 0 else 0
    return fasteners

def bearing_stress_v(fasteners, t_lug, t_sc, D_fo):
    for i in range(len(fasteners)):
        fasteners[i][3] = fasteners[i][2] / (D_fo * t_lug)
        fasteners[i][4] = fasteners[i][2] / (D_fo * t_sc)
    return fasteners

def WorstShearMargin_v(fasteners, yield_lug, yield_sc):
    margins = []
    for f in fasteners:
        margins.append(yield_lug / abs(f[3]) - 1 if f[3] != 0 else 10)
        margins.append(yield_sc / abs(f[4]) - 1 if f[4] != 0 else 10)
    return min(margins)

def tearout_vertical(fasteners, l, w_1, t_2):
    min_dist = min(min(w_1/2 + f[0], w_1/2 - f[0], l/2 + f[1], l/2 - f[1]) for f in fasteners)
    F_max = max(abs(f[2]) for f in fasteners)
    return (2 * min_dist * t_2 * sigma_y_lug) / F_max - 1