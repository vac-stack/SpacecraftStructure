import numpy as np
from scipy.optimize import differential_evolution
import Lugs as lg

# --- 1. Fixed Variables & Materials  ---
TAU_MAX_F = 480e6
TAU_MAX_L = 300e6
SIGMA_Y_LUG = 300e6
ULTIMATE_BENDING_LUG = 450e6
RHO_LUG = 2700

# Push variables to Lug module 
lg.tau_max_f = TAU_MAX_F
lg.tau_max_l = TAU_MAX_L
lg.sigma_y_lug = SIGMA_Y_LUG
lg.ultimate_bending_lug = ULTIMATE_BENDING_LUG

F_X, F_Y, F_Z = 2740.6, 2740.6, 8222
T_X = 0.0002

param_bounds = {
    'n_l': (4, 10), 'n_f': (1, 4), 'w_1': (0.018, 0.3), 
    'w_2': (0.018, 0.15), 't_1': (0.0005, 0.025), 
    't_2': (0.0005, 0.15), 'l': (0.018, 0.2)
}
keys = list(param_bounds.keys())
bounds = [param_bounds[k] for k in keys]

def objective(x, penalize_factor=1e12):
    n_l_val, n_f_val, w1, w2, t1, t2, l_dim = x
    n_l, n_f = int(round(n_l_val)), int(np.ceil(n_f_val / 2) * 2)

    try:
        sm_f, d_out, _ = lg.fastener_diameter(F_X, F_Z, n_l, n_f, TAU_MAX_F)
        Ixx = lg.Ixx(w1, t1)
        _, _, sm_b, M_val = lg.bending(Ixx, w2, F_Z, t2, n_l)
        f_list = lg.CreateFastenerList_v(n_f, l_dim, w1, d_out)
        f_list = lg.ForceAtFastener_v(f_list, w2, n_f, n_l, F_Z, F_Y)
        f_list = lg.bearing_stress_v(f_list, t1, T_X, d_out)
        sm_p = lg.WorstShearMargin_v(f_list, TAU_MAX_L, 400e6)
        sm_c, _ = lg.length_lug(F_Z, t2, TAU_MAX_L, SIGMA_Y_LUG, n_l, w1, M_val, w2)
        sm_t = lg.tearout_vertical(f_list, l_dim, w1, t2)
        margins = [sm_f, sm_b, sm_p, sm_t, sm_c]
    except: return 1e20

    penalty = sum(np.square(max(0, -m)) for m in margins)
    vol_per_lug = ((w1 * l_dim - n_f * (d_out/2)**2 * np.pi) * t1) + \
                  ((w1 * w2 - n_f * (d_out/2)**2 * np.pi) * t2)
    return (vol_per_lug * n_l * RHO_LUG) + (penalty * penalize_factor)

if __name__ == '__main__':
    res = differential_evolution(objective, bounds, strategy='best1bin', popsize=15, seed=42)
    p = dict(zip(keys, res.x))
    
    n_l, n_f = int(round(p['n_l'])), int(np.ceil(p['n_f'] / 2) * 2)
    sm_f, d_out, m_size = lg.fastener_diameter(F_X, F_Z, n_l, n_f, TAU_MAX_F)
    Ixx = lg.Ixx(p['w_1'], p['t_1'])
    _, _, sm_b, M_f = lg.bending(Ixx, p['w_2'], F_Z, p['t_2'], n_l)
    f_list = lg.CreateFastenerList_v(n_f, p['l'], p['w_1'], d_out)
    f_list = lg.ForceAtFastener_v(f_list, p['w_2'], n_f, n_l, F_Z, F_Y)
    f_list = lg.bearing_stress_v(f_list, p['t_1'], T_X, d_out)
    sm_p = lg.WorstShearMargin_v(f_list, TAU_MAX_L, 400e6)
    sm_c, _ = lg.length_lug(F_Z, p['t_2'], TAU_MAX_L, SIGMA_Y_LUG, n_l, p['w_1'], M_f, p['w_2'])
    sm_t = lg.tearout_vertical(f_list, p['l'], p['w_1'], p['t_2'])

    mass = (((p['w_1']*p['l'] - n_f*(d_out/2)**2*np.pi)*p['t_1']) + \
            ((p['w_1']*p['w_2'] - n_f*(d_out/2)**2*np.pi)*p['t_2'])) * n_l * RHO_LUG

    # --- Print Converged Results ---
    print("-" * 40)
    print("CONVERGED STRUCTURAL DESIGN")
    print("-" * 40)
    print(f"Final Mass: {mass:.6f} kg")
    print(f"Fastener: {m_size}")
    print(f"Bolt d_out: {d_out*1000:.2f} mm")
    print(f"Margins: Fastener({sm_f:.2f}), Bending({sm_b:.2f}), Combined({sm_c:.2f}), Tearout({sm_t:.2f})")
    
    print("\nCONVERGED GEOMETRIC VALUES (mm):")
    print(f"w1 (Width 1):  {p['w_1']*1000:.2f} mm")
    print(f"w2 (Width 2):  {p['w_2']*1000:.2f} mm")
    print(f"t1 (Thick 1):  {p['t_1']*1000:.2f} mm")
    print(f"t2 (Thick 2):  {p['t_2']*1000:.2f} mm")
    print(f"l  (Length):   {p['l']*1000:.2f} mm")
    print(f"n_l (Lugs):    {n_l}")
    print(f"n_f (Fasten):  {n_f}")
    print("-" * 40)