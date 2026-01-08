import numpy as np
from scipy.optimize import differential_evolution
import Lugs as lg

# --- 1. Fixed Variables & Materials ---
TAU_MAX_F = 480e6             # shear yield stress fastener [Pa] 
TAU_MAX_L = 300e6             # shear yield stress lug [Pa] 
ULTIMATE_BENDING_LUG = 450e6  # ultimate stress of lug [Pa] 
RHO_LUG = 2700                # Lug density (kg/m^3) 
KT = 3                        # Stress concentration factor

# --- 2. Determining Loads ---
F_X = 400 
F_Y = 400 
F_Z = 1200
T_X = 0.0002   # Thickness of the wall the lugs attach to 

# --- 3. Optimization Search Space ---
param_bounds = {
    'n_l': (2, 10),         # number of lugs 
    'n_f': (1, 4),          # number of fasteners 
    'w_1': (0.015, 0.15),   # width of lug 
    'w_2': (0.015, 0.15),   # length of lug horizontal plate 
    't_1': (0.002, 0.025),  # thickness lug bottom 
    't_2': (0.002, 0.025),  # thickness lug top 
    'l':   (0.02, 0.2)      # length of lug horizontal plate 
}

keys = list(param_bounds.keys())
# SciPy differential_evolution expects a list of (min, max) tuples
scipy_bounds = [param_bounds[k] for k in keys]

def objective(x, penalize_factor=1e10):
    """Objective function tailored for SciPy optimizer."""
    # SciPy passes an array 'x' of current parameter values
    n_l_val, n_f_val, w1, w2, t1, t2, l_dim = x
    
    n_l = int(round(n_l_val))
    # Enforce even number of fasteners
    n_f = int(np.ceil(n_f_val / 2) * 2)

    try:
        # 1. Fastener Diameter
        sm_fastener, d_out, m_size = lg.fastener_diameter(F_Z, n_l, n_f, TAU_MAX_F)

        # 2. Bending Margin
        Ixx = lg.Ixx(w1, t1)
        _, _, sm_bending = lg.bending(Ixx, KT, w2, F_Z, t2, n_l)

        # 3. Pull-through and Plate Shear (Vertical as per Lug logic)
        f_list = lg.CreateFastenerList_v(n_f, l_dim, w1, d_out)
        f_list = lg.ForceAtFastener_v(f_list, w2, n_f, n_l, F_Z, F_Y)
        f_list = lg.bearing_stress_v(f_list, t1, T_X, d_out)
        sm_plates = lg.WorstShearMargin_v(f_list, TAU_MAX_L, 400e6)

        # 4. Tear-out Margin
        sm_tearout = lg.tearout_vertical(f_list, l_dim, w1, t2)
        
        margins = [sm_fastener, sm_bending, sm_plates, sm_tearout]
        
    except Exception:
        return 1e20 

    # Quadratic penalty for any margin < 0 
    penalty = sum(np.square(max(0, -m)) for m in margins)

    # Mass Calculation (Preserved exactly as requested)
    vol_per_lug = ((w1 * l_dim - n_f * ((d_out)/2)**2 * np.pi) * t1) + \
                  ((w1 * w2 - n_f * 2 * ((d_out)/2)**2 * np.pi) * t2) 
    total_mass = (vol_per_lug * n_l * RHO_LUG)
    
    # Returning total mass with potential penalty
    return total_mass + (penalty * penalize_factor)

if __name__ == '__main__':
    print("Searching for Global Optimum using SciPy Differential Evolution...")
    
    # Using Differential Evolution for robust global searching 
    result = differential_evolution(
        objective, 
        scipy_bounds, 
        strategy='best1bin', 
        maxiter=1000, 
        popsize=15, 
        tol=0.01, 
        mutation=(0.5, 1), 
        recombination=0.7, 
        seed=42,
        polish=True # Polishes the result with a local minimizer
    )

    # Final best parameters found by SciPy
    best_params_array = result.x
    best_params = dict(zip(keys, best_params_array))
    
    # --- FINAL CLEAN CALCULATION ---
    n_l = int(round(best_params['n_l']))
    n_f = int(np.ceil(best_params['n_f'] / 2) * 2)
    w1, w2, t1, t2, l_dim = best_params['w_1'], best_params['w_2'], best_params['t_1'], best_params['t_2'], best_params['l']

    # Final margin evaluation
    sm_f, d_out, m_size = lg.fastener_diameter(F_Z, n_l, n_f, TAU_MAX_F)
    Ixx = lg.Ixx(w1, t1)
    _, _, sm_b = lg.bending(Ixx, w2, F_Z, t2, n_l)
    f_list = lg.CreateFastenerList_v(n_f, l_dim, w1, d_out)
    f_list = lg.ForceAtFastener_v(f_list, w2, n_f, n_l, F_Z, F_Y)
    f_list = lg.bearing_stress_v(f_list, t1, T_X, d_out)
    sm_p = lg.WorstShearMargin_v(f_list, TAU_MAX_L, 400e6)
    sm_t = lg.tearout_vertical(f_list, l_dim, w1, t2)

    # Physical Mass Calculation (User specified logic)
    vol_single = ((w1 * l_dim - n_f * (d_out/2)**2 * np.pi) * t1) + \
                 ((w1 * w2 - n_f * 2 * (d_out/2)**2 * np.pi) * t2)
    true_mass_total = vol_single * n_l * RHO_LUG

    print("\n" + "="*30)
    print("GLOBAL OPTIMUM FOUND")
    print("="*30)
    for k, v in best_params.items():
        if k == 'n_f':
            print(f"{k:4}: {n_f} (Even enforced)")
        elif k == 'n_l':
            print(f"{k:4}: {n_l}")
        else:
            print(f"{k:4}: {v:.6f}")
    
    print("\n" + "="*30)
    print("FINAL DESIGN SPECS")
    print("="*30)
    print(f"True Total Mass : {true_mass_total:.6f} kg")
    print(f"Mass per Lug    : {true_mass_total/n_l:.6f} kg")
    print(f"Fastener Size   : {m_size}")
    print(f"Fastener d_out  : {d_out*1000:.4f} mm")
    
    print("\n" + "="*30)
    print("FINAL MARGINS")
    print("="*30)
    print(f"Shear Fastener  : {sm_f:10.4f}")
    print(f"Bending Lug     : {sm_b:10.4f}")
    print(f"Shear Plates    : {sm_p:10.4f}")
    print(f"Tear-out        : {sm_t:10.4f}")
    print("="*30)