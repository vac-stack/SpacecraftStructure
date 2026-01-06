import numpy as np
import cma
import Lugs as lg

# --- 1. Fixed Variables & Materials  ---
TAU_MAX_F = 480e6             # shear yield stress fastener [Pa] 
TAU_MAX_L = 300e6             # shear yield stress lug [Pa] 
ULTIMATE_BENDING_LUG = 450e6  # ultimate stress of lug [Pa] 
RHO_LUG = 2700                # Lug density (kg/m^3) 

# --- 2. Determining Loads  ---
F_X = 100.0  
F_Y = 100.0  
F_Z = 300.0 
T_X = 0.0002   # Thickness of the wall the lugs attach to 

# --- 3. Optimization Search Space (R Removed)  ---
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
scales = [(param_bounds[k][0], param_bounds[k][1] - param_bounds[k][0]) for k in keys]

def unscale(x_scaled):
    return [low + x * rng for x, (low, rng) in zip(x_scaled, scales)]

def objective(x_scaled, penalize_factor=1e4):
    # Unpack variables (Removed R) 
    n_l_val, n_f_val, w1, w2, t1, t2, l_dim = unscale(x_scaled)
    
    n_l = int(round(n_l_val))
    # 1. Enforce even number of fasteners
    n_f = int(np.ceil(n_f_val / 2) * 2)

    try:
        # 1. Fastener Diameter 
        sm_fastener, d_out, m_size = lg.fastener_diameter(F_Z, n_l, n_f, TAU_MAX_F)

        # 2. Bending Margin 
        Ixx, _ = lg.Ixx_zz(w1, t1)
        _, _, sm_bending = lg.bending(Ixx, w2, F_Z, t2, n_l)

        # 3. Pull-through and Plate Shear 
        f_list = lg.CreateFastenerList(n_f * 2, l_dim, w1, d_out)
        f_list = lg.ForceAtFastener(f_list, w2, n_l)
        f_list = lg.ShearInPlates(f_list, t1, T_X, d_out)
        sm_plates = lg.WorstShearMargin(f_list, TAU_MAX_L, 400e6)

        # 4. Tear-out Margin 
        sm_tearout = lg.tearout(f_list, w1, l_dim, t2)
        
        margins = [sm_fastener, sm_bending, sm_plates, sm_tearout]
        
    except Exception:
        return 1e20 

    # Quadratic penalty for any margin < 0 
    penalty = sum(np.square(max(0, -m)) for m in margins)

    # Mass Calculation (True Mass) [cite: 1, 2]
    vol_per_lug = ((w1 * l_dim - n_f * ((d_out*10**(-3))/2)**2 * np.pi)  * t1) + ((w1 * w2 - n_f * 2 * ((d_out*10**(-3))/2)**2 * np.pi) * t2) 
    total_mass = (vol_per_lug * n_l * RHO_LUG)
    
    return total_mass + (penalty * penalize_factor)

def run_optimization(popsize=20, maxiter=300):
    dim = len(keys)
    x0 = [0.5] * dim 
    sigma0 = 0.3 
    
    opts = {'popsize': popsize, 'bounds': [[0.0]*dim, [1.0]*dim], 'maxiter': maxiter, 'verb_disp': 0}
    
    es = cma.CMAEvolutionStrategy(x0, sigma0, opts)
    es.optimize(objective)
    
    res_scaled = es.result.xbest
    return dict(zip(keys, unscale(res_scaled)))

if __name__ == '__main__':
    print("Searching for Global Optimum...")
    
    best_overall_params = None
    min_score = float('inf')
    
    num_runs = 5
    for i in range(num_runs):
        params = run_optimization(popsize=24)
        score = objective([(params[k]-param_bounds[k][0])/(param_bounds[k][1]-param_bounds[k][0]) for k in keys])
        if score < min_score:
            min_score = score
            best_overall_params = params
        print(f"Run {i+1}/{num_runs} complete.")

    # --- FINAL CLEAN CALCULATION (No Penalties) ---
    p = best_overall_params
    n_l = int(round(p['n_l']))
    # Enforce even number for final results
    n_f = int(np.ceil(p['n_f'] / 2) * 2)

    # Get final fastener details and margins 
    sm_f, d_out, m_size = lg.fastener_diameter(F_Z, n_l, n_f, TAU_MAX_F)
    """
    Ixx, _ = lg.Ixx_zz(p['w_1'], p['t_1'])
    _, _, sm_b = lg.bending(Ixx, p['w_2'], F_Z, p['t_2'], n_l)
    f_list = lg.CreateFastenerList(n_f * 2, p['l'], p['w_1'], d_out)
    f_list = lg.ForceAtFastener(f_list, p['w_2'], n_l)
    f_list = lg.ShearInPlates(f_list, p['t_1'], T_X, d_out)
    sm_p = lg.WorstShearMargin(f_list, TAU_MAX_L, 400e6)
    sm_t = lg.tearout(f_list, p['w_1'], p['l'], p['t_2'])"""

    # Calculate physical mass [cite: 1, 2]
    vol_single = ((p['w_1'] * p['l'] - n_f * ((d_out*10**(-3))/2)**2 * np.pi) * p['t_1']) + ((p['w_1'] * p['w_2'] - n_f * ((d_out*10**(-3))/2)**2 * np.pi) * p['t_2'])
    true_mass_total = vol_single * n_l * RHO_LUG

    print("\n" + "="*30)
    print("GLOBAL OPTIMUM FOUND")
    print("="*30)
    for k, v in best_overall_params.items():
        if k == 'n_f':
            print(f"{k:4}: {n_f}") # Print the even number
        else:
            fmt = ".0f" if k == 'n_l' else ".6f"
            print(f"{k:4}: {v:{fmt}}")
    
    print("\n" + "="*30)
    print("FINAL DESIGN SPECS")
    print("="*30)
    print(f"True Total Mass : {true_mass_total:.6f} kg")
    print(f"Mass per Lug    : {true_mass_total/n_l:.6f} kg")
    print(f"Fastener Size   : {m_size}")
    print(f"Fastener d_out  : {d_out*1000:.4f} mm")
    """
    print("\n" + "="*30)
    print("FINAL MARGINS")
    print("="*30)
    print(f"Shear Fastener  : {sm_f:10.4f}")
    print(f"Bending Lug     : {sm_b:10.4f}")
    print(f"Shear Plates    : {sm_p:10.4f}")
    print(f"Tear-out        : {sm_t:10.4f}")
    print("="*30)"""