import numpy as np
import cma
import Lugs as lg

# --- 1. Fixed Variables & Materials  ---
TAU_MAX_F = 480e6             # shear yield stress fastener [Pa]
TAU_MAX_L = 300e6             # shear yield stress lug [Pa]
ULTIMATE_BENDING_LUG = 450e6  # ultimate stress of lug [Pa]
RHO_LUG = 2700                # Lug density (kg/m^3)

# --- 2. Determining Loads  ---
# These should be determined before optimization as per your instructions
F_X = 100.0  
F_Y = 100.0  
F_Z = 300.0 
T_X = 0.0002   # Thickness of the wall the lugs attach to 

# --- 3. Optimization Search Space  ---
param_bounds = {
    'R':   (0.1, 0.5),      # Radius [m]
    'n_l': (2, 12),         # number of lugs
    'n_f': (1, 8),          # number of fasteners
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

def objective(x_scaled, penalize_factor=1e10):
    # Unpack variables based on Source 2 
    R, n_l_val, n_f_val, w1, w2, t1, t2, l_dim = unscale(x_scaled)
    
    n_l = int(round(n_l_val))
    n_f = int(round(n_f_val))

    try:
        # 1. Fastener Diameter 
        sm_fastener, d_out, m_size = lg.fastener_diameter(F_Z, n_l, n_f, TAU_MAX_F)

        # 2. Bending Margin 
        # Calling function from Lugs.py [cite: 1]
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
        return 1e20 # Geometry is mathematically impossible

    # Quadratic penalty for any margin < 0 
    penalty = sum(np.square(max(0, -m)) for m in margins)

    # Mass Calculation 
    # L-shape approximation: Bottom plate + Vertical plate
    vol_per_lug = (w1 * l_dim * t2) + (w1 * w2 * t1)
    total_mass = (vol_per_lug * n_l * RHO_LUG)
    
    # If any margin is failed, the penalty will dominate the mass 
    return total_mass + (penalty * penalize_factor)

def run_optimization(popsize=20, maxiter=300):
    """Executes CMA-ES and returns the best found dictionary."""
    dim = len(keys)
    x0 = [0.5] * dim 
    sigma0 = 0.3 # Higher initial sigma for better global exploration
    
    opts = {
        'popsize': popsize,
        'bounds': [[0.0]*dim, [1.0]*dim],
        'maxiter': maxiter,
        'verb_disp': 0 # Set to 1 for live updates
    }
    
    es = cma.CMAEvolutionStrategy(x0, sigma0, opts)
    es.optimize(objective)
    
    res_scaled = es.result.xbest
    fitness = es.result.fbest
    res_final = unscale(res_scaled)
    
    return dict(zip(keys, res_final)), fitness

if __name__ == '__main__':
    print("Searching for Global Optimum...")
    
    best_overall_params = None
    best_overall_mass = float('inf')
    
    # Multirun strategy to avoid local optima
    num_runs = 5
    for i in range(num_runs):
        params, final_mass = run_optimization(popsize=24)
        if final_mass < best_overall_mass:
            best_overall_mass = final_mass
            best_overall_params = params
        print(f"Run {i+1}/{num_runs} complete. Current best mass: {best_overall_mass:.4f} kg")

    print("\n--- GLOBAL OPTIMUM FOUND ---")
    for k, v in best_overall_params.items():
        fmt = ".0f" if k in ['n_l', 'n_f'] else ".6f"
        print(f"{k}: {v:{fmt}}")