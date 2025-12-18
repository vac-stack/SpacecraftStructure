import numpy as np
import cma

param_bounds = {
    'example_constraint':  (0.003, 0.080)
}

keys = ['example_constraint']
scales = [(param_bounds[k][0], param_bounds[k][1] - param_bounds[k][0]) for k in keys]

def unscale(x_scaled):
    vals = [low + x*rng for x,(low,rng) in zip(x_scaled, scales)]
    
    return vals

# ---------------------------
# Objective function (mass + pin mass) with constraints via penalty
# ---------------------------
def objective(x_scaled, penalize_factor=2e4, debug=False):
    
    # Get an educated guess for all needed variables 
    example_constraint = unscale(x_scaled)

    example_margin = max(0.0, -(0.8 - 1)) #0.8 represents division factor of some variable

    # Calculate total penalty from all margins
    penalty = 1.0 + penalize_factor * (1e4 * example_margin**2)

    total_mass = penalty * 30 # random value

    return total_mass

# ---------------------------
# CMA-ES wrapper
# ---------------------------
def optimize_cma(popsize=12, seed=0, maxiter=1000):
    rng = np.random.RandomState(seed)
    dim = len(keys)
    x0 = rng.rand(dim).tolist()
    sigma0 = 0.2
    opts = {'popsize':popsize,'bounds':[[0.0]*dim,[1.0]*dim],'CMA_active':True,'verb_disp':1,'maxiter':maxiter,'maxfevals':200000}
    es = cma.CMAEvolutionStrategy(x0, sigma0, opts)

    best = {'obj':float('inf'),'x':None}
    while True:
        sols = es.ask()
        fitnesses=[]
        for s in sols:
            val = objective(s)
            fitnesses.append(val)
            if val < best['obj']:
                best = {'obj':val,'x':np.array(s).copy()}
        es.tell(sols, fitnesses)
        es.disp()
        if es.stop(): break

    best_scaled = best['x']
    example_constraint = unscale(best_scaled)

    return {
        'scaled':best_scaled,
    }

# Main function
if __name__ == '__main__':
    print('Starting optimization (independent lug+pin materials)...')
    res = optimize_cma(popsize=14, seed=42, maxiter=200)
    print('\nBest result:')
    print(f" w = {res['w']:.6f} m")