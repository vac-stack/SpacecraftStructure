import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple, Callable, List
from scipy.optimize import minimize

# ----------------------------
# Materials
# ----------------------------
@dataclass(frozen=True)
class Material:
    name: str
    E: float      # Young's modulus [Pa]
    nu: float     # Poisson's ratio [-]
    rho: float    # density [kg/m^3]
    sigma_y: float  # yield [Pa]

MATERIALS: Dict[str, Material] = {
    "Al7075-T6": Material("Al7075-T6", E=70e9,  nu=0.33, rho=2810, sigma_y=5.03e8),
    "Ti-6Al-4V": Material("Ti-6Al-4V", E=114e9,  nu=0.34	, rho=4420, sigma_y=8.28e8),
    "CFRP": Material("CFRP", E=200e9,  nu=0.9	, rho=1600, sigma_y=1.2e9),
    # have to change to the correct material
}
# ----------------------------
# the geometry and the mass models
# ----------------------------
            #mass
def shell_surface_area(R: float, L: float) -> float:
    # cylinder lateral area [m^2]
    return 2.0 * np.pi * R * L

def shell_mass(R: float, t1: float, L: float, mat: Material) -> float:
    # thin-walled cylinder approximation [kg]
    return shell_surface_area(R, L) * t1 * mat.rho


            #geometry
def ring_cross_section_area(R: float, t1: float) -> float:
    # Thin-wall cross-section area for axial load: A ≈ 2πR t
    return 2.0 * np.pi * R * t1

def second_moment_area_thin_ring(R: float, t1: float) -> float:
    # Thin-walled ring bending stiffness about centroidal axis:
    # I ≈ π R^3 t
    return np.pi * (R ** 3) * t1

# ----------------------------
# Buckling models
# ----------------------------
def sigma_cr_euler(R: float, t1: float, L: float, E: float) -> float:
    """
    Eq (4.1): sigma_cr = (pi^2 E I) / (A L^2)
    using thin ring approximations for I and A.
    """
    A = ring_cross_section_area(R, t1)
    I = second_moment_area_thin_ring(R, t1)
    return (np.pi**2 * E * I) / (A * L**2)

def Q_param(p: float, R: float, t1: float, E: float) -> float:
    # Eq (4.3): Q = (p/E) * (R/t1)^2
    return (p / E) * (R / t1)**2

def k_param(L: float, R: float, t1: float, nu: float, lam: float) -> float:
    # Eq (4.4): k = lam + (12/pi^4) * (L^4/(R^2 t^2))*(1-nu^2)*(1/lam)
    return lam + (12.0 / np.pi**4) * ((L**4) / (R**2 * t1**2)) * (1.0 - nu**2) * (1.0 / lam)

def sigma_cr_shell(R: float, t1: float, L: float, mat: Material, p: float,
                   lam_grid: np.ndarray = None) -> float:
    """
    Eq (4.2) with minimization over lambda (lam) as requested in the text.
    """
    if lam_grid is None:
        # Reasonable search grid (you can refine)
        lam_grid = np.linspace(0.2, 20.0, 400)

    E, nu = mat.E, mat.nu
    Q = Q_param(p, R, t1, E)
    # bracket term: [1.983 - 0.983*exp(-23.14 Q)]
    bracket = 1.983 - 0.983 * np.exp(-23.14 * Q)

    # compute k(lam) and pick min k (or rather, use k that minimizes eq 4.4 contribution)
    k_vals = np.array([k_param(L, R, t1, nu, lam) for lam in lam_grid])
    k_min = np.min(k_vals)

    return bracket * k_min * (np.pi**2 * E) / (12.0 * (1.0 - nu**2)) * (t1 / L)**2

# ----------------------------
# Loads + constraints
# ----------------------------
@dataclass
class DesignInputs:
    m_fixed: float          # mass of everything except main shell [kg]
    p_internal: float       # pressure difference in shell [Pa]
    g_load: float           # launch axial acceleration [m/s^2] (or g*n)
    fos_buckling: float     # factor of safety for buckling
    fos_yield: float        # factor of safety for yield (optional)
    t_min: float            # min thickness [m]
    R_bounds: Tuple[float, float]
    L_bounds: Tuple[float, float]
    # constraints from tanks (example)
    R_required_min: float   # minimum radius to fit tanks [m]
    L_required_min: float   # minimum length to fit tanks [m]

def axial_compressive_stress(m_total: float, R: float, t1: float, g_load: float) -> float:
    # Conservative: total compressive force F = m_total * g_load
    F = m_total * g_load
    A = ring_cross_section_area(R, t1)
    return F / A

def hoop_stress_thin_cyl(p: float, R: float, t1: float) -> float:
    # sigma_hoop = pR/t
    return p * R / t1

def optimize_shell_for_mass_guess(m_guess: float,
                                 inputs: DesignInputs,
                                 mat: Material) -> Dict:
    """
    For a fixed m_guess (used in loads), find R, t1, L that minimize shell mass
    while satisfying buckling + (optional) yield + packaging constraints.
    """

    # Decision variables: x = [R, t1, L]
    # Initial guess: mid-bounds
    x0 = np.array([
        max(inputs.R_required_min, 0.5*(inputs.R_bounds[0] + inputs.R_bounds[1])),
        max(inputs.t_min, 1e-3),
        max(inputs.L_required_min, 0.5*(inputs.L_bounds[0] + inputs.L_bounds[1])),
    ])

    bounds = [
        (max(inputs.R_bounds[0], inputs.R_required_min), inputs.R_bounds[1]),
        (inputs.t_min, None),  # upper thickness not bounded here; you can add if needed
        (max(inputs.L_bounds[0], inputs.L_required_min), inputs.L_bounds[1]),
    ]

    def objective(x):
        R, t1, L = x
        return shell_mass(R, t1, L, mat)

    def constraints_fun(x):
        R, t1, L = x
        # Applied axial stress from m_guess
        sigma_app = axial_compressive_stress(m_guess, R, t1, inputs.g_load)

        # Buckling allowables
        sigma_eu = sigma_cr_euler(R, t1, L, mat.E) / inputs.fos_buckling
        sigma_sh = sigma_cr_shell(R, t1, L, mat, inputs.p_internal) / inputs.fos_buckling

        # Optional yield (axial + hoop combined crudely)
        sigma_hoop = hoop_stress_thin_cyl(inputs.p_internal, R, t1)
        sigma_allow = mat.sigma_y / inputs.fos_yield

        # Inequality constraints are defined as >= 0 for SLSQP
        return np.array([
            sigma_eu - sigma_app,        # Euler buckling safe
            sigma_sh - sigma_app,        # shell buckling safe
            sigma_allow - sigma_app,     # axial yield safe (optional)
            sigma_allow - sigma_hoop,    # hoop yield safe (optional)
        ])

    cons = [{"type": "ineq", "fun": lambda x, i=i: constraints_fun(x)[i]} for i in range(4)]

    res = minimize(objective, x0, method="SLSQP", bounds=bounds, constraints=cons,
                   options={"ftol": 1e-9, "maxiter": 300, "disp": False})

    if not res.success:
        return {
            "success": False,
            "message": res.message,
            "x": res.x,
            "m_shell": np.nan,
        }

    R, t1, L = res.x
    return {
        "success": True,
        "x": res.x,
        "m_shell": shell_mass(R, t1, L, mat),
        "R": R, "t1": t1, "L": L,
        "mat": mat.name,
    }

def mass_converging_optimizer(inputs: DesignInputs,
                             material_names: List[str],
                             m0: float,
                             eps: float = 0.01,
                             max_outer_iter: int = 50) -> Dict:
    """
    Outer loop: update total mass until convergence.
    Inner loop: optimize shell size for that mass guess.
    Material is handled by trying each candidate and picking the lowest total mass.
    """
    m_prev = float(m0)
    history = []

    for it in range(max_outer_iter):
        best = None

        # Try each material (discrete choice)
        for name in material_names:
            mat = MATERIALS[name]
            cand = optimize_shell_for_mass_guess(m_prev, inputs, mat)
            if not cand["success"]:
                continue

            m_total = inputs.m_fixed + cand["m_shell"]
            cand["m_total"] = m_total

            if (best is None) or (m_total < best["m_total"]):
                best = cand

        if best is None:
            return {
                "success": False,
                "message": "No feasible design found for any material with given bounds/constraints.",
                "history": history
            }

        m_new = best["m_total"]
        history.append({"iter": it, "m_prev": m_prev, "m_new": m_new, **best})

        if abs(m_new - m_prev) < eps:
            return {
                "success": True,
                "message": "Converged.",
                "solution": best,
                "history": history
            }

        m_prev = m_new

    return {
        "success": False,
        "message": "Did not converge within max_outer_iter.",
        "history": history
    }

# ----------------------------
# TO BE FILLED WITH CORRECT VALUES
# ----------------------------
if __name__ == "__main__":
    inputs = DesignInputs(
        m_fixed=801.32,        # kg (everything except main structural shell so payload and structure)
        p_internal=5.0e5,       # Pa (5 bar)
        g_load=9.0*9.80665,     # m/s^2 ( 6g axial with margin of 1.5 so 9g)
        fos_buckling=1.5,
        fos_yield=1.5,
        t_min=0.5e-3,           # 0.5 mm guess
        R_bounds=(0.2152, 0.25),  # m
        L_bounds=(1.35, 2.5),  # m
        R_required_min=0.10,    # m (from tank packaging)
        L_required_min=1.35,    # m (from tank packaging)
    )

    result = mass_converging_optimizer(
        inputs=inputs,
        material_names=["Al7075-T6", "Ti-6Al-4V", "CFRP"],
        m0=15.0,      # initial total mass guess [kg]
        eps=0.02,     # convergence tolerance [kg]
        max_outer_iter=30
    )


    print(result["message"])
    if result["success"]:
        sol = result["solution"]
        print(f"Material: {sol['mat']}")
        print(f"R={sol['R']:.4f} m, t1={sol['t1']*1e3:.3f} mm, L={sol['L']:.4f} m")
        print(f"Shell mass: {sol['m_shell']:.3f} kg")
        print(f"Total mass: {sol['m_total']:.3f} kg")
    else:
        print("No solution.")