from tank_redesign import tank_radius
from dimension import compare_geom
from spacecraft_mass import *
from L_shaped_opt import *


# Define overall height of the satellite (and closing panels)
h = 1.35 # m

# Max mass of components on a panel
m_components_max = 4 # kg

# Main function
if __name__ == '__main__':
    
    # Tank Calculation
    # ---------------------------------------------------------------------------------
    # Get tank properties
    R_tank_prop, t_Ti_prop, M_Ti_prop = tank_radius(h)
    
    # Add margin to structure radius so the propelant tank can fit
    R_struct = R_tank_prop + 0.01   # Add 1 cm

    # Sandwich panel geometry calculation
    # ---------------------------------------------------------------------------------
    # Calculate properties of the sandwich panels
    # geometry type, mass, parameter length (hexagon side, circle radius or rectangle side) and weight saved
    geom, min_mass, param_length, weight_saved = compare_geom(h)

    # Structure optimization
    # ---------------------------------------------------------------------------------
    inputs = DesignInputs(
        m_fixed=801.32,        # kg (everything except main structural shell so payload and structure)
        p_internal=5.0e5,       # Pa (5 bar)
        g_load=9*9.80665,     # m/s^2 ( 6g axial with margin of 1.5 so 9g)
        fos_buckling=1.5,
        fos_yield=1.5,
        t_min=0.35e-3,           # 0.35 mm guess
        R_bounds=(R_struct, 2 * R_struct),  # m +100mm for attachment
        L_bounds=(h, 2.5),  # m
        R_required_min= R_struct,    # m (from tank packaging)
        L_required_min= h,    # m (from tank packaging)
    )

    result = mass_converging_optimizer(
        inputs=inputs,
        material_names=["Al7075-T6", "Ti-6Al-4V", "CFRP"],
        m0=801.32,      # initial total mass guess [kg]
        eps=0.02,     # convergence tolerance [kg]
        max_outer_iter=30
    )

    if result["success"]:
        sol = result["solution"]
        print(f"Material: {sol['mat']}")
        print(f"R={sol['R']:.4f} m, t1={sol['t1']*1e3:.3f} mm, L={sol['L']:.4f} m")
        print(f"Shell mass: {sol['m_shell']:.3f} kg")
        print(f"Total mass: {sol['m_total']:.3f} kg")

    if result["success"]:
        sol = result["solution"]
        mat = MATERIALS[sol["mat"]]

        R = sol["R"]
        t1 = sol["t1"]
        L = sol["L"]
        m_total = sol["m_total"]

        # ----------------------------
        # Recompute stresses
        # ----------------------------
        sigma_app = axial_compressive_stress(
            m_total, R, t1, inputs.g_load
        )

        sigma_hoop = hoop_stress_thin_cyl(
            inputs.p_internal, R, t1
        )

        # ----------------------------
        # Buckling allowables
        # ----------------------------
        sigma_eu = sigma_cr_euler(R, t1, L, mat.E) / inputs.fos_buckling
        sigma_sh = sigma_cr_shell(R, t1, L, mat, inputs.p_internal) / inputs.fos_buckling

        sigma_buckling_allow = min(sigma_eu, sigma_sh)

        # ----------------------------
        # Yield allowable
        # ----------------------------
        sigma_yield_allow = mat.sigma_y / inputs.fos_yield

        # ----------------------------
        # Margins of safety
        # ----------------------------
        mos_buckling = sigma_buckling_allow / sigma_app - 1.0
        mos_yield_axial = sigma_yield_allow / sigma_app - 1.0
        mos_yield_hoop = sigma_yield_allow / sigma_hoop - 1.0

        mos_yield = min(mos_yield_axial, mos_yield_hoop)

        # ----------------------------
        # Prints
        # ----------------------------
        print(result["message"])
        print(f"Material: {sol['mat']}")
        print(f"R = {R:.4f} m")
        print(f"t = {t1*1e3:.3f} mm")
        print(f"L = {L:.4f} m")
        print(f"Shell mass: {sol['m_shell']:.3f} kg")
        print(f"Total mass: {m_total:.3f} kg")

        print("\n--- Stress summary ---")
        print(f"Applied axial stress: {sigma_app/1e6:.2f} MPa")
        print(f"Hoop stress: {sigma_hoop/1e6:.2f} MPa")

        print("\n--- Allowables (with FoS = 1.5) ---")
        print(f"Euler buckling allowable: {sigma_eu/1e6:.2f} MPa")
        print(f"Shell buckling allowable: {sigma_sh/1e6:.2f} MPa")
        print(f"Yield allowable: {sigma_yield_allow/1e6:.2f} MPa")

        print("\n--- Margins of Safety ---")
        print(f"Buckling MoS: {mos_buckling:.3f}")
        print(f"Yield MoS (axial): {mos_yield_axial:.3f}")
        print(f"Yield MoS (hoop): {mos_yield_hoop:.3f}")
        print(f"Governing Yield MoS: {mos_yield:.3f}")

    else:
        print("No solution.")

    # Lug optimization
    # ---------------------------------------------------------------------------------
    OptimizeLug()