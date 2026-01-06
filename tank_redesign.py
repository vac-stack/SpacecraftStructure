import numpy as np

# --------------------------
# Constants & vehicle data
# --------------------------

T = 25.0*2  # [N] two Aerojet GR-22 engines (used in the final design)
I_sp = 257  # [s] specific impulse (final monoprop design)
v_exhaust = I_sp * 9.81  # [m/s]
m_dot = T / v_exhaust    # [kg/s]

c_res = 0.2 # residual propellant/total propellant fraction
c_t = 0.2 # tank structure fraction from previous report

# orbiter
m_po = 304 # [kg] empty mass orbiter excluding tank 512.32s
c_tO = 0.05 # tank mass/propellant mass fraction

# orbit geometry
R_m = 2240e3   # central body radius (m)
h = 750e3      # altitude (m)
mu = 2.20319e13 # gravitational paramter Mercury
a0 = 91615e3   # initial semi-major axis (m)

# tank pressurerization
rho = 1470 # kg/m^3  density AF-M315E
P_0 = 310 * 10**5 # [Pa] initial tank pressure 
P_f = 27.6 * 10**5 # [Pa] min prop fed pressure 
R = 8.31446 # universal gas constant
M_N2 = 28.0134e-3 # molar mass nitrogen

# tank structure 
sigma_ult_carbon = 5.49e9 # ultimate tensile strenght T800 carbon fiber 
SF = 1.5 # safety factor 
V_f = 0.6 # fiber volume fraction
eta = 0.7 # translation efficiency
rho_carbon = 1810 # [kg/m^3] density T800 carbon fibers
rho_resin = 1250 # [kg/m^3] density resin
rho_Ti = 4470 # [kg/m^3] Ti6AlV titanium alloy density
sigma_yield_Ti = 1170e6 # [Pa] Ti6AlV max yield strengh
t_Ti = 2.6e-3 # [m] titanium liner thickness (taken from space shuttle gas tanks)

# maximum allowable stresses
sigma_allow = (sigma_ult_carbon * V_f * eta) / SF 
sigma_allow_Ti = sigma_yield_Ti / SF 

# calculated pressure vessel radius
R_gas = 0.203 + 9.202e-3

# -------------------------------------------------------
# Compute Δv requirements
# -------------------------------------------------------
r_init = R_m + h # pericenter/final orbit radius
V_elliptical = np.sqrt(2 * (mu / r_init - mu / (2 * a0)))
V_circular = np.sqrt(mu / r_init)
v_disturbance = 90 # [m/s] gravity, solar radiation pressure etc.
v_req_total = abs(V_elliptical - V_circular)


# Mass budget calculations (rewritten rocket equation)
m_fo = (m_po*((np.e**((v_req_total+v_disturbance)/v_exhaust))-1)) / (
    1 + c_res + c_tO - (c_res+c_tO) * np.e**((v_req_total+v_disturbance)/v_exhaust))  # [kg] propellant mass orbiter
mo = m_fo*(1 + c_tO + c_res) + m_po   # [kg] total mass orbiter
m_fo_2 = m_fo * (1 + c_res) # [kg] total orbiter fuel mass (with reserve)

# calculating required pressurization gas tank volume 
V_o = m_fo_2 / rho # [m^3] propellant tank volume
V_s = V_o / (P_0 / P_f - 1) # [m^3] gas tank volume
m_N2 = (P_0*V_s*M_N2)/(R*(273.15 + 27)) # [kg] nitrogen gas mass (assumed to be non compressible)


def find_vessel_radius(volume, height):
    """
    Calculates the radius of a cylindrical pressure vessel with hemispherical 
    endcaps given a fixed volume and total height.
    
    Equation: V = pi * r^2 * (h - (2/3) * r)
    Rearranged: (2/3)*pi*r^3 - pi*h*r^2 + V = 0
    
    Args:
        volume (float): The target internal volume.
        height (float): The total fixed height of the vessel.
        
    Returns:
        float: The radius r.
        None: If no valid physical solution exists (e.g., Volume is too large for the height).
    """
    if volume <= 0 or height <= 0:
        raise ValueError("Volume and height must be positive.")

    # Coefficients for the cubic equation: ar^3 + br^2 + cr + d = 0
    a = (2/3) * np.pi
    b = -np.pi * height
    c = 0
    d = volume

    # Find roots of the polynomial
    roots = np.roots([a, b, c, d])

    # Filter for physically valid roots:
    # 1. Must be a real number
    # 2. Must be positive
    # 3. Must satisfy 2*r <= h (The two hemispherical caps cannot exceed total height)
    valid_roots = []
    for r in roots:
        if np.isreal(r):
            r_real = r.real
            # Use a small epsilon for float comparison tolerance
            if r_real > 0 and (2 * r_real <= height + 1e-9):
                valid_roots.append(r_real)

    if not valid_roots:
        print("No valid physical dimensions found for this volume/height.")
        return None
        
    # There should be exactly one valid solution in the physical range
    return valid_roots[0]



# required tank (cylindrical with spherical endcaps) inner radius (r_int) for fixed volume and total height
# returns the outer radius of the tank and empty mass
def tank_radius(h): 
    r_int = find_vessel_radius(V_o, h)
    print(f"For V={V_o:.3f}, h={h}: Calculated r = {r_int:.3f}")

    # monoprop system
    R_tank_prop = r_int # [m] fuel tank radius
    t_Ti_prop = (P_f*R_tank_prop)/(sigma_allow_Ti)
    V_Ti_prop = (4/3) * np.pi * ((R_tank_prop + t_Ti_prop)**3 - R_tank_prop**3) + (h-2*r_int)*np.pi*((r_int+t_Ti_prop)**2-r_int**2)
    M_Ti_prop = V_Ti_prop * rho_Ti 


    return r_int , t_Ti_prop, M_Ti_prop


R_tank = ((3/4)*(V_s/np.pi))**(1/3) # [m] tank radius from volume
t_comp = (P_0*(R_tank + t_Ti))/(sigma_allow*V_f) # [m] required wall thickness using hoop stress relation
V_comp = (4/3) * np.pi * ((R_tank + t_comp + t_Ti)**3 - (R_tank + t_Ti)**3) # [m^3] composite volume
V_Ti = (4/3) * np.pi * ((R_tank + t_Ti)**3 - R_tank**3) # [m^3] titanium volume from tank geometry and fixed liner thickness
M_comp = V_comp * V_f * rho_carbon + V_comp * (1 - V_f) * rho_resin # [kg] mass composite
M_Ti = V_Ti * rho_Ti # [kg] mass titanium liner
M_tot = M_comp + M_Ti # [kg] total mass

h=1.3445
R_tank_prop,t_Ti_prop , M_Ti_prop = tank_radius(h)

# -------------------------------------------------------
# Outputs
# -------------------------------------------------------

print(f"\n--- Mass breakdown ---")
print(f"Propellant mass orbiter: {m_fo_2:.2f} kg for Δv: {(v_req_total+v_disturbance):.2f} m/s")

print(f"\n--- Characteristics propellant tank ---")
print(f"\n---V1---")
print(f"Tank volume: {V_o:.4f} m^3")
print(f"Radius prop tank: {R_tank_prop:.4f} m")
print(f"Tank wall thickness: {t_Ti_prop*10**3:.3f} mm")
print(f"Tank mass: {M_Ti_prop:.3f} kg")

print(f"\n--- Characteristics pressure vessel ---")
print(f"Pressurization gas vessel volume: {V_s:.4f} m^3")
print(f"Gas vessel radius: {R_tank + t_Ti + t_comp:.4f} m")
print(f"Composite wall thickness: {(t_comp*10**3):.3f} mm")
print(f"Total gas vessel wall thickness: {((t_comp+t_Ti)*10**3):.3f} mm")
print(f"Total gas vessel mass: {M_tot:.3f} kg")
print(f"Gas mass: {m_N2:.3f} kg")
print(f"Total gas vessel mass: {(M_tot+m_N2):.3f} kg")




