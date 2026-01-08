### finding the optimum shape, dimensions and weight of the sandwich plate

import math as math

R = 0.2152 #[m] radius of the load bearing structure
W = 0.55 #[m] length of one side of the reaction wheel (limiting dimension)

rho_nomex = 48.2 #[kg/m^3] density of nomex
rho_fabric = 1611 #[kg/m^3] density of fabric 
rho_panel = 1156 #[kg/m^3] density of the closing panels

n_s = 2 #[-] number of sandwich panel
n_cT = 1 #[-] number of closing panel top and/or bottom
h = 1.35 #[-] height of closing panel

#mass of sandwich per sqm
m_sand = 2*rho_fabric*0.19805/1000+rho_nomex*15/1000

#calcualtions of different geometries
def hex(R, W, n_s, n_cT, m_sand, rho_panel):

    h_hex = 2/math.sqrt(3)*(R+W)
    perimeter = 6*h_hex

    area_sandwich = (3*math.sqrt(3)/2*h_hex) - (math.pi*R**2)
    mass_sandwich = n_s * m_sand * area_sandwich

    area_closing_top = 3*math.sqrt(3)/2*h_hex
    area_closing_side = perimeter * h

    mass_top_bottom = area_closing_top * n_cT * rho_panel * 0.004
    mass_side = area_closing_side * rho_panel * 0.004

    total_mass = mass_sandwich + mass_top_bottom + mass_side

    return h_hex, area_sandwich, mass_sandwich, mass_top_bottom, mass_side, total_mass


def circ(R, W, n_s, n_cT, m_sand, rho_panel):

    r_circ = math.sqrt((W/2)**2 + (W+R)**2)
    perimeter = 2*math.pi*r_circ

    area_sandwich = math.pi*r_circ**2 - math.pi*R**2
    mass_sandwich = n_s * m_sand * area_sandwich

    area_closing_top = math.pi*r_circ**2
    area_closing_side = perimeter * h

    mass_top_bottom = area_closing_top * n_cT * rho_panel * 0.004
    mass_side = area_closing_side * rho_panel * 0.004

    total_mass = mass_sandwich + mass_top_bottom + mass_side

    return r_circ, area_sandwich, mass_sandwich, mass_top_bottom, mass_side, total_mass


def square(R, W, n_s, n_cT, m_sand, rho_panel):

    w_sq = (W*math.sqrt(2) + R) * math.sqrt(2)
    perimeter = 4*w_sq

    area_sandwich = w_sq**2 - math.pi*R**2
    mass_sandwich = n_s * m_sand * area_sandwich

    area_closing_top = w_sq**2
    area_closing_side = perimeter * h

    mass_top_bottom = area_closing_top * n_cT * rho_panel * 0.004
    mass_side = area_closing_side * rho_panel * 0.004

    total_mass = mass_sandwich + mass_top_bottom + mass_side

    return w_sq, area_sandwich, mass_sandwich, mass_top_bottom, mass_side, total_mass

h_hex, area_hex, m_s_hex, m_tb_hex, m_side_hex, m_hex = hex(R, W, n_s, n_cT, m_sand, rho_panel)
r_circ, area_circ, m_s_circ, m_tb_circ, m_side_circ, m_circ = circ(R, W, n_s, n_cT, m_sand, rho_panel)
w_sq, area_sq, m_s_sq, m_tb_sq, m_side_sq, m_sq = square(R, W, n_s, n_cT, m_sand, rho_panel)

def compare_geom():
    geometries = [
        ("hex", m_hex, h_hex),
        ("circle", m_circ, r_circ),
        ("square", m_sq, w_sq)
    ]
    
    geom, min_mass, length = min(geometries, key=lambda x: x[1])
    max_mass = max(geometries, key=lambda x: x[1])[1]
    weight_saved = max_mass - min_mass
    
    return geom, min_mass, length, weight_saved

geom, min_mass, length, weight_saved = compare_geom()

print(f"Minimum mass geometry: {geom}")
print(f"Total mass: {min_mass:.3f} kg")
print(f"Corresponding length: {length:.3f} m")
print(f"Weight saved: {weight_saved:.3f} kg\n")

print("---- Mass breakdown ----")
print(f"Hex sandwich mass: {m_s_hex:.3f} kg")
print(f"Hex top/bottom panel mass: {m_tb_hex:.3f} kg")
print(f"Hex side panel mass: {m_side_hex:.3f} kg\n")

print(f"Circle sandwich mass: {m_s_circ:.3f} kg")
print(f"Circle top/bottom panel mass: {m_tb_circ:.3f} kg")
print(f"Circle side panel mass: {m_side_circ:.3f} kg\n")

print(f"Square sandwich mass: {m_s_sq:.3f} kg")
print(f"Square top/bottom panel mass: {m_tb_sq:.3f} kg")
print(f"Square side panel mass: {m_side_sq:.3f} kg")

print("\nSandwich mass per m²:", m_sand, "kg/m²")
