import numpy as np
import math

E = 5  # [Pa] shell Young's modulus
R = 5  # [m] shell cross-section diameter
L = 5  # [] shell length
t1 = 5  # [] shell thickness
p = 5  # [Pa] pressure diff. inside / outside shell
v = 5  # [-] Poisson's ratio
A = np.pi*(R**2 - (R-t1)**2)  # [m^2] shell cross-section area
I = np.pi/4 *(R**4 - (R-t1)**4)  # [m^4] shell cross-section mom. of inertia
shell_density = 5  # [kg/m^3]

#mass
mass_shell = A * L * shell_density

#column buckling euler
col_buck = np.pi**2 * E * I / A / L**2

#shell buckling
#lambda optimization
halfwaves_unrounded = (12 / np.pi**4 * L**4 / R**2 / t1**2 * (1-v**2))**0.5
halfwaves_low = math.floor(halfwaves_unrounded)
halfwaves_high = math.ceil(halfwaves_unrounded)
k_low = halfwaves_low + 12 / np.pi**4 * L**4 / R**2 / t1**2 * (1-v**2) / halfwaves_low
k_high = halfwaves_high + 12 / np.pi**4 * L**4 / R**2 / t1**2 * (1-v**2) / halfwaves_high
if k_low < k_high: halfwaves = halfwaves_low
elif k_low >= k_high: halfwaves = halfwaves_high

k = halfwaves + 12 / np.pi**4 * L**4 / R**2 / t1**2 * (1-v**2) / halfwaves
Q = p * R**2 / E / t1**2
shell_buck = (1.983 - 0.983 * np.e**(-23.14*Q)) * k * np.pi**2 * E * t1**2 / L**2 /12 / (1-v**2)    

print("German fella69")
print("Mass of shell [kg]: ", mass_shell)

