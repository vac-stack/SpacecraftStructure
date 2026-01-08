import math
from Lugs import *

shear_margin, D_fo, M_size = fastener_diameter(F_z, n_l, n_f, tau_max_f)

def mass_lug(D_fo, t_1, t_2, t_x, w_1, l, w_2, n_f):
    mass_head = (0.75*D_fo/0.866)**2*3*math.sqrt(3)/2*0.8*D_fo #assuming nut and head have the same thickness
    mass_shaft = math.pi*(D_fo/2)**2*(max(t_1,t_2)+D_fo+t_x)
    mass_fastener_total = n_f*(2*mass_head + mass_shaft - math.pi*(D_fo/2)**2*0.8*(D_fo))

    mass_total_lug = w_1*l*t_2+w_2*w_1*t_1 - 2*mass_fastener_total

    return mass_total_lug

