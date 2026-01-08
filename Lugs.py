import math
import numpy as np

### shit to optimise: n_l, n_f, w_2, t_2, w_1, t_1

# ---------------------------------------------------------------------------------------------------------------
# dont forget constrians: add minimum distance of holes from edge, minimum l, minimum w_2 - as function of d_fo
# -> minimum t_1, t_2
# no the code is not duplicated, they're for horizontal and vertical part, im too lazy to unify them
# i know its unnecessary long im sorry
# ---------------------------------------------------------------------------------------------------------------


tau_max_f = 1 #yield shear stress of fastener
tau_max_l = 1 #yield shear stress of lug

#stress concentration factor

ultimate_bending_lug = 1


#total forces like in 5.5 diagram **NOT from wp4** [N]
F_x = 1 
F_y = 1
F_z = 1

### fasteners calculations

#list for metric bolt [mm]
metric_bolt_d_i = [0.729, 0.829, 0.929, 1.075, 1.221, 1.421, 1.567, 1.713, 2.013, 2.459, 2.850, 3.242, 3.688, 4.134, 4.917, 5.917, 6.647] 
metric_bolt_name = ['M1', 'M1.1', 'M1.2', 'M1.4', 'M1.6', 'M1.8', 'M2', 'M2.2', 'M2.5', 'M3', 'M3.5', 'M4', 'M4.5', 'M5', 'M6', 'M7', 'M8']
metric_bolt_d_o = [1, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8]

def fastener_diameter(F_x, F_z, n_l, n_f, tau_max_f):

    F_s = F_z / (n_l * n_f)   # shear
    F_n = F_x / (n_l * n_f)   # axial

    D = math.sqrt(
        (4 / (math.pi * tau_max_f)) *
        math.sqrt(F_n**2 + 3*F_s**2)
    )

    
    # Pick the smallest bolt >= D
    for i, d in enumerate(metric_bolt_d_i):
        if d >= D*1000:
            D_metric_i = metric_bolt_d_i[i]*0.001       # inner diameter in meters
            D_metric_o = metric_bolt_d_o[i]*0.001
            M_size = metric_bolt_name[i]
            break
    else:
        raise ValueError("No metric bolt large enough")
    
    #actual stress
    A = math.pi*D_metric_i**2/4
    sigma = F_n / A
    tau = F_s/A
    sigma_vm = math.sqrt(sigma**2+3*tau**2)

    shear_margin = tau_max_f / sigma_vm - 1
    
    return shear_margin, D_metric_o, M_size


# pull through from wp4, moment and force modified, fastener constrians modified
# pull through vertical and horizontal is done in a separate part

### vertical pull through
def ForceAtFastener_v(fasteners, w_2,n_f,n_l,F_z, F_y):
    
    # Calculate sum of distances of fasteners to cg
    ri_sq_tot = 0
    for j in range(len(fasteners)):
        ri_sq_tot += fasteners[j][0]**2 + fasteners[j][1]**2

    # Iterate through each fastener
    for i in range(len(fasteners)):
        fasteners[i][2] = F_y/n_l/n_f # Add force in y direction
        ri = np.sqrt(fasteners[i][0]**2 + fasteners[i][1]**2) # Calculate distance to cg of this fastener
        
        M_x = F_z * w_2/2 # Moment
        F_due_moment = M_x*ri/ri_sq_tot # Force due to moment

        if M_x > 0:
            if fasteners[i][0] > 0:
                fasteners[i][2] += F_due_moment
            else:
                fasteners[i][2] -= F_due_moment
        else:
            if fasteners[i][0] > 0:
                fasteners[i][2] -= F_due_moment
            else:
                fasteners[i][2] += F_due_moment
        
    return(fasteners)

# Calculate shear stress in each fastener
def bearing_stress_v(fasteners, t_lug, t_sc, D_fo):
    
    # Contact area with fastener per surface
    area_lug = np.pi*D_fo*t_lug
    area_sc = np.pi*D_fo*t_sc

    # Calculate stress per fastener per surface
    for i in range(len(fasteners)):
        Fy = fasteners[i][2]
        fasteners[i][3] = Fy / area_lug
        fasteners[i][4] = Fy / area_sc
    
    return fasteners

def WorstShearMargin_v(fasteners, yield_lug, yield_sc):
    margins = []
    for i in range(len(fasteners)):
        margins.append(yield_lug/np.abs(fasteners[i][3])-1)
        margins.append(yield_sc/np.abs(fasteners[i][4])-1)
        
    return min(margins)

def CreateFastenerList_v(n_f, l, w_1, D_fo):
    array = []
    n_side = n_f // 2

    if n_side == 1:
        y_positions = [0.0]
    else:
        d = (l - 3*D_fo) / (n_side - 1)
        y_positions = [-(2*D_fo + i*d) + l/2 for i in range(n_side)]

    for y in y_positions:
        array.append([+w_1/2 - 2*D_fo, y, 0, 0, 0])
        array.append([-w_1/2 + 2*D_fo, y, 0, 0, 0])

    return array


def FastenersRanked_v(fasteners):
    ranked1 = sorted(fasteners, key=lambda bolt: abs(bolt[2]), reverse=True)
    ranked2 = [bolt[:3] for bolt in ranked1]
    return ranked2



###pull through horizontal

def ForceAtFastener_h(fasteners,n_f, n_l,F_z, F_y,l):
    
    # Calculate sum of distances of fasteners to cg
    ri_sq_tot = 0
    for j in range(len(fasteners)):
        ri_sq_tot += fasteners[j][0]**2 + fasteners[j][1]**2

    # Iterate through each fastener
    for i in range(len(fasteners)):
        fasteners[i][2] = F_z/n_l/n_f # Add force in y direction
        ri = np.sqrt(fasteners[i][0]**2 + fasteners[i][1]**2) # Calculate distance to cg of this fastener
        
        M_x = F_y * l/2 # Moment
        F_due_moment = M_x*ri/ri_sq_tot # Force due to moment

        if M_x > 0:
            if fasteners[i][0] > 0:
                fasteners[i][2] += F_due_moment
            else:
                fasteners[i][2] -= F_due_moment
        else:
            if fasteners[i][0] > 0:
                fasteners[i][2] -= F_due_moment
            else:
                fasteners[i][2] += F_due_moment
        
    return(fasteners)

# Calculate shear stress in each fastener
def bearing_stress_h(fasteners, t_lug, t_sc, D_fo):
    
    # Contact area with fastener per surface
    area_lug = np.pi*D_fo*t_lug
    area_sc = np.pi*D_fo*t_sc

    # Calculate stress per fastener per surface
    for i in range(len(fasteners)):
        Fy = fasteners[i][2]
        fasteners[i][3] = Fy / area_lug
        fasteners[i][4] = Fy / area_sc
    
    return fasteners

def WorstShearMargin_h(fasteners, yield_lug, yield_sc):
    margins = []
    for i in range(len(fasteners)):
        margins.append(yield_lug/np.abs(fasteners[i][3])-1)
        margins.append(yield_sc/np.abs(fasteners[i][4])-1)
        
    return min(margins)

def CreateFastenerList_h(n_f, w_2, w_1, D_fo):
    array = []
    n_side = n_f // 2

    if n_side == 1:
        y_positions = [0.0]
    else:
        d = (w_2 - 3*D_fo) / (n_side - 1)
        y_positions = [-(2*D_fo + i*d) + w_2/2 for i in range(n_side)]

    for y in y_positions:
        array.append([+w_1/2 - 2*D_fo, y, 0, 0, 0])
        array.append([-w_1/2 + 2*D_fo, y, 0, 0, 0])

    return array


def FastenersRanked_h(fasteners):
    ranked1 = sorted(fasteners, key=lambda bolt: abs(bolt[2]), reverse=True)
    ranked2 = [bolt[:3] for bolt in ranked1]
    return ranked2


### bending and shear of the lug

def Ixx(w_1,t_1):
    I_xx = max((t_1**3*w_1)/12.0, 1e-15) #Ixx of the horizontal polate, same for Iyy
    return I_xx

#bending horizontal plate
def bending(Ixx, w_2, F_z, t_2, n_l):
    M = F_z/n_l * (w_2/2)
    bending_applied = (M * t_2/2)/Ixx
    bending_margin = ultimate_bending_lug/ bending_applied-1
    return t_2, w_2, bending_margin

#shear vertical
def length_lug(F_z, kt,t_2, tau_max_l, n_l):
    length = (3*F_z*kt/n_l)/(8*t_2*tau_max_l)
    l = max(length, 0.01)
    shear_margin = (8*t_2*tau_max_l/3/l)/(F_z/n_l) - 1
    return l, shear_margin

### tear-outs 

def tearout_vertical(fasteners, l, w_1, t_2):
    min_distance = []
    for f in fasteners:
        x, y = f[0], f[1]
        min_distance.extend([
            w_1/2 + x,
            w_1/2 - x,
            l/2 + y,
            l/2 - y

        ])
    edge_distance = min(min_distance)

    F_max = max(abs(f[2]) for f in fasteners) 
    F_tear = 2* edge_distance*t_2*tau_max_l
    tearout_margin = F_tear/F_max-1

    return tearout_margin

def tearout_horizontal(fasteners, w_1, w_2, t_2):
    min_distance = []
    for f in fasteners:
        x, y = f[0], f[1]
        min_distance.extend([
            w_1/2 + x,
            w_1/2 - x,
            w_2/2 + y,
            w_2/2 - y
        ])
    edge_distance = min(min_distance)

    F_max = max(abs(f[2]) for f in fasteners) 
    F_tear = 2* edge_distance*t_2*tau_max_l
    tearout_margin = F_tear/F_max-1

    return tearout_margin
