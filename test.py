import numpy as np
# Example parameters
n_f = 4         # number of fasteners
n_l = 2         # number of attachments
F_z = 40000     # [N]
F_y = 5000     # [N]
w_1 = 0.04      # horizontal plate width [m]
w_2 = 0.03      # vertical plate width [m]
l = 0.04         # vertical plate length [m]
t_1 = 0.01      # horizontal plate thickness [m]
t_2 = 0.015     # vertical plate thickness [m]
D_fo = 0.003    # fastener outer diameter [m]
yield_lug = 250e6   # [Pa]
yield_sc = 250e6    # [Pa]
tau_max_l = 480e6

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
        margins.append(yield_lug/np.abs(fasteners[i][3]))
        margins.append(yield_sc/np.abs(fasteners[i][4]))
        
    return min(margins)

def CreateFastenerList_v(n, l, w_1, D_fo):
    array = []
    n_side = n // 2

    if n_side == 1:
        y_positions = [0.0]
    else:
        d = (l - 2.5*D_fo) / (n_side - 1) 
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
        margins.append(yield_lug/np.abs(fasteners[i][3]))
        margins.append(yield_sc/np.abs(fasteners[i][4]))
        
    return min(margins)

def CreateFastenerList_h(n, w_2, w_1, D_fo):
    array = []
    n_side = n // 2

    if n_side == 1:
        y_positions = [0.0]
    else:
        d = (w_2 - 2*D_fo) / (n_side - 1) 
        y_positions = [-(2*D_fo + i*d) + w_2/2 for i in range(n_side)]

    for y in y_positions:
        array.append([+w_1/2 - 2*D_fo, y, 0, 0, 0])
        array.append([-w_1/2 + 2*D_fo, y, 0, 0, 0])

    return array


def FastenersRanked_h(fasteners):
    ranked1 = sorted(fasteners, key=lambda bolt: abs(bolt[2]), reverse=True)
    ranked2 = [bolt[:3] for bolt in ranked1]
    return ranked2


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

fasteners_v = CreateFastenerList_v(n_f, l, w_1, D_fo)
fasteners_v = ForceAtFastener_v(fasteners_v, w_2, n_f, n_l, F_z, F_y)
fasteners_v = bearing_stress_v(fasteners_v, t_2, t_1, D_fo)
vertical_margin = WorstShearMargin_v(fasteners_v, yield_lug, yield_sc)
tearout_margin_v = tearout_vertical(fasteners_v, l, w_1, t_2)

print("Vertical Pull-through:")
print("Fastener forces and stresses:", fasteners_v)
print("Worst shear margin:", vertical_margin)
print("Tear-out margin:", tearout_margin_v)


fasteners_h = CreateFastenerList_h(n_f, w_2, w_1, D_fo)
fasteners_h = ForceAtFastener_h(fasteners_h, n_f, n_l, F_z, F_y, l)
fasteners_h = bearing_stress_h(fasteners_h, t_2, t_1, D_fo)
horizontal_margin = WorstShearMargin_h(fasteners_h, yield_lug, yield_sc)
tearout_margin_h = tearout_horizontal(fasteners_h, w_1, w_2, t_2)

print("\nHorizontal Pull-through:")
print("Fastener forces and stresses:", fasteners_h)
print("Worst shear margin:", horizontal_margin)
print("Tear-out margin:", tearout_margin_h)