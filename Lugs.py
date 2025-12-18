import math
import numpy as np

n_l = 1 #number of attatchemnts
n_f = 1 #numebr of lugs
tau_max = 1 #ultimate shear stress of fastener

metric_bolt_d_i = [0.729, 0.829, 0.929, 1.075, 1.221, 1.421, 1.567, 1.713, 2.013, 2.459, 2.850, 3.242, 3.688, 4.134, 4.917, 5.917, 6.647]
metric_bolt_name = ['M1', 'M1.1', 'M1.2', 'M1.4', 'M1.6', 'M1.8', 'M2', 'M2.2', 'M2.5', 'M3', 'M3.5', 'M4', 'M4.5', 'M5', 'M6', 'M7', 'M8']
metric_bolt_d_o = [1, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8]

F_x = 1
F_y = 1
F_z = 1

shear_yield_lug = 1
shear_yield_sc = 1 


def fastener_diameter(F_z, n_l, n_f, shear_yield_lug):
    D = math.sqrt(32*F_z / n_l / n_f / (3*math.pi*shear_yield_lug))
    
    # Pick the smallest bolt >= D
    for i, d in enumerate(metric_bolt_d_i):
        if d >= D*1000:
            D_metric_i = metric_bolt_d_i[i]*0.001       # inner diameter in meters
            D_metric_o = metric_bolt_d_o[i]*0.001
            M_size = metric_bolt_name[i]
            break
    else:
        raise ValueError("No metric bolt large enough")
    
    F_applied = F_z / n_l / n_f
    F_allowable = 3 * shear_yield_lug * math.pi * D_metric_i**2 / 32
    shear_margin = F_allowable / F_applied - 1
    
    return shear_margin, D_metric_o, M_size





### pull through from wp4, moment and force modified

def ForceAtFastener(fasteners,w_2,n_f):
    
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
def ShearInPlates(fasteners, t_lug, t_sc, D_fi):
    
    # Contact area with fastener per surface
    area_lug = np.pi*D_fi*t_lug
    area_sc = np.pi*D_fi*t_sc

    # Calculate stress per fastener per surface
    for i in range(len(fasteners)):
        Fy = fasteners[i][2]
        fasteners[i][3] = Fy / area_lug
        fasteners[i][4] = Fy / area_sc
    
    return fasteners

def WorstShearMargin(fasteners, yield_lug, yield_sc):
    margins = []
    for i in range(len(fasteners)):
        margins.append(yield_lug/np.abs(fasteners[i][3]))
        margins.append(yield_sc/np.abs(fasteners[i][4]))
        
    return min(margins)

def CreateFastenerList(n: int, height, width, D2):
    array = []
    n_side = n/2
    d = (height - 4*D2)/(n_side - 1)

    for i in range(int(n_side)):
        if i != (n_side - 1):
            array.append([-2*D2 + width/2, -(2*D2 + i*d) + height/2, 0,0,0])
            array.append([-(width - 2*D2) + width/2, -(2*D2 + i*d) + height/2, 0,0,0])
        else:
            array.append([-2*D2 + width/2, -(height - 2*D2) + height/2, 0,0,0])
            array.append([-(width - 2*D2) + width/2, -(height - 2*D2) + height/2, 0,0,0])
    
    return array

def FastenersRanked(fasteners):
    ranked1 = sorted(fasteners, key=lambda bolt: abs(bolt[2]), reverse=True)
    ranked2 = [bolt[:3] for bolt in ranked1]
    return ranked2


shear_margin, d_metric_o, M_size = fastener_diameter(30000,4,2,480*10**6)

print(shear_margin, d_metric_o, M_size)