def Ixx_zz(w_1,t):
    I_xx = max((t**3*w_1)/12.0, 1e-15) #Ixx of the horizontal polate, same for Iyy
    I_zz = max((w_1**3*t)/12.0, 1e-15)
    return I_xx, I_zz
def bending(Ixx, M, w_2, Fz, t_2):
    M = Fz * (w_2/2)
    Bending = (M * t_2)/Ixx
    return Bending
def length(Fz,t_2, tau):
    L = (3*Fz)/(8*t_2*tau)
    