
def Ixx_zz(w,t):
    I_xx = max((w**3*t)/12.0, 1e-15)
    I_zz = max((w*t**3)/12.0, 1e-15)
    return I_xx, I_zz
