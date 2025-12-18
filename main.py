param_bounds = {
    'w':  (0.003, 0.080),
    't':  (0.001, 0.010),
    'D':  (0.0001, 0.020),
    't2': (0.0001, 0.020),
    'D2': (0.005, 0.012),
    'n_f': (4.0, 8.0),       # integer valued; will round
    'x':  (0.004, 0.012*2),
    'mat_lug': (0.0, 1.0),
    'mat_pin': (0.0, 1.0),
    'phi': (0.1,0.45)
}

keys = ['w','t','D','t2','D2','n_f','x','mat_lug','mat_pin', 'phi']
scales = [(param_bounds[k][0], param_bounds[k][1] - param_bounds[k][0]) for k in keys]

def unscale(x_scaled):
    vals = [low + x*rng for x,(low,rng) in zip(x_scaled, scales)]
    # round n_f to integer 1..8
    vals[5] = int(round(vals[5]))  # n_f index 5 in keys
    if vals[5] < int(param_bounds['n_f'][0]): vals[5] = int(param_bounds['n_f'][0])
    if vals[5] > int(param_bounds['n_f'][1]): vals[5] = int(param_bounds['n_f'][1])
    
    # Constraint the number of fasteners to even number -> creates some converegnce problems
    if int(vals[5]) % 2 == 1:
        if int(vals[5]) > 4:
            vals[5] -= 1
        else:
            vals[5] += 1

    return vals

# ---------------------------
# Objective function (mass + pin mass) with constraints via penalty
# ---------------------------

def objective(x_scaled, penalize_factor=2e4, debug=False):
    w,t,D,t2,D2,n_f,x,mat_lug_cont,mat_pin_cont,phi = unscale(x_scaled)
    mat_lug = "8630" if mat_lug_cont >= 0.5 else "7075-T6"
    mat_pin = "8630" if mat_pin_cont >= 0.5 else "7075-T6"

    # early sanity
    if any(val <= 0 for val in (w,t,D,t2,D2)) or n_f < 1:
        return 1e12
    
    Ftu,Fty,rho,E_lug = material_props(mat_lug)
    tau_y = Fty * 1/np.sqrt(3) # van Mises
    h = clearance_spacing(w, tau_y)
    
    # backplate
    marginbackplate, m_backplate = backplate(D2,x,h,t,t2, rho, Fty, n_f) # n_f from 2 to eight 
    
    MS_out, f_val, rho_lug = lug_limits(Fx/2,Fy/2,Fz/2,D,w,t,mat_lug)
    invalid_lug = (MS_out < -1e5) or (f_val < -1e5)

    m_lug, y_cg = mass_lug(D,w,t,rho_lug if rho_lug is not None else 3000.0)

    margin_bending = WorstBackplateMargin(D,t2,D2,Fty,w) # bending and torsion load on backplate
    
    # Calculate total y-coordinate of cg of the whole part
    part_cg_y = (2*m_lug*(y_cg + t2) + m_backplate*0.5*t2)/(m_backplate + 2*m_lug)
    
    # Fasteners
    # Position of each fastener in x,z and force acting on the fastener
    # [x-coord, z-coord, Fy, shear in flange, shear in S/C]
    plate_width = 2*(h/2 + t + x + 2*D2)
    
    fasteners = CreateFastenerList(n_f, w, plate_width, D2)
    fasteners = ForceAtFastener(fasteners,part_cg_y,n_f) # Calculate force at each fastener ------ took part_cg_y here do not know wether that is correct
    fasteners = ShearInPlates(fasteners, t2, t3, D2) # Calculate shear stresss at s/c wall and lug backplate
    
    margin_shear = WorstShearMargin(fasteners, Fty) # Find the worst margin for shear stress
    m_fastener = mass_fastner(t2,E_lug,D2,n_f,phi)
    margin_fastener_shear = fastener_shear(D2,n_f,x,t,w,h)
    
    for i in range(len(fasteners)):
        """temp = (plate_width/2 - np.abs(fasteners[i][0]))
        if temp - 2*D2 <= 1e-5 and temp - 2*D2 > 0:
            print("True: " + str(plate_width/2 - np.abs(fasteners[i][0])))
        else:
            print("False: " + str(plate_width/2 - np.abs(fasteners[i][0])))"""
        
        """temp = (w/2 - np.abs(fasteners[i][1]))
        if temp - 2*D2 <= 1e-5 and temp - 2*D2 > 0:
            print("True: " + str(w/2 - np.abs(fasteners[i][1])))
        else:
            print("False: " + str(w/2 - np.abs(fasteners[i][1])))"""
        
    """for i in range(len(fasteners)):
        for j in range(len(fasteners)):
            if i != j:
                dx = np.abs(fasteners[i][0] - fasteners[j][0])
                dy = np.abs(fasteners[i][1] - fasteners[j][1])

                print(i,j)
                if dx - 2*D2 < 0:
                    print("Horizontal violation: " + str(dx - 2*D2))
                elif dy - 2*D2 < 0:
                    print("Vertical violation: " + str(dy - 2*D2))"""

    # Check vertical distance to top and bottom edges
    penalty_w = 0.0
    for i in range(len(fasteners)):
        temp = (w/2 - np.abs(fasteners[i][1]))
        violation = max(0.0, 2*D2 - temp)  # >0 if too close to the edge
        penalty_w += violation**2  # quadratic penalty

    penalize_factor_w = 2e5
    
    # Check horizontal distance to left and right edges
    penalty_x = 0.0
    for i in range(len(fasteners)):
        temp = (plate_width/2 - np.abs(fasteners[i][0]))
        violation = max(0.0, 2*D2 - temp)  # >0 if too close to the edge
        penalty_x += violation**2  # quadratic penalty

    penalize_factor_x = 2e5

    # Check vertical distance between each fastener
    for i in range(len(fasteners)):
        for j in range(len(fasteners)):
            temp = np.abs(fasteners[i][1] - fasteners[j][1])
            violation = max(0.0, 2*D2 - temp)  # >0 if too close to the edge
            penalty_w += violation**2  # quadratic penalty

    # Check horizontal distance between each fastener
    for i in range(len(fasteners)):
            for j in range(len(fasteners)):
                temp = np.abs(fasteners[i][0] - fasteners[j][0])
                violation = max(0.0, 2*D2 - temp)  # >0 if too close to the edge
                penalty_x += violation**2  # quadratic penalty

    P = math.hypot(Fy, Fz)
    b = 0.5 * t + 0.25*h + PIN_GAP

    d_req, rho_pin = required_pin_diameter(P, t, b, mat_pin, sf, d_min=0.002, d_max=0.050)
    min_D_allowed = d_req + PIN_CLEARANCE

    m_pin = mass_pin(d_req, h + 2*PIN_GAP, rho_pin)

    if invalid_lug:
        pseudo = 2 * m_lug + m_pin + m_backplate + n_f * m_fastener
        return pseudo * (1.0 + penalize_factor * 10.0)

    deficit = max(0.0, (min_D_allowed - D) / max(min_D_allowed, 1e-12))
    pin_penalty = 1.0 + 1e5 * deficit**2

    # new constraints: marginbackplate and margin_shear must be > 0
    v_back = max(0.0, -marginbackplate)  # if margin < 0, penalize
    v_shear = max(0.0, -(margin_shear - 1.0))  # margin_shear is ratio yield/shear; require >1 -> margin>1 -> convert so >0 means violation
    # (we convert so that when margin_shear < 1 it becomes positive violation)

    v_bending = max(0.0, -(margin_bending - 1.0))
    v_shear_fastener = max(0.0, -(margin_fastener_shear - 1.0))

    v_MS = max(0.0, (MS_g - MS_out) / max(MS_g, 1e-12))
    v_f = max(0.0, (sf - f_val) / max(sf, 1e-12))

    design_pen = 1.0 + penalize_factor * (v_MS**2 + v_f**2 + 1e4 * v_back**2 + 1e4 * v_shear**2 + 1e4 * v_bending**2 + 1e4 * v_shear_fastener)

    r = w/D
    geo = 1.0 + 1e-3 * ((max(0.0, 0.8 - r))**2 + (max(0.0, r - 6.0))**2)

    # total mass: m_pin + 2*m_lug + m_backplate + n_f * m_fastener
    total_mass = (m_pin + 2.0 * m_lug + m_backplate + float(n_f) * m_fastener) * pin_penalty * design_pen * geo * (1.0 + penalize_factor_w * penalty_w) * (1.0 + penalize_factor_x * penalty_x)

    if debug:
        print("DBG:", {"D":D,"w":w,"t":t,"t2":t2,"D2":D2,"n_f":n_f,"x":x,"mat_lug":mat_lug,"mat_pin":mat_pin})
        print(" masses: lug=",m_lug,"pin=",m_pin,"back=",m_backplate,"fast=",m_fastener,"n_f=",n_f)
        print(" margins: MS_out=",MS_out,"f=",f_val,"m_back=",marginbackplate,"m_shear=",margin_shear,"d_req=",d_req)
        print(" objective mass", total_mass)

    return total_mass

# ---------------------------
# CMA-ES wrapper
# ---------------------------

def optimize_cma(popsize=12, seed=0, maxiter=1000):
    rng = np.random.RandomState(seed)
    dim = len(keys)
    x0 = rng.rand(dim).tolist()
    sigma0 = 0.2
    opts = {'popsize':popsize,'bounds':[[0.0]*dim,[1.0]*dim],'CMA_active':True,'verb_disp':1,'maxiter':maxiter,'maxfevals':200000}
    es = cma.CMAEvolutionStrategy(x0, sigma0, opts)

    best = {'obj':float('inf'),'x':None}
    while True:
        sols = es.ask()
        fitnesses=[]
        for s in sols:
            val = objective(s)
            fitnesses.append(val)
            if val < best['obj']:
                best = {'obj':val,'x':np.array(s).copy()}
        es.tell(sols, fitnesses)
        es.disp()
        if es.stop(): break

    best_scaled = best['x']
    w,t,D,t2,D2,n_f,x,mat_lug_cont,mat_pin_cont,phi = unscale(best_scaled)
    mat_lug = "8630" if mat_lug_cont>=0.5 else "7075-T6"
    mat_pin = "8630" if mat_pin_cont>=0.5 else "7075-T6"
    e = w/2.0
    MS,f_val,rho = lug_limits(Fx/2,Fy/2,Fz/2,D,w,t,mat_lug)
    Ftu, Fty, rho, E = material_props(mat_lug)
    P = math.hypot(Fy,Fz)
    b = 0.5 * t + PIN_GAP
    d_req, rho_pin = required_pin_diameter(P,t,b,mat_pin,sf)
    m_lug,y_cg = mass_lug(D,w,t,rho)
    h = clearance_spacing(w,Fty/np.sqrt(3))
    m_pin = mass_pin(d_req,2*h+2*PIN_GAP,rho_pin)
    marginbackplate, m_backplate = backplate(D2,x,h,t,t2,rho,Fty,n_f)
    m_fastener = mass_fastner(D2, t2, t3,n_f, phi)
    total_mass = m_pin + 2.0*m_lug + m_backplate + n_f * m_fastener

    plate_width = 2*(h/2 + t + x + 2*D2)
    fasteners = CreateFastenerList(n_f, w, plate_width, D2)

    return {
        'scaled':best_scaled,
        'w':w,'t':t,'D':D,'d_req':d_req,'t2':t2,'D2':D2,'n_f':n_f,'x':x,'h':h,
        'mat_lug':mat_lug,'mat_pin':mat_pin,
        'mass_lug':m_lug,'mass_pin':m_pin,'mass_backplate':m_backplate,'margin_backplate':marginbackplate,'mass_fastener':m_fastener,
        'mass_total':total_mass,'MS_out':MS,'f':f_val,'obj':best['obj'], 'fasteners':fasteners
    }

# Main function
if __name__ == '__main__':
    
    # Create excel sheet
    wb = Workbook()
    ws = wb.active

    print('Starting optimization (independent lug+pin materials)...')
    res = optimize_cma(popsize=14, seed=42, maxiter=200)
    print('\nBest result:')
    print(f" w = {res['w']:.6f} m")
    print(f" h = {res['h']:.6f} m")
    print(f" t = {res['t']:.6f} m")
    print(f" D = {res['D']:.6f} m")
    print(f" t2 = {res['t2']:.6f} m")
    print(f" D2 = {res['D2']:.6f} m")
    print(f" n_f = {res['n_f']}")
    print(f" x = {res['x']:.6f} m")
    print(f" mat_lug = {res['mat_lug']}")
    print(f" mat_pin = {res['mat_pin']}")
    print(f" mass_flange (each) = {res['mass_lug']:.6f} kg")
    print(f" mass_pin = {res['mass_pin']:.6f} kg")
    print(f" mass_backplate = {res['mass_backplate']:.6f} kg")
    print(f" mass_fastener (each) = {res['mass_fastener']:.6f} kg")
    print(f" total mass (approx) = {res['mass_total']:.6f} kg")
    print(f" MS_out (lug) = {res['MS_out']}")
    print(f" Margin backplate = {res['margin_backplate']:.6f}")
    print(f" f (lug) = {res['f']:.6f}")
    print(f" objective = {res['obj']:.6f}")
    print(f" fasteners = {res['fasteners']}")

    """fasteners = CreateFastenerList(4, 0.5, 1, 0.003)
    print(fasteners)"""
    #plot_lug_with_pin(res['D'], res['w'], res['t'], res['d_req'])

    file_path = r"C:\Users\Kristian Nemcek\Desktop\Optimization.xlsx"
    
    # Load or create workbook
    if os.path.exists(file_path):
        wb = load_workbook(file_path)
        ws = wb.active
    else:
        wb = Workbook()
        ws = wb.active

    # Append headers only if the sheet is empty
    if ws.max_row == 1 and ws.cell(row=1, column=1).value is None:
        headers = [
            'w', 't', 'D', 't2', 'D2', 'n_f', 'x', 'mat_lug', 'mat_pin',
            'mass_lug', 'mass_pin', 'mass_backplate', 'mass_fastener',
            'mass_total', 'MS_out', 'margin_backplate', 'f', 'objective'
        ]
        ws.append(headers)

    # Prepare data row
    row = [
        res['w'], res['t'], res['D'], res['t2'], res['D2'], res['n_f'],
        res['x'], res['mat_lug'], res['mat_pin'], res['mass_lug'],
        res['mass_pin'], res['mass_backplate'], res['mass_fastener'],
        res['mass_total'], res['MS_out'], res['margin_backplate'],
        res['f'], res['obj']
    ]

    # Append the row
    ws.append(row)
    wb.save(file_path)