from tank_redesign import tank_radius
from dimension import compare_geom


# Define overall height of the satellite (and closing panels)
h = 1.35 # m

# Main function
if __name__ == '__main__':
    
    # Get tank properties
    R_tank_prop, t_Ti_prop, M_Ti_prop = tank_radius(h)
    
    # Add margin to structure radius so the propelant tank can fit
    R_struct = R_tank_prop + 0.01   # Add 1 cm

    # Calculate properties of the sandwich panels
    # geometry type, mass, parameter length (hexagon side, circle radius or rectangle side) and weight saved
    geom, min_mass, param_length, weight_saved = compare_geom(h)

    

    

    

