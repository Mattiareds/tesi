import numpy as np

def plane_from_points(P1, P2, P3):
    # Convert points to numpy arrays for easy vector calculations
    P1, P2, P3 = np.array(P1), np.array(P2), np.array(P3)
    
    # Calculate vectors
    v1 = P2 - P1
    v2 = P3 - P1
    
    # Cross product to get the normal vector
    normal_vector = np.cross(v1, v2)
    A, B, C = normal_vector
    
    # Calculate D using point P1
    D = -np.dot(normal_vector, P1)
    
    # Return the equation coefficients and the normal vector
    return A, B, C, D, normal_vector

# Read points from the terminal
def get_point_input(point_name):
    x = float(input(f"Enter x-coordinate for {point_name}: "))
    y = float(input(f"Enter y-coordinate for {point_name}: "))
    z = float(input(f"Enter z-coordinate for {point_name}: "))
    return (x, y, z)

# Get the points from user input
P1 = get_point_input("P1")
P2 = get_point_input("P2")
P3 = get_point_input("P3")

# Calculate the plane equation and the normal vector
A, B, C, D, normal_vector = plane_from_points(P1, P2, P3)

# Display the equation of the plane and the normal vector
print(f"The equation of the plane is: {A}x + {B}y + {C}z + {D} = 0")
print(f"The perpendicular (normal) vector to the plane is: {normal_vector}")
