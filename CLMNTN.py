'''
CLMNTN (Citrus Logistics Models for Nautical Trade Novelties AKA Clementine) is a suite of functions designed to allow 
the user to model various hull shapes at different dimensions and assess their aptitude based on weight, cargo packing, and drag. 

Glossary:
    - Hull: outer shell of the boat
    - Keel: lowest point along bottom of hull
    - Gunwale: upper edge of the hull
    - Draught: submerged part of the hull

Inputs:
    - dim = [x,y,z]: width, length, height 
    - hull_f = function describing the shape of the hull in the xy plane (must be even, differentiable, and have f(0)= 0)
'''

import math
import numpy as np
import scipy as sc
import sympy as sp


# This function will calculate the coefficient needed for the characterstic function (eg. x**2 or cosh(x)-1) to adhere to the 
# height constraints of the selected dimensions. This is a necessary step before proceeding as otherwise functions will output
# unreliable values.
def hull_f_adj_coeff(dim,ihull_f):

    # Find and return normalisation coefficient.
    coeff = (dim[2])/ihull_f(dim[0]/2)

    return coeff
# from this point on we will take hull_f to be a properly adjusted hull function (coeff * ihull_f).


# This function models the shape of the hull and calculates it's weight according to industry standard hull thickness
# and building materials.
def hull_weight(dim,hull_f):

    # Calculating area to be covered by steel, assuming no top cover.

    # Ship Face Area: by taking rectangle of dim[0] and dim[2], removing integral of hull_f.
    ship_face = (dim[0] * dim[2]) - sc.integrate.quad(hull_f, -0.5 * dim[0], 0.5 * dim[0])[0]

    # Length of the curve from keel to gunwale (along z axis) using standard formula.
    x = sp.Symbol('x')
    hull_f = sp.sympify(x,hull_f)

    deriv = sp.Derivative(hull_f).doit()
    curve_length = float(sp.integrate(sp.sqrt(deriv**2 + 1),(x, 0, 0.5 * dim[0])).doit())

    # Putting them together.
    total_surface = 2 * ship_face + 2 * dim[1] * curve_length

    # Finding volume of hull.
    hull_volume = total_surface * 0.014 # where 0.014 is average hull thickness in meters

    # finding weight of a steel hull.
    hull_weight = hull_volume * 7970 # where 7950 is average density of steel in kg/m**3 


    return hull_weight


# Using Citrus LTD's patented orange stacking method, this function efficiently packs a container of given 
# dimensions with Italy's third finest export, after Tommaso Tufarelli and Gerardo Adesso.
def orange_weight(cont_dim):

    # ORCs (Orange Related Constants)
    orange_volume_density = 988.07 # obtained through real world testing (yum)
    max_packing_efficiency = 0.74 # assuming max efficiency stacking

    # Volume and weight calculations
    cont_V = cont_dim[0] * cont_dim[1] * cont_dim[2]
    orange_weight = orange_volume_density * max_packing_efficiency * cont_V
    
    return orange_weight 


# Our packing solution calculates the maximum number of standard shipping containers that can fit
# in a given volume of hull. for ease of embark and disembark, the containers are all oriented in the
# direction of the ship. The most common sizes of shipping containers are 6.06 and 12.2 meters long.
# This function will fill the boat with long containers, fill the remaining space with short containers,
# and return the respective number of each. We favour long containers for weight reasons; they are lighter for
# the same volume of cargo. 
def n_containers(dim,hull_f):

    # These are the dimensions of industrial shipping containers.
    long_cont_dim = [2.4,12.2,2.59]
    short_cont_dim = [2.4,6.06,2.59]

    # Setting up our outputs.
    max_n_long_cont = 0
    max_n_short_cont = 0

    # Subdividing the hull into planes, iterating from top to bottom until we reach the maximum depth of the hull.
    max_floors = math.floor(dim[2] / long_cont_dim[2])
    for i in range(1, max_floors + 1):

        # Finding the x values for each of the "floors" by subtracting height and taking roots.
        x = sp.Symbol('x')
        temp_hull_f = hull_f(x)
        temp_hull_f = temp_hull_f - (dim[2] - i * (long_cont_dim[2]))
        temp_hull_f = sp.lambdify(x,temp_hull_f)
        width_vals = sc.optimize.fsolve(temp_hull_f, [dim[0],dim[0]])

        # if the plane is wide enough to accept cargo containers:
        if 2 * width_vals[-1] > long_cont_dim[0]:

            # determining maximum number of rows of containers a plane can house.
            temp_max_cont_rows = math.floor(2 * width_vals[-1] / long_cont_dim[0])
            # determining maximum number of long containers on the plane with respect to boat length.
            temp_max_n_long_cont = temp_max_cont_rows * math.floor(dim[1] / long_cont_dim[1])
            # appending to our long container output value.
            max_n_long_cont += temp_max_n_long_cont 


            # if plane is long enough to also include a short container:
            if dim[1] % long_cont_dim[1] >= 6.02:

                # add column of short containers to the plane.
                max_n_short_cont += temp_max_cont_rows

        # as soon as the ship hull is too narrow to include containers, we stop the iteration.
        else:
            break

    return (max_n_long_cont,max_n_short_cont)


# Now that we have the maximum number of containers that can be put placed into the hull, we fill them with ideally 
# stacked oranges (using previous 2 functions) and add the weight of the containers to find the weight of the cargo.
def cargo_weight(dim, hull_f):
    max_long_cont, max_short_cont = n_containers(dim, hull_f)
    
    # Calculate weight of cargo containers.
    long_cont_weight = max_long_cont * 3750 # where 3750 is the weight of a long shipping container in kg
    short_cont_weight = max_short_cont * 2300 # where 2300 is the weight of a short shipping container in kg
    total_cont_weight = long_cont_weight + short_cont_weight

    # Calculate weight of orange cargo with respect to size and number of containers.
    long_cont_orange_weight = orange_weight([2.4,12.2,2.59])
    short_cont_orange_weight = orange_weight([2.4,6.06,2.59])
    total_orange_weight = long_cont_orange_weight * max_long_cont + short_cont_orange_weight * max_short_cont

    # Add them together for total cargo and container weight.
    total_cargo_weight = total_cont_weight + total_orange_weight
    
    return total_cargo_weight


def boat_weight(dim, hull_f):

    boat_weight = hull_weight(dim, hull_f) + cargo_weight(dim, hull_f) + 2300000 
    # where 2300000 is the weight of an industry standard cargo ship engine.

    return boat_weight


# This function calculates the distance from the lowest y value to the waterline. to be used to determine cross section.
def draught(dim, hull_f, weight):

    # Establish rootfinding protocol (Newton-Raphson).
    
    def newton_stop(f,dfdx,p0,Nmax,TOL): 
    
        # Initialise the array with zeros and set first approximation.
        p_array = np.zeros(Nmax)
        p = p0 - (f(p0)/dfdx(p0))
        p_array[0] = p
    
        # Start loop
        for i in range(1,Nmax):
            # Store (to be compared to p for tolerance comparison)
            p1 = p
            # Begin iteration 
            p = p1 - (f(p1)/dfdx(p1))
            # Store 
            p_array[i] = p
            # Check against the tolerance
            if abs(p-p1)<=TOL and p1>p:
                # In the case where the tolerance is crossed, break the loop
                break
             
        # Trim the unecessary zeros from the array of approximations
        p_array = np.trim_zeros(p_array,'b')
    
        # Method finished 
        return p_array
    
    # Making Sympy conversions.
    x = sp.Symbol('x')
    fs = sp.sympify(hull_f(x))
    y0 = sp.Symbol('r')
    
    # Find quantity of water to be displaced.
    p = 1025 # where 1025 is the density of salt water in kg/m^3.
    Vs = weight/p

    # Check to see if the boat can physically float.
    max_boat_volume = np.abs(float(sp.integrate(fs-dim[2],( x, -0.5 * dim[0], 0.5 * dim[0])).doit()) * dim[1])
    if max_boat_volume < Vs:

        # If not, politely inform user:
        return ("vaffanculo e riprova")

    # Create function of x where x is the width at which the boat sits in the water.
    NRFs = 2 * x * dim[1] * fs - 2 * dim[1] * sp.Integral(fs,(x,0,x)).doit() - Vs
    dNRFs = sp.Derivative(NRFs,x,evaluate=True)
    
    # Further conversions.
    NRF = sp.lambdify(x,NRFs)
    dNRF = sp.lambdify(x,dNRFs)
    
    # Implement Newton-Raphson rootfinding solution.
    x0 = newton_stop(NRF,dNRF,20,100,0.0000000000000000001)[-1]
    y0 = hull_f(x0)
    submerged_volume =  np.abs(sp.integrate(fs-y0,(x, -x0, x0)).doit()) * dim[1]
    draft_surface = np.abs(float(sp.integrate(fs-y0,(x, -x0, x0)).doit()))

    return y0, draft_surface


# This function calculates the cross sectional drag with respect to the speed of the boat in knots. 
def cross_section_drag(dim, hull_f, weight, speed):

    adj_speed = 0.514444 * speed
    A = draught(dim, hull_f, weight)[1]
    p = 1025
    drag = 0.5 * (1.12) * p * A * adj_speed**2

    return drag

# This function normalises the cosh curve and then dampens it by a factor of i.
# (Due to time limitations, we would kindy ask you to only dampen to a fact of 6.)
def cosh_dampening(dim,i):

    norm_coeff = hull_f_adj_coeff(dim,lambda x:sp.cosh(x)-1)

    def f(x):
        return norm_coeff * (sp.cosh(x)-1)
    
    damp_coeff = lambda x:((2*x)/(dim[0]))**i
    def f2(x):
        return damp_coeff(x) * f(x)
    
    return f2

#--------------------------------
# Input your test dimensions and function here (replace dim and placeholder cosh curve)
dim = [38,230,58]

norm_coeff = hull_f_adj_coeff(dim,lambda x:sp.cosh(x)-1)
def f(x):
    return norm_coeff * (sp.cosh(x)-1)

# Example base testing suite
print("----------------")
print("The weight of the hull is " + str(hull_weight(dim,f)) + "kg.")
print("The number of long and short containers is " + str(n_containers(dim,f)) + " respectively.")
print("The weight of the cargo in their containers is " + str(cargo_weight(dim,f)) + "kg.")
print("The total weight of the ship is " + str(boat_weight(dim,f)) + "kg.")
print("The draught and drag surface are " + str(draught(dim, f, boat_weight(dim,f))))
print("----------------")
print("")

# The following is the cosh dampening function, deconstructed for ease of modification

i = 3 # dampening exponent
damp_coeff = lambda x:((2*x)/(dim[0]))**i

def f2(x):
    return damp_coeff(x) * f(x)

lc = 4200 # number of long containers filled with oranges, so as to easily vary capacity and check draught and drag values.
adjusted_weight = hull_weight(dim,f2) + lc * orange_weight([2.4,12.2,2.59]) + lc * 3750


#cross_section_drag(dim, hull_f, weight, speed)
# This is the basic cosh curve, simply fitted to the given dimensions.
print("Basic Cosh:")
d1 = draught(dim, f, boat_weight(dim, f))
dr1 = cross_section_drag(dim, f, boat_weight(dim, f),1)
cp1 = n_containers(dim,f)
print("The draft, drag, and capacity are "+str(d1[0])+"m, "+str(dr1)+"*speed^2 N, and "+str(cp1)+".")
print("")

# This demonstrates the same cosh curve, but with variable weight via lc
print("Adjusted Cosh:")
d2 = draught(dim, f, adjusted_weight)
dr2 = cross_section_drag(dim, f, adjusted_weight,1)
cp2 = n_containers(dim,f)
print("The draft, drag, and capacity are "+str(d2[0])+"m, "+str(dr2)+"*speed^2 N, and "+str(lc)+" long containers.")
print("")

# This demonstrates the last cosh curve, with a dampening exponent of i.
print("Dampened Cosh:")
d3 = draught(dim, f2, adjusted_weight)
dr3 = cross_section_drag(dim, f2, adjusted_weight,1)
cp3 = n_containers(dim,f2)
print("The draft, drag, and capacity are "+str(d3[0])+"m, "+str(dr3)+"*speed^2 N, and "+str(lc)+" long containers.")