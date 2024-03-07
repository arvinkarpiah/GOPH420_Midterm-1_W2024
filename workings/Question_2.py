###### Question 2a

import numpy as np
import math
from math import sqrt

#Define constants

g = 3.718
c = 47.0
m = 93
vf = g*m/c
vo = 100
b = 20
a = 0

#Generate t
dt = 2   # Change del_t accordingly
t = np.linspace(0, 20, num=int(((b-a)/dt)+1))
N = len(t)-1
print(t)

#Define function to integrate

def v(t):
    return vf + (vo - vf)*math.exp((-c/m)*t)

#Generate v for t
V = np.zeros(len(t), dtype=float)

for i in range(len(t)):
    V[i] = v(t[i])
#print(V)
#print(V[1:-1])

#Compute integration using trapezium method

total_sum = sum(V[1:-1])  # sum of all elements except the first and last

I_trap = ((b - a) / (2 * N)) * (V[0] + 2 * total_sum + V[-1])

#Compute integration using Simpsons 1/3

total_sum_1 = sum(V[1:-1:2])  # sum of odd elements except the first and last
total_sum_2 = sum(V[2:-1:2])  # sum of even elements except the first and last

I_simps = ((b - a) / (3 * N)) * (V[0] + 4 * total_sum_1 + 2 * total_sum_2 + V[-1])
print('for delta_t=',dt,'trapezoid integration=',I_trap)
print('for delta_t=',dt,'simpsons integration=',I_simps)

#Compute integration using 5-point Gauss Legendre Quadrature

def integrate_gausslegendre(f, lims, npts):
    if len(lims) != 2:
        raise ValueError("lims must contain exactly two values")
    
    try:
        a, b = map(float, lims)
    except (ValueError, TypeError):
        raise ValueError("lims must be convertible to float")

    valid_npts = [1, 2, 3, 4, 5]
    if npts not in valid_npts:
        raise ValueError(f"Invalid value for npts. Must be one of the following values: {valid_npts}")

    if callable(f):
        # Change limits
        A = (b + a) / 2
        B = (b - a) / 2

        if npts == 1:
            c_0 = 2.0
            x_0 = 0.0

            f_0 = f(A + B * x_0)
                
            I = B * c_0 * f_0
       

        elif npts == 2:
            c_0 = 1
            c_1 = 1
            x_0 = -1/sqrt(3)
            x_1 = 1/sqrt(3)
    

            f_0=f(A + B * x_0)
            f_1=f(A + B * x_1)

            I = B*(c_0 * f_0 + c_1 * f_1)
                  

        elif npts == 3:

            c_0 = 5/9
            c_1 = 8/9
            c_2 = 5/9
            x_0 = -sqrt(3/5)
            x_1 = 0
            x_2 = sqrt(3/5)

            f_0=f(A + B * x_0)
            f_1=f(A + B * x_1)
            f_2=f(A + B * x_2)
                
            I = B*(c_0 * f_0 + c_1 * f_1 + c_2 * f_2)

        if npts == 4:

            c_0 = (18-sqrt(30))/36
            c_1 = (18+sqrt(30))/36
            c_2 = (18+sqrt(30))/36
            c_3 = (18-sqrt(30))/36
            x_0 = -sqrt(3/7+(2/7)*sqrt(6/5))
            x_1 = -sqrt(3/7-(2/7)*sqrt(6/5))
            x_2 = sqrt(3/7-(2/7)*sqrt(6/5))
            x_3 = sqrt(3/7+(2/7)*sqrt(6/5))

            f_0=f(A + B * x_0)
            f_1=f(A + B * x_1)
            f_2=f(A + B * x_2)
            f_3=f(A + B * x_3)
    
            I = B*(c_0 * f_0 + c_1 * f_1 + c_2 * f_2 + c_3 * f_3)


        if npts == 5:


            c_0 = (322+13*sqrt(70))/900
            c_1 = (322-13*sqrt(70))/900
            c_2 = 128/225
            c_3 = (322-13*sqrt(70))/900
            c_4 = (322+13*sqrt(70))/900
            x_0 = -(1/3)*sqrt(5-2*sqrt(10/7))
            x_1 = -(1/3)*sqrt(5+2*sqrt(10/7))
            x_2 = 0.0
            x_3 = (1/3)*sqrt(5+2*sqrt(10/7))
            x_4 = (1/3)*sqrt(5-2*sqrt(10/7))

            f_0=f(A + B * x_0)
            f_1=f(A + B * x_1)
            f_2=f(A + B * x_2)
            f_3=f(A + B * x_3)
            f_4=f(A + B * x_4)
    
            I = B*(c_0 * f_0 + c_1 * f_1 + c_2 * f_2 + c_3 * f_3 + c_4 * f_4)

        return float(I)
    else:
        raise TypeError("f is not a callable")

I_gauss = integrate_gausslegendre(v,[a,b],5)
print('using 5 point GLQ',I_gauss)








