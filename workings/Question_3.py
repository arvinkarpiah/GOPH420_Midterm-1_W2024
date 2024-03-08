###### Solutions for Question 3.

###### By default the code will print roots (in exact form and to 4sf) for fixed point iteration and Newton Raphson 
###### User will be able to choose other variables/plots to display at the bottom of this script 

import numpy as np
import math
from math import sqrt
import matplotlib.pyplot as plt

#Define constants

V = 1.38 * 10**6
W = 2.15 * 10**6
k = 0.825
Q = 1.29 * 10**5

## Part (a)

# Define g1 and g2

def g1(c): 
    return ((W - Q * c)/(k * V)) ** 2

def g2(c): 
    return float(((W - (k * V) * math.sqrt(c)) / Q))

#Checking for convergence

c_true = 2.56
c_ini = 4

g_prime_ksi = (g2(c_true) - g2(c_ini)) / (c_true - c_ini)


# Initialize c and g1(c) and g2(c)

c = np.linspace(1,10,num=100)
#print("c=",c)
G1 = np.zeros(len(c), dtype=float)
G2 = np.zeros(len(c), dtype=float)

for j in range(len(c)):
    G1[j] = g1(c[j])
    G2[j] = g2(c[j])

# #Plot
plt.figure()
fig, axes = plt.subplots(2)  # 2 rows, 1 column

# Plot data on the first subplot
axes[0].plot(c,c,label="y_1=c")
axes[0].plot(c,G1,label="y_2=g_1")
axes[0].set_title('y1=c and y1=g1')
axes[0].set_xlim(1, 10)  # Set x-axis limits
axes[0].set_ylim(1, 10)  # Set y-axis limits
axes[0].legend()  # Display legend

# Plot data on the second subplot
axes[1].plot(c,c,label="y_1=c")
axes[1].plot(c,G2,label="y_2=g_2")
axes[1].set_title('y1=c and y1=g2')
axes[1].set_xlim(1, 10)  # Set x-axis limits
axes[1].set_ylim(1, 10)  # Set y-axis limits

plt.xlabel('x')
plt.ylabel('y')
plt.legend()

## Part (b)

#Fixed Point Iteration

c_r = np.zeros(1000, dtype=float)
c_r[0] = 4
epsilon_a = np.zeros(1000, dtype=float)
epsilon_a[0] = 1
epsilon_s = 0.5*10**-4

i = 0  
for i in range(15):
    c_r[i+1] = g1(c_r[i])
    epsilon_a[i] = abs(c_r[i+1] - c_r[i]) / abs(c_r[i+1])

    result_str = "{:.4g}".format(c_r[i])

    print("Root with FPI=",result_str,"iteration=",i); # View exact root and rounded to 4sf for fixed point iteration method
    print("Root with FPI=",c_r[i],"iteration=",i)

#Netwon-Raphson
    
# Functions
def f(c):
    return (W - Q*c - k*V*math.sqrt(c))

def fprime(c):
    return (-Q - 0.5 * k * V * c**(-0.5))

for i in range(15):
    c_r[i+1] = c_r[i] - f(c_r[i])/fprime(c_r[i])
    epsilon_a[i+1] = abs(c_r[i+1] - c_r[i]) / abs(c_r[i+1])

    result_str = "{:.4g}".format(c_r[i])

    print("Root with NR=",result_str,"iteration=",i) # View exact root and rounded to 4sf for Newton-Raphson method
    print("Root with NR=",c_r[i],"iteration=",i)
    

#Compute f''(x)'f'(x) near root

def fdoubleprime(c):
    return (0.25 * k * V * c**(-1.5))

A = fdoubleprime(2.55) / fprime(2.55)


#print("g_prime_ksi=",g_prime_ksi)                     # View convergence criterion

#plt.show() # View plot used for discussion in Q3a

#print("A=",A) #View additional criteria to check for Newton-Raphson convergence speed







