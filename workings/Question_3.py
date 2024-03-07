###### Question 3

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

# Initialize c and g1(c) and g2(c)

c = np.linspace(1,10,num=100)
print(c)
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
axes[0].set_xlim(0, 10)  # Set x-axis limits
axes[0].set_ylim(0, 10)  # Set y-axis limits
axes[0].legend()  # Display legend


# Plot data on the second subplot
axes[1].plot(c,c,label="y_1=c")
axes[1].plot(c,G2,label="y_2=g_2")
axes[1].set_title('y1=c and y1=g2')
axes[1].set_xlim(0, 10)  # Set x-axis limits
axes[1].set_ylim(0, 10)  # Set y-axis limits

# Add labels and title2
plt.xlabel('x')
plt.ylabel('y')

plt.legend()

#Show plot
#plt.show()

## Part (b)

#Fixed Point Iteration

c_r = np.zeros(1000, dtype=float)
c_r[0] = 4
epsilon_a = np.zeros(1000, dtype=float)
epsilon_a[0] = 1
epsilon_s = 0.5*10**-4


# # Perform iteration
i = 0  # Counter initialization
for i in range(20):
    c_r[i+1] = g1(c_r[i])
    epsilon_a[i] = abs(c_r[i+1] - c_r[i]) / abs(c_r[i+1])

    result_str = "{:.4g}".format(c_r[i])
    print("Root with FPI=",result_str,"iteration=",i)

    print("Root with FPI=",c_r[i],"iteration=",i)


#Netwon-Raphson

# Functions
def f(c):
    return (W - Q*c - k*V*math.sqrt(c))

def fprime(c):
    return (-Q - 0.5 * k * V * c**(-0.5))

# Iteration
for i in range(20):
    c_r[i+1] = c_r[i] - f(c_r[i])/fprime(c_r[i])
    epsilon_a[i+1] = abs(c_r[i+1] - c_r[i]) / abs(c_r[i+1])

    result_str = "{:.4g}".format(c_r[i])
    print("Root with NR=",result_str,"iteration=",i)

    print("Root with NR=",c_r[i],"iteration=",i)


def fdoubleprime(c):
    return (0.25 * k * V * c**(-1.5))

#Compute f''(x)'f'(x) near root

A = fdoubleprime(2.55) / fprime(2.55)
print("A=",A)











