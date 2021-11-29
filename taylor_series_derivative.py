from scipy.misc import derivative
import math

i = 0
def f(x):
    global i
    p = 1
    for k in range(i + 1):
        p *= (x + k) 
    return p/(math.factorial(i + 1))

arr = []
t = 0
for j in range(20):
    arr.append(derivative(f, 0, dx=1e-6, n = 2))
    i += 1

print(arr)