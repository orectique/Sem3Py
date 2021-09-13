from src.Sem3Py.Numerical import roots

k = roots('x**2')
arr = k.bisection(-1.0, 2.0, 10)
print(arr)

