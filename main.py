from src.Sem3Py.Numerical import roots

func = lambda x : 3*x**3 + x*5 - 9
dx = lambda x: 9*x**2 + 5
k = roots(func, dx)
arr = k.newton_raphson(1, 10)
print(arr)