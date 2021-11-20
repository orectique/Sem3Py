import numpy as np

class Roots:

    """
    Instantiate a function.

    :param func: The lambda function to evaluate.
    :type func: Lambda function object
    
    :param dx: The derivative of the function passed to 'func'. Optional - Used only in Newton-Raphson method.
    :type dx: Lambda function object
    """
    def __init__(self, func, dx = None):
        self.func = func
        self.dx = dx

    def bisection(self, lbound: float, ubound: float, max_itn = 100, delta = 1e-6):
        """
        Finding the root of a function using bisection method.

        :param lbound: The lower bound of the range to evaluate.
        :type lbound: int

        :param ubound: The upper bound of the range to evaluate.
        :type ubound: int

        :param max_itn: The maximum number of iterations to evaluate.
        :type itn: int
        
        :param delta: The error threshhold.
        :type itn: float

        :return: Numpy array of iteration steps.
        :rtype: ndarray
        """
  
        a = lbound
        b = ubound
        itn = 0
        
        if self.func(a)*self.func(b) > 0:
            return 'Invalid upper and lower limits'
        
        else:
            while itn <= max_itn and abs(a-b) >= delta:
                x = (a + b)/2
                if itn == 0:
                    arr = [a, b, x, self.func(x)]
                else:
                    arr = np.vstack([arr, [a, b, x, self.func(x)]])
                if self.func(x) > 0:
                    b = x
                elif self.func(x) < 0:
                    a = x
                else:
                    return arr
                itn += 1

            return arr
        
    def incremental():
        return
    
    def false_position(self, x0: float, x1: float, max_itn = 100, delta = 1e-6):
        """
        Finding the root of a function using false position method.

        :param x0: The lower bound of the range to evaluate.
        :type lbound: int

        :param x1: The upper bound of the range to evaluate.
        :type ubound: int

        :param max_itn: The maximum number of iterations to evaluate.
        :type itn: int
        
        :param delta: The error threshhold.
        :type itn: float

        :return: Numpy array of iteration steps.
        :rtype: ndarray
        """

        a = x0
        b = x1
        lst = [a, b]
        i = 2
        itn = 0
        
        if self.func(a)*self.func(b) > 0:
            return 'Invalid upper and lower limits'
        
        else:
            while itn <= max_itn and abs(a - b) >= delta:
                lst.append(lst[0] - (lst[i - 1] - lst[0])*(self.func(lst[0]))/(self.func(lst[i-1]) - self.func(lst[0])))
                if itn == 0:
                    arr = [itn, i, lst[i-1], lst[i], self.func(lst[i])]
                else:
                    arr = np.vstack([arr, [itn, i, lst[i-1], lst[i], self.func(lst[i])]])
                a = lst[i - 2]
                b = lst[i - 1]
                itn += 1
                i += 1
            return arr
    
    def secant(self, x0: float, x1: float, max_itn = 100, delta = 1e-6):
        """
        Finding the root of a function using secant method.

        :param x0: The lower bound of the range to evaluate.
        :type lbound: int

        :param x1: The upper bound of the range to evaluate.
        :type ubound: int

        :param max_itn: The maximum number of iterations to evaluate.
        :type itn: int
        
        :param delta: The error threshhold.
        :type itn: float

        :return: Numpy array of iteration steps.
        :rtype: ndarray
        """

        a = x0
        b = x1
        lst = [a, b]
        i = 2
        itn = 0
        
        if self.func(a)*self.func(b) > 0:
            return 'Invalid upper and lower limits'
        
        else:
            while itn <= max_itn and abs(a - b) >= delta:
                lst.append(lst[i-1] - (self.func(lst[i - 1])*(lst[i-1] - lst[i-2]))/(self.func(lst[i-1]) - self.func(lst[i-2])))
                if itn == 0:
                    arr = [itn, i, lst[i-1], lst[i], self.func(lst[i])]
                else:
                    arr = np.vstack([arr, [itn, i, lst[i-1], lst[i], self.func(lst[i])]])
                a = lst[i - 2]
                b = lst[i - 1]
                itn += 1
                i += 1
            return arr
        return
    
    def newton_raphson(self, x0: float, max_itn = 100, delta = 0.0001):
        """
        Finding the root of a function using Newton-Raphson method.

        :param x0: The initial value to evaluate.
        :type x0: int

        :param max_itn: The maximum number of iterations to evaluate.
        :type itn: int
        
        :param delta: The error threshhold.
        :type itn: float

        :return: Numpy array of iteration steps.
        :rtype: ndarray
        """
        
        x = x0
        itn = 0
        
        while itn <= max_itn and abs(self.func(x)/self.dx(x)) >= delta:
            x = x - (self.func(x)/self.dx(x))
            
            if itn == 0:
                arr = [itn, x, self.func(x), self.dx(x), abs(self.func(x)/self.dx(x))]
            
            else: 
                arr = np.vstack([arr, [itn, x, self.func(x), self.dx(x), abs(self.func(x)/self.dx(x))]])
            itn += 1
            
        return arr
    
    
class Solutions:
    """
        Instantiate a function.
        
        :param arr_a: The coefficient matrix.
        :type arr_a: Numpy ndarray
           
        :param arr_b: The constant matrix.
        :type arr_b: Numpy ndarray
        
    """  
    
    def __init__(self, arr_a: np.ndarray, arr_b: np.ndarray):
        self.arr_a = np.array(arr_a).astype(float)
        self.arr_b = np.array(arr_b).astype(float)
                
    def gauss(self):
        """
            Finding the solution to a system of equations using Gauss Elimination.
            
            :return: The solution vector.
            :rtype: Numpy ndarray
            
        """
        a = self.arr_a
        b = self.arr_b
        n = len(b)
        x = [0]*n
        for k in range(0,n-1):
            for i in range(k+1,n):
                if a[i,k] != 0.0:
                    lam = a[i,k]/a[k,k]
                    a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                    b[i] = b[i] - lam*b[k]
                        
        for k in range(n-1,-1,-1):
            b[k] = (b[k] - np.dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
                        
        return b
    

    def gauss_seidel(self, funcs: list, guess: list, delta = 0.01):
        
        """
            Gauss Seidel approximation of the solution to a system of equations.
            
            :param funcs: List of lambda functions of lists.
            :type funcs: list
            
            :param guess: Vector of guessed values.
            :type guess: list
            
            :param delta: Tolerance value.
            :type delta: float
           
            :return: Iterated approximations of solution.
            :rtype: Numpy ndarray
            
        """
        
        funcs = funcs
        arr1 = guess
        arr2 = guess
        flag = False
        itn = 1
        
        while flag == False:
                
            temp_arr2 = []
            for i in range(len(funcs)):
                temp_arr = []
                for j in range(len(funcs)):
                    if i != j:
                        temp_arr.append(arr1[j])
                
                temp_arr2.append(funcs[i](temp_arr))
                               
            arr2 = np.vstack([arr2, temp_arr2])
            arr1 = temp_arr2
                        
            ls = [abs((arr2[itn, i] - arr2[itn-1, i])/arr2[itn, i]) for i in range(len(arr1))]

            for val in ls:
                if val < delta:
                    flag = True
                    
            itn += 1
            
        return arr2   
    
    def LUdecomp(self):
        """
            LU Decomposition of symmetric matrix.
            
            :param arr_a: The coefficient matrix.
            :type arr_a: Numpy ndarray
           
            :return: Decomposed L and U matrices.
            :rtype: Numpy ndarray, Numpy ndarray
            
        """
               
        a = self.arr_a
        n = len(a)
        for k in range(0,n-1):
            for i in range(k+1,n):
                if a[i,k] != 0.0:
                    lam = a [i,k]/a[k,k]
                    a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                    a[i,k] = lam
                    
        L = np.array([[0.0]*3]*3)
        U = np.array([[0.0]*3]*3)
        
        for i in range(n):
            for j in range(n):
                if i > j:
                    L[i, j] = a[i, j]
                elif i == j:
                    L[i, j] = 1
                    U[i, j] = a[i, j]
                else:
                    U[i, j] = a[i, j]
                
        return L, U 
                
    def LUsolve(self):
        """
            Finding the solution to a system of equations using LU Decomposition.
            
            :param arr_a: The coefficient matrix.
            :type arr_a: Numpy ndarray
           
            :param arr_b: The constant matrix.
            :type arr_b: Numpy ndarray
            
            :return: The solution vector.
            :rtype: Numpy ndarray
            
        """
        a = self.arr_a
        b = self.arr_b
        n = len(a)
        for k in range(1,n):
            b[k] = b[k] - np.dot(a[k,0:k],b[0:k])
        b[n-1] = b[n-1]/a[n-1,n-1]
        
        for k in range(n-2,-1,-1):
            b[k] = (b[k] - np.dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
        
        return b

    
class Eigenvalue:
    def __init__(self, matrix, init_vector):
        self.matrix = matrix
        self.init_vector = init_vector
        
        return
    
    def jacobi(self):
        return
    
    def power(self, max_itn = 100, delta = 0.0001):
        eig_val = 0.0000
        eig_prev = 1.0000
        eig_vect = self.init_vector
        matrix = self.matrix
        itn = 1
        
        table = np.array(['itn', 'eig_val', 'eig_vect'], dtype = object)
        
        while abs(eig_val - eig_prev) >= delta and itn <= max_itn:
            eig_vect = np.dot(matrix, eig_vect)
            
            eig_prev = eig_val
            
            eig_val = abs(eig_vect).max()
            eig_vect = eig_vect/eig_vect.max()
            
            table = np.vstack([table, [itn, eig_val, eig_vect]])
                
            itn += 1
                 
        return table
    
class Interpolation:
    
    """
    Instantiate a function.
    
    :param arr_a: Known 'x' values.
    :type arr_a: Numpy ndarray or list or series
           
    :param arr_b: Known 'y' values.
    :type arr_b: Numpy ndarray or list or series
    
    """
    def __init__(self, x_val, y_val):
        self.x_val = x_val
        self.y_val = y_val
            
    def lagrange(self, x: float) -> float:
        
        """
            Implementation of Lagrange's interpolation formula.
            
            :param x: The 'x' value in question.
            :type x: float
            
            :return: Interpolated 'y' value.
            :rtype: float
            
        """
        
        result = 0.0
        for i in range(len(self.x_val)):
            temp = self.y_val[i]
            
            for j in range(len(self.x_val)):
                if i != j:
                    temp = temp* (x - self.x_val[j])/(self.x_val[i] - self.x_val[j])
                    
            result += temp
        
        return result
    
    def newton_dd(self, x_new: float, forward : bool = True):
        """
            Implementation of Newton's divided difference algorithm.
            
            :param x_new: The 'x' value in question.
            :type x_new: float
            
            :param forward: Selection of forward or backward algorithm
            :type forward: bool
            
            :return: Interpolated 'y' value.
            :rtype: float

        """
    
        x = self.x_val
        y = self.y_val 
    
        n = len(y)
        coef = np.zeros([n, n])
        coef[:,0] = y
        for j in range(1,n):
            for i in range(n-j):
                coef[i][j] = (coef[i+1][j-1]-coef[i][j-1])/(x[i+j]-x[i])
            
            
        if forward:
            cf = coef[0, :]
        else:
            cf = np.flipud(coef).diagonal()
        
            
        k = 1
        p = 0
        for i in range(len(cf)):
            p += cf[i]*k
            k *= x_new - x[i]
                
    
        return coef, p

    def differentiate(self, x_new: float, h : float, second :bool = False, forward: bool = True):

        """
            Implementation of Newton's forward and backward difference differentiation algorithms.
            
            :param x_new: The 'x' value in question.
            :type x_new: float

            :param h: The h value - the step size
            :type h: float

            :param second: Selection of order of derivation
            :type forward: bool
            
            :param forward: Selection of forward or backward algorithm
            :type forward: bool
            
            :return: Interpolated 'y' value.
            :rtype: float

        """
        
        x = self.x_val
        y = self.y_val


        n = len(y)

        arr = [0, 0.0, 0.9999999999217246, 0.9999999999746642, 0.9166666666919701, 0.8333333333563364, 0.7611111111427086, 0.7000000000510863, 0.6482142857779337, 0.6039682540329495, 0.5657936508246969, 0.5325396825532482, 0.5033128908190265, 0.47741702740370173, 0.45430482214953094, 0.4335416435253855, 0.4147786241242193, 0.39773282270807536, 0.3821725024909123, 0.36790611347003926, 0.35477396569519837]
        
        if forward:

            x_index = x.index(x_new)

            coef = np.zeros([n, n])
            coef[:,0] = y
            for j in range(1,n):
                for i in range(n-j):
                    coef[i][j] = (coef[i+1][j-1]-coef[i][j-1])

            
            delt = coef[x_index,:]
            val = 0

            if second:
                for i in range(2, n):
                    if i%2 == 0:
                        val += arr[i]*delt[i]  
                    else:
                        val -= arr[i]*delt[i]
                dy = (1/h**2)*val
                
                
            else:

                i = 1
                while i<n:
                    if i%2 == 0:
                        val -= (1/i)*delt[i]  
                    else:
                        val += (1/i)*delt[i]
                    i += 1
                dy = (1/h)*val
                

        else:

            x.reverse()
            y.reverse()

            x_index = x.index(x_new)

            coef = np.zeros([n, n])
            coef[:,0] = y
            for j in range(1,n):
                for i in range(n-j):
                    coef[i][j] = -(coef[i+1][j-1]-coef[i][j-1])

            delt = coef[x_index,:]
            val = 0


            if second:
                for i in range(2, n):
                    val += arr[i]*delt[i]
                dy = (1/h**2)*val

            
            else:
                               
                i = 1
                while i<n:
                   val += (1/i)*delt[i]
                   i += 1
                dy = (1/h)*val
        
        return dy

    def trapezoidal(self, h : float):
        """
        Implementation of Trapezoidal method of Integration

        :param h: The 'h' value - step size.
        :type h: float
        
        """
        y = self.y_val

        return (h/2)*(2*sum(y[1: -1]) + y[0] + y[-1])

    def simpsons(self, h: float, three_eighths : bool= False ):
        """
        Implementation of Simpson's methods of Integration

        :param h: The 'h' value - step size.
        :type h: float

        :param three_eights: True - Simpson's 3/8 method, False(Default) - Simspson's 1/3 method
        :type three_eights: bool
        
        """
        y = self.y_val
        f = y.pop(0)
        l = y.pop(-1)

        if ~three_eighths:
            
            return (h/3)*(f + l + 4*sum([y[i] for i in range(len(y)) if i%2 == 0]) + 2*sum([y[i] for i in range(len(y)) if i%2 == 1]))
        
        else:
            return (3*h/8)*(f + l + 3*sum([y[i] for i in range(len(y)) if (i + 1)%3 != 0]) + 2*sum([y[i] for i in range(len(y)) if (i + 1)%3 == 0]))



        
    
class O_D_E:
    """
    Instantiate a function.

    :param func: Function of x and y
    :type func: Lambda function  
    
    :param x_init: Initial x value
    :type x_init: float
           
    :param y_init: Initial y value
    :type y_init: float
    
    """
    def __init__(self, func, x_init: float, y_init: float):
        self.func = func
        self.x_init = x_init
        self.y_init = y_init

    def euler(self, x_new : float, h: float):

        """
    Implementation of Euler's method for solving Initial Value Problem.
    
    :param x_new: The x value as input
    :type x_new: float
           
    :param h: The h value - step size.
    :type h: float
    
    """

        fun = self.func
        x = []
        y = []

        x_val = self.x_init
        y_val = self.y_init

        while x_val < x_new:
            y_val += h*fun(x_val, y_val)
            x_val += h
            x.append(x_val)
            y.append(y_val)
            
        
        return list(zip(x, y))

    def mod_euler(self, x_new : float, h: float):

        """
    Implementation of Modified Euler's method for solving Initial Value Problem.
    
    :param x_new: The x value as input
    :type x_new: float
           
    :param h: The h value - step size.
    :type h: float
    
    """

        fun = self.func
        x = []
        y = []

        x_val = self.x_init
        y_val = self.y_init

        while x_val < x_new:
            y_val += h*fun(x_val + h/2, y_val + (h/2)*fun(x_val, y_val))
            x_val += h
            x.append(x_val)
            y.append(y_val)
        
        return list(zip(x, y))

    def runge_kutta2(self, x_new : float, h: float):

        """
    Implementation of Second Order Runge-Kutta method for solving Initial Value Problem.
    
    :param x_new: The x value as input
    :type x_new: float
           
    :param h: The h value - step size.
    :type h: float
    
    """

        fun = self.func
        x = []
        y = []

        x_val = self.x_init
        y_val = self.y_init

        while x_val < x_new:

            k1 = h*fun(x_val, y_val)
            k2 = h*fun(x_val + h, y_val + k1)

            y_val += (k1 + k2)/2
            x_val += h
            x.append(x_val)
            y.append(y_val)
        
        return list(zip(x, y))

    def runge_kutta4(self, x_new : float, h: float):

        """
    Implementation of Fourth Order Runge-Kutta method for solving Initial Value Problem.
    
    :param x_new: The x value as input
    :type x_new: float
           
    :param h: The h value - step size.
    :type h: float
    
    """

        fun = self.func
        x = []
        y = []

        x_val = self.x_init
        y_val = self.y_init

        while x_val < x_new:

            k1 = h*fun(x_val, y_val)
            k2 = h*fun(x_val + h/2, y_val + k1/2)
            k3 = h*fun(x_val + h/2, y_val + k2/2)
            k4 = h*fun(x_val + h, y_val + k3)

            y_val += (1/6)*(k1 + 2*k2 + 2*k3 + k4)
            x_val += h
            x.append(x_val)
            y.append(y_val)
        
        return list(zip(x, y))