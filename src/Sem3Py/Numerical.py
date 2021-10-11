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
    
    def newton_dd(self):
        
        
        return