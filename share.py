import numpy as np
def gauss_seidel(funcs: list, guess: list, delta = 0.01):
        
        funcs = funcs
        arr1 = guess
        arr2 = guess
        flag = False
        itn = 1
        
        while flag == False:
                
            for i in range(len(funcs)):
                temp_arr = []
                for j in range(len(funcs)):
                    if i != j:
                        temp_arr.append(arr1[j])
                
                arr1[i] = funcs[i](temp_arr)
                               
            arr2 = np.vstack([arr2, arr1])
                        
            ls = [abs((arr2[itn, i] - arr2[itn-1, i])/arr2[itn, i]) for i in range(len(arr1))]
            flag = True
            for val in ls:
                if val > delta:
                    flag = False
                    
            itn += 1
        
        return arr2
    

a = lambda l: (1 - 3*l[0] + 5*l[1])/12
b = lambda l: (28 - l[0] - 3*l[1])/5
c = lambda l: (76 - 3*l[0] -7*l[1])/13

ls = [a, b, c]

print(gauss_seidel(ls, [ 0, 0, 0]))