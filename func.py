import numpy as np



def init_images(n, m, origK, origB, origDx, origDy):
    
    lImg = np.zeros( (n, m) )
    rImg = np.zeros( (n, m) )
    
    
# choose intervals:[-1; 1] for example
    absX = 1
    absY = 1
    lX = np.linspace(-absX, absX, n)
    lY = np.linspace(-absY, absY, m)
    print 'lX:\n', lX
    print 'lY:\n', lY, '\n'
    
    
    rX = np.linspace(-absX+origDx, absX+origDx, n)
    rY = np.linspace(-absY+origDy, absY+origDy, m)
    print 'rX:\n', rX
    print 'rY:\n', rY, '\n'
        
    
    # creating meshgrid
    lXX, lYY = np.meshgrid(lX, lY)
    rXX, rYY = np.meshgrid(rX, rY)
    
    
    lImg = np.e ** (- lXX**2 - lYY**2)
    # rImg = k * lImg(x+dx, y+dy) + b
    rImg = origK * np.e ** (- rXX**2 - rYY**2) + origB
    
    
    print 'lImg:\n', lImg, '\n'
    print 'rImg:\n', rImg, '\n'
    
    
    return (lImg, rImg)




def kramer_solve_SLE_2(a11, a12, a21, a22, b1, b2):
    """Solves System of Linear Equations rank==2 using Kramer's rule"""
    _delta = a11 * a22 - a12 * a21
    _delta1 = b1 * a22 - b2 * a12
    _delta2 = b2 * a11 - b1 * a21
    
    x = _delta1 / _delta
    y = _delta2 / _delta
    return (x, y)




def calc_derivatives_Ax_Ay(A):
    n = A.shape[0]
    m = A.shape[1]
    
    Ax = np.zeros( (n, m) )
    Ay = np.zeros( (n, m) )
    
    for i in range(n):
        for j in range(m):
            
            if i == 0:
    #            only Ay_+
                Ay[i, j] = A[i+1, j] - A[i, j]
            
            elif i == n-1:
    #            only Ay_-
                Ay[i, j] = A[i, j] - A[i-1, j]
            
            else:
    #            Ay = 0.5 * { (Ay_-) + (Ay_+) }
                Ay[i, j] = 0.5*( A[i+1, j] - A[i-1, j] )
            
            
            if j == 0:
    #            only Ax_+
                Ax[i, j] = A[i, j+1] - A[i, j]
            
            elif j == m-1:
    #            only Ax_-
                Ax[i, j] = A[i, j] - A[i, j-1]
            
            else:
    #            Ax = 0.5 * { (Ax_-) + (Ax_+) }
                Ax[i, j] = 0.5*( A[i, j+1] - A[i, j-1] )
            
    
    return Ax, Ay

