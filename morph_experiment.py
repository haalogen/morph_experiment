#-*-coding: utf-8-*-
import sys
import numpy as np
import matplotlib.pyplot as plt
import func as fu


if len(sys.argv) == 1:
    n = 5 # default value
    print "Default: N = %r " % n
    
elif len(sys.argv) == 2:
    n = int(sys.argv[1])
else:
    print """
Usage:
    python morph_experiment.py [side_size]
    
    side_size == N -- length of image side (square image)
"""
    sys.exit(-1)
    



# create 2 matrices NxN: lImg, rImg
lImg = np.zeros( (n, n) )
rImg = np.zeros( (n, n) )


# define k, b, dx, dy
origK = 1.5
origB = 0
origDx = 0.0
origDy = 0.6


# choose intervals:[-1; 1] for example
absX = 1
absY = 1
lX = np.linspace(-absX, absX, n)
lY = np.linspace(-absY, absY, n)
print 'lX:\n', lX
print 'lY:\n', lY, '\n'


rX = np.linspace(-absX+origDx, absX+origDx, n)
rY = np.linspace(-absY+origDy, absY+origDy, n)
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



# lam ("lambda" letter) = 
#       sum (lImg[i,j] - rImg[i,j])**2 , 
# for all i,j = 0 .. N-1
#lam = 0
#for i in range(n):
#    for j in range(n):
#        lam += (lImg[i, j] - rImg[i, j])**2

#print 'lam: %r\n' % lam


A = lImg
B = rImg
# Let lImg == A(x, y), rImg == B(x, y)
# [!] Suppose: 
# B(x, y) = k*A(x+dx, y+dy) + b + dx*k*(dA/dx) + dy*k*(dA/dy)
# [2D Tailor's Row Decomposition till O(dx) + O(dy) including ]


Ax = np.zeros( (n, n) )
Ay = np.zeros( (n, n) )
#Calculating partial derivatives Ax = dA/dx & Ay = dA/dy

for i in range(n):
    for j in range(n):
        
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
        
        elif j == n-1:
#            only Ax_-
            Ax[i, j] = A[i, j] - A[i, j-1]
        
        else:
#            Ax = 0.5 * { (Ax_-) + (Ax_+) }
            Ax[i, j] = 0.5*( A[i, j+1] - A[i, j-1] )
        
        

#print "Ax: \n", Ax, "\n"
#print "Ay: \n", Ay, "\n"


# Calculating REVERSE (!!!) k, b 
# from 2x2 System of Linear Equations with Kramer's Rule

sA = sum( sum(A) )
sB = sum( sum(B) )
sAB = sum( sum( A*B ) )
sSqrB = sum( sum( B*B ) )

sAx = sum( sum(Ax) )
sAy = sum( sum(Ay) )
sSqrAx = sum( sum(Ax*Ax) )
sSqrAy = sum( sum(Ay*Ay) )

sAAx = sum( sum(A*Ax) )
sAAy = sum( sum(A*Ay) )
sBAx = sum( sum(B*Ax) )
sBAy = sum( sum(B*Ay) )
sAxAy = sum( sum(Ax*Ay) )

#print A*A, '\n'
#print 'sA: ', sA, '\n'
#print 'sB: ', sB, '\n'
#print 'sSqrB: ', sSqrB, '\n'
#print 'sAB: ', sAB, '\n'

#delta = sSqrB * (n*n) - sB * sB


#delta1 = sAB * (n*n) - sA * sB

#delta2 = sSqrB * sA - sB * sAB


#print 'delta: ', delta
#print 'delta1: ', delta1
#print 'delta2: ', delta2

#Actually, we calculated k^-1 and -b/k
#My fault, now too lazy to rewrite
#So
#k = delta / delta1
#b = -k * delta2 / delta

a11 = sSqrB
a12 = sB
a21 = sB
a22 = n*n
b1 = sAB
b2 = sA


k, b = fu.kramer_solve_SLE_2(a11, a12, a21, a22, b1, b2)
k = 1/k
b = -k*b

# Calculating dx, dy

#delta = sSqrAx * sSqrAy - sAxAy * sAxAy

#delta1 = ( (sBAx - k*sAAx - b*sAx) * sSqrAy - 
#            (sBAy - k*sAAy - b*sAy) * sAxAy      ) / k

#delta2 = ( (sBAy - k*sAAy - b*sAy) * sSqrAx - 
#            (sBAx - k*sAAx - b*sAx) * sAxAy     ) / k


#dx = delta1 / delta
#dy = delta2 / delta

a11 = sSqrAx
a12 = sAxAy
a21 = sAxAy
a22 = sSqrAy
b1 = (sBAx - k * sAAx - b * sAx) / k
b2 = (sBAy - k * sAAy - b * sAy) / k


dx, dy = fu.kramer_solve_SLE_2(a11, a12, a21, a22, b1, b2)

print 'Original k: ', origK
print 'Original b: ', origB
print 'Original dx: ', origDx
print 'Original dy: ', origDy


print 'Calculated k: ', k, '\n'
print 'Calculated b: ', b, '\n'
print 'Calculated dx: ', dx, '\n'
print 'Calculated dy: ', dy, '\n'






# --- PLOTTING ---
#pcolormesh -- colormap
#cmap = plt.cm.Greys -- for grayscale
colorMap = plt.cm.gray
plt.subplot(1, 2, 1)
plt.imshow(lImg, cmap=colorMap, interpolation='none')

plt.subplot(1, 2, 2)
plt.imshow(rImg, cmap=colorMap, interpolation='none')

plt.show()
