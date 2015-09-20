#-*-coding: utf-8-*-
import sys
import numpy as np
import matplotlib.pyplot as plt
import func as fu


if len(sys.argv) == 1:
    n = 5 # default value
    m = n
    print "Default: N = M = %r " % n
    
elif len(sys.argv) == 2:
    n = int(sys.argv[1])
    m = int(sys.argv[1])
else:
    print """
Usage:
    python morph_experiment.py [side_size]
    
    side_size == N = M -- length of image side (square image)
"""
    sys.exit(-1)
    



# create 2 matrices NxM: lImg, rImg
lImg = np.zeros( (n, m) )
rImg = np.zeros( (n, m) )

# define k, b, dx, dy
origK = 1.00
origB = 5
origDx = 0.2
origDy = 0.5

lImg, rImg = fu.init_images(n, m, origK, origB, origDx, origDy)



A = lImg
B = rImg
# Let lImg == A(x, y), rImg == B(x, y)
# [!] Suppose: 
# B(x, y) = k*A(x+dx, y+dy) + b + dx*k*(dA/dx) + dy*k*(dA/dy)
# [ 2D Tailor's Row Decomposition including linear terms ]


Ax = np.zeros( (n, m) )
Ay = np.zeros( (n, m) )
#Calculating partial derivatives Ax = dA/dx & Ay = dA/dy

Ax, Ay = fu.calc_derivatives_Ax_Ay(A)



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


#Actually, we calculated k^-1 and -b/k
#My fault, now too lazy to rewrite

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
