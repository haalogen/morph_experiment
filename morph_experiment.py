#-*-coding: utf-8-*-
import sys
import numpy as np
import matplotlib.pyplot as plt


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

# choose intervals:[-1; 1] for example
lX = np.linspace(-1, 1, n)
lY = np.linspace(-1, 1, n)
print lX, '\n'

rX = np.linspace(-0.8, 1.2, n)
rY = np.linspace(-1, 1, n)
print rX, '\n'



# creating meshgrid
lXX, lYY = np.meshgrid(lX, lY)
rXX, rYY = np.meshgrid(rX, rY)

lImg = np.e ** (- lXX**2 - lYY**2)
rImg = np.e ** (- rXX**2 - rYY**2)


print lImg, '\n'
print rImg,'\n'

#pcolormesh -- colormap
#cmap = plt.cm.Greys -- for grayscale
colorMap = plt.cm.gray
plt.subplot(1, 2, 1)
plt.imshow(lImg, cmap=colorMap, interpolation='none')

plt.subplot(1, 2, 2)
plt.imshow(rImg, cmap=colorMap, interpolation='none')


plt.show()
