from pulp import *
from math import *
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
# this is the circle partition plot code for algorithm 8.1



# initialize the constant for this algorithm
alpha = 2
b1 = 1
b2 = 0.5
rho = 1
epsilon = 0.2

# initialize the node coordinates
XY = [(0.1, 0.5), (1.1, 0.7), (0.4, 0.1)]
(x1, y1) = XY[0]
(x2, y2) = XY[1]
(x3, y3) = XY[2]
# D = [sqrt((x-OA[0])*(x-OA[0])+(y-OA[1])*(y-OA[1])) for (x, y) in XY]


# write a function to calculate the origin and radius to identify SED
def circle(x1, y1, x2, y2, x3, y3):
    a = x1 - x2
    b = y1 - y2
    c = x1 - x3
    d = y1 - y3
    a1 = ((x1 * x1 - x2 * x2) + (y1 * y1 - y2 * y2)) / 2.0
    a2 = ((x1 * x1 - x3 * x3) + (y1 * y1 - y3 * y3)) / 2.0
    theta = b * c - a * d
    x0 = (b * a2 - d * a1) / theta
    y0 = (c * a1 - a * a2) / theta
    r = sqrt(pow((x1 - x0), 2) + pow((y1 - y0), 2))
    return x0, y0, r


# get the origin and radius of SED
x0, y0, R = circle(x1, y1, x2, y2, x3, y3)
# calculate the lower and upper bound on CiB for each node i
D = [sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)) for (x, y) in XY]
Cmin = [b1 for xy in XY]
Cmax = [(b1 + b2 * pow((d + R), alpha))for d in D]

# plot the circle of SED with origin OA####################################################
fig = plt.figure()
ax = fig.add_subplot(111)
cir = Circle(xy = (x0, y0), radius = R, color='r',linestyle= '-',alpha=0.4, fill=False)
ax.add_patch(cir)
plt.axis('scaled')
plt.axis('equal')
plt.plot(x1,y1,'r*',x2,y2,'r*',x3,y3,'r*')

##############################################################################################

# draw circles centered at each node i, with cost C[i]###################################

# find Hi
h = int(round((log(1 + (b2/b1 )* pow((R + R), 2)))/(log(1+epsilon)))+1)
# find C[iFor]
C = [0 for i in range(h)]
for i in range(h):
    C[i] = b1 * pow((b1 + epsilon),i+1)
# calculate the corresponding radius for each cost C[i]
Rs = [sqrt((c - b1)/b2) for c in C]
for r in Rs:
    for (x, y) in XY:
        c = Circle(xy = (x, y), radius=r, linestyle = ':', alpha=0.4, fill = False)
        ax.add_patch(c)
        plt.axis('scaled')
        plt.axis('equal')
ax.set_xlim(0,1.5)
ax.set_ylim(0,1.5)
plt.show()
####################################################################################

print(C)
print(h)
print(x0, y0, R)
#print(x0, y0, R)


