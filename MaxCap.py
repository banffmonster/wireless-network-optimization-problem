from pulp import *
from math import *
import numpy as np
import matplotlib.pyplot as plt

# initialize the value of the constant needed in this algorithm
rho = 25*pow(10, -9)
e = 25000
T = 50*24*3600
beta1 = 50*pow(10, -9)
beta2 = 0.0013*pow(10, -12)
alpha = 4

# initialize the variables fij, fiB and the objective ri
f = [[0 for i in range(10)] for n in range(10)]
fB = [0 for i in range(10)]
ri = [0 for i in range(10)]

# node index, in python index start from 0
node = [0,1,2,3,4,5,6,7,8,9]

# initialize note coordinates
XY = [(400, -320), (300, 440), (-300, -420), (320, -100), (-120, 340), (-500, 100), (-400, 0), (420, 120), (200, 140), (220, -340)]

# calculate the distance between these two nodes, and nodes to the base staion
diB = [sqrt(x*x+y*y) for (x,y) in XY]
dij = [[sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)) for (xi, yi) in XY ] for (xj, yj) in XY]
# calculate the corresponding power consumption
Cij = [[(beta1 + beta2*pow(d, 4)) for d in di] for di in dij]
CiB = [beta1+beta2*pow(d, 4) for d in diB]

# implement the function of MaxCap algorithm using pulp linear programme package
def LpMAX(node):
    # define the problem
    prob = LpProblem("MaxCap", LpMaximize)

    # create problem variables
    for i in node:
        ri[i] = LpVariable("r%d"%(i), 0)
        fB[i] = LpVariable("flow%dB"%(i), 0)
        for j in node:
            if j != i:
                f[i][j] = LpVariable("flow%d%d"%(i, j), 0)
    # objective function is sum(ri)
    prob += sum(ri)

    # constraints for the objective
    for i in node:
        powerj = sum(f[i][j] for j in (y for y in range(10) if y != i))
        powerk = sum(f[k][i] for k in (y for y in range(10) if y != i))
        prob += fB[i] + powerj - powerk - ri[i] == 0

    for i in node:
        powerj = sum(Cij[i][j]*f[i][j]* T for j in (y for y in range(10) if y != i))
        powerk = sum(rho*T*f[k][i] for k in (y for y in range(10) if  y != i))
        prob += CiB[i]*T*fB[i] + powerj + powerk <= e

    # save the problem date into a file
    prob.writeLP("MaxCap.lp")
    # solve the problem
    prob.solve()
    # get the value of the objective
    for i in node:
        for v in prob.variables():
            if v.name == "r%d"%i:
                ri[i] = v.varValue * 0.001
                print("the rate of node%d:" % (i+1), v.varValue)
    return ri                       # return the objective
    # for v in prob.variables():
    #     print(v.name, "=", v.varValue)
    # for i in node:
    #     for v in prob.variables():
    #         if v.name == "flow%dB" % i:
    #             fB[i] = v.varValue
    #     for j in (y for y in node if y != i):
    #         for v in prob.variables():
    #             if v.name == "flow%d%d" % (i, j):
    #                 f[i][j] = v.varValue
    #                 # print("f%d%d"%(i,j), f[i][j])
    # return fB, f, value(prob.objective)




# solve the problem
unsorted_rate = LpMAX(node)


sorted_index = [1,2,3,4,5,6,7,8,9,10]        # the sorted index from 1 to 10
sorted_rate = sorted(unsorted_rate)            # sort the rate level of these nodes
np.savetxt('MaxCap_plot.dat', sorted_rate)      # save data for plot
# f = open("plotdata.txt", "w+")
# f.write(str(sorted_rate)+"\n")

#plot MaxCap
plt.figure()
plt.title("MaxCap")
plt.xlabel("sorted node index")
plt.ylabel("rate level")
plt.plot(sorted_index, sorted_rate, '-s')
plt.show()







# for v in prob.variables():
#
#     print(v.name, "=", v.varValue)










