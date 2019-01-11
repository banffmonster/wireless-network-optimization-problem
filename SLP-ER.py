from pulp import *
from math import *
import numpy as np

import matplotlib.pyplot as plt

# initialize the constants needed in this problem
rho = 50*pow(10, -9)
e = [25000 for i in range(10)]
T = 50*24*3600
beta1 = 50*pow(10, -9)
beta2 = 0.0013*pow(10, -12)
alpha = 4

#initialize the variables
f = [[0 for i in range(10)] for n in range(10)]
fB = [0 for i in range(10)]
fk = [[0 for i in range(10)] for n in range(10)]

# node index and coordinates
node = [0,1,2,3,4,5,6,7,8,9]
XY = [(400, -320), (300, 440), (-300, -420), (320, -100), (-120, 340), (-500, 100), (-400, 0), (420, 120), (200, 140), (220, -340)]

# calculate the distance in this two nodes and the corresponding power consumption
diB = [sqrt(x*x+y*y) for (x,y) in XY]
dij = [[sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)) for (xi, yi) in XY ] for (xj, yj) in XY]
Cij = [[(beta1 + beta2*pow(d, 4)) for d in di] for di in dij]
CiB = [beta1+beta2*pow(d, 4) for d in diB]

# implement the function of SLP-ER algorithm
def LpMAX(node):
    # define the problem
    prob = LpProblem("SLP-ER", LpMaximize)

    # objective of this problem is lambda denoted as l
    l = LpVariable("lambda", 0)
    prob += l

    # create flow variables for this problem
    for i in node:
        fB[i] = LpVariable("flow%dB"%(i), 0)
        for j in node:
            if j != i:
                f[i][j] = LpVariable("flow%d%d"%(i, j), 0)

    # create constraints
    for i in node:
        powerj = sum(f[i][j] for j in (y for y in range(10) if y != i))
        powerk = sum(f[k][i] for k in (y for y in range(10) if y != i))
        prob += fB[i] + powerj - powerk - l == 0

    for i in node:
        powerj = sum(Cij[i][j]*f[i][j]* T for j in (y for y in range(10) if y != i))
        powerk = sum(rho*T*f[k][i] for k in (y for y in range(10) if y != i))
        prob += CiB[i]*T*fB[i] + powerj + powerk <= e[i]

    # write the problem data into file and solve it
    prob.writeLP("SLP-ER.lp")
    prob.solve()

    # get the value of the objective and flow
    for i in node:
        for v in prob.variables():
            if v.name == "flow%dB" % i:
                fB[i] = v.varValue
        for j in (y for y in node if y != i):
            for v in prob.variables():
                if v.name == "flow%d%d" % (i, j):
                    f[i][j] = v.varValue
                    # print("f%d%d"%(i,j), f[i][j])
    return fB, f, value(prob.objective)


# initialize the rate level
rate = [0 for i in range(10)]

# perform iteration
while node != []:
    fB, f, maxlabda= LpMAX(node)        # get the result for the problem
    print("\nmax lambda: ", maxlabda)
    nodetemp = node
    for i in node:
        rate[i] = maxlabda * 0.001  # b/s to kb/s
    node =  []

    #calculate the remaining energy for every node
    for i in nodetemp:
        powerj = sum(Cij[i][j] * f[i][j] * T for j in (y for y in range(10) if y != i))
        powerk = sum(rho * T * f[k][i] for k in (y for y in range(10) if y != i))
        energy = CiB[i] * T * fB[i] + powerj + powerk
        powerj = sum(f[i][j] for j in (y for y in range(10) if y != i))
        powerk = sum(f[k][i] for k in (y for y in range(10) if y != i))
        fh = fB[i] + powerj - powerk
        #print("energy consumed in node%d:" % (i + 1), energy)
        #print("optimal rate of node%d:"%(i+1), fh)

        # if there are nodes that still have remaining energy, add them to the next iteration set
        if energy < 0.9999*e[i]:
             node.append(i)
             #e[i] = e[i]- energy
             #print(e)
    set = [item for item in nodetemp if item not in node]
    print("optimal set: ", [i+1 for i in set])              # print the node set that reach their optimal rate

# print the optimal rate for every node
for i, r  in enumerate(rate):
    print("optimal rate of node%d i is:"%i, r)


# get the sorted rate and save it to datafile
sorted_index = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
sorted_rate = sorted(rate)
np.savetxt('ER_plot.dat', sorted_rate)
# f = open("plotdata.txt", "a")
# f.write(str(sorted_rate)+"\n")
# plot the Rate allocation of SLP-ER
plt.figure()
plt.title("SLP-ER")
plt.xlabel("sorted node index")
plt.ylabel("rate level")
plt.plot(sorted_index, sorted_rate, '-D')
plt.show()



# for v in prob.variables():
#
#     print(v.name, "=", v.varValue)










