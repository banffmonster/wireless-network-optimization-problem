from pulp import *
from math import *
import numpy as np

from matplotlib import pyplot as plt

#initialize the constant
rho = 50*pow(10, -9)
e = 25000
T = 50*24*3600
beta1 = 50*pow(10, -9)
beta2 = 0.0013*pow(10, -12)
alpha = 4

# initialize the variables
f = [[0 for i in range(10)] for n in range(10)]
fB = [0 for i in range(10)]
lamb = [0 for i in range(10)]
lamb[0] = 0
delta = [LpVariable("delta%d" % i, 0) for i in range(11)]

# initialize the node index and coordinates
node = [0,1,2,3,4,5,6,7,8,9]
nodetemp = [0,1,2,3,4,5,6,7,8,9]
L = [1,2,3,4,5,6,7,8,9,10]
nodeh = [[]]
XY = [(400, -320), (300, 440), (-300, -420), (320, -100), (-120, 340), (-500, 100), (-400, 0), (420, 120), (200, 140), (220, -340)]

# calculate the distance and the corresponding power consumption
diB = [sqrt(x*x+y*y) for (x,y) in XY]
dij = [[sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)) for (xi, yi) in XY ] for (xj, yj) in XY]
Cij = [[(beta1 + beta2*pow(d, 4)) for d in di] for di in dij]
CiB = [beta1+beta2*pow(d, 4) for d in diB]

# implement the function of SLP-PA
def LpMAX(S, l):
    # define the problem
    prob = LpProblem("SLP-PA", LpMaximize)
    # lamb[l] = LpVariable("lambda%d"%l, 0)
    # delta = [LpVariable("delta%d" % i, 0) for i in range(11)]

    # objective is delta[l]
    prob += delta[l]
    # fB = [LpVariable("flow%dB"%i, 0) for i in node]
    # f = [[LpVariable("flow%d%d"%(i, j),0) for j in range(10) if j != i] for i in range(10)]

    # create variables for this problem
    for i in node:
        fB[i] = LpVariable("flow%dB"%(i), 0)
        for j in node:
            if j != i:
                f[i][j] = LpVariable("flow%d%d"%(i, j), 0)

    # create constraints
    # nodes that has not reached their optimal rate
    for i in node:
        t = 1
        for n in S:
            if i in n:
                t = 0
                break
        if t == 1:
            powerj = sum(f[i][j] for j in (x for x in range(10) if x != i))
            powerk = sum(f[k][i] for k in (x for x in range(10) if x != i))
            prob += fB[i] + powerj - powerk - delta[l] == lamb[l - 1]

    for i in node:
        t = 1
        for n in S:
            if i in n:
                t = 0
                break
        if t == 1:
            powerj = sum(Cij[i][j]*f[i][j] for j in (y for y in range(10) if y != i))
            powerk = sum(rho*f[k][i] for k in (y for y in range(10) if y != i))
            prob += T*(CiB[i]*fB[i] + powerj + powerk) <= e

    # nodes that have reached their optimal rate are in the corresponding set
    for index, S in enumerate(S):
        for i in node:
            if i in S:
                prob += fB[i] + sum(f[i][j] for j in (y for y in range(10) if y != i)) - sum(f[k][i] for k in (y for y in range(10) if y != i)) == lamb[index]
                prob += T*(sum(Cij[i][j]*f[i][j] for j in (y for y in range(10) if y != i)) + sum(rho*f[k][i] for k in (y for y in range(10) if y != i)) + CiB[i]*fB[i]) == e

    # solve and save the problem
    prob.writeLP("SLP-PA.lp")
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


# implement the function to determine the minimum set
# increase the rate node ln by arg and solve the problem again to calculate the new delta
def GetSet(S, l, ln, arg):
    prob = LpProblem("SLP-PA", LpMaximize)
    # lamb[l] = LpVariable("lambda%d"%l, 0)
    # delta = [LpVariable("delta%d" % i, 0) for i in range(11)]
    prob += delta[l]
    # fB = [LpVariable("flow%dB"%i, 0) for i in node]
    # f = [[LpVariable("flow%d%d"%(i, j),0) for j in range(10) if j != i] for i in range(10)]
    for i in node:
        fB[i] = LpVariable("flow%dB" % (i), 0)
        for j in node:
            if j != i:
                f[i][j] = LpVariable("flow%d%d" % (i, j), 0)

    for i in node:
        t = 1
        for n in S:
            if i in n:
                t = 0
                break
        if t == 1:
            if i == ln:
                powerj = sum(f[i][j] for j in (x for x in range(10) if x != i))
                powerk = sum(f[k][i] for k in (x for x in range(10) if x != i))
                prob += fB[i] + powerj - powerk  == lamb[l - 1] + delta[l] + arg        # for node index ln, increase the rate by arg and calculate delta again
            else:
                powerj = sum(f[i][j] for j in (x for x in range(10) if x != i))
                powerk = sum(f[k][i] for k in (x for x in range(10) if x != i))
                prob += fB[i] + powerj - powerk - delta[l] == lamb[l-1]

    for i in node:
        t = 1
        for n in S:
            if i in n:
                t = 0
                break
        if t == 1:
            powerj = sum(Cij[i][j] * f[i][j] for j in (y for y in range(10) if y != i))
            powerk = sum(rho * f[k][i] for k in (y for y in range(10) if y != i))
            prob += T * (CiB[i] * fB[i] + powerj + powerk) <= e

    for index, S in enumerate(S):
        for i in node:
            if i in S:
                if i == ln:
                    prob += fB[i] + sum(f[i][j] for j in (y for y in range(10) if y != i)) - sum(f[k][i] for k in (y for y in range(10) if y != i)) == lamb[index]
                prob += T * (sum(Cij[i][j] * f[i][j] for j in (y for y in range(10) if y != i)) + sum(rho * f[k][i] for k in (y for y in range(10) if y != i)) + CiB[i] * fB[i]) == e

    prob.writeLP("SLP-PA.lp")
    prob.solve()
    # for i in node:
    #     for v in prob.variables():
    #         if v.name == "flow%dB" % i:
    #             fB[i] = v.varValue
    #     for j in (y for y in node if y != i):
    #         for v in prob.variables():
    #             if v.name == "flow%d%d" % (i, j):
    #                 f[i][j] = v.varValue
    #                 # print("f%d%d"%(i,j), f[i][j])
    return value(prob.objective)


# for i in node:
#     powerj = sum(Cij[i][j]*fj[i][j]* T for j in (y for y in range(10) if y != i))
#     powerk = sum(rho*T*fk[i][k] for k in (y for y in range(10) if  y != i))
#     print(powerj)
    # if CiB[i]*T*fB[i] + powerj + powerk > e:
    #    print(i)
# for i in node:
#     for j in node:

# print("status:", LpStatus[prob.status])
# while node != []:
Rate = []
unsorted_rate = [0 for i in range(10)]
rate = 0

# begin iteration for l = 1, 2,..., n,
for l in L:
    fB, f, maxdelta= LpMAX(nodeh, l)        # get the result for the problem
    rate += maxdelta                        # rate[l] = rate[l-1] + delta[l]
    Rate.append(rate * 0.001)                       # add rate to the rate level set
    nodetemp = []
    for i in node:
        newdelta1 = GetSet(nodeh, l, i, 10.0)   # increase the rate by 10 of node i and get the new delta
        newdelta2 = GetSet(nodeh, l, i, 20.0)   # increase the rate by 20 of node i and get the new delta
        sub = newdelta1 - newdelta2             # calculate the difference between newdelta1 and newdelata2
        #print(sub)
        if sub > 0:
            nodetemp.append(i)                  # if sub > 0 , which means delta increase when rate decrease in node i,
                                                # add i to the optimal set
    nodeh.append(nodetemp)

    # print the results
    print("max additional rate: ", maxdelta)
    print("optimal rate level: ", Rate)
    print("corresponding sets: ", [[i + 1 for i in S] for S in nodeh], "\n")

    lamb[l] = lamb[l-1] + maxdelta          # update lambda
    if sum([len(s) for s in nodeh]) == 10:      # if all 10 nodes have reached their optimal rate, terminate the process
        break

# print the minimum set and the rate level
for num, s in enumerate(nodeh):
    print("Set%d: "%num, [i+1 for i in s])
    if s != []:
        for i in s:
            unsorted_rate[i] = Rate[num-1]
print("unsorted rate level: ", unsorted_rate)

# save the data into datafile and plot it
sorted_rate = sorted(unsorted_rate)
np.savetxt('PA_plot.dat', sorted_rate)
# f = open("plotdata.txt", "a")
# f.write(str(sorted_rate)+"\n")
sorted_index = [i+1 for i in range(10)]

plt.figure()
plt.title("SLP-PA")
plt.xlabel("sorted node index")
plt.ylabel("rate level")
plt.plot(sorted_index, sorted_rate, '-o')
plt.show()








