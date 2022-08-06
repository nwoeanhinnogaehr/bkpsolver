"""
author: Margarida Carvalho
License: MIT
PYTHON 2
Run CCLW approach

Requires to install knapsack generator: http://hjemmesider.diku.dk/~pisinger/codes.html
In this code the knapsack generator is ./gen2
"""


from time import time
import os
from random import randint,seed
from gurobipy import *

env = Env(empty=True)
env.setParam("OutputFlag", 0)
env.setParam("LogFile", "gurobi.log")
env.setParam("LogToConsole", 0)
env.start()

# make xOpt from the MIP maximal (with respect to the MIP)
# ratio = Pi/wi by decreasing order
def MakeMaximal(ratio,F,c_F,xOpt,n,R,L):
	i = 0
	while c_F>0 and i<=n-1:
		lixo, item = ratio[i]
		if xOpt[item] <0.5: # xOpt=0
			c_F = c_F - F[item]
		i = i+1
	while i<=n-1:
		lixo, item = ratio[i]
		if xOpt[item]==0 and R-L[item]>=0:
			R = R - L[item]
			xOpt[item] = 1
		i = i+1
	return xOpt





# ADD CUT: Follower does not improve constraint
# INPUT
# model - name of the model
# yOpt - we dont want solutions that let these items available
# n - number of items
# k - if k>2 then we can erase the last constraint added
# Profit - best (minor) value for the objective function found
# wmax - item with the largest weight
# OUTPUT
# xnew -  return new solution of the leader
# m  - new model
def NoImpFollower4(m,x,z,yOpt,n,Profit,wmax,k,P,upRow,oldProfit):
	if k>2:
		m.remove(m.getConstrs()[-1])
		m.update()
	# row generation cut
	rowGen = LinExpr(0)
	# add no improvement constraint for the follower
	for i in range(1,n+1):
		if yOpt[i-1]>0.5: # sum_ {i: yi=1} xi >=1
			rowGen = rowGen + P[i-1]*(1-x[i])
	m.addConstr(rowGen <= Profit-1)
	m.update()
	# it remains to update the rhs of rowGen constraints
	if upRow == 1:
		for j in range(k-2): # k-2 cutting plane constraints to update
			m.getConstrs()[-j-2].setAttr("rhs",m.getConstrs()[-j-2].getAttr("rhs")+Profit-oldProfit) # remove the old profit to the new one
			#m.getConstrs()[-j-2].setAttr("rhs",Profit-oldProfit) # remove the old profit to the new one
			m.update()
	m.update()
	# add improve4 idea
	m.addConstr( m.getObjective()-z[1]*wmax<= Profit-1)
	m.update()
#	# save problem - OPTIONAL
#	m.write('current.lp')
	# solve
	m.optimize()
	# show solution x and
	xnew=[] # when it is empty we know that the problem is infeasible
	value = 0
	if m.Status == GRB.OPTIMAL:
		value = m.ObjVal
		for i in range(1,n+1):
			xnew.append(x[i].x)
	if m.status == 9:
		xnew = 0 # time limit reached
	return xnew,m,value



# Just stops when it has found a solution (i.e., the problem becames infeasible)
# using bounds for the followers profit (using Nogood function)
# INPUT
# m - model
# xOpt - current best solution
# OUTPUT
# iteraOpt - iteration in which the optimal solution is found
# iteraProve - number of iteration needed to prove optimality
# PLUS XOPT ALWAYS MAXIMAL WITH RESPECT TO THE MIP SOLUTION
# PLUS cut that avoid computing bi-level feasible solutions uninteresting (difference from improve2)
def EnumStopImproved4(m,x,z,xOpt,n,P,F,c_F,L,c_L,ti):
	first_MIP_value = m.ObjVal
	MIP_value_anterior = first_MIP_value
	MIP_time = m.Runtime # running time of all solved MIPs
	first_MIP_time = MIP_time
	worst_MIP_time = MIP_time
	first_MIP_nodes = m.NodeCount
	worst_MIP_nodes = first_MIP_nodes
	KP_time = 0 # running time of all solved KP
	ratio =[((P[i]*1.)/F[i],i) for i in range(n)]
	ratio.sort()
	ratio.reverse()
	# make xOpt maximal with respect to the mip solution
	value = WeightInterdict(F,L,c_L,n)
	xOpt = MakeMaximal(ratio,F,c_F,xOpt,n,m.getConstrs()[n*2+1].getAttr(GRB.Attr.Slack),L)
	wmax,pmax = stopCritImproved(F,c_F,n,ratio,value,P)
	#pmax,wmax = stopCrit(P,F,c_F,n,ratio)
	#pmax = max(P)
	#wmax = max(F)
	yOpt,Profit,FollowerTime = ReactionFollower(n,P,F,c_F,xOpt,1) #maybe here we can opt for a differente algorithm
	KP_time = KP_time+FollowerTime
	x_new = xOpt
	y_new = yOpt
	s = 1
	i=2
	iteraOpt = 1 # step in which we find the best solution
	upRow = 1 # do not update rowGen constraints rhs, if 1 update
	oldProfit = Profit
	while s==1:
		x_new,m,LastMIPValue = NoImpFollower4(m,x,z,y_new,n,Profit,wmax,i,P,upRow,oldProfit)	#to add the constraint sum_ {i: yi=1} xi >=1
		upRow = 0
		if m.Runtime>worst_MIP_time:
			worst_MIP_time=m.Runtime
			worst_MIP_nodes = m.NodeCount
		MIP_time = MIP_time + m.Runtime
#		fnovo = open('Evolucao4.txt','a')
#		fnovo.write(str(m.ObjVal)+' ##########   '+str(Profit)+'\n')
#		fnovo.close()
		if x_new ==[] or x_new==0 or time()-ti>3600: # the last model m is an infeasible problem
			if x_new ==[]:
				#print('\nStopped in iteration ', i,'\nTHE LAST MODEL BUILT IS INFEASIBLE')
				s = 0
				LastMIPValue = MIP_value_anterior
			else:
				#print('\nStopped in iteration ', i,'\nTIME LIMIT REACHED')
				s = 0
				LastMIPValue = MIP_value_anterior
				xOpt = []
		else:
			if Profit + pmax <= LastMIPValue:
				#print('\nStopped in iteration ', i,'\nIT IS NOT POSSIBLE TO IMPROVE')
				s = 0
			else:
				x_new = MakeMaximal(ratio,F,c_F,x_new,n,m.getConstrs()[n*2+1].getAttr(GRB.Attr.Slack),L)
				y_new,profit_new,FollowerTime = ReactionFollower(n,P,F,c_F,x_new,1)
				KP_time = KP_time+FollowerTime
				if profit_new<Profit:
					upRow = 1
					oldProfit = Profit
					Profit = profit_new
					xOpt = x_new
					yOpt = y_new
					iteraOpt = i
		MIP_value_anterior = LastMIPValue
		i = i+1
	# iteraProve = i-1
	return xOpt,yOpt,Profit,iteraOpt,i-1,LastMIPValue,MIP_time,KP_time,first_MIP_time,first_MIP_nodes,worst_MIP_time,worst_MIP_nodes,first_MIP_value

# INPUT
# n - number of items
# P - list with the items profit
# F - list with the weights of the items
# c_F - follower capacity
# xOpt - strategy of the leader
# Bi - 1 if the variables are binary and 0 otherwise
# OUTPUT
# yOpt -  best reaction of the follower
# profit
def ReactionFollower(n,P,F,c_F,xOpt,Bi):
	Reaction = Model("ReactionF", env=env)
	y = {}
	for i in range(1,n+1):
		# create decision variables and objective function
		if Bi == 1:
			y[i] = Reaction.addVar(obj=P[i-1],vtype="B",name="yf"+str(i))
		else:
			y[i] = Reaction.addVar(obj=P[i-1],vtype="C",name="yf"+str(i))
		Reaction.update()
		Reaction.addConstr(y[i]>=0)
		Reaction.addConstr(xOpt[i-1]+y[i]<=1)
	Reaction.update()
	# constrain
	Reaction.addConstr(quicksum(F[i]*y[i+1] for i in range(n))<=c_F)
	Reaction.ModelSense = -1 # maximize
	Reaction.update()
	# solve
	Reaction.optimize()
	# view solution
	profit = Reaction.ObjVal
	yOpt=[]
	for i in range(1,n+1):
		yOpt.append(y[i].x)
	return yOpt, profit,Reaction.Runtime

############## PRE-PROCESSING

def stopCrit(F,c_F,n,ratio,P):
	C = 0
	i = 0
	while C <= c_F:
		a,b = ratio[i]
		C = C+ F[b]
		i = i+1
	wmax = 0
	pmax = 0
	for j in range(i-1,n):	# items that may be critical
		a,b = ratio[j]
		if F[b]> wmax:
			wmax = F[b]
		if P[b] >pmax:
			pmax = P[b]
	return wmax,pmax

############## PRE-PROCESSING Improved
# ratio (of the Follower) - sorted by decreasing order

def stopCritImproved(F,c_F,n,ratio,value,P):
	aux1,aux2 = ratio[-1]
	if sum(F)>c_F+value+F[aux2]:# F[aux2]=wn - item with smaller ratio
		C = 0
		i = 0
		while C <= c_F:
			a,b = ratio[i]
			C = C+ F[b]
			i = i+1
		wmax = F[b]
		pmax = P[b]
		a,b = ratio[i]
		C = C+F[b]
		while C-F[b]<c_F+value:
			if F[b]> wmax:
				wmax = F[b]
			if P[b]>pmax:
				pmax = P[b]
			i = i+1
			a,b = ratio[i]
			C = C+F[b]
	else:
		wmax,pmax = stopCrit(F,c_F,n,ratio,P) #finalitem is the one above which the items are never critical
	return wmax,pmax

### computation of the maximum weight of the follower that the leader can interdict (relaxed version)

def WeightInterdict(F,L,c_L,n):
	value = 0
	# sort the items by their value wi/vi
	ratio1 = list(range(n))
	ratio1.sort(key=lambda i: (-F[i]/float(L[i])))
	R = c_L
	aux3 = 0
	while R>=0 and aux3<=n-1:
		item = ratio1[aux3]
		R = R-L[item]
		value = value+F[item]
		aux3 = aux3+1
	if R < 0:
		R = R +L[item]
		value = value-F[item]+ F[item]*(c_L-R)/float(L[item])
	return int(value)


def Run_CCLW(n,P,F,c_F,L,c_L):
    ti = time()
####################################### FIRST ITERATIONS
    setParam("TimeLimit", 3600)


    # BUILT the first MIP
    LeaderRelax = Model("Leader_Relax", env=env)

    x,z = {},{}
    z[1] = LeaderRelax.addVar(obj = c_F,vtype="C",name="z1")
    LeaderRelax.update()
    LeaderRelax.addConstr(z[1]>=0)
    for j in range(1,n+1):
    	x[j] = LeaderRelax.addVar(vtype="B",name="x"+str(j))
    	z[j+1] = LeaderRelax.addVar(obj = 1,vtype="C",name="z"+str(j+1))
    	LeaderRelax.update()
    	LeaderRelax.addConstr(z[j+1]>=0)
    	# these dual constraints were updated to the new formulation of the followers primal problem
    	LeaderRelax.addConstr(F[j-1] * z[1] + z[j+1]>=P[j-1]*(1-x[j]))
    	LeaderRelax.update()
    LeaderRelax.update()

    # knapsasck constraint
    LeaderRelax.addConstr(quicksum(L[i]*x[i+1] for i in range(n))<=c_L)
    LeaderRelax.update()

    # OBJECTIVE FUNCTION : add xi to the objective function in order to have a maximal solution
    #LeaderRelax.setObjective(z[1]*c_F+sum(u[i] for i in range(1,n+1)),GRB.MINIMIZE)
    LeaderRelax.ModelSense = 1 # to minimize
    LeaderRelax.update()

    # save problem - OPTIONAL
    # LeaderRelax.write('Experiencia2.lp')

    # solve
    LeaderRelax.optimize()
    # save solution file
    # LeaderRelax.write('Experiencia2.sol')

    if LeaderRelax.status != 9: # Optimization terminated because the time exceeded the value specified in the TimeLimit parameter
    	# x optimal:
    	xOpt=[]
    	for j in range(1,n+1):
    			xOpt.append(x[j].x)
################################# end of first iteration
    if LeaderRelax.status != 9:
        xOpt,yOpt,Profit,iteraOpt,iteraProve,LastMIPValue,MIP_time,KP_time,first_MIP_time,first_MIP_nodes,worst_MIP_time,worst_MIP_nodes,first_MIP_value = EnumStopImproved4(LeaderRelax,x,z,xOpt,n,P,F,c_F,L,c_L,ti)
    else:
        xOpt = []
        timecpu = time()-ti
    return xOpt,yOpt,Profit,iteraOpt,iteraProve,LastMIPValue,MIP_time,KP_time,first_MIP_time,first_MIP_nodes,worst_MIP_time,worst_MIP_nodes,first_MIP_value


if __name__ == "__main__":

    seed(1)

    ################## Run a single instance #####################################
    n = int(input())
    c_F = int(input())
    c_L = int(input())
    F = list(map(int, input().split()))
    L = list(map(int, input().split()))
    P = list(map(int, input().split()))
    xOpt,yOpt,Profit,iteraOpt,iteraProve,LastMIPValue,MIP_time,KP_time,first_MIP_time,first_MIP_nodes,worst_MIP_time,worst_MIP_nodes,first_MIP_value = Run_CCLW(n,P,F,c_F,L,c_L)
    #print("yOpt: ", yOpt)
    print(int(Profit))

    ################## Run many instances #########################################

    # fResults = open('Results.txt','a')
    # fResults.write(' ins & FirstMIP.Obj & LastMIP.Obj & OPT & cpu_time\n')
    # fResults.close()
    #
    # #for n in [10,15,20,25,30,35,40]:
    # for n in [10,15,20]:
    #     for t in [1,2,3,4,5]:
    #         for i in range(1,101):
    #             # generate instance
    #             f = open('trash.txt','w')
    #             f.write(str(n)+'\n') # define n
    #             r = 100 # coeficients
    #             f.write(str(r)+'\n')
    #             f.write(str(t)+'\n')
    #             f.write(str(i)+'\n') # instance number
    #             f.write(str(1000))	# number of tests
    #             f.close()
    #             # Built instance using gen2.c for the follower
    #             os.system("./gen2<trash.txt")
    #             # read input file with the followers instances: test.in
    #             f = open('test.in','r')
    #             n = int(f.readline()) # number of items
    #             P_ord = [] # benefits - items values
    #             L_ord = [] # leader weights
    #             F_ord = [] # followers weights
    #             # formulate the knapsack problem for F
    #             for j in range(n):
    #                 a,p,wf = map(int,f.readline().split())
    #                 P_ord.append(p)
    #                 F_ord.append(wf)
    #                 wl = randint(1,r) # r is given in trash.txt
    #                 L_ord.append(wl)
    #             # capacity
    #             c_F = int(f.readline())
    #             f.close()
    #             c_L = randint(c_F-10,c_F+10)
    #
    #             xOpt,yOpt,Profit,iteraOpt,iteraProve,LastMIPValue,MIP_time,KP_time,first_MIP_time,first_MIP_nodes,worst_MIP_time,worst_MIP_nodes,first_MIP_value = Run_CCLW(n,P,F,c_F,L,c_L)
    #             fResults = open('Results.txt','a')
    #             fResults.write(str(i)+' & '+ str(first_MIP_value)+' & '+str(LastMIPValue)+' & '+str(Profit)+' & '+str(timecpu)+'\\\ \n\hline \n')
    #             fResults.close()
