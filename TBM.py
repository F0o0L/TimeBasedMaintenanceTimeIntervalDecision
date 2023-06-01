from gurobipy import *
import math

def wx(n=1):
    cy = 1
    cs = 20
    cz = 3
    Tt = 1 / 170
    Rmin = 0.8

    a = 4
    mpowa = 0.4  # m^a

    alpha = []
    u = []
    uproduct = []  # Πui
    for i in range(n):
        # alpha.append(i / (15 * i + 5))
        # u.append((17 * i + 1) / (16 * i + 1))
        # alpha.append(0)
        # u.append(1)
        alpha.append(i / (50 * i + 5))
        u.append((50 * i + 1) / (49 * i + 1))
        if i == 0:
            uproduct.append(u[i])
        else:
            uproduct.append(uproduct[i - 1] * u[i])

    # print(alpha)
    # print(uproduct)

    model = Model("model")

    Tk = model.addVars(n, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="Tk")
    T1powa = model.addVar(vtype=GRB.CONTINUOUS, name="T1powa")  # T1^a
    C = model.addVar(vtype=GRB.CONTINUOUS, name="C")
    T = model.addVar(vtype=GRB.CONTINUOUS, name="T")
    Ca = model.addVar(vtype=GRB.CONTINUOUS, name="Ca")
    sumnk = model.addVar(vtype=GRB.CONTINUOUS, name="sumnk")

    if n == 1:
        model.addGenConstrPow(Tk[0], T1powa, a)
        model.addConstr(T1powa / mpowa == sumnk)
        model.addConstr(-T1powa / mpowa >= math.log(Rmin, math.e))
    else:
        model.addGenConstrPow(Tk[0], T1powa, a)
        model.addConstr(-T1powa / mpowa >= math.log(Rmin, math.e))

        temp1 = model.addVars(n - 1, vtype=GRB.CONTINUOUS, name="temp1")  # Tk+ΣαiTi
        temp1powa = model.addVars(n - 1, vtype=GRB.CONTINUOUS, name="temp1powa")  # (Tk+ΣαiTi)^a
        temp2 = model.addVars(n - 1, vtype=GRB.CONTINUOUS, name="temp2")  # ΣαiTi
        temp2powa = model.addVars(n - 1, vtype=GRB.CONTINUOUS, name="temp2powa")  # (ΣαiTi)^a
        temp3 = model.addVars(n - 1, vtype=GRB.CONTINUOUS, name="temp3")  # Tk+2ΣαiTi
        temp3powa = model.addVars(n - 1, vtype=GRB.CONTINUOUS, name="temp3powa")  # (Tk+2ΣαiTi)^a

        for i in range(1, n):
            model.addConstr(Tk[i] + quicksum((alpha[j] * Tk[j]) for j in range(i)) == temp1[i - 1])
            model.addConstr(quicksum(alpha[j] * Tk[j] for j in range(i)) == temp2[i - 1])
            model.addConstr(Tk[i] + 2 * quicksum((alpha[j] * Tk[j]) for j in range(i)) == temp3[i - 1])

            model.addGenConstrPow(temp1[i - 1], temp1powa[i - 1], a)
            model.addGenConstrPow(temp2[i - 1], temp2powa[i - 1], a)
            model.addGenConstrPow(temp3[i - 1], temp3powa[i - 1], a)

            model.addConstr(-1 / mpowa * (temp2powa[i - 1] + uproduct[i - 1] * (
                    temp3powa[i - 1] - math.pow(2, a) * temp2powa[i - 1])) >= math.log(Rmin, math.e))

        model.addConstr(
            1 / mpowa * (T1powa + quicksum(uproduct[i] * (temp1powa[i] - temp2powa[i]) for i in range(n - 1))) == sumnk)

    model.addConstr(cs + cy * (n - 1) + cz * sumnk == C)
    model.addConstr(Tt * sumnk + quicksum(Tk[i] for i in range(n)) == T)
    model.addConstr(T * Ca == C)

    model.setObjective(Ca)

    model.write("wx4.lp")
    model.Params.NonConvex = 2
    model.setParam('OutputFlag', 0)
    model.optimize()

    print(n)
    print('Ca:' + str(Ca.X))
    print('C:' + str(C.X))
    print('T:' + str(T.X))
    sumTk = 0
    Tk1=[]
    for i in range(n):
        # print(Tk[i].X)
        Tk1.append(Tk[i].X)
        sumTk += Tk[i].X
    print('sumTk:' + str(sumTk))
    print('sumnk' + str(sumnk.X))
    print()
    return [n, Ca.X, C.X, sumTk, sumnk.X, Tk1]


n = []
Ca = []
C = []
sumTk = []
sumnk = []
Tk=[]
for i in range(23, 24):
    para = wx(i)
    n.append(para[0])
    Ca.append(para[1])
    C.append(para[2])
    sumTk.append(para[3])
    sumnk.append(para[4])
    Tk.append(para[5])
print(n)
print('Ca'+str(Ca))
print('C'+str(C))
print('sumTk'+str(sumTk))
print('sumnk'+str(sumnk))
for i in range(len(Tk)):
    print('n:'+str(i+1)+'Tk:'+str(Tk[i]))
