import scipy.io as sci
import matplotlib.pyplot as plt
import os
import numpy as np

mn1 = 20
mx1 = 200
mn2 = 40
mx2 = 200
cut = 250

for num in range(15):
    win1 = np.arange(20, 530, 30)
    width1 = len(win1)
    win2 = np.arange(100, 550, 30)
    width2 = len(win2)
    win3 = np.arange(40, 1060, 60)
    width3 = len(win3)
    win4 = np.arange(200, 1100, 60)
    width4 = len(win4)
    win5 = np.arange(60, 1590, 90)
    width5 = len(win5)
    win6 = np.arange(300, 1650, 90)
    width6 = len(win6)

    print(num)              ## print processing number
    PATH = 'C:\\Users\\hirokane\\PycharmProjects\\Wheel_Analyze\\wateranalysis_Agra'
    os.chdir(PATH)
    matdata = sci.loadmat(str(num+1))    #create file name
    Data = matdata["data"]
    peg = matdata["RLPegTouchAll"]
    drink = matdata["DrinkOnArray"]
    oneturn = Data[:, 4][np.nonzero(Data[:, 4])]

    n = (len(peg))-1
    dr = (len(drink))-1
    s = np.zeros(n+1)

    A = np.zeros(width1)
    B = np.zeros(width1)
    C = np.zeros(width1)
    D = np.zeros(width1)
    X = np.zeros(width1)
    Y = np.zeros(width1)

    turn = 1
    time = oneturn[0]

    for i in range(1, n-30):  # counting putative rpeg touch number
        s[i] = peg[i] - peg[i-1]
        if peg[i] > time:
            time = time + oneturn[turn]
            turn = turn + 1
        if s[i] < 600:
            for k in range(0, width1-1):
                counter0 = 0
                counter1 = 0
                counter2 = 0
                if win1[k] < s[i] < win1[k + 1]:  # filtering interval
                    for j in range(2, dr):
                        if ((peg[i-1]) < drink[j] < (peg[i]))and(mn1 < drink[j+2]-drink[j+1] < mx1)and\
                                (mn1 < drink[j+1]-drink[j] < mx1)and(mn1 < drink[j]-drink[j-1] < mx1)and\
                                (mn1 < drink[j-1]-drink[j-2] < mx1):
                            Y[k] = Y[k] + 1
                            counter0 = counter0 + 1

                            if oneturn[turn] < oneturn[turn-1]: ###どんどん速くなっていている時
                                C[k] = C[k] + 1
                                counter1 = counter1 + 1
                            else:
                                D[k] = D[k] + 1
                                counter2 = counter2 + 1
                if 0.5 < counter1:
                    A[k] = A[k] + 1
                if 0.5 < counter2:
                    B[k] = B[k] + 1
                if 0.5 < counter0:
                    X[k] = X[k] + 1

    for i in range(width1):
        if not(A[i] == 0):
            C[i] = (C[i]/A[i])
        else:
            C[i] = 0
    for i in range(width1):
        if not(B[i] == 0):
            D[i] = (D[i]/B[i])
        else:
            D[i] = 0
    for i in range(width1):
        if not(X[i] == 0):
            Y[i] = (Y[i]/X[i])
        else:
            Y[i] = 0

    plt.figure(num=1)
    win1 = np.delete(win1, -1)
    C = np.delete(C, -1)
    D = np.delete(D, -1)
    Y = np.delete(Y, -1)
    x = win1
    y = Y
    y1 = C
    y2 = D

    plt.plot(x, y1, label="Speed up", color="r")
    plt.plot(x, y2, label="Speed down", color="b")
    plt.title("Drink_Histogram")
    plt.ylim(0, 6)
    plt.xlabel("interval(ms)")
    plt.ylabel("drink frequency (time/interval)")

    plt.figure(num=2)
    plt.plot(x, y, color="b")
    plt.title("Drink_Histogram_All")
    plt.ylim(0, 6)
    plt.xlabel("interval(ms)")
    plt.ylabel("drink frequency (time/interval)")

    #### ｎ歩ごとに区切って何回水飲みが入ったかを出す。n = 2,3 (5もいるか?)

    A = np.zeros(width3)
    B = np.zeros(width3)
    C = np.zeros(width3)
    D = np.zeros(width3)
    X = np.zeros(width3)
    Y = np.zeros(width3)

    turn = 1
    time = oneturn[0]

    for i in range(2, n - 30, 2):  # counting putative rpeg touch number
        s[i] = peg[i] - peg[i-2]
        if peg[i] > time:
            time = time + oneturn[turn]
            turn = turn + 1
        if s[i] < 1200:
            for k in range(0, width3-1):
                counter0 = 0
                counter1 = 0
                counter2 = 0
                if win3[k] < s[i] < win3[k+1]:  # filtering interval
                    for j in range(2, dr):
                        if ((peg[i - 2]) < drink[j] < (peg[i])) and (mn1 < drink[j + 2] - drink[j + 1] < mx1) and \
                                (mn1 < drink[j + 1] - drink[j] < mx1) and (mn1 < drink[j] - drink[j - 1] < mx1) and \
                                (mn1 < drink[j - 1] - drink[j - 2] < mx1):
                            Y[k] = Y[k] + 1
                            counter0 = counter0 + 1

                            if oneturn[turn] < oneturn[turn - 1]:  ###どんどん速くなっていている時
                                C[k] = C[k] + 1
                                counter1 = counter1 + 1
                            else:
                                D[k] = D[k] + 1
                                counter2 = counter2 + 1
                if 0.5 < counter1:
                    A[k] = A[k] + 1
                if 0.5 < counter2:
                    B[k] = B[k] + 1
                if 0.5 < counter0:
                    X[k] = X[k] + 1

    for i in range(width3):
        if not (A[i] == 0):
            C[i] = (C[i] / A[i])
        else:
            C[i] = 0
    for i in range(width3):
        if not (B[i] == 0):
            D[i] = (D[i] / B[i])
        else:
            D[i] = 0
    for i in range(width3):
        if not (X[i] == 0):
            Y[i] = (Y[i] / X[i])
        else:
            Y[i] = 0


    plt.figure(num=5)
    win3 = np.delete(win3, -1)
    C = np.delete(C, -1)
    D = np.delete(D, -1)
    Y = np.delete(Y, -1)
    x = win3
    y = Y
    y1 = C
    y2 = D

    plt.plot(x, y1, label="Speed up", color="r")
    plt.plot(x, y2, label="Speed down", color="b")
    plt.title("Drink_Histogram (2 Step)")
    plt.ylim(0, 10)
    plt.xlabel("interval(ms)")
    plt.ylabel("drink frequency (time/interval)")

    plt.figure(num=6)
    plt.plot(x, y, color="b")
    plt.title("Drink_Histogram_All(2 Step)")
    plt.ylim(0, 10)
    plt.xlabel("interval(ms)")
    plt.ylabel("drink frequency (time/interval)")

    A = np.zeros(width5)
    B = np.zeros(width5)
    C = np.zeros(width5)
    D = np.zeros(width5)
    X = np.zeros(width5)
    Y = np.zeros(width5)

    turn = 1
    time = oneturn[0]

    for i in range(3, n - 30, 3):  # counting putative rpeg touch number
        s[i] = peg[i] - peg[i - 3]
        if peg[i] > time:
            time = time + oneturn[turn]
            turn = turn + 1
        if s[i] < 1800:
            for k in range(0, width5 - 1):
                counter0 = 0
                counter1 = 0
                counter2 = 0
                if win5[k] < s[i] <win5[k + 1]:  # filtering interval
                    for j in range(2, dr):
                        if ((peg[i - 3]) < drink[j] < (peg[i])) and (mn1 < drink[j + 2] - drink[j + 1] < mx1) and \
                                (mn1 < drink[j + 1] - drink[j] < mx1) and (mn1 < drink[j] - drink[j - 1] < mx1) and \
                                (mn1 < drink[j - 1] - drink[j - 2] < mx1)and(mn1 < drink[j - 2] - drink[j - 3] < mx1):
                            Y[k] = Y[k] + 1
                            counter0 = counter0 + 1

                            if oneturn[turn] < oneturn[turn - 1]:  ###どんどん速くなっている時
                                C[k] = C[k] + 1
                                counter1 = counter1 + 1
                            else:
                                D[k] = D[k] + 1
                                counter2 = counter2 + 1
                if 0.5 < counter1:
                    A[k] = A[k] + 1
                if 0.5 < counter2:
                    B[k] = B[k] + 1
                if 0.5 < counter0:
                    X[k] = X[k] + 1

    for i in range(width5):
        if not (A[i] == 0):
            C[i] = (C[i] / A[i])
        else:
            C[i] = 0
    for i in range(width5):
        if not (B[i] == 0):
            D[i] = (D[i] / B[i])
        else:
            D[i] = 0
    for i in range(width5):
        if not (X[i] == 0):
            Y[i] = (Y[i] / X[i])
        else:
            Y[i] = 0

    plt.figure(num=9)
    win5 = np.delete(win5, -1)
    C = np.delete(C, -1)
    D = np.delete(D, -1)
    Y = np.delete (Y, -1)
    x = win5
    y = Y
    y1 = C
    y2 = D

    plt.plot(x, y1, label="Speed up", color="r")
    plt.plot(x, y2, label="Speed down", color="b")
    plt.title("Drink_Histogram (3 Step)")
    plt.ylim(0, 15)
    plt.xlabel("interval(ms)")
    plt.ylabel("drink frequency (time/interval)")

    plt.figure (num=10)
    plt.plot (x, y, color="b")
    plt.title ("Drink_Histogram_All(3 Step)")
    plt.ylim (0, 15)
    plt.xlabel ("interval(ms)")
    plt.ylabel ("drink frequency (time/interval)")

####ここまで水飲みの階段解析
########################################################################################################################

    E = np.zeros (width2)
    F = np.zeros (width2)
    Z = np.zeros (width2)
    drinkhist0 = np.zeros (width2)
    drinkhist1 = np.zeros (width2)
    drinkhist2 = np.zeros (width2)
    turn = 1
    time = oneturn[0]
    for i in range(1, n-30):  # counting putative rpeg touch number
        s[i] = peg[i] - peg[i-1]
        if peg[i] > time:
            time = time + oneturn[turn]
            turn = turn + 1
        for k in range(width2-1):
            if win2[k] < s[i] < win2[k + 1]:  # filtering interval
                for j in range(dr-1):
                    if ((peg[i-1]) < drink[j])and(drink[j] < peg[i]):
                        if (mn2 < drink[j+2]-drink[j+1] < mx2)and(mn2 < drink[j+1]-drink[j] < mx2)and\
                                (mn2 < drink[j]-drink[j-1] < mx2)and(mn2 < drink[j-1]-drink[j-2] < mx2):
                            Z[k] = Z[k] + 1
                            drinkhist0[k] = drinkhist0[k] + (drink[j+1]-drink[j])
                            if oneturn[turn] < oneturn[turn-1]:
                                E[k] = E[k] + 1
                                drinkhist1[k] = drinkhist1[k] + (drink[j+1]-drink[j])
                            else:
                                F[k] = F[k] + 1
                                drinkhist2[k] = drinkhist2[k] + (drink[j+1]-drink[j])

    for i in range(width2-1):
        if not(E[i] == 0):
            drinkhist1[i] = (drinkhist1[i]/E[i])
        else:
            drinkhist1[i] = 0
    for i in range(width2-1):
        if not(F[i] == 0):
            drinkhist2[i] = (drinkhist2[i]/F[i])
        else:
            drinkhist2[i] = 0
    for i in range(width2-1):
        if not(Z[i] == 0):
            drinkhist0[i] = (drinkhist0[i]/Z[i])
        else:
            drinkhist0[i] = 0

    plt.figure(num=3)
    x = win2
    y1 = drinkhist1
    y2 = drinkhist2
    plt.plot(x, y1, label="Speed up", color="r")
    plt.plot(x, y2, label="Speed down", color="b")
    plt.title("Drink interval histogram")
    plt.xlabel("Running interval (ms)")
    plt.ylabel("Averaged drink interval (ms)")
    plt.ylim(0, 400)

    plt.figure(num=4)
    x = win2
    y = drinkhist0
    plt.plot(x, y, label="Through trial")
    plt.title("Drink interval histogram_All")
    plt.xlabel("Running interval (ms)")
    plt.ylabel("Averaged drink interval (ms)")
    plt.ylim(0, 400)


    E = np.zeros (width4)
    F = np.zeros (width4)
    Z = np.zeros (width4)
    drinkhist0 = np.zeros (width4)
    drinkhist1 = np.zeros (width4)
    drinkhist2 = np.zeros (width4)
    turn = 1
    time = oneturn[0]
    for i in range(2, n - 30, 2):  # counting putative rpeg touch number
        s[i] = peg[i] - peg[i - 1]
        if peg[i] > time:
            time = time + oneturn[turn]
            turn = turn + 1
        for k in range (width4 - 1):
            if win4[k] < s[i] < win4[k + 1]:  # filtering interval
                for j in range (dr - 1):
                    if (peg[i - 2]) < drink[j] < (peg[i]):
                        if (mn2 < drink[j + 2] - drink[j + 1] < mx2) and (mn2 < drink[j + 1] - drink[j] < mx2) and \
                                (mn2 < drink[j] - drink[j - 1] < mx2) and (mn2 < drink[j - 1] - drink[j - 2] < mx2):
                            Z[k] = Z[k] + 1
                            drinkhist0[k] = drinkhist0[k] + (drink[j + 1] - drink[j])
                            if oneturn[turn] < oneturn[turn - 1]:
                                E[k] = E[k] + 1
                                drinkhist1[k] = drinkhist1[k] + (drink[j + 1] - drink[j])
                            else:
                                F[k] = F[k] + 1
                                drinkhist2[k] = drinkhist2[k] + (drink[j + 1] - drink[j])

    for i in range (width4 - 1):
        if not (E[i] == 0):
            drinkhist1[i] = (drinkhist1[i] / E[i])
        else:
            drinkhist1[i] = 0
    for i in range (width4 - 1):
        if not (F[i] == 0):
            drinkhist2[i] = (drinkhist2[i] / F[i])
        else:
            drinkhist2[i] = 0
    for i in range (width4 - 1):
        if not (Z[i] == 0):
            drinkhist0[i] = (drinkhist0[i] / Z[i])
        else:
            drinkhist0[i] = 0

    plt.figure (num=7)
    x = win4
    y1 = drinkhist1
    y2 = drinkhist2
    plt.plot (x, y1, label="Speed up", color="r")
    plt.plot (x, y2, label="Speed down", color="b")
    plt.title ("Drink interval histogram(2 Step)")
    plt.xlabel ("Running interval (ms)")
    plt.ylabel ("Averaged drink interval (ms)")
    plt.ylim (0, 400)

    plt.figure (num=8)
    x = win4
    y = drinkhist0
    plt.plot (x, y, label="Through trial")
    plt.title ("Drink interval histogram_All(2 Step)")
    plt.xlabel ("Running interval (ms)")
    plt.ylabel ("Averaged drink interval (ms)")
    plt.ylim (0, 400)

    E = np.zeros (width6)
    F = np.zeros (width6)
    Z = np.zeros (width6)
    drinkhist0 = np.zeros (width6)
    drinkhist1 = np.zeros (width6)
    drinkhist2 = np.zeros (width6)
    turn = 1
    time = oneturn[0]
    for i in range (3, n - 30, 3):  # counting putative rpeg touch number
        s[i] = peg[i] - peg[i - 1]
        if peg[i] > time:
            time = time + oneturn[turn]
            turn = turn + 1
        for k in range (width6 - 1):
            if win6[k] < s[i] < win6[k + 1]:  # filtering interval
                for j in range (dr - 1):
                    if (peg[i - 3]) < drink[j] < (peg[i]):
                        if (mn2 < drink[j + 2] - drink[j + 1] < mx2) and (mn2 < drink[j + 1] - drink[j] < mx2) and \
                                (mn2 < drink[j] - drink[j - 1] < mx2) and (mn2 < drink[j - 1] - drink[j - 2] < mx2):
                            Z[k] = Z[k] + 1
                            drinkhist0[k] = drinkhist0[k] + (drink[j + 1] - drink[j])
                            if oneturn[turn] < oneturn[turn - 1]:
                                E[k] = E[k] + 1
                                drinkhist1[k] = drinkhist1[k] + (drink[j + 1] - drink[j])
                            else:
                                F[k] = F[k] + 1
                                drinkhist2[k] = drinkhist2[k] + (drink[j + 1] - drink[j])

    for i in range (width6 - 1):
        if not (E[i] == 0):
            drinkhist1[i] = (drinkhist1[i] / E[i])
        else:
            drinkhist1[i] = 0
    for i in range (width6 - 1):
        if not (F[i] == 0):
            drinkhist2[i] = (drinkhist2[i] / F[i])
        else:
            drinkhist2[i] = 0
    for i in range (width6 - 1):
        if not (Z[i] == 0):
            drinkhist0[i] = (drinkhist0[i] / Z[i])
        else:
            drinkhist0[i] = 0

    plt.figure (num=11)
    x = win6
    y1 = drinkhist1
    y2 = drinkhist2
    plt.plot (x, y1, label="Speed up", color="r")
    plt.plot (x, y2, label="Speed down", color="b")
    plt.title ("Drink interval histogram(3 Step)")
    plt.xlabel ("Running interval (ms)")
    plt.ylabel ("Averaged drink interval (ms)")
    plt.ylim (0, 400)

    plt.figure (num=12)
    x = win6
    y = drinkhist0
    plt.plot (x, y, label="Through trial")
    plt.title ("Drink interval histogram_All(3 Step)")
    plt.xlabel ("Running interval (ms)")
    plt.ylabel ("Averaged drink interval (ms)")
    plt.ylim (0, 400)

####ここまで水飲みのインターバル解析
########################################################################################################################

    Z = np.zeros(width2)
    lastint = np.zeros(width2)
    turn = 1
    time = oneturn[0]
    for i in range(1, n - 30):  # counting putative rpeg touch number
        s[i] = peg[i] - peg[i - 1]
        if peg[i] > time:
            time = time + oneturn[turn]
            turn = turn + 1
        for k in range(width2 - 1):
            if win2[k] < s[i] < win2[k + 1]:  # filtering interval
                for j in range(dr - 1):
                    if ((peg[i - 1]) < drink[j])&(drink[j] < peg[i])&(drink[j+1]>peg[i]):
                        if (mn2 < drink[j + 2] - drink[j + 1])&(drink[j + 2] - drink[j + 1] < mx2)&\
                                (mn2 < drink[j + 1] - drink[j])&(drink[j + 1] - drink[j] < mx2)& \
                                (mn2 < drink[j] - drink[j - 1])&(drink[j] - drink[j - 1] < mx2)&\
                                (mn2 < drink[j - 1] - drink[j - 2])&(drink[j - 1] - drink[j - 2] < mx2):
                            Z[k] = 1
                            lastint[k] = lastint[k] + (peg[i] - drink[j])

    for i in range(width2 - 1):
        if not (Z[i] == 0):
            lastint[i] = (lastint[i] / Z[i])
        else:
            lastint[i] = 0

    plt.figure(num=13)
    x = win2
    y = lastint
    plt.plot(x, y, color="r")
    plt.title("Last interval")
    plt.xlabel("Running interval (ms)")
    plt.ylabel("Last interval (ms)")
    plt.ylim(0, 200)


    Z = np.zeros(width4)
    lastint = np.zeros(width4)
    turn = 1
    time = oneturn[0]
    for i in range(1, n - 30):  # counting putative rpeg touch number
        s[i] = peg[i] - peg[i - 1]
        if peg[i] > time:
            time = time + oneturn[turn]
            turn = turn + 1
        for k in range(width4 - 1):
            if win4[k] < s[i] < win4[k + 1]:  # filtering interval
                for j in range(dr - 1):
                    if ((peg[i - 1]) < drink[j]) & (drink[j] < peg[i]) & (drink[j + 1] > peg[i]):
                        if (mn2 < drink[j + 2] - drink[j + 1]) & (drink[j + 2] - drink[j + 1] < mx2) & \
                                (mn2 < drink[j + 1] - drink[j]) & (drink[j + 1] - drink[j] < mx2) & \
                                (mn2 < drink[j] - drink[j - 1]) & (drink[j] - drink[j - 1] < mx2) & \
                                (mn2 < drink[j - 1] - drink[j - 2]) & (drink[j - 1] - drink[j - 2] < mx2):
                            Z[k] = 1
                            lastint[k] = lastint[k] + (peg[i] - drink[j])

    for i in range(width4 - 1):
        if not (Z[i] == 0):
            lastint[i] = (lastint[i] / Z[i])
        else:
            lastint[i] = 0

    plt.figure(num=14)
    x = win4
    y = lastint
    plt.plot(x, y, color="r")
    plt.title("Last interval(2 Step)")
    plt.xlabel("Running interval (ms)")
    plt.ylabel("Last interval (ms)")
    plt.ylim(0, 200)


    Z = np.zeros (width6)
    lastint = np.zeros(width6)
    turn = 1
    time = oneturn[0]
    for i in range (1, n - 30):  # counting putative rpeg touch number
        s[i] = peg[i] - peg[i - 1]
        if peg[i] > time:
            time = time + oneturn[turn]
            turn = turn + 1
        for k in range (width6 - 1):
            if win6[k] < s[i] < win6[k + 1]:  # filtering interval
                for j in range(dr - 1):
                    if ((peg[i - 1]) < drink[j]) & (drink[j] < peg[i]) & (drink[j + 1] > peg[i]):
                        if (mn2 < drink[j + 2] - drink[j + 1]) & (drink[j + 2] - drink[j + 1] < mx2) & \
                                (mn2 < drink[j + 1] - drink[j]) & (drink[j + 1] - drink[j] < mx2) & \
                                (mn2 < drink[j] - drink[j - 1]) & (drink[j] - drink[j - 1] < mx2) & \
                                (mn2 < drink[j - 1] - drink[j - 2]) & (drink[j - 1] - drink[j - 2] < mx2):
                            Z[k] = 1
                            lastint[k] = lastint[k] + (peg[i] - drink[j])

    for i in range (width6 - 1):
        if not (Z[i] == 0):
            lastint[i] = (lastint[i] / Z[i])
        else:
            lastint[i] = 0

    plt.figure (num=15)
    x = win6
    y = lastint
    plt.plot(x, y, color="r")
    plt.title("Last interval(3 Step)")
    plt.xlabel("Running interval (ms)")
    plt.ylabel("Last interval (ms)")
    plt.ylim(0, 500)

####ここまでペグの手前最後の水飲みから次のペグまでのインターバルのヒストグラム
########################################################################################################################

    plt.show()

###### 5s Inter Peg Interval = about400ms with avarage 2 drink on.  So,10s IPI = about800ms with 4 drink on
###### therefore, need 5 samples each of 10s 8s 6s 5s 4s trials, for prediction drink on nember 4,3,2,2,2 of each trial.
###### マウス1匹づつ階段グラフを出したもの。