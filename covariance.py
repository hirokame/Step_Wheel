import scipy.io as sci
import matplotlib.pyplot as plt
import numpy as np
import os

f = [213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 301, 302, 303, 304, 305, 306, 307,
     308, 309, 310]
ID = [1101, 1201, 1301, 1701, 1801, 1901, 2001, 2101, 2201, 2301]
AllRCov = list(range(12))
AllLCov = list(range(12))

for number in range(10):                                 #put mouse number
    mouse = ID[number]
    print(mouse)
    RCovariance = list(range(12))
    LCovariance = list(range(12))
    for i in range(12):
        RCovariance[i] = 0
    for i in range(12):
        LCovariance[i] = 0

    for num in range(4, 14):                                #put file nuber
        PATH = 'C:\\Users\\hirokane\\PycharmProjects\\Wheel_Analyze\\LearnCurve_201701\\learncurve'
        os.chdir(PATH)
        matdata = sci.loadmat(str(mouse) + "_170" + str(f[num]))    #create file name
        RPegnumber = matdata["RpegTimeArray2D"]
        LPegnumber = matdata["LpegTimeArray2D"]       ##LPegnumber(a,b)でa行b列を抜き出せる

        ##print(RPegnumber)
        ##print(LPegnumber)
        Pegnumber = [[]]

        for i in range(len(RPegnumber)-1):
            C = 0
            for j in range(len(RPegnumber[1])):
                if 100 < RPegnumber[i,j] < 500000:
                    C = C + 1
            for j in range(len(LPegnumber[1])):
                if  100 < LPegnumber[i,j] <500000:
                    C = C + 1
            if C == 24:
                Pegnumber.extend(RPegnumber[i])
                Pegnumber.extend(LPegnumber[i])
                Pegnumber[i].sort()
                #####一周すべて足取りが取れている時だけ抜き出してタッチ時間順に並べ直す

        if len(RPegnumber) > len(LPegnumber):
            TurnNumber = len(RPegnumber)
        else:
            TurnNumber = len(LPegnumber)
        RSize = len(RPegnumber[1])
        LSize = len(LPegnumber[1])

        RAve = list(range(RSize))
        LAve = list(range(LSize))
        for i in range(RSize):
            RAve[i] = 0
        for i in range(LSize):
            LAve[i] = 0

        RCov = list(range(RSize))
        LCov = list(range(LSize))
        for i in range(RSize):
            RCov[i] = 0
        for i in range(LSize):
            LCov[i] = 0

        A = list(range(RSize))
        B = list(range(LSize))
        for i in range(RSize):
            A[i] = 0
        for i in range(LSize):
            B[i] = 0

        for i in range(TurnNumber-2):
            for j in range(RSize):
                if 0 < RPegnumber[i, j] < 500000:  ####whether pegtouch exists or not
                    if not ((j == 0) and (i == 0)):
                        if j == 0:   ####if pegnumber is 1
                            if (0 < RPegnumber[i-1, RSize-1] < 500000) and (0 < RPegnumber[i-1, RSize-2] < 500000) and (0 < RPegnumber[i, j+1] < 500000) and (0 < RPegnumber[i, j+2] < 500000):
                                if 100 < (RPegnumber[i, j] - RPegnumber[i-1, RSize-1]) < 600:
                                    RAve[j] = RAve[j] + (RPegnumber[i, j] - RPegnumber[i-1, RSize-1])
                                    A[j] = A[j] + 1     ####RAve is cumulative pegtouch interval
                        elif j == 1:
                            if (0 < RPegnumber[i, j-1] < 500000) and (0 < RPegnumber[i-1, RSize-1] < 500000) and (0 < RPegnumber[i, j+1] < 500000) and (0 < RPegnumber[i, j+2] < 500000):
                                if 100 < (RPegnumber[i, j] - RPegnumber[i, j-1]) < 600:
                                    RAve[j] = RAve[j] + (RPegnumber[i, j] - RPegnumber[i, j-1])
                                    A[j] = A[j] + 1
                        elif j == RSize-1:
                            if (0 < RPegnumber[i, j-1] < 500000) and (0 < RPegnumber[i, j-2] < 500000) and (0 < RPegnumber[i+1, 1] < 500000) and (0 < RPegnumber[i+1, 2] < 500000):
                                if 100 < (RPegnumber[i, j] - RPegnumber[i, j-1]) < 600:
                                    RAve[j] = RAve[j] + (RPegnumber[i, j] - RPegnumber[i, j-1])
                                    A[j] = A[j] + 1
                        elif j == RSize-2:
                            if (0 < RPegnumber[i, j-1] < 500000) and (0 < RPegnumber[i, j-2] < 500000) and (0 < RPegnumber[i, j+1] < 500000) and (0 < RPegnumber[i+1, 1] < 500000):
                                if 100 < (RPegnumber[i, j] - RPegnumber[i, j-1]) < 600:
                                    RAve[j] = RAve[j] + (RPegnumber[i, j] - RPegnumber[i, j-1])
                                    A[j] = A[j] + 1
                        else:
                            if (0 < RPegnumber[i, j-1] < 500000) and (0 < RPegnumber[i, j-2] < 500000) and (0 < RPegnumber[i, j+1] < 500000) and (0 < RPegnumber[i, j+2] < 500000):
                                if 100 < (RPegnumber[i, j] - RPegnumber[i, j-1]) < 600:
                                    RAve[j] = RAve[j] + (RPegnumber[i, j] - RPegnumber[i, j-1])
                                    A[j] = A[j] + 1
            for j in range(LSize):
                if 0 < LPegnumber[i, j] < 500000:
                    if not((j == 0) and (i == 0)):
                        if j == 0:
                            if (0 < LPegnumber[i-1, LSize-1] < 500000) and (0 < LPegnumber[i-1, LSize-2] < 500000) and (0 < LPegnumber[i, j+1] < 500000) and (0 < LPegnumber[i, j+2] < 500000):
                                if 100 < (LPegnumber[i, j] - LPegnumber[i - 1, LSize - 1]) < 600:
                                    LAve[j] = LAve[j] + (LPegnumber[i, j] - LPegnumber[i-1, LSize-1])
                                    B[j] = B[j] + 1
                        elif j == 1:
                            if (0 < LPegnumber[i, j-1] < 500000) and (0 < LPegnumber[i-1, LSize-1] < 500000) and (0 < LPegnumber[i, j+1] < 500000) and (0 < LPegnumber[i, j+2] < 500000):
                                if 100 < (LPegnumber[i, j] - LPegnumber[i, j - 1]) < 600:
                                    LAve[j] = LAve[j] + (LPegnumber[i, j] - LPegnumber[i, j-1])
                                    B[j] = B[j] + 1
                        elif j == LSize-1:
                            if (0 < LPegnumber[i, j-1] < 500000) and (0 < LPegnumber[i, j-2] < 500000) and (0 < LPegnumber[i+1, 1] < 500000) and (0 < LPegnumber[i+1, 2] < 500000):
                                if 100 < (LPegnumber[i, j] - LPegnumber[i, j - 1]) < 600:
                                    LAve[j] = LAve[j] + (LPegnumber[i, j] - LPegnumber[i, j-1])
                                    B[j] = B[j] + 1
                        elif j == LSize-2:
                            if (0 < LPegnumber[i, j-1] < 500000) and (0 < LPegnumber[i, j-2] < 500000) and (0 < LPegnumber[i, j+1] < 500000) and (0 < LPegnumber[i+1, 1] < 500000):
                                if 100 < (LPegnumber[i, j] - LPegnumber[i, j - 1]) < 600:
                                    LAve[j] = LAve[j] + (LPegnumber[i, j] - LPegnumber[i, j-1])
                                    B[j] = B[j] + 1
                        else:
                            if (0 < LPegnumber[i, j-1] < 500000) and (0 < LPegnumber[i, j-2] < 500000) and (0 < LPegnumber[i, j+1] < 500000) and (0 < LPegnumber[i, j+2] < 500000):
                                if 100 < (LPegnumber[i, j] - LPegnumber[i, j - 1]) < 600:
                                    LAve[j] = LAve[j] + (LPegnumber[i, j] - LPegnumber[i, j-1])
                                    B[j] = B[j] + 1
        for i in range(RSize):
            RAve[i] = RAve[i] / A[i]

        for i in range(LSize):
            LAve[i] = LAve[i] / B[i]

        #####RAve: average peg touch interval for each peg
        for i in range(RSize):
            A[i] = 0
        for i in range(LSize):
            B[i] = 0

        for i in range(TurnNumber-2):
            for j in range(RSize):
                if not((i == 0) and (j == 0)):
                    if j == (RSize-1):
                        if (0 < RPegnumber[i, j] < 500000) and (0 < RPegnumber[i, j-1] < 500000) and (0 < RPegnumber[i+1, 1] < 500000) and (0 < RPegnumber[i, j-2] < 500000) and (0 < RPegnumber[i+1, 2] < 500000):
                            if 100 < (RPegnumber[i, j] - RPegnumber[i, j - 1]) < 600:
                                RCov[j] = RCov[j]+((RPegnumber[i, j]-RPegnumber[i, j-1])-RAve[j])*((RPegnumber[i+1, 1]-RPegnumber[i, j])-RAve[1])
                                A[j] = A[j] + 1
                    elif j == (RSize-2):
                        if (0 < RPegnumber[i, j] < 500000) and (0 < RPegnumber[i, j-1] < 500000) and (0 < RPegnumber[i+1, 1] < 500000) and (0 < RPegnumber[i, j-2] < 500000) and (0 < RPegnumber[i, j+1] < 500000):
                            if 100 < (RPegnumber[i, j] - RPegnumber[i, j - 1]) < 600:
                                RCov[j] = RCov[j]+((RPegnumber[i, j]-RPegnumber[i, j-1])-RAve[j])*((RPegnumber[i, j+1]-RPegnumber[i, j])-RAve[j+1])
                                A[j] = A[j] + 1
                    elif j == 1:
                        if (0 < RPegnumber[i, j] < 500000) and (0 < RPegnumber[i, j+2] < 500000) and (0 < RPegnumber[i, j+1] < 500000) and (0 < RPegnumber[i, j-1] < 500000) and (0 < RPegnumber[i-1,RSize-1] < 500000):
                            if 100 < (RPegnumber[i, j] - RPegnumber[i, j - 1]) < 600:
                                RCov[j] = RCov[j]+((RPegnumber[i, j]-RPegnumber[i, j-1])-RAve[j])*((RPegnumber[i, j+1]-RPegnumber[i, j])-RAve[j+1])
                                A[j] = A[j] + 1
                    elif j == 0:
                        if (0 < RPegnumber[i, j] < 500000) and (0 < RPegnumber[i, j+1] < 500000) and (0 < RPegnumber[i-1, RSize-1] < 500000) and (0 < RPegnumber[i-1, RSize-2] < 500000) and (0 < RPegnumber[i, j+2] < 500000):
                            if 100 < (RPegnumber[i, j] - RPegnumber[i - 1, RSize - 1]) < 600:
                                RCov[j] = RCov[j]+((RPegnumber[i, j+1] - RPegnumber[i, j])-RAve[j+1])*((RPegnumber[i, j]-RPegnumber[i-1, RSize-1]) - RAve[1])
                                A[j] = A[j] + 1
                    else:
                        if (0 < RPegnumber[i, j] < 500000) and (0 < RPegnumber[i, j-1] < 500000) and (0 < RPegnumber[i, j+1] < 500000) and (0 < RPegnumber[i, j+2] < 500000) and (0 < RPegnumber[i, j-2] < 500000):
                            if 100 < (RPegnumber[i, j] - RPegnumber[i, j - 1]) < 600:
                                RCov[j] = RCov[j]+((RPegnumber[i, j]-RPegnumber[i, j-1])-RAve[j])*((RPegnumber[i, j+1]-RPegnumber[i, j])-RAve[j+1])
                                A[j] = A[j] + 1

        for i in range(TurnNumber-2):
            for j in range(LSize):
                if not((i == 0) and(j == 0)):
                    if j == (LSize - 1):
                        if (0 < LPegnumber[i, j] < 500000) and (0 < LPegnumber[i, j-1] < 500000) and (0 < LPegnumber[i+1, 1] < 500000) and (0 < LPegnumber[i+1, 2] < 500000) and (0 < LPegnumber[i, j-2] < 500000):
                            if 100 < (LPegnumber[i, j] - LPegnumber[i, j - 1]) < 600:
                                LCov[j] = LCov[j]+((LPegnumber[i, j]-LPegnumber[i, j-1])-LAve[j])*((LPegnumber[i+1, 1]-LPegnumber[i, j])-LAve[1])
                                B[j] = B[j] + 1
                    elif j == (LSize-2):
                        if (0 < LPegnumber[i, j] < 500000) and (0 < LPegnumber[i, j-1] < 500000) and (0 < LPegnumber[i+1, 1] < 500000) and (0 < LPegnumber[i, j-2] < 500000) and (0 < LPegnumber[i, j+1] < 500000):
                            if 100 < (LPegnumber[i, j] - LPegnumber[i, j - 1]) < 600:
                                LCov[j] = LCov[j]+((LPegnumber[i, j]-LPegnumber[i, j-1])-LAve[j])*((LPegnumber[i, j+1]-LPegnumber[i, j])-LAve[j+1])
                                B[j] = B[j] + 1
                    elif j == 1:
                        if (0 < LPegnumber[i, j] < 500000) and (0 < LPegnumber[i, j+2] < 500000) and (0 < LPegnumber[i, j+1] < 500000) and (0 < LPegnumber[i, j-1] < 500000) and (0 < LPegnumber[i-1,LSize-1] < 500000):
                            if 100 < (LPegnumber[i, j] - LPegnumber[i, j - 1]) < 600:
                                LCov[j] = LCov[j]+((LPegnumber[i, j]-LPegnumber[i, j-1])-LAve[j])*((LPegnumber[i, j+1]-LPegnumber[i, j])-LAve[j+1])
                                B[j] = B[j] + 1
                    elif j == 0:
                        if (0 < LPegnumber[i, j] < 500000) and (0 < LPegnumber[i, j+1] < 500000) and (0 < LPegnumber[i-1, LSize-1] < 500000) and (0 < LPegnumber[i-1, LSize-2] < 500000) and (0 < LPegnumber[i, j+2] < 500000):
                            if 100 < (LPegnumber[i, j] - LPegnumber[i - 1, LSize - 1]) < 600:
                                LCov[j] = LCov[j]+((LPegnumber[i, j+1] - LPegnumber[i, j])-LAve[j+1])*((LPegnumber[i, j]-LPegnumber[i-1, LSize-1]) - LAve[1])
                                B[j] = B[j] + 1
                    else:
                        if (0 < LPegnumber[i, j] < 500000) and (0 < LPegnumber[i, j-1] < 500000) and (0 < LPegnumber[i, j+1] < 500000) and (0 < LPegnumber[i, j+2] < 500000) and (0 < LPegnumber[i, j-2] < 500000):
                            if 100 < (LPegnumber[i, j] - LPegnumber[i, j - 1]) < 600:
                                LCov[j] = LCov[j]+((LPegnumber[i, j]-LPegnumber[i, j-1])-LAve[j])*((LPegnumber[i, j+1]-LPegnumber[i, j])-LAve[j+1])
                                B[j] = B[j] + 1
        for i in range(RSize):
            RCov[i] = RCov[i] / A[i]

        for i in range(LSize):
            LCov[i] = LCov[i] / B[i]

        for i in range(RSize):
            RCovariance[i] = RCovariance[i] + RCov[i]

        for i in range(LSize):
            LCovariance[i] = LCovariance[i] + LCov[i]

    for i in range(RSize):
        RCovariance[i] = RCovariance[i] / 10

    for i in range(LSize):
        LCovariance[i] = LCovariance[i] / 10

    #print(RCovariance)
    #print(LCovariance)
    ##同じグラフしか出力されない
    CoIntmouse = []
    for i in range(12):
        CoIntmouse.append(RCovariance[i])
        CoIntmouse.append(LCovariance[i])

    plt.figure()
    x = np.arange(24)
    y = CoIntmouse
    plt.bar(x, y, tick_label=x+1)
    plt.title("Peg_Covariance_eachmouse")
    plt.xlabel("Peg_number")
    plt.ylabel("Covariance")
    plt.ylim(-500, 500)
    plt.show()

    for i in range(RSize):
        AllRCov[i] = AllRCov[i] + RCovariance[i]/10

    for i in range(LSize):
        AllLCov[i] = AllLCov[i] + LCovariance[i]/10

CoInt = []
for i in range(12):
    CoInt.append(AllRCov[i])
    CoInt.append(AllLCov[i])

plt.figure()
x = np.arange(24)
y = CoInt
plt.bar(x, y, tick_label=x+1)
plt.xlabel("Peg_number")
plt.ylabel("Covariance")
plt.title("Peg_Covariance_All")
plt.ylim(-200, 200)
plt.show()


##ペグ番号(1~12)のそれを挟むペグインターバルのcovaarianceを計算
##マウス10匹それぞれに対して26日分のデータ
##10匹の平均のfigureを作成
##covarianceの計算を右足だけ、左足だけでやってしまった
##前後2歩(合計5歩)連続で歩いている時のみ解析