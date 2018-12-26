import scipy.io as sci
import os
import matplotlib.pyplot as plt
import numpy as np


f = [213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 301, 302, 303, 304, 305, 306, 307,
     308, 309, 310]
ID = [1101, 1201, 1301, 1701, 1801, 1901, 2001, 2101, 2201, 2301]
start = 20
finish = 21
high = 500
low = 0
day = finish - start

for number in range(10):                                 #put mouse number
    mouse = ID[number]
    print(mouse)
    PSave = list(range(24))   #####crate prg touch array of each mouse.
    for i in range(24):
        PSave[i] = 0

    for num in range(start, finish):                                #put file nuber
        PATH = 'C:\\Users\\hirokane\\PycharmProjects\\Wheel_Analyze\\LearnCurve_201701\\learncurve'
        os.chdir(PATH)
        matdata = sci.loadmat(str(mouse) + "_170" + str(f[num]))    #create file name
        RPegnumber = matdata["RpegTimeArray2D"]
        LPegnumber = matdata["LpegTimeArray2D"]
        Pegnumber_sort = []

        for i in range(len(RPegnumber)-1):
            Pegnumber = []
            C = 0
            for j in range(12):
                if (RPegnumber[i, j] == RPegnumber[i, j])and(LPegnumber[i, j] == LPegnumber[i, j]):
                    C = C + 1
            if C == 12:
                Pegnumber.extend(RPegnumber[i])
                Pegnumber.extend(LPegnumber[i])
                Pegnumber.sort()
                Pegnumber_sort.append(Pegnumber)
                #####align only clear peg touch
        TurnNumber = len(Pegnumber_sort)
        Size = 24

        Ave = list(range(24))
        A = list(range(24))
        for i in range(24):
            Ave[i] = 0
            A[i] = 0

        for i in range(1, TurnNumber-1):
            for j in range(3, Size-1):
                if (low < (Pegnumber_sort[i][j+1] - Pegnumber_sort[i][j]) < high)and(low < (Pegnumber_sort[i][j] - Pegnumber_sort[i][j-1]) < high)and(low < (Pegnumber_sort[i][j-1] - Pegnumber_sort[i][j-2]) < high)and(low < (Pegnumber_sort[i][j-2] - Pegnumber_sort[i][j-3]) < high):
                    Phasepost = (Pegnumber_sort[i][j]-Pegnumber_sort[i][j-1])/(Pegnumber_sort[i][j+1]-Pegnumber_sort[i][j-1])
                    Phasepre = (Pegnumber_sort[i][j-2]-Pegnumber_sort[i][j-3])/(Pegnumber_sort[i][j-1]-Pegnumber_sort[i][j-3])
                    Phaseshift = Phasepost - Phasepre
                    Ave[j] = Ave[j] + Phaseshift
                    A[j] = A[j] + 1
                    ### analize only good walking zone ex.)before and after 2steps are natural walking interval.

            ##j == 0
            if (low < (Pegnumber_sort[i][1] - Pegnumber_sort[i][0]) < high)and(low < (Pegnumber_sort[i][0] - Pegnumber_sort[i-1][23]))and(low < (Pegnumber_sort[i-1][23] - Pegnumber_sort[i-1][22]))and(low < (Pegnumber_sort[i-1][22] - Pegnumber_sort[i-1][21])):
                Phasepost = (Pegnumber_sort[i][0] - Pegnumber_sort[i-1][23]) / (Pegnumber_sort[i][1] - Pegnumber_sort[i-1][23])
                Phasepre = (Pegnumber_sort[i-1][22] - Pegnumber_sort[i-1][21]) / (Pegnumber_sort[i-1][23] - Pegnumber_sort[i-1][21])
                Phaseshift = Phasepost - Phasepre
                Ave[0] = Ave[0] + Phaseshift
                A[0] = A[0] + 1
            ##j == 1
            if (low < (Pegnumber_sort[i][2] - Pegnumber_sort[i][1]) < high)and(low < (Pegnumber_sort[i][1] - Pegnumber_sort[i][0]))and(low < (Pegnumber_sort[i][0] - Pegnumber_sort[i-1][23]))and(low < (Pegnumber_sort[i-1][23] - Pegnumber_sort[i-1][22])):
                Phasepost = (Pegnumber_sort[i][1] - Pegnumber_sort[i][0]) / (Pegnumber_sort[i][2] - Pegnumber_sort[i][0])
                Phasepre = (Pegnumber_sort[i-1][23] - Pegnumber_sort[i-1][22]) / (Pegnumber_sort[i][0] - Pegnumber_sort[i-1][22])
                Phaseshift = Phasepost - Phasepre
                Ave[1] = Ave[1] + Phaseshift
                A[1] = A[1] + 1
            ##j == 2
            if (low < (Pegnumber_sort[i][3] - Pegnumber_sort[i][2]) < high) and (low < (Pegnumber_sort[i][2] - Pegnumber_sort[i][1])) and (low < (Pegnumber_sort[i][0] - Pegnumber_sort[i - 1][23])) and (low < (Pegnumber_sort[i - 1][23] - Pegnumber_sort[i - 1][22])):
                Phasepost = (Pegnumber_sort[i][2] - Pegnumber_sort[i][1]) / (Pegnumber_sort[i][3] - Pegnumber_sort[i][1])
                Phasepre = (Pegnumber_sort[i][0]-Pegnumber_sort[i-1][23])/(Pegnumber_sort[i][1] - Pegnumber_sort[i-1][23])
                Phaseshift = Phasepost - Phasepre
                Ave[2] = Ave[2] + Phaseshift
                A[2] = A[2] + 1
            ##j == 23
            if (low < (Pegnumber_sort[i+1][0] - Pegnumber_sort[i][23]) < high) and (low < (Pegnumber_sort[i][23] - Pegnumber_sort[i][22])) and (low < (Pegnumber_sort[i][22] - Pegnumber_sort[i][21])) and (low < (Pegnumber_sort[i][21] - Pegnumber_sort[i][20])):
                Phasepost = (Pegnumber_sort[i][23] - Pegnumber_sort[i][22]) / (Pegnumber_sort[i+1][0] - Pegnumber_sort[i][22])
                Phasepre = (Pegnumber_sort[i][21] - Pegnumber_sort[i][20]) / (Pegnumber_sort[i][22] - Pegnumber_sort[i][20])
                Phaseshift = Phasepost - Phasepre
                Ave[23] = Ave[23] + Phaseshift
                A[23] = A[23] + 1
            print(Ave)

        for i in range(Size):
            if not A[i] == 0:
                Ave[i] = Ave[i] / A[i]

        for i in range(24):
            PSave[i] = PSave[i] + abs(Ave[i] / day)/2

    x = np.arange(1, 25)
    y = PSave
    plt.bar(x, y)
    plt.title("Phaseshift of Pegtouch")
    plt.xlabel("PegNumber")
    plt.ylabel("Phase shift")
    os.chdir('C:\\Users\\hirokane\\PycharmProjects\\Wheel_Analyze\\Phase_shift_figure')
    plt.savefig(str(mouse) + ".png")
    plt.show()