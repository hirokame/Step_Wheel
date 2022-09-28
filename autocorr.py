import scipy.io as sci
import os
import seaborn as sns
import scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import heapq
import numpy.ma as ma
from scipy.signal import savgol_filter
from scipy.fftpack import fft
from scipy import stats
import scipy

np.set_printoptions (threshold=np.inf, linewidth=1000, suppress=True)

def tukey_hsd( ind, *args):
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    data_arr = np.hstack(args)

    ind_arr = np.array([])
    for x in range(len(args)):
        ind_arr = np.append(ind_arr, np.repeat(ind[x], len(args[x])))
    print(pairwise_tukeyhsd(data_arr, ind_arr))

f2 = [1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1211, 1212, 1213] ##day 0~8, 9~11
f3 = [1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1210] ##day0~8
ID2 = [1, 2, 4, 5, 7, 8, 9] ##mouse 0~6
ID3 = [2, 5, 7, 9]  ##mouse 0~4 #8

start = 0
finish = 9
day = finish - start
Inthigh =800
high =500
low = 0
cut = 600
thr = 15000
All_Variance_ave_learn = np.zeros (day)
All_Variance_peak_learn = np.zeros (day)  # 上から3つ
All_Variance_base_learn = np.zeros (day)  # 下から3つ
All_Variance_peak2_learn = np.zeros (day)  # １SD超えたもの
High_Variance_Phase_first = np.zeros(15)
High_Variance_Int_first = np.zeros(15)
High_Variance_Phase_last = np.zeros(15)
High_Variance_Int_last = np.zeros(15)
drinkst = [[] for i in range(9)]
drinkend = [[] for i in range(9)]
drinkstpower = [[] for i in range(9)]
drinkendpower =[[] for i in range(9)]
mn1 = 20
mx1 = 200
mn2 = 40
mx2 = 200
N = 600
dt = 0.01
Autocorr = np.zeros(200)
drinkhist = np.zeros(200)
for number in range (2,3):  # put mouse number
    mouse = ID2[number]
    #mouse = ID3[number]

    Variance = np.zeros ((24, 24))
    RVariance = np.zeros ((12, 12))
    LVariance = np.zeros ((12, 12))

    Variance_mean = np.zeros (24)
    Var_ave_learn = np.zeros (day)
    Var_peak_learn = np.zeros (day)
    Var_base_learn = np.zeros (day)
    Var_peak2_learn = np.zeros (day)

    for num in range (start, finish):  # put file nuber
        print (mouse)
        corrwin = np.arange (-300, 300, 2)
        corr = np.zeros (len (corrwin))
        PATH = 'C:\\Users\\hirokane\\PycharmProjects\\Wheel_Analyze\\mat2\\mat2learn'
        os.chdir (PATH)
        matdata = sci.loadmat (str (mouse) + "08" + str (f2[num]))
        # matdata = sci.loadmat (str (mouse) + "08" + str (f3[num]))
        # print(matdata.keys())
        Data = matdata["data"]
        peg = matdata["RLPegTouchAll"]
        drink = np.array(matdata["DrinkOnArray"])
        for j in range(5):
            delete = []
            eachdrinkdiff = np.diff(drink)
            for i in range(len(eachdrinkdiff)):
                if eachdrinkdiff[i] < 50:
                    delete.append(i)
            drink = np.delete(drink, delete)
        RPegTimeArray_turn = matdata["RpegTimeArray2D_turn"]
        LPegTimeArray_turn = matdata["LpegTimeArray2D_turn"]
        RpegTouchCell = matdata["RpegTouchCell"]
        LpegTouchCell = matdata["LpegTouchCell"]
        RpegTouchTurn = matdata["RpegTouchTurn"]
        LpegTouchTurn = matdata["LpegTouchTurn"]
        wateroffarray = matdata["WaterOffArray"]
        wateronarray = matdata["WaterOnArray"]
        oneturn = Data[:, 4][np.nonzero (Data[:, 4])]

        drinkdiff = np.diff(drink, axis=0)
        print(drinkdiff)
        for i in range (len (drinkdiff)):
            if drinkdiff[i] <= 200:
                for j in range (200):
                    if j  <= drinkdiff[i] < j + 1:
                        drinkhist[j] = drinkhist[j] + 1
                        continue
        for i in range (len (drink)):
            for k in [j for j, x in enumerate (drink) if drink[i] - 500 < x < drink[i] + 500]:
                for l in range (200):
                    if l * 5 - 500 <= drink[k] - drink[i] < (l + 1) * 5 - 500:
                        Autocorr[l] = Autocorr[l] + 1
                    continue

plt.figure(num=33)
Autocorr[100] = 0
for i in range(190):
    Autocorr[i+5] = np.average(Autocorr[i:i+10])
plt.title("Auto correlogram")
plt.xlabel("Lag")
plt.ylabel("Auto correlation")
x = np.arange(-500, 500, 5)
y = Autocorr
plt.plot(x, y)

plt.figure(num=34)
for i in range(190):
    drinkhist[i+5] = np.average(drinkhist[i:i+10])
plt.title("Drink histogram")
plt.xlabel("interval")
plt.ylabel("count")
x = np.arange(0, 200)
y = drinkhist
plt.plot(x, y)

plt.show()