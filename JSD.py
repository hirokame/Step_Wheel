import scipy.io as sci
from scipy import stats
import scipy
import os
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.stats import entropy
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
from sklearn import linear_model

def tukey_hsd( ind, *args):
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    data_arr = np.hstack(args)

    ind_arr = np.array([])
    for x in range(len(args)):
        ind_arr = np.append(ind_arr, np.repeat(ind[x], len(args[x])))
    print(pairwise_tukeyhsd(data_arr, ind_arr))

f = [213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 301, 302, 303, 304, 305, 306, 307,
     308, 309, 310]  ## day 0~3, 4~13, 14~20, 21~25
f2 = [1127, 1128, 1129, 1130, 1202, 1203, 1204, 1205, 1206, 1207]  #1201==6s
ID = [1101, 1201, 1301, 1701, 1801, 1901, 2101, 2201, 2301]  # n=9
ID2 = [198, 199, 200, 201, 204, 205, 207] # n=7

start = 1
finish = 10
day = finish - start
high = 800
low = 0
Inthigh = 700
cut = 250
thr = 20000
n = 7    # sample number
size = day*n

KL1 = np.zeros([n*day,n*day])
KL2 = np.zeros([n*day,n*day])
JSD = np.zeros([n*day,n*day])
Tensor = np.zeros([n, day, 200])
KLShift1 = np.zeros([n * day, n * day])
KLInt1 = np.zeros([n * day, n * day])
KLVar1 = np.zeros([n * day, n * day])
KLCov1 = np.zeros([n * day, n * day])
KLVarFlat1 = np.zeros([n * day, n * day])
KLShift2 = np.zeros([n * day, n * day])
KLInt2 = np.zeros([n * day, n * day])
KLVar2 = np.zeros([n * day, n * day])
KLCov2 = np.zeros([n * day, n * day])
KLVarFlat2 = np.zeros([n * day, n * day])

JSDShift = np.zeros([n * day, n * day])
JSDInt = np.zeros([n * day, n * day])
JSDVar = np.zeros([n * day, n * day])
JSDCov = np.zeros([n * day, n * day])
JSDVarFlat = np.zeros([n * day, n * day])
PhaseTensor = np.zeros([n, day, 24])
IntTensor = np.zeros([n, day, 24])
VarTensor = np.zeros([n, day, 24])
CovTensor = np.zeros([n, day, 24])
VarTensorFlat = np.zeros([n, day, 24*24])

for number in range(n):                                 #put mouse number
    mouse = ID2[number]
    print(mouse)

    for num in range (start, finish):  # put file nuber
        print(num)
        #PATH = 'C:\\Users\\hirokane\\PycharmProjects\\Wheel_Analyze\\LearnCurve_201701\\learncurve'
        PATH = 'C:\\Users\\hirokane\\PycharmProjects\\Wheel_Analyze\\analyze_new\\behav_data\\All_in_one'
        os.chdir (PATH)
        #matdata = sci.loadmat (str (mouse) + "_170" + str (f[num]))  # create file name
        matdata = sci.loadmat (str (mouse) + "01_19" + str (f2[num]))  # new data C9 10days
        # print(matdata.keys())
        Data = matdata["data"]
        peg = matdata["RLPegTouchAll"]
        RPegTimeArray_turn = matdata["RpegTimeArray2D_turn"]
        LPegTimeArray_turn = matdata["LpegTimeArray2D_turn"]
        RpegTouchCell = matdata["RpegTouchCell"]
        LpegTouchCell = matdata["LpegTouchCell"]
        wateroffarray = matdata["WaterOffArray"]
        wateronarray = matdata["WaterOnArray"]
        ########################################################################################################################
        stop = Data[:, 2][0]
        WaterOnArray = np.zeros (stop)
        for i in range (len (wateronarray)):
            for j in range (int (wateronarray[i]), int (wateroffarray[i])):
                WaterOnArray[j] = 1
        WaterOnArray[0] = 0

        #### WaterOnArray にwateronなら1をwateroffなら0をあてがった配列
        ########################################################################################################################

        oneturn = Data[:, 4][np.nonzero (Data[:, 4])]
        Rtouchall = Data[:, 13]
        turnnum = len (oneturn)
        RPegTouchTimeArray2D = np.zeros ([turnnum, 12])
        RPegTouchTimeArray2D_turn = np.zeros ([turnnum, 12])
        TurnmarkerTime = 0
        for i in range (turnnum - 1):
            for j in range (12):
                R = Data[:, j + 13][np.nonzero (Data[:, j + 13])]
                TurnPegAll = peg[TurnmarkerTime < peg]
                TurnPegAll = TurnPegAll[TurnPegAll < TurnmarkerTime + oneturn[i]]
                if np.any ((TurnmarkerTime < R) & (R < TurnmarkerTime + oneturn[i])) == False:
                    RPegTouchTimeArray2D[i, j] = 0
                    RPegTouchTimeArray2D_turn[i, j] = 0
                    continue
                else:
                    Rnow = R[(TurnmarkerTime < R) & (R < TurnmarkerTime + oneturn[i])]
                    for now in Rnow:
                        for touch in TurnPegAll:
                            if now == touch:
                                RPegTouchTimeArray2D[i, j] = now
                                RPegTouchTimeArray2D_turn[i, j] = now - TurnmarkerTime
                                break
                        else:
                            continue
                        break

            TurnmarkerTime = TurnmarkerTime + oneturn[i]

        #### RPegTouchTimeArray2D という二次元配列に右足の各ペグのタッチタイミングを1周ごとに入れた
        #########################################################################################################################
        TurnmarkerTime = 0
        Ltouchall = Data[:, 41]
        LPegTouchTimeArray2D = np.zeros ([turnnum, 12])
        LPegTouchTimeArray2D_turn = np.zeros ([turnnum, 12])
        for i in range (turnnum - 1):
            for j in range (12):
                L = Data[:, j + 41][np.nonzero (Data[:, j + 41])]
                TurnPegAll = peg[TurnmarkerTime < peg]
                TurnPegAll = TurnPegAll[TurnPegAll < TurnmarkerTime + oneturn[i]]
                if np.any ((TurnmarkerTime < L) & (L < TurnmarkerTime + oneturn[i])) == False:
                    LPegTouchTimeArray2D[i, j] = 0
                    LPegTouchTimeArray2D_turn[i, j] = 0
                    continue
                else:
                    Lnow = L[(TurnmarkerTime < L) & (L < TurnmarkerTime + oneturn[i])]
                    for now in Lnow:
                        for touch in TurnPegAll:
                            if np.any (now == touch):
                                LPegTouchTimeArray2D[i, j] = now
                                LPegTouchTimeArray2D_turn[i, j] = now - TurnmarkerTime
                                break
                        else:
                            continue
                        break

            TurnmarkerTime = TurnmarkerTime + oneturn[i]
        # 左足も同じ操作をする

        TouchHist_turn = np.zeros ([turnnum, 200])
        for i in range (turnnum):
            for j in range (12):
                for k in range (100):
                    if (k * 41 < RPegTouchTimeArray2D_turn[i, j]) and (RPegTouchTimeArray2D_turn[i, j] < (k + 1) * 41):
                        TouchHist_turn[i, k] = 1
                        break

        for i in range (turnnum):
            for j in range (12):
                for k in range (100):
                    if (k * 41 < LPegTouchTimeArray2D_turn[i, j]) and (LPegTouchTimeArray2D_turn[i, j] < (k + 1) * 41):
                        TouchHist_turn[i, k + 100] = 1
                        break

        Touchvector = np.zeros (200)
        for i in range (200):
            for j in range (turnnum):
                Touchvector[i] = Touchvector[i] + TouchHist_turn[j, i]
        Touchvector = (Touchvector * 100) / np.sum (Touchvector)
        ##正規化
        for i in range (200):
            Tensor[number, num - start, i] = Touchvector[i]

        PegTouch2DArray = np.zeros ([turnnum, 24])
        PegTouch2DArray_sort = np.zeros ([turnnum, 24])  # タッチ時間でsortした方
        Peg2DArray = np.zeros ([turnnum, 24])
        PegTouch2DArray_turn = np.zeros ([turnnum, 24])

        for i in range (turnnum - 1):
            for j in range (12):
                PegTouch2DArray[i, 2 * j + 1] = RPegTouchTimeArray2D[i, j]
                PegTouch2DArray[i, 2 * j] = LPegTouchTimeArray2D[i, j]
                PegTouch2DArray_turn[i, 2 * j + 1] = RPegTouchTimeArray2D_turn[i, j]
                PegTouch2DArray_turn[i, 2 * j] = LPegTouchTimeArray2D_turn[i, j]
                if RPegTimeArray_turn[i, j] > LPegTimeArray_turn[i, j]:
                    PegTouch2DArray_sort[i, 2 * j + 1] = RPegTouchTimeArray2D[i, j]
                    PegTouch2DArray_sort[i, 2 * j] = LPegTouchTimeArray2D[i, j]
                    Peg2DArray[i, 2 * j + 1] = RPegTimeArray_turn[i, j]
                    Peg2DArray[i, 2 * j] = LPegTimeArray_turn[i, j]
                else:
                    PegTouch2DArray_sort[i, 2 * j] = RPegTouchTimeArray2D[i, j]
                    PegTouch2DArray_sort[i, 2 * j + 1] = LPegTouchTimeArray2D[i, j]
                    Peg2DArray[i, 2 * j] = RPegTimeArray_turn[i, j]
                    Peg2DArray[i, 2 * j + 1] = LPegTimeArray_turn[i, j]
        # PegTouch2DArray: L→R→L......の順に入れたもの
        # PegTouch2DArray_sort:(L, R)の組で早い方を先に入れた配列。

        PegHist = np.zeros ([turnnum, 24])
        PegHist_sort = np.zeros ([turnnum, 24])
        for i in range (turnnum):
            for j in range (24):
                if PegTouch2DArray[i, j] > 0:
                    PegHist[i, j] = 1
                if PegTouch2DArray_sort[i, j] > 0:
                    PegHist_sort[i, j] = 1
        ########################################################################################################################

        PHIS = np.zeros (24)
        Rpeg_number = np.arange (1, 13)
        Lpeg_number = np.arange (1, 13)
        peg_number = np.arange (1, 25)
        delete = []
        Rdelete = []
        Ldelete = []
        PHIS = np.sum (PegHist, axis=0)
        plt.imshow (PegHist)
        for i in range (24):
            if (PHIS[i] < turnnum / 3):
                delete.append (i)
                if i % 2 == 0:
                    Ldelete.append (i // 2)
                else:
                    Rdelete.append (i // 2)
        length = 24 - len (delete)
        Rlength = 12 - len (Rdelete)
        Llength = 12 - len (Ldelete)
        Rpeg_number = np.delete (Rpeg_number, Rdelete)
        Lpeg_number = np.delete (Lpeg_number, Ldelete)
        peg_number = np.delete (peg_number, delete)
        # 　取り敢えず消す配列を決めてみる

        nearpeg = []
        PegTouch2DArray_turn25 = np.insert (PegTouch2DArray_turn, 24, 0, axis=1)
        for i in range (24):
            PegTouch2DArray_turn25[i, 24] = PegTouch2DArray_turn[i + 1, 0]
        for i in range (24):
            if abs (np.average (PegTouch2DArray_turn25[:i + 1]) - np.average (PegTouch2DArray_turn25[:i])) < 80:
                nearpeg = np.append (nearpeg, i)
        for i in range (len (delete)):
            if (PHIS[delete[i]] < turnnum / 3):
                for j in range (turnnum):
                    if PegTouch2DArray_turn25[j, i] != PegTouch2DArray_turn25[j, i]:
                        PegTouch2DArray_turn25[j, i] = PegTouch2DArray_turn25[j, i + 1]

        # 前後で接近しているペグを選んでそのどちらかが抜けていれば片方にコピーする
        # 一周ごとに接近している違うペグを交互に使っているようなマウスがいる気がするのでこのような処理をする。
        # 果たしてこの処理が適当なのか？

        # 削る配列決め
        PegTouch2DArray_del = np.delete (PegTouch2DArray, delete, axis=1)
        PegTouch2DArray_sort_del = np.delete (PegTouch2DArray_sort, delete, axis=1)
        Peg2DArray_del = np.delete (PegTouch2DArray, delete, axis=1)
        PegTouch2DArray_turn_del = np.delete (PegTouch2DArray_turn, delete, axis=1)
        MedPegTouch = np.median (PegTouch2DArray_turn_del, axis=1)
        RPegTouchTimeArray2D_del = np.delete (RPegTouchTimeArray2D, Rdelete, axis=1)
        LPegTouchTimeArray2D_del = np.delete (LPegTouchTimeArray2D, Ldelete, axis=1)
        RPegTouchTimeArray2D_turn_del = np.delete (RPegTouchTimeArray2D_turn, Rdelete, axis=1)
        LPegTouchTimeArray2D_turn_del = np.delete (LPegTouchTimeArray2D_turn, Ldelete, axis=1)
        MedRPegTouch = np.median (RPegTouchTimeArray2D_turn_del, axis=0)
        MedLPegTouch = np.median (LPegTouchTimeArray2D_turn_del, axis=0)
        for i in range (len (RPegTouchTimeArray2D_turn)):
            for j in range (12):
                if RPegTouchTimeArray2D_turn[i, j] < 0:
                    RPegTouchTimeArray2D_turn[i, j] = RPegTouchTimeArray2D_turn[i, j] + 4000
        for i in range (len (LPegTouchTimeArray2D_turn)):
            for j in range (12):
                if LPegTouchTimeArray2D_turn[i, j] < 0:
                    LPegTouchTimeArray2D_turn[i, j] = LPegTouchTimeArray2D_turn[i, j] + 4000
        for i in range (len (MedRPegTouch)):
            if MedRPegTouch[i] < 0:
                MedRPegTouch[i] = MedRPegTouch[i] + 4000
                MedRPegTouch.sort ()
        for i in range (len (MedLPegTouch)):
            if MedLPegTouch[i] < 0:
                MedLPegTouch[i] = MedLPegTouch[i] + 4000
                MedLPegTouch.sort ()
        for i in range (len (MedPegTouch)):
            if MedPegTouch[i] < 0:
                MedPegTouch[i] = MedPegTouch[i] + 4000
                MedPegTouch.sort ()
        ######################################################################################################################

        PegTouch2DArray_del_2cycle = PegTouch2DArray_del
        for i in range (length):
            PegTouch2DArray_del_2cycle = np.insert (PegTouch2DArray_del_2cycle, len (PegTouch2DArray_del_2cycle[0]), 0,
                                                    axis=1)
            PegTouch2DArray_del_2cycle = np.insert (PegTouch2DArray_del_2cycle, 0, 0, axis=1)
        for i in range (1, turnnum - 1):
            for j in range (length):
                PegTouch2DArray_del_2cycle[i, j] = PegTouch2DArray_del[i - 1, j]
                PegTouch2DArray_del_2cycle[i, length * 2 + j] = PegTouch2DArray_del[i + 1, j]

        ########################################################################################################################
        Ave = np.zeros ((length, length))
        Var = np.zeros ((length, length))
        A = np.zeros ((length, length))

        for i in range (1, turnnum - 1):
            for j in range (length, length * 2):
                for k in range (length):
                    if (WaterOnArray[int (PegTouch2DArray_del_2cycle[i, j + 1 + k])] == 1) and (
                            WaterOnArray[int (PegTouch2DArray_del_2cycle[i, j])] == 1) and (
                            (PegTouch2DArray_del_2cycle[i, j + 1 + k] - PegTouch2DArray_del_2cycle[i, j]) > 0):
                        if (PegTouch2DArray_del_2cycle[i, j + 1 + k] != 0) and (PegTouch2DArray_del_2cycle[i, j] != 0):
                            Ave[j - length, k] = Ave[j - length, k] + (
                                        PegTouch2DArray_del_2cycle[i, j + 1 + k] - PegTouch2DArray_del_2cycle[i, j])
                            A[j - length, k] = A[j - length, k] + 1

        for i in range (length):
            for j in range (length):
                if not A[i, j] == 0:
                    Ave[i, j] = Ave[i, j] / A[i, j]
        for i in range (length):
            for j in range (length):
                A[i, j] = 0

        for i in range (turnnum-1):
            for j in range (length, length*2):
                for k in range (length):
                    if (WaterOnArray[int (PegTouch2DArray_del_2cycle[i, j+1+k])] == 1) and (
                            WaterOnArray[int (PegTouch2DArray_del_2cycle[i, j])] == 1):
                        if (Ave[j-length, k] - (PegTouch2DArray_del_2cycle[i, j+1+k] - PegTouch2DArray_del_2cycle[i, j])) ** 2 < thr:
                            Var[j-length, k] = Var[j-length, k] + (
                                        Ave[j-length, k] - (PegTouch2DArray_del_2cycle[i, j+1+k] - PegTouch2DArray_del_2cycle[i, j])) ** 2
                            A[j-length, k] = A[j-length, k] + 1

        for i in range (length):
            A[i, i] = 1
            Var[i, i] = 1
        for i in range (length):
            for j in range (length):
                if A[i, j] > (turnnum / 5):
                    Var[i, j] = Var[i, j] / A[i, j]
                else:
                    Var[i, j] = 0

        Varsum = np.sum (Var, axis=0)
        Varsum24 = Varsum
        for i in range (len (delete)):
             Varsum24 = np.insert (Varsum24, delete[i], 0)

        # ここまでVariance Heatmap
        ########################################################################################################################

        PhAve = np.zeros (length)
        E = np.zeros (length)

        for i in range (1, turnnum - 1):
            for j in range (5, length + 5):
                iprc = 2  # ipsi-precounter 同側の足で何歩もどるか(標準だと2歩前)
                ipoc = 2  # ipsi-postcounter　同側の足で何歩進むか(標準だと2歩後)
                cprc = 1  # contra-precounter 反対側の足で何歩戻るか(標準だと1歩前)
                cpoc = 1  # contra-postcounter 反対側の足で何歩進むか(標準だと1歩後)

                if WaterOnArray[int (PegTouch2DArray_del_2cycle[i][j])] == 0:
                    continue

                if PegTouch2DArray_del_2cycle[i][j - 2] != PegTouch2DArray_del_2cycle[i][j - 2]:
                    iprc = 4
                if PegTouch2DArray_del_2cycle[i][j + 2] != PegTouch2DArray_del_2cycle[i][j + 2]:
                    ipoc = 4
                if ipoc == 4:
                    cpoc = 3
                    if PegTouch2DArray_del_2cycle[i][j + 3] != PegTouch2DArray_del_2cycle[i][j + 3]:
                        cpoc = 1
                if iprc == 4:
                    if PegTouch2DArray_del_2cycle[i][j - 1] != PegTouch2DArray_del_2cycle[i][j - 1]:
                        cprc = 3
                ## target pegの2(1)歩前（後）がnanの場合4(3)歩前（後）を参照するようにする。ipsi(contra)
                ## 基本的には反対側の足の１歩の中で、２歩ついていた場合、遅い方の足を採用する
                ## よってpreではiprc=4の時cprc=1がデフォルトだが、ipoc=4の場合ではcpoc=3がデフォルトとなる

                if (low < (
                        PegTouch2DArray_del_2cycle[i][j + ipoc] - PegTouch2DArray_del_2cycle[i][j + cpoc]) < high) & (
                        low < (PegTouch2DArray_del_2cycle[i][j + cpoc] - PegTouch2DArray_del_2cycle[i][j]) < high) & (
                        low < (PegTouch2DArray_del_2cycle[i][j] - PegTouch2DArray_del_2cycle[i][j - cprc]) < high) & (
                        low < (
                        PegTouch2DArray_del_2cycle[i][j - cprc] - PegTouch2DArray_del_2cycle[i][j - iprc]) < high):
                    Phasepost = (PegTouch2DArray_del_2cycle[i][j + cpoc] - PegTouch2DArray_del_2cycle[i][j]) / (
                            PegTouch2DArray_del_2cycle[i][j + ipoc] - PegTouch2DArray_del_2cycle[i][j])
                    Phasepre = (PegTouch2DArray_del_2cycle[i][j - cprc] - PegTouch2DArray_del_2cycle[i][j - iprc]) / (
                            PegTouch2DArray_del_2cycle[i][j] - PegTouch2DArray_del_2cycle[i][j - iprc])
                    Phaseshift = Phasepost - Phasepre
                    if Phaseshift < 0.8:
                        PhAve[j - 5] = PhAve[j - 5] + Phaseshift
                        E[j - 5] = E[j - 5] + 1  ## jの値を 5～length+5 で動かしてるので j-5 にする
                        continue

        for i in range (length):
            if E[i] > (turnnum / 3):
                PhAve[i] = abs (PhAve[i] / E[i])
            else:
                PhAve[i] = 0
        PhAve24 = PhAve
        for i in range (len (delete)):
             PhAve24 = np.insert (PhAve24, delete[i], 0)  ##24行配列に直すために解析しなかったペグに0を入れた


        # ここまでPhase Shift
        ########################################################################################################################

        IntAve = np.zeros (length)
        F = np.zeros (length)
        for i in range (1, turnnum - 1):
            for j in range (length, length * 2):
                if WaterOnArray[int (PegTouch2DArray_del_2cycle[i][j])] == 0:
                    continue
                if (low < (PegTouch2DArray_del_2cycle[i, j] - PegTouch2DArray_del_2cycle[i, j - 2]) < Inthigh) & (
                        low < (PegTouch2DArray_del_2cycle[i, j + 2] - PegTouch2DArray_del_2cycle[i, j]) < Inthigh):
                    Intervalpost = (PegTouch2DArray_del_2cycle[i, j + 2] - PegTouch2DArray_del_2cycle[i, j])
                    Intervalpre = (PegTouch2DArray_del_2cycle[i, j] - PegTouch2DArray_del_2cycle[i, j - 2])
                    Intervalshift = abs (Intervalpost - Intervalpre)
                    IntAve[j - length] = IntAve[j - length] + Intervalshift
                    F[j - length] = F[j - length] + 1

                ### analize only good walking zone ex.)before and after 2steps are natural walking interval.
        for i in range (length):
            if F[i] > (turnnum / 4):
                IntAve[i] = abs (IntAve[i] / F[i])
            else:
                IntAve[i] = 0
        IntAve24 = IntAve
        for i in range (len (delete)):
            IntAve24 = np.insert (IntAve24, delete[i], 0)

        # ここまで interval shift
        ########################################################################################################################

        if np.sum (PhAve) > 0:
            PhAve = PhAve / np.sum (PhAve)
        if np.sum (IntAve) > 0:
            IntAve = IntAve / np.sum (IntAve)
        Var24 = Var
        for i in range(len(delete)):
            Var24 = np.insert(Var24, delete[i], 0, axis=0)
            Var24 = np.insert(Var24, delete[i], 0, axis=1)
        for i in range (24):
            PhaseTensor[number, num - start, i] = PhAve24[i]
            IntTensor[number, num - start, i] = IntAve24[i]
            VarTensor[number, num - start, i] = Varsum24[i]
            for j in range (length):
                VarTensorFlat[number, num - start, i * 24 + j] = Var24[i, j]

Tensor[Tensor == 0] = 0.0001
for i in range (n * day):
    for j in range (n * day):
        KL1[i, j] = entropy (Tensor[i // day, i % day], (Tensor[i // day, i % day] + Tensor[j // day, j % day]) / 2)
        KL2[i, j] = entropy (Tensor[j // day, j % day], (Tensor[i // day, i % day] + Tensor[j // day, j % day]) / 2)
JSD = (KL1 + KL2)

JSD_sqr = squareform (JSD)
Z = linkage (JSD, 'ward')

PhaseTensor[PhaseTensor == 0] = 0.01
IntTensor[IntTensor == 0] = 0.01
VarTensor[VarTensor == 0] = 0.01
VarTensorFlat[VarTensorFlat == 0] = 0.01
# ma.masked_less(PhaseTensor, 0.00001)
# ma.masked_less(IntTensor, 0.00001)
for i in range (n * day):
    for j in range (n * day):
        KLShift1[i, j] = entropy (PhaseTensor[i // day, i % day],
                                  (PhaseTensor[i // day, i % day] + PhaseTensor[j // day, j % day]) / 2)
        KLShift2[i, j] = entropy (PhaseTensor[j // day, j % day],
                                  (PhaseTensor[i // day, i % day] + PhaseTensor[j // day, j % day]) / 2)
        KLInt1[i, j] = entropy (IntTensor[i // day, i % day],
                                (IntTensor[i // day, i % day] + IntTensor[j // day, j % day]) / 2)
        KLInt2[i, j] = entropy (IntTensor[j // day, j % day],
                                (IntTensor[i // day, i % day] + IntTensor[j // day, j % day]) / 2)
        KLVar1[i, j] = entropy (VarTensor[i // day, i % day],
                                (VarTensor[i // day, i % day] + VarTensor[j // day, j % day]) / 2)
        KLVar2[i, j] = entropy (VarTensor[j // day, j % day],
                                (VarTensor[i // day, i % day] + VarTensor[j // day, j % day]) / 2)
        KLVarFlat1[i, j] = entropy (VarTensorFlat[i // day, i % day],
                                    (VarTensorFlat[i // day, i % day] + VarTensorFlat[j // day, j % day]) / 2)
        KLVarFlat2[i, j] = entropy (VarTensorFlat[j // day, j % day],
                                    (VarTensorFlat[i // day, i % day] + VarTensorFlat[j // day, j % day]) / 2)

JSDShift = (KLShift1 + KLShift2)
JSDInt = (KLInt1 + KLInt2)
JSDVar = (KLVar1 + KLVar2)
JSDVarFlat = (KLVarFlat1 + KLVarFlat2)

JSDShift_sqr = squareform(JSDShift)
ZShift = linkage(JSDShift, 'ward')
JSDInt_sqr = squareform(JSDInt)
ZInt = linkage(JSDInt, 'ward')
JSDVar_sqr = squareform(JSDVar)
ZVar = linkage(JSDVar, 'ward')
JSDVarFlat_sqr = squareform(JSDVarFlat)
ZVarFlat = linkage(JSDVarFlat, "ward")

JSD_flat = []
JSDPhase_flat = []
JSDInt_flat = []
JSDVar_flat = []

JSDlen = len(JSD)
for i in range(JSDlen):
    for j in range(JSDlen):
        if (JSD[i,j]<0.95)and(JSDInt[i,j]<0.95)and(JSDShift[i,j]<0.95)and(JSDVar[i,j]<0.95):
            JSD_flat.append(JSD[i,j])
            JSDInt_flat.append(JSDInt[i,j])
            JSDPhase_flat.append(JSDShift[i,j])
            JSDVar_flat.append(JSDVar[i,j])

# JSDShift = ma.masked_equal(JSDShift, 0)
# JSDShift = ma.masked_greater(JSDShift, 1)
# JSDInt = ma.masked_equal(JSDInt, 0)
# JSDInt = ma.masked_greater(JSDInt, 1)
for i in range(63):
    for j in range(63):
        if JSD[i,j] > 1:
            JSD[i,j] = 0.5
        if JSDVar[i,j] > 1:
            JSDVar[i,j] = 0.5
        if JSDInt[i,j] > 1:
            JSDInt[i,j] = 0.5
        if JSDShift[i,j] > 1:
            JSDShift[i,j] = 0.5
        if JSDVarFlat[i,j] > 1:
            JSDVarFlat[i,j] = 0.5
plt.figure(num=1)
plt.title("Jnesen Shannon divergence")
x = np.arange(1, n*day+1)
df = pd.DataFrame(data=JSD, index=x, columns=x)
sns.heatmap(df, cmap="coolwarm",robust=True, center=0.3)
plt.gca().xaxis.tick_top()
plt.tight_layout()

plt.figure(num=2)
plt.title("Jensen Shannon Divergence Phase")
x = np.arange(1, n*day+1)
df = pd.DataFrame(data=JSDShift, index=x, columns=x)
sns.heatmap(df, robust=True, cmap="coolwarm",center=0.5)
plt.gca().xaxis.tick_top()
plt.tight_layout()

plt.figure(num=3)
plt.title("Jensen Shannon Divergence Int")
x = np.arange(1, n*day+1)
df = pd.DataFrame(data=JSDInt, index=x, columns=x)
sns.heatmap(df, robust=True, cmap="coolwarm",center=0.3)
plt.gca().xaxis.tick_top()
plt.tight_layout()

plt.figure(num=4)
plt.title("Jensen Shannon Divergence Variance")
x = np.arange(1, n*day+1)
df = pd.DataFrame(data=JSDVar, index=x, columns=x)
sns.heatmap(df, robust=True, cmap="coolwarm",center=0.3)
plt.gca().xaxis.tick_top()
plt.tight_layout()

plt.figure (num=39)
plt.title ("Jensen Shannon Divergence Variance map")
x = np.arange (1, n * day + 1)
df = pd.DataFrame (data=JSDVarFlat, index=x, columns=x)
sns.heatmap (df, robust=True, cmap="coolwarm",center=0.3)
plt.gca ().xaxis.tick_top ()
plt.tight_layout ()

plt.figure(num=5)
dendrogram(Z, color_threshold=100)
plt.title("Clustering histogram")

plt.figure(num=6)
dendrogram(ZShift, color_threshold=100)
plt.title("Clustering Phase Shift")

plt.figure(num=7)
dendrogram(ZInt, color_threshold=100)
plt.title("Clustering Interval Shift")

plt.figure(num=8)
dendrogram(ZVar, color_threshold=100)
plt.title("Clustering Variance")

plt.figure(num=9)
x = JSD_flat
y = JSDPhase_flat
a, b = np.polyfit(x, y, 1)
plt.scatter(x, y)
x1 = np.linspace(0, 1)
y1 = a * x1 + b
plt.plot(x1, y1)
# data = pd.DataFrame([x, y])
# el = ConfidenceEllipse(data, p=0.95)
# plt.plot(el.get_patch(face_color="blue", alpha=0.5))
plt.title("Correlation Phase to touch")
plt.xlabel("JSD of Peg touch")
plt.ylabel("JSD of Phase Shift")

plt.figure(num=10)
y = JSDInt_flat
a, b = np.polyfit(x, y, 1)
plt.scatter(x, y)
x1 = np.linspace(0, 1)
y1 = a * x1 + b
plt.plot(x1, y1)
# data = pd.DataFrame([x, y])
# el = ConfidenceEllipse(data, p=0.95)
# plt.plot(el.get_patch(face_color="blue", alpha=0.5))
plt.title("Correlation Interval to touch")
plt.xlabel("JSD of Peg touch")
plt.ylabel("JSD of Interval Shift")

plt.figure(num=11)
y = JSDVar_flat
a, b = np.polyfit(x, y, 1)
plt.scatter(x, y)
x1 = np.linspace(0, 1)
y1 = a * x1 + b
plt.plot(x1, y1)
# data = pd.DataFrame([x, y])
# el = ConfidenceEllipse(data, p=0.95)
# plt.plot(el.get_patch(face_color="blue", alpha=0.5))
plt.title("Correlation Variance to touch")
plt.xlabel("JSD of Peg touch")
plt.ylabel("JSD of Variance")

plt.figure(num=21)
y = JSDVar_flat
x = JSDPhase_flat

clf = linear_model.LinearRegression()
x2 = [[X] for X in x]
clf.fit(x2, y)

plt.scatter(x, y)
plt.plot(x2, clf.predict(x2))
# data = pd.DataFrame([x, y])
# el = ConfidenceEllipse(data, p=0.95)
# plt.plot(el.get_patch(face_color="blue", alpha=0.5))
plt.title("Correlation Phase to variance")
plt.xlabel("JSD (Variance)")
plt.ylabel("JSD (Phase Shift)")
print("回帰係数",clf.coef_)
print("切片", clf.intercept_)
print("決定係数", clf.score(x2, y))


plt.figure(num=22)
y = JSDVar_flat
x = JSDInt_flat

clf = linear_model.LinearRegression()
x2 = [[X] for X in x]
clf.fit(x2, y)

plt.scatter(x, y)
plt.plot(x2, clf.predict(x2))
# data = pd.DataFrame([x, y])
# el = ConfidenceEllipse(data, p=0.95)
# plt.plot(el.get_patch(face_color="blue", alpha=0.5))
plt.title("Correlation Interval to Variance")
plt.xlabel("JSD (Variance)")
plt.ylabel("JSD (Phase Shift)")
print("回帰係数",clf.coef_)
print("切片", clf.intercept_)
print("決定係数", clf.score(x2, y))

plt.figure(num=12)
Kotainai = []
Kotaikan = []
for i in range(size):
    for j in range(i, size):
        if i//day == j//day:
            Kotainai.append(JSD[i, j])
        else:
            Kotaikan.append(JSD[i, j])
Kotaikan = random.sample(Kotaikan, len(Kotainai))
x = ["Intra", "Inter"]
y = [np.mean(Kotainai), np.mean(Kotaikan)]
sem = [np.std(Kotainai)/(np.sqrt(len(Kotainai))), np.std(Kotaikan)/(np.sqrt(len(Kotaikan)))]
plt.title("JSD-PegTouch")
plt.bar(x, y, yerr=sem, width= 0.4)

result = stats.ttest_ind(Kotainai, Kotaikan)
print("JSD peg touch (t-value, p-value)", result)

plt.figure(num=13)
Kotainai = []
Kotaikan = []
for i in range(size):
    for j in range(i, size):
        if i//day == j//day:
            Kotainai.append(JSDShift[i, j])
        else:
            Kotaikan.append(JSDShift[i, j])
Kotaikan = random.sample(Kotaikan, len(Kotainai))
x = ["Intra", "Inter"]
y = [np.mean(Kotainai), np.mean(Kotaikan)]
sem = [np.std(Kotainai)/(np.sqrt(len(Kotainai))), np.std(Kotaikan)/(np.sqrt(len(Kotaikan)))]
plt.title("JSD-PhaseShift")
plt.bar(x, y, yerr=sem, width=0.4)

result = stats.ttest_ind(Kotainai, Kotaikan)
print("JSD Phase shift (t-value, p-value)", result)

plt.figure(num=14)
Kotainai = []
Kotaikan = []
for i in range(size):
    for j in range(i, size):
        if i//day == j//day:
            Kotainai.append(JSDInt[i, j])
        else:
            Kotaikan.append(JSDInt[i, j])
Kotaikan = random.sample(Kotaikan, len(Kotainai))
x = ["Intra", "Inter"]
y = [np.mean(Kotainai), np.mean(Kotaikan)]
sem = [np.std(Kotainai)/(np.sqrt(len(Kotainai))), np.std(Kotaikan)/(np.sqrt(len(Kotaikan)))]
plt.title("JSD-IntervalShift")
plt.bar(x, y, yerr=sem, width = 0.4)

result = stats.ttest_ind(Kotainai, Kotaikan)
print("JSD Interval Shift (t-value, p-value)", result)

plt.figure(num=15)
Kotainai = []
Kotaikan = []
for i in range(size):
    for j in range(i, size):
        if i//day == j//day:
            Kotainai.append(JSDVar[i, j])
        else:
            Kotaikan.append(JSDVar[i, j])
Kotaikan = random.sample(Kotaikan, len(Kotainai))
x = ["Intra", "Inter"]
y = [np.mean(Kotainai), np.mean(Kotaikan)]
sem = [np.std(Kotainai)/(np.sqrt(len(Kotainai))), np.std(Kotaikan)/(np.sqrt(len(Kotaikan)))]
plt.title("JSD-Variance")
plt.bar(x, y, yerr=sem, width = 0.4)

result = stats.ttest_ind(Kotainai, Kotaikan)
print("JSD variance (t-value, p-value)", result)

plt.figure (num=40)
Kotainai = []
Kotaikan = []
for i in range (size):
    for j in range (i, size):
        if i // day == j // day:
            Kotainai.append (JSDVarFlat[i, j])
        else:
            Kotaikan.append (JSDVarFlat[i, j])
Kotaikan = random.sample (Kotaikan, len (Kotainai))
x = ["Intra", "Inter"]
y = [np.mean (Kotainai), np.mean (Kotaikan)]
sem = [np.std (Kotainai) / (np.sqrt (len (Kotainai))), np.std (Kotaikan) / (np.sqrt (len (Kotaikan)))]
plt.title ("JSD-Variancemap")
plt.bar (x, y, yerr=sem, width=0.4)

result = stats.ttest_ind (Kotainai, Kotaikan)
print ("JSD variance (t-value, p-value)", result)

plt.figure(num=16)
Kotainai_day = [[] for i in range(day-1)]
Kotaikan_day = [[] for i in range(day-1)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                Kotainai_day[i].append(JSD[day*j+i, day*k+i+1])
            else:
                Kotaikan_day[i].append(JSD[day*j+i, day*k+i+1])

x = np.arange(1,day)
Intra = np.mean(Kotainai_day, axis=1)
Inter = np.mean(Kotaikan_day, axis=1)
sem1 = np.zeros(day-1)
sem2 = np.zeros(day-1)
for i in range(day-1):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra")
plt.errorbar(x, Inter, yerr=sem2, label="Inter")
plt.title("Histogram JSD transition")
plt.legend()

result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotainai_day[2], Kotainai_day[3], Kotainai_day[4],
                              Kotainai_day[5], Kotainai_day[6], Kotainai_day[7])
print("Kotainai", result.pvalue)
result = scipy.stats.f_oneway(Kotaikan_day[0], Kotaikan_day[1], Kotaikan_day[2], Kotaikan_day[3], Kotaikan_day[4],
                              Kotaikan_day[5], Kotaikan_day[6], Kotaikan_day[7])
print("Kotaikan", result.pvalue)
for i in range(day-1):
    result = stats.ttest_ind(Kotainai_day[i], Kotaikan_day[i])
    print("JSD transition Histogram (t-value, p-value)", result)

plt.figure(num=17)
Kotainai_day = [[] for i in range(day-1)]
Kotaikan_day = [[] for i in range(day-1)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                Kotainai_day[i].append(JSDShift[day*j+i, day*k+i+1])
            else:
                Kotaikan_day[i].append(JSDShift[day*j+i, day*k+i+1])

x = np.arange(1,day)
Intra = np.mean(Kotainai_day, axis=1)
Inter = np.mean(Kotaikan_day, axis=1)
sem1 = np.zeros(day-1)
sem2 = np.zeros(day-1)
for i in range(day-1):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra")
plt.errorbar(x, Inter, yerr=sem2, label="Inter")
plt.title("Phase shift JSD transition")
plt.legend()

result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotainai_day[2], Kotainai_day[3], Kotainai_day[4],
                              Kotainai_day[5], Kotainai_day[6], Kotainai_day[7])
print("Kotainai", result.pvalue)
result = scipy.stats.f_oneway(Kotaikan_day[0], Kotaikan_day[1], Kotaikan_day[2], Kotaikan_day[3], Kotaikan_day[4],
                              Kotaikan_day[5], Kotaikan_day[6], Kotaikan_day[7])
print("Kotaikan", result.pvalue)
for i in range(day-1):
    result = stats.ttest_ind(Kotainai_day[i], Kotaikan_day[i])
    print("JSD transition Phase shift (t-value, p-value)", result)

plt.figure(num=18)
Kotainai_day = [[] for i in range(day-1)]
Kotaikan_day = [[] for i in range(day-1)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                Kotainai_day[i].append(JSDInt[day*j+i, day*k+i+1])
            else:
                Kotaikan_day[i].append(JSDInt[day*j+i, day*k+i+1])

x = np.arange(1,day)
Intra = np.mean(Kotainai_day, axis=1)
Inter = np.mean(Kotaikan_day, axis=1)
sem1 = np.zeros(day-1)
sem2 = np.zeros(day-1)
for i in range(day-1):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra")
plt.errorbar(x, Inter, yerr=sem2, label="Inter")
plt.title("Interval shift JSD transition")
plt.legend()

result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotainai_day[2], Kotainai_day[3], Kotainai_day[4],
                              Kotainai_day[5], Kotainai_day[6], Kotainai_day[7])
print("Kotainai", result.pvalue)
result = scipy.stats.f_oneway(Kotaikan_day[0], Kotaikan_day[1], Kotaikan_day[2], Kotaikan_day[3], Kotaikan_day[4],
                              Kotaikan_day[5], Kotaikan_day[6], Kotaikan_day[7])
print("Kotaikan", result.pvalue)
for i in range(day-1):
    result = stats.ttest_ind(Kotainai_day[i], Kotaikan_day[i])
    print("JSD transition Interval shift (t-value, p-value)", result)

plt.figure(num=19)
Kotainai_day = [[] for i in range(day-1)]
Kotaikan_day = [[] for i in range(day-1)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                Kotainai_day[i].append(JSDVar[day*j+i, day*k+i+1])
            else:
                Kotaikan_day[i].append(JSDVar[day*j+i, day*k+i+1])

x = np.arange(1,day)
Intra = np.mean(Kotainai_day, axis=1)
Inter = np.mean(Kotaikan_day, axis=1)
sem1 = np.zeros(day-1)
sem2 = np.zeros(day-1)
for i in range(day-1):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra")
plt.errorbar(x, Inter, yerr=sem2, label="Inter")
plt.title("Variance JSD transition")
plt.legend()

result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotainai_day[2], Kotainai_day[3], Kotainai_day[4],
                              Kotainai_day[5], Kotainai_day[6], Kotainai_day[7])
print("Kotainai", result.pvalue)
result = scipy.stats.f_oneway(Kotaikan_day[0], Kotaikan_day[1], Kotaikan_day[2], Kotaikan_day[3], Kotaikan_day[4],
                              Kotaikan_day[5], Kotaikan_day[6], Kotaikan_day[7])
print("Kotaikan", result.pvalue)
for i in range(day-1):
    result = stats.ttest_ind(Kotainai_day[i], Kotaikan_day[i])
    print("JSD transition Variance (t-value, p-value)", result)


plt.figure(num=23)
Kotainai_day = [[] for i in range(day-1)]
Kotaikan_day = [[] for i in range(day-1)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                Kotainai_day[i].append(JSD[day*j+i, day*k+day-1])
            else:
                Kotaikan_day[i].append(JSD[day*j+i, day*k+day-1])

x = np.arange(1,day)
Intra = np.mean(Kotainai_day, axis=1)
Inter = np.mean(Kotaikan_day, axis=1)
sem1 = np.zeros(day-1)
sem2 = np.zeros(day-1)
for i in range(day-1):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra")
plt.errorbar(x, Inter, yerr=sem2, label="Inter")
plt.title("Histogram JSD transition from lastday")
plt.legend()

result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotainai_day[2], Kotainai_day[3], Kotainai_day[4],
                              Kotainai_day[5], Kotainai_day[6], Kotainai_day[7])
print("Kotainai", result.pvalue)
result = scipy.stats.f_oneway(Kotaikan_day[0], Kotaikan_day[1], Kotaikan_day[2], Kotaikan_day[3], Kotaikan_day[4],
                              Kotaikan_day[5], Kotaikan_day[6], Kotaikan_day[7])
print("Kotaikan", result.pvalue)
for i in range(day-1):
    result = stats.ttest_ind(Kotainai_day[i], Kotaikan_day[i])
    print("JSD transition Histogram from lastday (t-value, p-value)", result)

plt.figure(num=24)
Kotainai_day = [[] for i in range(day-1)]
Kotaikan_day = [[] for i in range(day-1)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                Kotainai_day[i].append(JSDShift[day*j+i, day*k+day-1])
            else:
                Kotaikan_day[i].append(JSDShift[day*j+i, day*k+day-1])

x = np.arange(1,day)
Intra = np.mean(Kotainai_day, axis=1)
Inter = np.mean(Kotaikan_day, axis=1)
sem1 = np.zeros(day-1)
sem2 = np.zeros(day-1)
for i in range(day-1):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra")
plt.errorbar(x, Inter, yerr=sem2, label="Inter")
plt.title("Phase shift JSD transition from lastday")
plt.legend()

result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotainai_day[2], Kotainai_day[3], Kotainai_day[4],
                              Kotainai_day[5], Kotainai_day[6], Kotainai_day[7])
print("Kotainai", result.pvalue)
result = scipy.stats.f_oneway(Kotaikan_day[0], Kotaikan_day[1], Kotaikan_day[2], Kotaikan_day[3], Kotaikan_day[4],
                              Kotaikan_day[5], Kotaikan_day[6], Kotaikan_day[7])
print("Kotaikan", result.pvalue)
for i in range(day-1):
    result = stats.ttest_ind(Kotainai_day[i], Kotaikan_day[i])
    print("JSD transition Phase shift from lastday (t-value, p-value)", result)

plt.figure(num=25)
Kotainai_day = [[] for i in range(day-1)]
Kotaikan_day = [[] for i in range(day-1)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                Kotainai_day[i].append(JSDInt[day*j+i, day*k+day-1])
            else:
                Kotaikan_day[i].append(JSDInt[day*j+i, day*k+day-1])

x = np.arange(1,day)
Intra = np.mean(Kotainai_day, axis=1)
Inter = np.mean(Kotaikan_day, axis=1)
sem1 = np.zeros(day-1)
sem2 = np.zeros(day-1)
for i in range(day-1):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra")
plt.errorbar(x, Inter, yerr=sem2, label="Inter")
plt.title("Interval shift JSD transition from lastday")
plt.legend()

result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotainai_day[2], Kotainai_day[3], Kotainai_day[4],
                              Kotainai_day[5], Kotainai_day[6], Kotainai_day[7])
print("Kotainai", result.pvalue)
result = scipy.stats.f_oneway(Kotaikan_day[0], Kotaikan_day[1], Kotaikan_day[2], Kotaikan_day[3], Kotaikan_day[4],
                              Kotaikan_day[5], Kotaikan_day[6], Kotaikan_day[7])
print("Kotaikan", result.pvalue)
for i in range(day-1):
    result = stats.ttest_ind(Kotainai_day[i], Kotaikan_day[i])
    print("JSD transition Interval shift from lastday (t-value, p-value)", result)

plt.figure(num=26)
Kotainai_day = [[] for i in range(day-1)]
Kotaikan_day = [[] for i in range(day-1)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                Kotainai_day[i].append(JSDVar[day*j+i, day*k+day-1])
            else:
                Kotaikan_day[i].append(JSDVar[day*j+i, day*k+day-1])

x = np.arange(1,day)
print(Kotainai_day)
Intra = np.mean(Kotainai_day, axis=1)
Inter = np.mean(Kotaikan_day, axis=1)
sem1 = np.zeros(day-1)
sem2 = np.zeros(day-1)
for i in range(day-1):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra")
plt.errorbar(x, Inter, yerr=sem2, label="Inter")
plt.title("Variance JSD transition from lastday")
plt.legend()

result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotainai_day[2], Kotainai_day[3], Kotainai_day[4],
                              Kotainai_day[5], Kotainai_day[6], Kotainai_day[7])
print("Kotainai", result.pvalue)
result = scipy.stats.f_oneway(Kotaikan_day[0], Kotaikan_day[1], Kotaikan_day[2], Kotaikan_day[3], Kotaikan_day[4],
                              Kotaikan_day[5], Kotaikan_day[6], Kotaikan_day[7])
print("Kotaikan", result.pvalue)
for i in range(day-1):
    result = stats.ttest_ind(Kotainai_day[i], Kotaikan_day[i])
    print("JSD transition Variance from lastday (t-value, p-value)", result)

plt.figure(num=27)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 3.5:
                    Kotainai_day[0].append(JSD[day*j+i, day*k+day-1])
                else:
                    Kotainai_day[1].append(JSD[day*j+i, day*k+day-1])
            else:
                if i < 3.5:
                    Kotaikan_day[0].append(JSD[day*j+i, day*k+day-1])
                else:
                    Kotaikan_day[1].append(JSD[day*j+i, day*k+day-1])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 4 day","After 4 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Histogram JSD transition from lastday split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition Histogram from lastday First 4 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition Histogram from lastday Last 4 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])


plt.figure(num=28)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 3.5:
                    Kotainai_day[0].append(JSDShift[day*j+i, day*k+day-1])
                else:
                    Kotainai_day[1].append(JSDShift[day*j+i, day*k+day-1])
            else:
                if i < 3.5:
                    Kotaikan_day[0].append(JSDShift[day*j+i, day*k+day-1])
                else:
                    Kotaikan_day[1].append(JSDShift[day*j+i, day*k+day-1])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 4 day","After 4 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Phase Shift JSD transition from lastday split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition Phase from lastday First 4 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition Phase from lastday Last 4 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])


plt.figure(num=29)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 3.5:
                    Kotainai_day[0].append(JSDInt[day*j+i, day*k+day-1])
                else:
                    Kotainai_day[1].append(JSDInt[day*j+i, day*k+day-1])
            else:
                if i < 3.5:
                    Kotaikan_day[0].append(JSDInt[day*j+i, day*k+day-1])
                else:
                    Kotaikan_day[1].append(JSDInt[day*j+i, day*k+day-1])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 4 day","After 4 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Interval Shift JSD transition from lastday split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition interval from lastday First 4 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition interval from lastday Last 4 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])

plt.figure(num=30)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 3.5:
                    Kotainai_day[0].append(JSDVar[day*j+i, day*k+day-1])
                else:
                    Kotainai_day[1].append(JSDVar[day*j+i, day*k+day-1])
            else:
                if i < 3.5:
                    Kotaikan_day[0].append(JSDVar[day*j+i, day*k+day-1])
                else:
                    Kotaikan_day[1].append(JSDVar[day*j+i, day*k+day-1])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 4 day","After 4 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Variance JSD transition from lastday split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition variance from lastday First 4 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition variance from lastday Last 4 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])

plt.figure(num=31)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 2.5:
                    Kotainai_day[0].append(JSD[day*j+i, day*k+day-1])
                elif i > 4.5:
                    Kotainai_day[1].append(JSD[day*j+i, day*k+day-1])
            else:
                if i < 2.5:
                    Kotaikan_day[0].append(JSD[day*j+i, day*k+day-1])
                elif i > 4.5:
                    Kotaikan_day[1].append(JSD[day*j+i, day*k+day-1])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 3 day","After 3 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Histogram JSD transition from last3day split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition Histogram from lastday First 3 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition Histogram from lastday Last 3 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])


plt.figure(num=32)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 2.5:
                    Kotainai_day[0].append(JSDShift[day*j+i, day*k+day-1])
                elif i > 4.5:
                    Kotainai_day[1].append(JSDShift[day*j+i, day*k+day-1])
            else:
                if i < 2.5:
                    Kotaikan_day[0].append(JSDShift[day*j+i, day*k+day-1])
                elif i > 4.5:
                    Kotaikan_day[1].append(JSDShift[day*j+i, day*k+day-1])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 3 day","After 3 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Phase Shift JSD transition from last3day split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition Phase from lastday First 3 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition Phase from lastday Last 3 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])

plt.figure(num=33)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 2.5:
                    Kotainai_day[0].append(JSDInt[day*j+i, day*k+day-1])
                elif i > 4.5:
                    Kotainai_day[1].append(JSDInt[day*j+i, day*k+day-1])
            else:
                if i < 2.5:
                    Kotaikan_day[0].append(JSDInt[day*j+i, day*k+day-1])
                elif i > 4.5:
                    Kotaikan_day[1].append(JSDInt[day*j+i, day*k+day-1])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 3 day","After 3 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Interval Shift JSD transition from last3day split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition interval from lastday First 3 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition interval from lastday Last 3 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])

plt.figure(num=34)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 2.5:
                    Kotainai_day[0].append(JSDVar[day*j+i, day*k+day-1])
                elif i > 4.5:
                    Kotainai_day[1].append(JSDVar[day*j+i, day*k+day-1])
            else:
                if i < 2.5:
                    Kotaikan_day[0].append(JSDVar[day*j+i, day*k+day-1])
                elif i > 4.5:
                    Kotaikan_day[1].append(JSDVar[day*j+i, day*k+day-1])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 3 day","After 3 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Variance JSD transition from last3day split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition variance from lastday First 3 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition variance from lastday Last 3 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])

plt.figure(num=35)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-2):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 2.5:
                    Kotainai_day[0].append(JSD[day*j+i, day*k+day-1])
                    Kotainai_day[0].append(JSD[day*j+i, day*k+day-2])
                    Kotaikan_day[0].append(JSD[day*j+i, day*k+day-1])
                    Kotaikan_day[0].append(JSD[day*j+i, day*k+day-2])
                elif i > 3.5:
                    Kotainai_day[1].append(JSD[day*j+i, day*k+day-1])
                    Kotainai_day[1].append(JSD[day*j+i, day*k+day-2])
                    Kotaikan_day[1].append(JSD[day*j+i, day*k+day-1])
                    Kotaikan_day[1].append(JSD[day*j+i, day*k+day-2])
            else:
                if i < 2.5:
                    Kotaikan_day[0].append(JSD[day*j+i, day*k+day-1])
                    Kotaikan_day[0].append(JSD[day*j+i, day*k+day-2])
                elif i > 3.5:
                    Kotaikan_day[1].append(JSD[day*j+i, day*k+day-1])
                    Kotaikan_day[1].append(JSD[day*j+i, day*k+day-2])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 3 day","After 3 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Histogram JSD transition from last2day split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition Histogram from last2day First 3 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition Histogram from last2day Last 3 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])


plt.figure(num=36)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-2):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 2.5:
                    Kotainai_day[0].append (JSDShift[day * j + i, day * k + day - 1])
                    Kotainai_day[0].append (JSDShift[day * j + i, day * k + day - 2])
                    Kotaikan_day[0].append (JSDShift[day * j + i, day * k + day - 1])
                    Kotaikan_day[0].append (JSDShift[day * j + i, day * k + day - 2])
                elif i > 3.5:
                    Kotainai_day[1].append (JSDShift[day * j + i, day * k + day - 1])
                    Kotainai_day[1].append (JSDShift[day * j + i, day * k + day - 2])
                    Kotaikan_day[1].append (JSDShift[day * j + i, day * k + day - 1])
                    Kotaikan_day[1].append (JSDShift[day * j + i, day * k + day - 2])
            else:
                if i < 2.5:
                    Kotaikan_day[0].append (JSDShift[day * j + i, day * k + day - 1])
                    Kotaikan_day[0].append (JSDShift[day * j + i, day * k + day - 2])
                elif i > 3.5:
                    Kotaikan_day[1].append (JSDShift[day * j + i, day * k + day - 1])
                    Kotaikan_day[1].append (JSDShift[day * j + i, day * k + day - 2])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 3 day","After 3 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Phase Shift JSD transition from last2day split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition Phase from last2day First 3 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition Phase from last2day Last 3 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])

plt.figure(num=37)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-2):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 2.5:
                    Kotainai_day[0].append (JSDInt[day * j + i, day * k + day - 1])
                    Kotainai_day[0].append (JSDInt[day * j + i, day * k + day - 2])
                    Kotaikan_day[0].append (JSDInt[day * j + i, day * k + day - 1])
                    Kotaikan_day[0].append (JSDInt[day * j + i, day * k + day - 2])
                elif i > 3.5:
                    Kotainai_day[1].append (JSDInt[day * j + i, day * k + day - 1])
                    Kotainai_day[1].append (JSDInt[day * j + i, day * k + day - 2])
                    Kotaikan_day[1].append (JSDInt[day * j + i, day * k + day - 1])
                    Kotaikan_day[1].append (JSDInt[day * j + i, day * k + day - 2])
            else:
                if i < 2.5:
                    Kotaikan_day[0].append (JSDInt[day * j + i, day * k + day - 1])
                    Kotaikan_day[0].append (JSDInt[day * j + i, day * k + day - 2])
                elif i > 3.5:
                    Kotaikan_day[1].append (JSDInt[day * j + i, day * k + day - 1])
                    Kotaikan_day[1].append (JSDInt[day * j + i, day * k + day - 2])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 3 day","After 3 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Interval Shift JSD transition from last2day split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition interval from last2day First 3 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition interval from last2day Last 3 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])

plt.figure(num=38)
Kotainai_day = [[] for i in range(2)]
Kotaikan_day = [[] for i in range(2)]
for i in range(day-2):
    for j in range(7):
        for k in range(7):
            if j == k:
                if i < 2.5:
                    Kotainai_day[0].append (JSDVar[day * j + i, day * k + day - 1])
                    Kotainai_day[0].append (JSDVar[day * j + i, day * k + day - 2])
                    Kotaikan_day[0].append (JSDVar[day * j + i, day * k + day - 1])
                    Kotaikan_day[0].append (JSDVar[day * j + i, day * k + day - 2])
                elif i > 3.5:
                    Kotainai_day[1].append (JSDVar[day * j + i, day * k + day - 1])
                    Kotainai_day[1].append (JSDVar[day * j + i, day * k + day - 2])
                    Kotaikan_day[1].append (JSDVar[day * j + i, day * k + day - 1])
                    Kotaikan_day[1].append (JSDVar[day * j + i, day * k + day - 2])
            else:
                if i < 2.5:
                    Kotaikan_day[0].append (JSDVar[day * j + i, day * k + day - 1])
                    Kotaikan_day[0].append (JSDVar[day * j + i, day * k + day - 2])
                elif i > 3.5:
                    Kotaikan_day[1].append (JSDVar[day * j + i, day * k + day - 1])
                    Kotaikan_day[1].append (JSDVar[day * j + i, day * k + day - 2])

for i in range(2):
    Kotaikan_day[i] = random.sample(Kotaikan_day[i], len(Kotainai_day[i]))

x = ["First 3 day","After 3 day" ]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros(2)
sem2 = np.zeros(2)
for i in range(2):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar(x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim(0,0.8)
plt.title("Variance JSD transition from last2day split")
plt.legend()

result = stats.ttest_ind(Kotainai_day[0], Kotainai_day[1])
print("Kotainai_t", result.pvalue)
result = stats.ttest_ind(Kotaikan_day[0], Kotaikan_day[1])
print("Kotaikan_t", result.pvalue)
result = stats.ttest_ind(Kotainai_day[0], Kotaikan_day[0])
print("JSD transition variance from last2day First 3 day (t-value, p-value)", result)
result = stats.ttest_ind(Kotainai_day[1], Kotaikan_day[1])
print("JSD transition variance from last2day Last 3 day (t-value, p-value)", result)
result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])

plt.figure (num=41)
Kotainai_day = [[] for i in range (2)]
Kotaikan_day = [[] for i in range (2)]
for i in range (day - 2):
    for j in range (7):
        for k in range (7):
            if j == k:
                if i < 2.5:
                    Kotainai_day[0].append (JSDVarFlat[day * j + i, day * k + day - 1])
                    Kotainai_day[0].append (JSDVarFlat[day * j + i, day * k + day - 2])
                    Kotaikan_day[0].append (JSDVarFlat[day * j + i, day * k + day - 1])
                    Kotaikan_day[0].append (JSDVarFlat[day * j + i, day * k + day - 2])
                elif i > 3.5:
                    Kotainai_day[1].append (JSDVarFlat[day * j + i, day * k + day - 1])
                    Kotainai_day[1].append (JSDVarFlat[day * j + i, day * k + day - 2])
                    Kotaikan_day[1].append (JSDVarFlat[day * j + i, day * k + day - 1])
                    Kotaikan_day[1].append (JSDVarFlat[day * j + i, day * k + day - 2])
            else:
                if i < 2.5:
                    Kotaikan_day[0].append (JSDVarFlat[day * j + i, day * k + day - 1])
                    Kotaikan_day[0].append (JSDVarFlat[day * j + i, day * k + day - 2])
                elif i > 3.5:
                    Kotaikan_day[1].append (JSDVarFlat[day * j + i, day * k + day - 1])
                    Kotaikan_day[1].append (JSDVarFlat[day * j + i, day * k + day - 2])

for i in range (2):
    Kotaikan_day[i] = random.sample (Kotaikan_day[i], len (Kotainai_day[i]))

x = ["First 3 day", "After 3 day"]
Intra = np.zeros(2)
Intra[0] = np.mean(Kotainai_day[0])
Intra[1] = np.mean(Kotainai_day[1])
Inter = np.zeros(2)
Inter[0] = np.mean(Kotaikan_day[0])
Inter[1] = np.mean(Kotaikan_day[1])
sem1 = np.zeros (2)
sem2 = np.zeros (2)
for i in range (2):
    sem1[i] = np.std (Kotainai_day[i]) / (np.sqrt (len (Kotainai_day[i])))
    sem2[i] = np.std (Kotaikan_day[i]) / (np.sqrt (len (Kotaikan_day[i])))
plt.errorbar (x, Intra, yerr=sem1, label="Intra", color="k")
plt.errorbar (x, Inter, yerr=sem2, label="Inter", color="k")
plt.ylim (0, 0.8)
plt.title ("Variance JSD transition from last2day split")
plt.legend ()

result = stats.ttest_ind (Kotainai_day[0], Kotainai_day[1])
print ("Kotainai_t", result.pvalue)
result = stats.ttest_ind (Kotaikan_day[0], Kotaikan_day[1])
print ("Kotaikan_t", result.pvalue)
result = stats.ttest_ind (Kotainai_day[0], Kotaikan_day[0])
print ("JSD transition variancemap from last2day First 3 day (t-value, p-value)", result)
result = stats.ttest_ind (Kotainai_day[1], Kotaikan_day[1])
print ("JSD transition variancemap from last2day Last 3 day (t-value, p-value)", result)
result = scipy.stats.f_oneway (Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])
print ("Kotainai[0], Kotainai[1], Kotaikan[0], Kotaikan[1],ANOVA", result)
if result.pvalue < 0.05:
    tukey_hsd (list ("ABCD"), Kotainai_day[0], Kotainai_day[1], Kotaikan_day[0], Kotaikan_day[1])

# plt.figure(num=20)
# Kotainai = []
# Kotaikan = []
# Right = []
# Switch = []
# Null = []
# Modulation = []
# for i in range(size):
#     for j in range(i, size):
#         if i//day == j//day:
#             Kotainai.append(JSD[i, j])
#         else:
#             Kotaikan.append(JSD[i, j])
#             if (i//day in right)and(j//day in right):
#                 Right.append(JSD[i, j])
#             elif (i//day in switch)and(j//day in switch):
#                 Switch.append(JSD[i, j])
#             elif (i//day in null)and(j//day in null):
#                 Null.append(JSD[i, j])
#             elif (i // day in modulation) and (j // day in modulation):
#                 Modulation.append (JSD[i, j])
#
# x = ["Intra", "Inter", "Right", "Switch", "Null", "Modulation"]
# y = [np.mean(Kotainai), np.mean(Kotaikan), np.mean(Right), np.mean(Switch), np.mean(Null), np.mean(Modulation)]
# sem = [np.std(Kotainai)/(np.sqrt(len(Kotainai))), np.std(Kotaikan)/(np.sqrt(len(Kotaikan))), np.std(Right)/(np.sqrt(len(Right))), np.std(Switch)/(np.sqrt(len(Switch))), np.std(Null)/(np.sqrt(len(Null))), np.std(Modulation)/(np.sqrt(len(Modulation)))]
# plt.title("JSD-PegTouch")
# plt.bar(x, y, yerr=sem, width= 0.4)
#
# result = scipy.stats.f_oneway(Kotainai, Kotaikan, Right, Switch, Null, Modulation)
# print("IroIro comparison", result.pvalue)
# if result.pvalue < 0.05:
#     tukey_hsd(list("ABCDEF"), Kotainai, Kotaikan, Right, Switch, Null, Modulation)


plt.figure(num=16)
Kotainai_day = [[] for i in range(day-1)]
Kotaikan_day = [[] for i in range(day-1)]
for i in range(day-1):
    for j in range(7):
        for k in range(7):
            if j == k:
                Kotainai_day[i].append(JSD[day*j+i, day*k+i+1])
            else:
                Kotaikan_day[i].append(JSD[day*j+i, day*k+i+1])

x = np.arange(1,day)
Intra = np.mean(Kotainai_day, axis=1)
Inter = np.mean(Kotaikan_day, axis=1)
sem1 = np.zeros(day-1)
sem2 = np.zeros(day-1)
for i in range(day-1):
    sem1[i] = np.std(Kotainai_day[i])/(np.sqrt(len(Kotainai_day[i])))
    sem2[i] = np.std(Kotaikan_day[i])/(np.sqrt(len(Kotaikan_day[i])))
plt.errorbar(x, Intra, yerr=sem1, label="Intra")
plt.errorbar(x, Inter, yerr=sem2, label="Inter")
plt.title("Histogram JSD transition")
plt.legend()

result = scipy.stats.f_oneway(Kotainai_day[0], Kotainai_day[1], Kotainai_day[2], Kotainai_day[3], Kotainai_day[4],
                              Kotainai_day[5], Kotainai_day[6], Kotainai_day[7])
print("Kotainai", result.pvalue)
result = scipy.stats.f_oneway(Kotaikan_day[0], Kotaikan_day[1], Kotaikan_day[2], Kotaikan_day[3], Kotaikan_day[4],
                              Kotaikan_day[5], Kotaikan_day[6], Kotaikan_day[7])
print("Kotaikan", result.pvalue)
for i in range(day-1):
    result = stats.ttest_ind(Kotainai_day[i], Kotaikan_day[i])
    print("JSD transition Histogram (t-value, p-value)", result)

plt.show ()
