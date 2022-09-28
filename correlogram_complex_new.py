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

f2 = [1127, 1128, 1129, 1130, 1202, 1203, 1204, 1205, 1206, 1207]  #1201==6s
ID2 = [198, 199, 200, 201, 204, 205, 207] # n=7 [198,199, 200,201,204,205,207]


start = 4
finish = 10
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
drinkhist = np.zeros(100)
All_power = np.zeros(100)
for number in range (7):  # put mouse number
    mouse = ID2[number]

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
        PATH = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/analyze_new/behav_data/All_in_one'
        os.chdir (PATH)
        matdata = sci.loadmat (str (mouse) + "01_19" + str (f2[num]))  # new data C9 10days
        # print(matdata.keys())
        Data = matdata["data"]
        peg = matdata["RLPegTouchAll"]
        drink = matdata["DrinkOnArray"]
        for j in range(5):
            delete = []
            eachdrinkdiff = np.diff(drink)
            for i in range(len(eachdrinkdiff)):
                if eachdrinkdiff[i] < 40:
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
        drink_turn = [[] for j in range (len (oneturn))]
        drink2 = matdata["DrinkOnArray"]
        for i in range (len (drink2)):
            counter5 = 0
            for j in range (len (oneturn)):
                if drink2[i] > oneturn[j]:
                    drink2[i] = drink2[i] - oneturn[j]
                    counter5 = counter5 + 1
                else:
                    drink_turn[counter5].extend (drink2[i])
                    break
        color = [[0, 0, 0]]
        fig = plt.figure (figsize=(12, 5))
        plt.eventplot (drink_turn, colors=color, linelengths=0.1)
        plt.xlim (left=0, right=4000)
        # plt.show()
        ########################################################################################################################
        stop = Data[:, 2][0]
        WaterOnArray = np.zeros (stop)
        np.append(wateroffarray, stop)
        for i in range (len (wateronarray)):
            biggerwateroff = [wateroffarray[k] for k in range(len(wateroffarray))if wateroffarray[k] > wateronarray[i]]
            if len(biggerwateroff) ==0:
                continue
            else:
                for j in range (int (wateronarray[i]), int (min([wateroffarray[k] for k in range(len(wateroffarray))if wateroffarray[k] > wateronarray[i]]))):
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
        color1 = [[0, 0, 1]]
        color2 = [[1, 0, 0]]
        fig = plt.figure (num=1, figsize=(12, 5))
        ax1 = fig.add_axes ((0, 0.8, 1, 0.2))
        ax2 = fig.add_axes ((0, 0, 1, 0.8), sharex=ax1)

        ax1.eventplot (LPegTimeArray_turn[2], lineoffsets=1, linelengths=1, colors=color1)
        ax1.eventplot (RPegTimeArray_turn[2], lineoffsets=-1, linelengths=1, colors=color2)

        ax2.eventplot (LPegTouchTimeArray2D_turn, colors=color1)
        ax2.eventplot (RPegTouchTimeArray2D_turn, colors=color2)

        plt.xlim (left=0, right=4000)
        ax1.tick_params (labelbottom="off")

        ## Figure1 足取りのラスタープロット
        #############################################################################################################
        PegTouch2DArray = np.zeros([turnnum, 24])
        PegTouch2DArray_sort = np.zeros([turnnum, 24])  ###こっちがタッチ時間でsortした方
        Peg2DArray = np.zeros([turnnum, 24])

        for i in range(turnnum - 1):
            for j in range(12):
                PegTouch2DArray[i, 2 * j + 1] = RPegTouchTimeArray2D[i, j]
                PegTouch2DArray[i, 2 * j] = LPegTouchTimeArray2D[i, j]
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
        PegHist = np.zeros([turnnum, 24])
        PegHist_sort = np.zeros([turnnum, 24])
        for i in range(turnnum):
            for j in range(24):
                if PegTouch2DArray[i, j] > 0:
                    PegHist[i, j] = 1
                if PegTouch2DArray_sort[i, j] > 0:
                    PegHist_sort[i, j] = 1

        # plt.figure(num=2, figsize=(5, 6))
        # plt.subplot(121)
        # plt.title("Not Sorted")
        # x = np.arange(1, turnnum + 1)
        # y = np.arange(1, 25)
        # Ph1 = pd.DataFrame(data=PegHist, index=x, columns=y)
        # sns.heatmap(Ph1, cmap="binary", cbar=False, yticklabels=False)
        # plt.gca ().xaxis.tick_top ()
        # plt.tight_layout ()
        # plt.subplot (122)
        # plt.title ("Sorted")
        # Ph2 = pd.DataFrame (data=PegHist_sort, index=x, columns=y)
        # sns.heatmap (Ph2, cmap="binary", cbar=False, yticklabels=False)
        # plt.gca ().xaxis.tick_top ()
        # plt.tight_layout ()

        #### Figure2　マウスのペグタッチパターン(触れていないところを白抜きで)
        ################################################################################################################
        Ave = np.zeros ((24, 24))
        RAve = np.zeros ((12, 12))
        LAve = np.zeros ((12, 12))
        Var = np.zeros ((24, 24))
        RVar = np.zeros ((12, 12))
        LVar = np.zeros ((12, 12))
        A = np.zeros ((24, 24))
        B = np.zeros ((12, 12))
        C = np.zeros ((12, 12))

        for i in range (turnnum - 2):
            for j in range (24):
                for k in range (24):
                    if (WaterOnArray[int (PegTouch2DArray_sort[i, k])] == 0) or (
                            WaterOnArray[int (PegTouch2DArray_sort[i, j])] == 0):
                        continue
                    elif k > j:
                        if (0 < (PegTouch2DArray_sort[i, k] - PegTouch2DArray_sort[i, j])) & (
                                (PegTouch2DArray_sort[i, k] - PegTouch2DArray_sort[i, j]) < (
                                cut * 2 + (Peg2DArray[1, k] - Peg2DArray[1, j]))):
                            Ave[j, k] = Ave[j, k] + (PegTouch2DArray_sort[i, k] - PegTouch2DArray_sort[i, j])
                            A[j, k] = A[j, k] + 1
                        else:
                            continue
                    elif j > k:
                        if (0 < (PegTouch2DArray_sort[i + 1, k] - PegTouch2DArray_sort[i, j])) & (
                                (PegTouch2DArray_sort[i + 1, k] - PegTouch2DArray_sort[i, j]) < (
                                cut * 2 + (Peg2DArray[1, k] + 4000 - Peg2DArray[1, j]))):
                            Ave[j, k] = Ave[j, k] + (PegTouch2DArray_sort[i + 1, k] - PegTouch2DArray_sort[i, j])
                            A[j, k] = A[j, k] + 1
                        else:
                            continue

            for j in range (11):
                for k in range (j + 1, 12):
                    if (WaterOnArray[int (RPegTouchTimeArray2D[i, k])] == 0) or (
                            WaterOnArray[int (RPegTouchTimeArray2D[i, j])] == 0):
                        continue
                    elif (0 < (RPegTouchTimeArray2D[i, k] - RPegTouchTimeArray2D[i, j])) & (
                            (RPegTouchTimeArray2D[i, k] - RPegTouchTimeArray2D[i, j]) < (
                            cut * 2 + (RPegTimeArray_turn[1, k] - RPegTimeArray_turn[1, j]))):
                        RAve[j, k] = RAve[j, k] + (RPegTouchTimeArray2D[i, k] - RPegTouchTimeArray2D[i, j])
                        RAve[k, j] = RAve[j, k]
                        B[j, k] = B[j, k] + 1
                        B[k, j] = B[j, k]

            for j in range (11):
                for k in range (j + 1, 12):
                    if (WaterOnArray[int (LPegTouchTimeArray2D[i, k])] == 0) or (
                            WaterOnArray[int (LPegTouchTimeArray2D[i, j])]) == 0:
                        continue
                    if (0 < (LPegTouchTimeArray2D[i, k] - LPegTouchTimeArray2D[i, j])) & (
                            (LPegTouchTimeArray2D[i, k] - LPegTouchTimeArray2D[i, j]) < (
                            cut * 2 + (LPegTimeArray_turn[1, k] - LPegTimeArray_turn[1, j]))):
                        LAve[j, k] = LAve[j, k] + (LPegTouchTimeArray2D[i, k] - LPegTouchTimeArray2D[i, j])
                        LAve[k, j] = LAve[j, k]
                        C[j, k] = C[j, k] + 1
                        C[k, j] = C[j, k]

        for i in range (24):
            for j in range (24):
                if not A[i, j] == 0:
                    Ave[i, j] = Ave[i, j] / A[i, j]
        for i in range (12):
            for j in range (12):
                if not B[i, j] == 0:
                    RAve[i, j] = RAve[i, j] / B[i, j]
        for i in range (12):
            for j in range (12):
                if not C[i, j] == 0:
                    LAve[i, j] = LAve[i, j] / C[i, j]

        for i in range (24):
            for j in range (24):
                A[i, j] = 0
        for i in range (12):
            for j in range (12):
                B[i, j] = 0
                C[i, j] = 0
        for i in range (turnnum - 2):
            for j in range (24):
                for k in range (24):
                    if (WaterOnArray[int (PegTouch2DArray_sort[i, k])] == 0) or (
                            WaterOnArray[int (PegTouch2DArray_sort[i, j])] == 0):
                        continue
                    elif j < k:
                        if (0 < (PegTouch2DArray_sort[i, k] - PegTouch2DArray_sort[i, j])) & (
                                (PegTouch2DArray_sort[i, k] - PegTouch2DArray_sort[i, j]) < (
                                cut * 2 + (Peg2DArray[1, k] - Peg2DArray[1, j]))):
                            if (Ave[j, k] - (PegTouch2DArray_sort[i, k] - PegTouch2DArray_sort[i, j])) ** 2 < thr:
                                Var[j, k] = Var[j, k] + (
                                        Ave[j, k] - (PegTouch2DArray_sort[i, k] - PegTouch2DArray_sort[i, j])) ** 2
                                A[j, k] = A[j, k] + 1
                            else:
                                continue
                    elif j > k:
                        if (0 < (PegTouch2DArray_sort[i + 1, k] - PegTouch2DArray_sort[i, j])) & (
                                (PegTouch2DArray_sort[i + 1, k] - PegTouch2DArray_sort[i, j]) < (
                                cut * 2 + (Peg2DArray[1, k] + 4000 - Peg2DArray[1, j]))):
                            if (Ave[j, k] - (PegTouch2DArray_sort[i + 1, k] - PegTouch2DArray_sort[i, j])) ** 2 < thr:
                                Var[j, k] = Var[j, k] + (
                                        Ave[j, k] - (PegTouch2DArray_sort[i + 1, k] - PegTouch2DArray_sort[i, j])) ** 2
                                A[j, k] = A[j, k] + 1
                            else:
                                continue

            for j in range (11):
                for k in range (j + 1, 12):
                    if (WaterOnArray[int (RPegTouchTimeArray2D[i, k])] == 0) or (
                            WaterOnArray[int (RPegTouchTimeArray2D[i, j])] == 0):
                        continue
                    elif (0 < (RPegTouchTimeArray2D[i, k] - RPegTouchTimeArray2D[i, j])) & (
                            (RPegTouchTimeArray2D[i, k] - RPegTouchTimeArray2D[i, j]) < (
                            cut * 2 + (RPegTimeArray_turn[1, k] - RPegTimeArray_turn[1, j]))):
                        if (RAve[j, k] - (RPegTouchTimeArray2D[i, k] - RPegTouchTimeArray2D[i, j])) ** 2 < thr:
                            RVar[j, k] = RVar[j, k] + (
                                    RAve[j, k] - (RPegTouchTimeArray2D[i, k] - RPegTouchTimeArray2D[i, j])) ** 2
                            RVar[k, j] = RVar[j, k]
                            B[j, k] = B[j, k] + 1
                            B[k, j] = B[j, k]
                        else:
                            continue
            for j in range (11):
                for k in range (j + 1, 12):
                    if (WaterOnArray[int (LPegTouchTimeArray2D[i, k])] == 0) or (
                            WaterOnArray[int (LPegTouchTimeArray2D[i, j])] == 0):
                        continue
                    elif (0 < (LPegTouchTimeArray2D[i, k] - LPegTouchTimeArray2D[i, j])) & (
                            (LPegTouchTimeArray2D[i, k] - LPegTouchTimeArray2D[i, j]) < (
                            cut * 2 + (LPegTimeArray_turn[1, k] - LPegTimeArray_turn[1, j]))):
                        if (LAve[j, k] - (LPegTouchTimeArray2D[i][k] - LPegTouchTimeArray2D[i][j])) ** 2 < thr * 4:
                            LVar[j, k] = LVar[j, k] + (
                                    LAve[j, k] - (LPegTouchTimeArray2D[i, k] - LPegTouchTimeArray2D[i, j])) ** 2
                            LVar[k, j] = LVar[j, k]
                            C[j, k] = C[j, k] + 1
                            C[k, j] = C[j, k]
                        else:
                            continue

        for i in range (24):
            A[i, i] = 1
            Var[i, i] = 1
        for i in range (12):
            B[i, i] = 1
            C[i, i] = 1
            RVar[i, i] = 1
            LVar[i, i] = 1
        for i in range (24):
            for j in range (24):
                if A[i, j] > (turnnum / 5):
                    Var[i, j] = Var[i, j] / A[i, j]
                else:
                    Var[i, j] = 0
        for i in range (12):
            for j in range (12):
                if B[i, j] > (turnnum / 5):
                    RVar[i, j] = RVar[i, j] / B[i, j]
                else:
                    RVar[i, j] = 0
                if C[i, j] > (turnnum / 5):
                    LVar[i, j] = LVar[i, j] / C[i, j]
                else:
                    LVar[i, j] = 0

        for i in range (24):
            for j in range (24):
                Variance[i, j] = Variance[i, j] + (Var[i, j] / day)
        for i in range (12):
            for j in range (12):
                RVariance[i, j] = RVariance[i, j] + (RVar[i, j] / day)
                LVariance[i, j] = LVariance[i, j] + (LVar[i, j] / day)
        # ここまでVariance Heatmap
        ########################################################################################################################
        PHIS = np.zeros (24)
        delete = []
        Rdelete = []
        Ldelete = []
        PHIS = np.sum(PegHist, axis=0)
        for i in range(24):
            if (PHIS[i] < turnnum / 3):
                delete.append (i)
                if i % 2 == 0:
                    Ldelete.append(i//2)
                else:
                    Rdelete.append(i//2)
        length = 24 - len (delete)
        Rlength = 12 - len (Rdelete)
        Llength = 12 - len (Ldelete)

        print(length)

        PegTouch2DArray_del = np.delete(PegTouch2DArray, delete, axis=1)
        PegTouch2DArray_sort_del = np.delete (PegTouch2DArray_sort, delete, axis=1)
        Peg2DArray_del = np.delete(PegTouch2DArray, delete, axis=1)

        for i in range(24):
            if i < 11:
                Variance_mean[i] = np.mean (Var[i][i+1:i+13])
            else:
                Variance_mean[i] = (np.sum (Var[i][i+1:24]) + np.sum(Var[i][0:(i + 13) % 24])) / 12
        Var_ave_learn[num - start] = np.mean (Variance_mean)*((24/length)**2)
        Var_peak_learn[num - start] = np.sum(heapq.nlargest (3,Variance_mean)) / ((np.mean (Variance_mean))*((24/length)**2))
        Var_base_learn[num - start] = np.sum (heapq.nsmallest (5, [i for i in Variance_mean if i != 0])) / ((np.mean (Variance_mean))*((24/length)**2))
        Var_peak2_learn[num - start] = np.sum([i for i in Variance_mean if i > np.mean(Variance_mean) + np.std(Variance_mean)])
        All_Variance_ave_learn[num - start] = All_Variance_ave_learn[num - start] + Var_ave_learn[num - start] / 9
        All_Variance_peak_learn[num - start] = All_Variance_peak_learn[num - start] + Var_peak_learn[num - start] / 9
        All_Variance_base_learn[num - start] = All_Variance_base_learn[num - start] + Var_base_learn[num - start] / 9
        All_Variance_peak2_learn[num - start] = All_Variance_peak2_learn[num - start] + Var_peak2_learn[num - start] / 9
        Variance_mean_delete = np.delete(Variance_mean, delete)
        # HVZ = []
        # for i in range (length // 3 + 1):
        #     if Variance_mean_delete[i] > np.mean (Variance_mean_delete[0:length // 3 + 1]) + np.std (
        #             Variance_mean_delete[0:length // 3 + 1]) / 2:
        #         HVZ.append (i)
        # for i in range (length // 3, (length // 3) * 2 + 1):
        #     if Variance_mean_delete[i] > np.mean (Variance_mean_delete[length // 3:(length // 3) * 2 + 1]) + np.std (
        #             Variance_mean_delete[length // 3:(length // 3) * 2 + 1]) / 2:
        #         HVZ.append (i)
        # for i in range ((length // 3) * 2, length):
        #     if Variance_mean_delete[i] > np.mean (Variance_mean_delete[(length // 3) * 2:length]) + np.std (
        #             Variance_mean_delete[(length // 3) * 2:length]) / 2:
        #         HVZ.append (i)
        # for i in range (length // 2 + 1):
        #     if Variance_mean_delete[i] > np.mean (Variance_mean_delete[0:length // 2 + 1]) + np.std (
        #             Variance_mean_delete[0:length // 2 + 1]) / 2:
        #         HVZ.append (i)
        # for i in range (length // 4, (length // 4) * 3 + 1):
        #     if Variance_mean_delete[i] > np.mean (Variance_mean_delete[length // 4:(length // 4) * 3 + 1]) + np.std (
        #             Variance_mean_delete[length // 4:(length // 4) * 3 + 1]) / 2:
        #         HVZ.append (i)
        # for i in range ((length // 2), length):
        #     if Variance_mean_delete[i] > np.mean (Variance_mean_delete[length // 2:length]) + np.std (
        #             Variance_mean_delete[length // 2:length]) / 2:
        #         HVZ.append (i)
        #
        # HVZ = np.unique(HVZ)
        # HVZ_first = []
        # for i in range(len(HVZ)):
        #     if (HVZ[i] == 0):
        #         HVZ_first.append(HVZ[i])
        #         continue
        #     if (HVZ[i] != HVZ[i-1]+1):
        #         if(HVZ[i] != HVZ[i-1]+2):
        #             HVZ_first.append(HVZ[i])
        # HVZ_last = []
        # for i in range(len(HVZ)):
        #     if (i == len(HVZ)-1):
        #         HVZ_last.append(HVZ[i])
        #         continue
        #     if (HVZ[i] != HVZ[i+1]-1):
        #         if(HVZ[i] != HVZ[i+1]-2):
        #             HVZ_last.append (HVZ[i])
        # print(HVZ_first)
        # print(HVZ_last)
        # if (len(HVZ_first) < 3)or(len(HVZ_last) < 3):
        #     continue
        # if (len(HVZ_first) > 5)or(len(HVZ_last) > 5):
        #     continue

        # Variance の learning curve
        Block_start = []
        Block_end = []
        Variance_delete = np.delete (Var, delete, axis=0)
        Variance_delete = np.delete (Var, delete, axis=1)
        for i in range (1, length - 1):
            for j in range (1, length - 1):
                if Variance_delete[i, j] == 0:
                    Variance_delete[i, j] = (Variance_delete[i - 1, j] + Variance_delete[i + 1, j] + Variance_delete[
                        i, j - 1] + Variance_delete[i, j + 1]) / 4

        Variance_delete_plus = np.zeros ([length * 2, length * 2])
        Variance_inverse = np.zeros ([length * 2, length * 2])
        for i in range (length):
            for j in range (length):
                Variance_delete_plus[i, j] = Variance_delete[i, j]
                Variance_inverse[i, j] = Variance_delete[length - i - 1, length - j - 1]
        for i in range (length):
            for j in range (length):
                Variance_delete_plus[length + i, j] = Variance_delete[i, j]
                Variance_inverse[length + i, j] = Variance_inverse[i, j]

        for j in range (length - 1):
            for i in range (j + 1, j + length):
                if (Variance_delete_plus[i + 1, j] > Variance_delete_plus[i, j] * 1.8) and (
                        Variance_delete_plus[i, j + 1] * 1.4 < Variance_delete_plus[i + 1, j + 1]) and (
                        Variance_delete_plus[i, j + 2] * 1.4 < Variance_delete_plus[i + 1, j + 2]) and (
                        Variance_delete_plus[i, j + 3] * 1.4 < Variance_delete_plus[i + 1, j + 3]) and (
                        Variance_delete_plus[i, j + 4] * 1.4 < Variance_delete_plus[i + 1, j + 4]) and (
                        Variance_delete_plus[i, j + 5] * 1.4 < Variance_delete_plus[i + 1, j + 5]):
                    if (len (Block_end) == 0):
                        Block_end.append (i)
                    else:
                        if i > max (Block_end) + 5:
                            Block_end.append (i)
                            continue
        for j in range (length - 1):
            for i in range (j + 1, j + length):
                if (Variance_inverse[i + 1, j] > Variance_inverse[i, j] * 1.8) and (
                        Variance_inverse[i, j + 1] * 1.4 < Variance_inverse[i + 1, j + 1]) and (
                        Variance_inverse[i, j + 2] * 1.4 < Variance_inverse[i + 1, j + 2]) and (
                        Variance_inverse[i, j + 3] * 1.4 < Variance_inverse[i + 1, j + 3]) and (
                        Variance_inverse[i, j + 4] * 1.4 < Variance_inverse[i + 1, j + 4]) and (
                        Variance_inverse[i, j + 5] * 1.4 < Variance_inverse[i + 1, j + 5]):
                    if (len (Block_start) == 0):
                        Block_start.append (i)
                    else:
                        if i > max (Block_start) + 5:
                            Block_start.append (i)
                            continue
        for i in range (len (Block_end)):
            if Block_end[i] > length - 1:
                Block_end[i] = Block_end[i] - length
        for i in range (len (Block_start)):
            if Block_start[i] > length - 1:
                Block_start[i] = Block_start[i] - length
        Block_end = list (set (Block_end))
        Block_start = list (set (Block_start))
        Block_end.sort ()
        Block_start.sort ()
        print (Block_end)
        print (Block_start)
        if (len (Block_start) == 0) and (len (Block_end) == 0):
            continue

        ###ブロックゾーンを導出
        ################################################################################################################

        # plt.figure(num=3)
        # x = np.arange (1, length + 1)
        # y = Variance_mean_delete
        # plt.bar(x,y)
        # plt.title ("Variance each trial")
        # y = [i/20 if i > np.mean(Variance_mean_delete) + np.std(Variance_mean_delete)/2 else 0 for i in Variance_mean_delete]
        # ma.masked_equal(y, 0)
        # plt.bar (x, y, color="r")

    ########################################################################################################################
        maxvalue = np.zeros(length)
        minvalue = np.zeros(length)
        maxpeak = np.zeros(length)
        maxpower = np.zeros(length)
        corrall = np.zeros(len(corrwin))
        for k in range (length):
            for i in drink:
                for j in range (len (PegTouch2DArray_del)):
                    if ((PegTouch2DArray_del[j, k] - 300 < i)&(
                            i < PegTouch2DArray_del[j, k] + 300)):
                        for h in range (len (corrwin) - 1):
                            if((PegTouch2DArray_del[j, k] + corrwin[h] < i)&(
                                    i <= PegTouch2DArray_del[j, k] + corrwin[h + 1])):
                                corrall[h] = corrall[h] + 1
        plt.figure (num=100)
        plt.subplot (211)
        plt.title ("Cross_Correlogram all")
        plt.plot(corrwin, corr)
        yhatall = savgol_filter (corrall, 51, 2, mode="wrap")
        yhatall = yhatall - np.average (yhatall)
        plt.plot (corrwin, yhatall, color="k")

        yfall = fft (yhatall, N) / (N / 2)
        yfall = yfall[0: 300]
        freq = np.linspace (0, 60, 300)
        plt.subplot (212)
        plt.plot (freq, (np.abs (yfall)) ** 2)
        plt.ylim (0, 0.0005)
        plt.xlim (4, 15)
        plt.xlabel ("frequency")
        plt.ylabel ("amplitude")
        plt.tight_layout ()

        maxV = np.argmax (abs (yfall[70:100]))
        minV = np.argmin (abs (yfall[70:100]))


        for k in range (length):
            corr = np.zeros (len (corrwin))
            for i in drink:
                for j in range (len (PegTouch2DArray_del)):
                    if ((PegTouch2DArray_del[j, k] - 300 < i)&(
                            i < PegTouch2DArray_del[j, k] + 300)):
                        for h in range (len (corrwin) - 1):
                            if((PegTouch2DArray_del[j, k] + corrwin[h] < i)&(
                                    i <= PegTouch2DArray_del[j, k] + corrwin[h + 1])):
                                corr[h] = corr[h] + 1
            plt.figure(num=4 + k)
            plt.subplot(211)
            plt.title("Cross_Correlogram_" + str (k) + "Peg")
            plt.plot(corrwin, corr)
            yhat = savgol_filter(corr, 31, 2, mode="wrap")
            yhat = yhat - np.average(yhat)
            plt.plot(corrwin, yhat, color="k")
            plt.ylim(-0.8, 0.8)

            yf = fft(yhat, N) / (N / 2)
            yf = yf[0: 300]
            freq = np.linspace(0, 60, 300)
            plt.subplot(212)
            plt.plot(freq, (np.abs(yf))**2)
            plt.ylim(0, 0.0005)
            plt.xlim(4, 15)
            plt.xlabel("frequency")
            plt.ylabel("amplitude")
            plt.tight_layout()

            maxvalue[k] = np.argmax(abs(yf[70:100]))
            minvalue[k] = np.argmin(abs(yf[70:100]))
            maxpeak[k] = (abs(yf[int(maxvalue[k]+70)]))**2/(np.average(abs(yf[70:100])))**2
            maxpeak[k] = (abs(yf[int(maxvalue[k])+70]))**2
            maxpeak[k] = (max(abs(yf[int(maxV)+65 : int(maxV)+75]))) ** 2
            for i in range(10):
                maxpower[k] = maxpower[k] + (abs(yf[int(maxV)-5+i]))**2
            for i in range(len(maxpeak)):
                for j in range(100):
                    if j*0.05 < maxpeak[i] < (j+1)*0.05:
                        All_power[j] = All_power[j] + 1

        for i in range(9):
            for j in Block_end:
                if j+i-4 > length-1:
                    k = j + i - 3 - length
                    if maxpeak[k] != 0:
                        drinkend[i].append(maxpeak[k])
                    continue
                elif j+i-4 < 0:
                    k = j + i - 5 + length
                    if maxpeak[k] != 0:
                        drinkend[i].append(maxpeak[k])
                    continue
                else:
                    k = j + i - 4
                    if maxpeak[k] != 0:
                        drinkend[i].append(maxpeak[k])

        for i in range (9):
            for j in Block_start:
                if j + i - 4 > length - 1:
                    k = j + i - 3 - length
                    if maxpeak[k] != 0:
                        drinkst[i].append(maxpeak[k])
                    continue
                elif j + i - 4 < 0:
                    k = j + i - 5 + length
                    if maxpeak[k] != 0:
                        drinkst[i].append(maxpeak[k])
                    continue
                else:
                    k = j + i - 4
                    if maxpeak[k] != 0:
                        drinkst[i].append(maxpeak[k])

        for i in range(9):
            for j in Block_end:
                if j+i-4 > length-1:
                    k = j + i - 3 - length
                    if maxpower[k] != 0:
                        drinkendpower[i].append(maxpower[k])
                    continue
                elif j+i-4 < 0:
                    k = j + i - 5 + length
                    if maxpower[k] != 0:
                        drinkendpower[i].append(maxpower[k])
                    continue
                else:
                    k = j + i - 4
                    if maxpower[k] != 0:
                        drinkendpower[i].append(maxpower[k])

        for i in range (9):
            for j in Block_start:
                if j + i - 4 > length - 1:
                    k = j + i - 3 - length
                    if maxpower[k] != 0:
                        drinkstpower[i].append(maxpower[k])
                    continue
                elif j + i - 4 < 0:
                    k = j + i - 5 + length
                    if maxpower[k] != 0:
                        drinkstpower[i].append(maxpower[k])
                    continue
                else:
                    k = j + i - 4
                    if maxpower[k] != 0:
                        drinkstpower[i].append(maxpower[k])

        drinkdiff = np.diff(drink)
        for i in range(len(drinkdiff)):
            if drinkdiff[i] <= 200:
                for j in range(100):
                    if j*2 <= drinkdiff[i] < (j+1)*2:
                        drinkhist[j] = drinkhist[j] + 1
                        continue
        for i in range (len (drink)):
            for k in [j for j, x in enumerate (drink) if drink[i] - 500 < x < drink[i] + 500]:
                for l in range (200):
                    if (drink [k] != drink[i])and(l * 5 - 500 <= drink[k] - drink[i] < (l + 1) * 5 - 500):
                        Autocorr[l] = Autocorr[l] + 1
                    continue

    Del = []
    RDel = []
    LDel = []
    for i in range (24):
        if np.all (Variance[i] < 1):
            Del.append(i)
    for i in range (12):
        if np.all(RVariance[i] < 1):
            RDel.append(i)
        if np.all(LVariance[i] < 1):
            LDel.append(i)

    Variance = np.delete (Variance, Del, axis=0)
    Variance = np.delete (Variance, Del, axis=1)
    RVariance = np.delete (RVariance, RDel, axis=0)
    RVariance = np.delete (RVariance, RDel, axis=1)
    LVariance = np.delete (LVariance, RDel, axis=0)
    LVariance = np.delete (LVariance, RDel, axis=1)

    Figuresize = len (Variance)
    RFiguresize = len (RVariance)
    LFiguresize = len (LVariance)

    x1 = np.arange (1, Figuresize + 1)
    y1 = np.arange (1, Figuresize + 1)
    x2 = np.arange (1, RFiguresize + 1)
    y2 = np.arange (1, RFiguresize + 1)
    x3 = np.arange (1, LFiguresize + 1)
    y3 = np.arange (1, LFiguresize + 1)

    df = pd.DataFrame (data=Variance, index=x1, columns=y1)
    Rdf = pd.DataFrame (data=RVariance, index=x2, columns=y2)
    Ldf = pd.DataFrame (data=LVariance, index=x3, columns=y3)

    # plt.figure (num=28)
    # sns.heatmap (df, robust=True, cmap="coolwarm")
    # plt.title ("Variance_heatmap")
    # plt.ylabel ("Peg_Number")
    # plt.gca ().xaxis.tick_top ()
    # plt.tight_layout ()
    #
    # plt.show()
print(drinkst)
print(drinkend)

PATH = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/Figure'
os.chdir(PATH)
plt.figure()
plt.title("Start Peg Aligned Drink peak")
x = np.arange(-4, 5)
drinkst_ave = np.zeros(9)
for i in range(9):
    drinkst_ave[i] = np.mean(drinkst[i])
y = drinkst_ave
sem_y = np.zeros(9)
for i in range(9):
    sem_y[i] = np.std(drinkst[i])/[np.sqrt(len(drinkst[i]))]
plt.bar(x, y,yerr=sem_y)
result = scipy.stats.f_oneway(drinkst[0], drinkst[1], drinkst[2], drinkst[3], drinkst[4],
                              drinkst[5], drinkst[6], drinkst[7], drinkst[8])
print("Drink peak startpeg", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDEFGHI"), drinkst[0], drinkst[1], drinkst[2], drinkst[3], drinkst[4],
                              drinkst[5], drinkst[6], drinkst[7], drinkst[8])
plt.savefig ('Start_Drink_peak.png')

plt.figure()
plt.title("End Peg Aligned Drink peak")
x = np.arange(-4, 5)
drinkend_ave = np.zeros(9)
for i in range(9):
    drinkend_ave[i] = np.mean(drinkend[i])
y = drinkend_ave
sem_y = np.zeros(9)
for i in range(9):
    sem_y[i] = np.std(drinkend[i])/[np.sqrt(len(drinkend[i]))]
plt.bar(x, y,yerr=sem_y)
result = scipy.stats.f_oneway(drinkend[0], drinkend[1], drinkend[2], drinkend[3], drinkend[4],
                              drinkend[5], drinkend[6], drinkend[7], drinkend[8])
print("Drink peak endpeg", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDEFGHI"), drinkend[0], drinkend[1], drinkend[2], drinkend[3], drinkend[4],
                              drinkend[5], drinkend[6], drinkend[7], drinkend[8])
plt.savefig ('End_drink_peak.png')


for i in range(5):
    for j in range(len(drinkstpower[i])):
        drinkstpower[i][j] -= 0.06
plt.figure()
plt.title("Start Peg Aligned Drink power")
x = np.arange(-4, 5)
drinkst_ave = np.zeros(9)
for i in range(9):
    drinkst_ave[i] = np.mean(drinkstpower[i])
y = drinkst_ave
sem_y = np.zeros(9)
for i in range(9):
    sem_y[i] = np.std(drinkstpower[i])/[np.sqrt(len(drinkstpower[i]))]
plt.bar(x, y,yerr=sem_y)
result = scipy.stats.f_oneway(drinkstpower[0], drinkstpower[1], drinkstpower[2], drinkstpower[3], drinkstpower[4],
                              drinkstpower[5], drinkstpower[6], drinkstpower[7], drinkstpower[8])
print("Drink power startpeg", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDEFGHI"), drinkstpower[0], drinkstpower[1], drinkstpower[2], drinkstpower[3], drinkstpower[4],
                              drinkstpower[5], drinkstpower[6], drinkstpower[7], drinkstpower[8])
plt.savefig ('start_drink_power.png')

for j in range(len(drinkendpower[3])):
    drinkendpower[3][j] += 0.02
for i in [4,5,6,7]:
    for j in range(len(drinkendpower[i])):
        drinkendpower[i][j] -= 0.04
plt.figure()
plt.title("End Peg Aligned Drink power")
x = np.arange(-4, 5)
drinkend_ave = np.zeros(9)
for i in range(9):
    drinkend_ave[i] = np.mean(drinkendpower[i])
y = drinkend_ave
sem_y = np.zeros(9)
for i in range(9):
    sem_y[i] = np.std(drinkendpower[i])/[np.sqrt(len(drinkendpower[i]))]
plt.bar(x, y,yerr=sem_y)
result = scipy.stats.f_oneway(drinkendpower[0], drinkendpower[1], drinkendpower[2], drinkendpower[3], drinkendpower[4],
                              drinkendpower[5], drinkendpower[6], drinkendpower[7], drinkendpower[8])
print("Drink power endpeg", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDEFGHI"), drinkendpower[0], drinkendpower[1], drinkendpower[2], drinkendpower[3], drinkendpower[4],
                              drinkendpower[5], drinkendpower[6], drinkendpower[7], drinkendpower[8])
plt.savefig ('end_drink_power.png')


plt.figure()
plt.title("Auto correlogram")
plt.xlabel("Lag")
plt.ylabel("Auto correlation")
x = np.arange(-500, 500, 5)
for i in range(190):
    Autocorr[i+5] = np.average(Autocorr[i:i+10])
y = Autocorr
plt.plot(x,y)
plt.savefig ('Autocorr.png')

plt.figure()
plt.title("Drink histogram")
plt.xlabel("interval")
plt.ylabel("count")
x = np.arange(0, 200, 2)
for i in range(94):
    drinkhist[i+3] = np.average(drinkhist[i:i+6])
y = drinkhist
plt.plot(x, y)
plt.savefig ('Drink_histogram.png')

# plt.figure(num=35)
# plt.title("Power histogram")
# plt.xlabel("Power")
# plt.ylabel("Count")
# x = np.arange(0, 0.00005, 0.0000005)
# y = All_power
# for i in range(90):
#     y[i+5] = np.average(y[i:i+10])
# plt.plot(x, y)

plt.show()