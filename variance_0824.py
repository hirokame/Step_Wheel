import csv
import scipy.io as sci
import scipy
import os
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import heapq
import math
import numpy.ma as ma
from scipy.stats import chi2
from scipy import signal
from scipy import stats

np.set_printoptions (threshold=np.inf, linewidth=1000, suppress=True)
def tukey_hsd( ind, *args):
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    data_arr = np.hstack(args)

    ind_arr = np.array([])
    for x in range(len(args)):
        ind_arr = np.append(ind_arr, np.repeat(ind[x], len(args[x])))
    res = pairwise_tukeyhsd(data_arr, ind_arr, alpha=0.05)
    print(pairwise_tukeyhsd(data_arr, ind_arr, alpha=0.05))
    print(vars(res))
    
def calcint(n):
    return int(n)

f = [213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 301, 302, 303, 304, 305, 306, 307,
     308, 309, 310]  ## day 0~3, 4~13, 14~20, 21~25
f2 = [1127, 1128, 1129, 1130, 1202, 1203, 1204, 1205, 1206, 1207]  #1201==6s
ID = [1101, 1201, 1301, 1701, 1801, 1901, 2101, 2201, 2301]  # n=9
ID2 = [198, 199, 200, 201, 204, 205, 207] # n=7 [198,199, 200,201,204,205,207]

start = 5
finish = 10
day = finish - start
Inthigh = 800
high = 600
low = 0
cut = 250
thr = 20000

Phst = [[] for i in range(11)]
Phend = [[] for i in range(11)]
Phend_R = [[] for i in range(11)]
Phend_L = [[] for i in range(11)]
Phinst_st = [[] for i in range(11)]
Phinst_end = [[] for i in range(11)]
Intst = [[] for i in range(15)]
Intend = [[] for i in range(15)]
Intend_R = [[] for i in range(15)]
Intend_L = [[] for i in range(15)]
Intst_R = [[] for i in range(15)]
Intst_L = [[] for i in range(15)]
Intend_ideal = [[] for i in range(15)]
Intst_ideal = [[] for i in range(15)]
Inst_st = [[] for i in range(11)]
Inst_end = [[] for i in range(11)]
Cost = [[] for i in range(11)]
Coend = [[] for i in range(11)]
whole_interval= []
whole_phaseR=[]
whole_phaseL=[]
whole_intervalshift = []
whole_phaseshift = []
modulation = [[] for i in range(7)]##[intend_R,intend_L,intst_R,intst_L]*7 の順に格納

Intervalshift_in = []
Phaseshift_in = []
Intervalshift_out = []
Phaseshift_out = []

Intervalshift_out_boundary = []
Phaseshift_out_boundary = []

for number in range (7):  # put mouse number
    # mouse = ID[number]
    mouse = ID2[number]  # new data C9 10days
    print (mouse)

    Variance = np.zeros ((24, 24))
    RVariance = np.zeros ((12, 12))
    LVariance = np.zeros ((12, 12))

    Covariance = np.zeros (24)
    PSave = np.zeros (24)
    ISave = np.zeros (24)
    RISave = np.zeros (12)
    LISave = np.zeros (12)

    Variance_mean = np.zeros (24)
    RVariance_mean = np.zeros (12)
    LVariance_mean = np.zeros (12)
    
    Intst_R_mouse = [[] for i in range(5)]
    Intst_L_mouse = [[] for i in range(5)]
    Intend_R_mouse = [[] for i in range(5)]
    Intend_L_mouse = [[] for i in range(5)]
    
    for num in range (start, finish):  # put file nuber
        print(num)
        PATH = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/analyze_new/behav_data/All_in_one'
        os.chdir (PATH)
        # matdata = sci.loadmat (str (mouse) + "_170" + str (f[num]))  # create file name day 0~3, 4~13, 14~20, 21~25
        matdata = sci.loadmat (str (mouse) + "01_19" + str (f2[num]))  # new data C9  day0~9
        
        Data = matdata["data"]
        '''
        Dataの内容
        ２：タスクの合計時間
        ４：一周にかかった時間(Turnmarker)
        １３〜２５：右足のペグタッチタイミング
        ４１〜５２：左足のペグタッチタイミング
        '''
        peg = matdata["RLPegTouchAll"]
        RPegTimeArray = matdata["RpegTimeArray2D"] ## ペグの時間が格納されていいる(turn数*12)
        LPegTimeArray = matdata["LpegTimeArray2D"]
        RPegTimeArray_turn = matdata["RpegTimeArray2D_turn"]##ペグの時間(ターンマーカーでアライン)
        LPegTimeArray_turn = matdata["LpegTimeArray2D_turn"]
        RpegTouchCell = matdata["RpegTouchCell"] ## 1*12*(Totalタッチ数)の配列
        LpegTouchCell = matdata["LpegTouchCell"] ## ペグごとに１周の中で最初の一発のみを配列にぶちこんだもの。
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
        turnnum = len (oneturn)
        RPegTouchTimeArray2D = np.zeros ([turnnum, 12])
        RPegTouchTimeArray2D_turn = np.zeros ([turnnum, 12])
        TurnmarkerTime = 0
        for i in range (turnnum - 1):
            for j in range (12):
#                 R = Data[:, j + 13][np.nonzero (Data[:, j + 13])]
                R = RpegTouchCell[0][j]
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
        LPegTouchTimeArray2D = np.zeros ([turnnum, 12])
        LPegTouchTimeArray2D_turn = np.zeros ([turnnum, 12])
        for i in range (turnnum - 1):
            for j in range (12):
#                 L = Data[:, j + 41][np.nonzero (Data[:, j + 41])]
                L = LpegTouchCell[0][j]
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

        PegTouch2DArray = np.zeros ([turnnum, 24])
        PegTouch2DArray_sort = np.zeros ([turnnum, 24])  #タッチ時間でsortした方
        Peg2DArray = np.zeros ([turnnum, 24])
        PegTouch2DArray_turn = np.zeros ([turnnum, 24])
        for j in range(turnnum-2):
            for i in range(3):
                if (0 < RPegTouchTimeArray2D_turn[j,9+i] < 500):
                    if (RPegTouchTimeArray2D_turn[j,9+i] < RPegTouchTimeArray2D_turn[j+1,0]):
                        RPegTouchTimeArray2D_turn[j,9+i] += oneturn[j]
                    else:
                        RPegTouchTimeArray2D_turn[j,i+9] = 0
                if (0 < LPegTouchTimeArray2D_turn[j, 9 + i] < 500):
                    if (LPegTouchTimeArray2D_turn[j, 9 + i] < LPegTouchTimeArray2D_turn[j + 1, 0]):
                        LPegTouchTimeArray2D_turn[j, 9 + i] += oneturn[j]
                    else:
                        LPegTouchTimeArray2D_turn[j, i + 9] = 0
                if (0 < RPegTimeArray_turn[j, 9 + i] < 500):
                    if (RPegTimeArray_turn[j, 9 + i] < RPegTimeArray_turn[j + 1, 0]):
                        RPegTimeArray_turn[j, 9 + i] += oneturn[j]
                    else:
                        RPegTimeArray_turn[j, i + 9] = 0
                if (0 < LPegTimeArray_turn[j, 9 + i] < 500):
                    if (LPegTimeArray_turn[j, 9 + i] < LPegTimeArray_turn[j + 1, 0]):
                        LPegTimeArray_turn[j, 9 + i] += oneturn[j]
                    else:
                        LPegTimeArray_turn[j, i + 9] = 0
        for j in range(turnnum-2):
            for i in range(2):
                if (0 < RPegTouchTimeArray2D[j,10+i] < RPegTouchTimeArray2D[j, 9])or(0 < RPegTouchTimeArray2D[j,10+i] < RPegTouchTimeArray2D[j, 1])or(0 < RPegTouchTimeArray2D[j,10+i] < RPegTouchTimeArray2D[j, 4]):
                    if (RPegTouchTimeArray2D[j,10+i] < RPegTouchTimeArray2D[j+1,0]):
                        RPegTouchTimeArray2D[j,10+i] += oneturn[j]
                    else:
                        RPegTouchTimeArray2D[j,i+10] = 0
                if (0 < LPegTouchTimeArray2D[j,10+i] < LPegTouchTimeArray2D[j, 9])or(0 < LPegTouchTimeArray2D[j,10+i] < LPegTouchTimeArray2D[j, 1])or(0 < LPegTouchTimeArray2D[j,10+i] < LPegTouchTimeArray2D[j, 4]):
                    if (LPegTouchTimeArray2D[j,10+i] < LPegTouchTimeArray2D[j+1,0]):
                        LPegTouchTimeArray2D[j,10+i] += oneturn[j]
                    else:
                        LPegTouchTimeArray2D[j,i+10] = 0
                if (0 < RPegTimeArray[j,10+i] < RPegTimeArray[j, 9])or(0 < RPegTimeArray[j,10+i] < RPegTimeArray[j, 1])or(0 < RPegTimeArray[j,10+i] < RPegTimeArray[j, 4]):
                    if (RPegTimeArray[j,10+i] < RPegTimeArray[j+1,0]):
                        RPegTimeArray[j,10+i] += oneturn[j]
                    else:
                        RPegTimeArray[j,i+10] = 0
                if (0 < LPegTimeArray[j,10+i] < LPegTimeArray[j, 9])or(0 < LPegTimeArray[j,10+i] < LPegTimeArray[j, 1])or(0 < LPegTimeArray[j,10+i] < LPegTimeArray[j, 4]):
                    if (LPegTimeArray[j,10+i] < LPegTimeArray[j+1,0]):
                        LPegTimeArray[j,10+i] += oneturn[j]
                    else:
                        LPegTimeArray[j,i+10] = 0


        for i in range(2,turnnum-1):
            for j in range(12):
                if (RPegTouchTimeArray2D[i, j] == 0)and(RPegTouchTimeArray2D[(i-1),j]!=0):
                    if np.random.random() < 0.7:
                        RPegTouchTimeArray2D[i, j] = RPegTouchTimeArray2D[i-1,j] + oneturn[i-1] + 100*(np.random.random()-0.5)
                elif (RPegTouchTimeArray2D[i, j] == 0)and(RPegTouchTimeArray2D[(i-2),j]!=0):
                    if np.random.random() < 0.4:
                        RPegTouchTimeArray2D[i, j] = RPegTouchTimeArray2D[i-2,j]+oneturn[i-1]+oneturn[i-2]+100*(np.random.random() - 0.5)
        for i in range(2,turnnum-1):
            for j in range(12):
                if (LPegTouchTimeArray2D[i, j] == 0)and(LPegTouchTimeArray2D[(i-1),j]!=0):
                    if np.random.random() < 0.7:
                        LPegTouchTimeArray2D[i, j] = LPegTouchTimeArray2D[i-1,j] + oneturn[i-1] + 100*(np.random.random()-0.5)
                elif (LPegTouchTimeArray2D[i, j] == 0)and(LPegTouchTimeArray2D[(i-2),j]!=0):
                    if np.random.random() < 0.4:
                        LPegTouchTimeArray2D[i, j] = LPegTouchTimeArray2D[i-2,j]+oneturn[i-1]+oneturn[i-2]+100*(np.random.random() - 0.5)

        for i in range(len(RPegTouchTimeArray2D)-1):
            for j in range(12):
                RPegTouchTimeArray2D[i,j] = int(RPegTouchTimeArray2D[i,j])
                LPegTouchTimeArray2D[i,j] = int(LPegTouchTimeArray2D[i,j])
                if RPegTimeArray[i,j] == RPegTimeArray[i,j]:
                    RPegTimeArray[i,j] = int(RPegTimeArray[i,j])
                if LPegTimeArray[i,j] == LPegTimeArray[i,j]:
                    LPegTimeArray[i,j] = int(LPegTimeArray[i,j])

        file_path = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/StepWheel_data_hirokane/RightTouchNo'+str(number)+'Day'+str(num)+'.txt'
        with open(file_path, 'w') as file:
            writer = csv.writer (file, lineterminator='\n')
            writer.writerows(RPegTouchTimeArray2D)

        file_path = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/StepWheel_data_hirokane/LeftTouchNo'+str(number)+'Day'+str(num)+'.txt'
        with open (file_path, 'w') as file:
            writer = csv.writer (file, lineterminator='\n')
            writer.writerows (LPegTouchTimeArray2D)

        file_path = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/StepWheel_data_hirokane/RightPegNo'+str(number)+'Day'+str(num)+'.txt'
        with open (file_path, 'w') as file:
            writer = csv.writer (file, lineterminator='\n')
            writer.writerows (RPegTimeArray)

        file_path = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/StepWheel_data_hirokane/LeftPegNo'+str(number)+'Day'+str(num)+'.txt'
        with open (file_path, 'w') as file:
            writer = csv.writer (file, lineterminator='\n')
            writer.writerows (LPegTimeArray)

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

        plt.figure (num=2, figsize=(5, 6))
        plt.subplot (121)
        plt.title ("Not Sorted")
        x = np.arange (1, turnnum + 1)
        y = np.arange (1, 25)
        Ph1 = pd.DataFrame (data=PegHist, index=x, columns=y)
        sns.heatmap (Ph1, cmap="binary", cbar=False, yticklabels=False)
        plt.gca ().xaxis.tick_top ()
        plt.tight_layout ()
        plt.subplot (122)
        plt.title ("Sorted")
        Ph2 = pd.DataFrame (data=PegHist_sort, index=x, columns=y)
        sns.heatmap (Ph2, cmap="binary", cbar=False, yticklabels=False)
        plt.gca ().xaxis.tick_top ()
        plt.tight_layout ()

        # Figure2　マウスのペグタッチパターン(触れていないところを白抜きで)
        # それぞれのトライアルごとに出力
        ################################################################################################################
        PHIS = np.zeros (24)
        Rpeg_number = np.arange (2, 26, 2)
        Lpeg_number = np.arange (1, 25, 2)
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
        #　取り敢えず消す配列を決めてみる

        nearpeg = []
        PegTouch2DArray_turn25 = np.insert(PegTouch2DArray_turn, 24, 0, axis=1)
        for i in range(24):
            PegTouch2DArray_turn25[i,24] = PegTouch2DArray_turn[i+1,0]
        for i in range(24):
            if abs(np.average(PegTouch2DArray_turn25[:i+1])-np.average(PegTouch2DArray_turn25[:i])) < 80:
                nearpeg = np.append(nearpeg, i)
        for i in range(len(delete)):
            if (PHIS[delete[i]] < turnnum / 3):
                for j in range(turnnum):
                    if PegTouch2DArray_turn25[j,i] != PegTouch2DArray_turn25[j,i]:
                        PegTouch2DArray_turn25[j,i] = PegTouch2DArray_turn25[j, i+1]

        #前後で接近しているペグを選んでそのどちらかが抜けていれば片方にコピーする
        #一周ごとに接近している違うペグを交互に使っているようなマウスがいる気がするのでこのような処理をする。
        #果たしてこの処理が適当なのか？
        for i in range (turnnum):
            for j in range (24):
                if PegTouch2DArray_turn25[i, j] > 0:
                    PegHist[i, j] = 1


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

        # 削る配列を決めた
        ########################################################################################################################

        plt.figure
        allpeghistR = np.zeros (4000)
        allpeghistL = np.zeros (4000)
        R_peghist = np.zeros (4000)
        L_peghist = np.zeros (4000)
        Rturn = len (RPegTouchTimeArray2D_turn)
        Lturn = len (LPegTouchTimeArray2D_turn)
        for i in range (Rturn - 1):
            for j in range (Rlength):
                for k in range (4000):
                    if k < RPegTouchTimeArray2D_turn_del[i, j] <= (k + 1):
                        allpeghistR[k] = allpeghistR[k] + 1
                        continue
        for i in range (Lturn - 1):
            for j in range (Llength):
                for k in range (4000):
                    if k < LPegTouchTimeArray2D_turn_del[i, j] <= (k + 1):
                        allpeghistL[k] = allpeghistL[k] + 1
                        continue

        x = np.arange (0, 4000)
        for i in range (3900):
            R_peghist[i + 50] = np.mean (allpeghistR[i:i + 100])
        for i in range (3900):
            L_peghist[i + 50] = np.mean (allpeghistL[i:i + 100])
        plt.subplot (211)
        plt.plot (x, R_peghist, color="b")
        plt.subplot (212)
        plt.plot (x, L_peghist, color="r")
        PATH = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/Figure'
        os.chdir (PATH)
        plt.savefig (str (number) + '_' + str (num) + '_histogram.png')
        
        
        color1 = [[0, 0, 1]]
        color2 = [[1, 0, 0]]
        color3 = [[0, 0.2, 0.8]]
        color4 = [[0.8, 0.2, 0]]
        #
        # fig = plt.figure (figsize=(12, 8))
        # ax1 = fig.add_axes ((0, 0.9, 1, 0.1))
        # ax2 = fig.add_axes ((0, 0.4, 1, 0.4), sharex=ax1)
        # ax3 = fig.add_axes ((0, 0, 1, 0.4), sharex=ax1)
        # ax4 = fig.add_axes ((0, 0.8, 1, 0.1), sharex=ax1)
        #
        # ax1.eventplot (LPegTimeArray_turn[2], lineoffsets=1, linelengths=1, colors="r")
        # ax1.eventplot (RPegTimeArray_turn[2], lineoffsets=-1, linelengths=1, colors="b")
        #
        # ax2.eventplot (LPegTouchTimeArray2D_turn, colors="r")
        # ax3.eventplot (RPegTouchTimeArray2D_turn, colors="b")
        #
        # # ax3.plot (grad_instantaneous_phase)
        #
        # ax4.eventplot (MedLPegTouch, lineoffsets=1, linelengths=1, colors="r")
        # ax4.eventplot (MedRPegTouch, lineoffsets=-1, linelengths=1, colors="b")
        #
        # ax1.tick_params (labelbottom="off")
        # PATH = 'C:\\Users\\hirokane\\Desktop\\Paper_Figure\\new_analysis'
        # os.chdir (PATH)
        # plt.savefig(str(number)+'_'+str(num)+'_raster.png')

        # Figure1 足取りのラスタープロット
        ################################################################################################################

        PegTouch2DArray_del_2cycle = PegTouch2DArray_del
        for i in range (length):
            PegTouch2DArray_del_2cycle = np.insert(PegTouch2DArray_del_2cycle, len(PegTouch2DArray_del_2cycle[0]), 0, axis=1)
            PegTouch2DArray_del_2cycle = np.insert(PegTouch2DArray_del_2cycle, 0, 0, axis=1)
        for i in range (1,turnnum-1):
            for j in range (length):
                PegTouch2DArray_del_2cycle[i, j] = PegTouch2DArray_del[i-1, j]
                PegTouch2DArray_del_2cycle[i, length*2+j] = PegTouch2DArray_del[i+1, j]

        # PegTouch2DArray_del_2cycle　という配列に前後1周ぶん配列をinsert、条件分けせずに処理
        ################################################################################################################
        # ここからVariance解析)

        Ave = np.zeros ((length, length))
        RAve = np.zeros ((12, 12))
        LAve = np.zeros ((12, 12))
        Var = np.zeros ((length, length))
        RVar = np.zeros ((12, 12))
        LVar = np.zeros ((12, 12))
        A = np.zeros ((length, length))
        B = np.zeros ((12, 12))
        C = np.zeros ((12, 12))

        for i in range (1,turnnum-1):
            for j in range (length, length*2):
                for k in range (length):
                    if (WaterOnArray[int (PegTouch2DArray_del_2cycle[i, j+1+k])] == 1) and (
                            WaterOnArray[int (PegTouch2DArray_del_2cycle[i, j])] == 1) and (
                            (PegTouch2DArray_del_2cycle[i, j + 1 + k] - PegTouch2DArray_del_2cycle[i, j]) > 0):
                        if (PegTouch2DArray_del_2cycle[i, j+1+k] != 0) and (PegTouch2DArray_del_2cycle[i, j] != 0):
                            Ave[j-length, k] = Ave[j-length, k]+(PegTouch2DArray_del_2cycle[i, j+1+k] - PegTouch2DArray_del_2cycle[i, j])
                            A[j-length, k] = A[j-length, k] + 1


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
            for j in range (length):
                if A[i, j] > (turnnum / 3):
                    Var[i, j] = Var[i, j] / A[i, j]
                else:
                    Var[i, j] = 1

        # ここまでVariance Heatmap
        ########################################################################################################################
        # Covariance 解析

        CovAve = np.zeros (length)
        Covar = np.zeros (length)
        D = np.zeros (length)
        for i in range (1, turnnum-1):
            for j in range (length, length*2):
                if (WaterOnArray[int (PegTouch2DArray_del_2cycle[i, j])] == 0):
                    continue
                if (low < (PegTouch2DArray_del_2cycle[i, j] - PegTouch2DArray_del_2cycle[i, j-1])) and (
                        (PegTouch2DArray_del_2cycle[i, j] - PegTouch2DArray_del_2cycle[i, j-1]) < high):
                    CovAve[j-length] = CovAve[j-length] + (PegTouch2DArray_del_2cycle[i, j] - PegTouch2DArray_del_2cycle[i, j-1])
                    D[j-length] = D[j-length] + 1

        for i in range (length):
            if D[i] > (turnnum / 4):
                CovAve[i] = CovAve[i] / D[i]
            else:
                CovAve[i] = 0

        D = np.zeros (length)
        for i in range (1,turnnum-1):
            for j in range (length, length*2):
                if WaterOnArray[int (PegTouch2DArray_del_2cycle[i,j])] == 0:
                    continue
                if (low < (PegTouch2DArray_del_2cycle[i,j+1] - PegTouch2DArray_del_2cycle[i,j]) < high) & (
                        low < (PegTouch2DArray_del_2cycle[i,j] - PegTouch2DArray_del_2cycle[i,j-1]) < high) & (
                        PegTouch2DArray_del_2cycle[i,j] != 0) & (PegTouch2DArray_del_2cycle[i,j-1] != 0) & (
                        PegTouch2DArray_del_2cycle[i,j+1] != 0):
                    if -30000 < (CovAve[j-length] - (PegTouch2DArray_del_2cycle[i,j] - PegTouch2DArray_del_2cycle[i,j-1])) * (
                            CovAve[(j+1-length)%length] - (PegTouch2DArray_del_2cycle[i,j+1] - PegTouch2DArray_del_2cycle[i,j])) < 0:
                        Covar[j-length] = Covar[j-length] + (
                                CovAve[j-length] - (PegTouch2DArray_del_2cycle[i,j] - PegTouch2DArray_del_2cycle[i,j-1])) * (
                                           CovAve[(j+1-length)%length]-(PegTouch2DArray_del_2cycle[i,j+1] - PegTouch2DArray_del_2cycle[i,j]))
                        D[j-length] = D[j-length] + 1
                        continue

        for i in range (length):
            if D[i] > (turnnum / 4):
                Covar[i] = Covar[i] / (D[i])
            else:
                Covar[i] = 0
        Covar24 = Covar
        for i in range (len (delete)):
            Covar24 = np.insert (Covar24, delete[i], 0)
        for i in range (24):
            Covariance[i] = Covariance[i] + Covar24[i] / day

        ####ここまでCovariance
        
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

                if PegTouch2DArray_del_2cycle[i][j-2] != PegTouch2DArray_del_2cycle[i][j-2]:
                    iprc = 4
                if PegTouch2DArray_del_2cycle[i][j+2] != PegTouch2DArray_del_2cycle[i][j+2]:
                    ipoc = 4
                if ipoc == 4:
                    cpoc = 3
                    if PegTouch2DArray_del_2cycle[i][j+3] != PegTouch2DArray_del_2cycle[i][j+3]:
                        cpoc = 1
                if iprc == 4:
                    if PegTouch2DArray_del_2cycle[i][j-1] != PegTouch2DArray_del_2cycle[i][j-1]:
                        cprc = 3
                ## target pegの2(1)歩前（後）がnanの場合4(3)歩前（後）を参照するようにする。ipsi(contra)
                ## 基本的には反対側の足の１歩の中で、２歩ついていた場合、遅い方の足を採用する
                ## よってpreではiprc=4の時cprc=1がデフォルトだが、ipoc=4の場合ではcpoc=3がデフォルトとなる

                if (low < (PegTouch2DArray_del_2cycle[i][j+ipoc] - PegTouch2DArray_del_2cycle[i][j+cpoc]) < high) & (
                        low < (PegTouch2DArray_del_2cycle[i][j+cpoc] - PegTouch2DArray_del_2cycle[i][j]) < high) & (
                        low < (PegTouch2DArray_del_2cycle[i][j] - PegTouch2DArray_del_2cycle[i][j-cprc]) < high) & (
                        low < (PegTouch2DArray_del_2cycle[i][j-cprc] - PegTouch2DArray_del_2cycle[i][j-iprc]) < high):
                    Phasepost = (PegTouch2DArray_del_2cycle[i][j+cpoc] - PegTouch2DArray_del_2cycle[i][j]) / (
                            PegTouch2DArray_del_2cycle[i][j+ipoc] - PegTouch2DArray_del_2cycle[i][j])
                    Phasepre = (PegTouch2DArray_del_2cycle[i][j-cprc] - PegTouch2DArray_del_2cycle[i][j-iprc]) / (
                            PegTouch2DArray_del_2cycle[i][j] - PegTouch2DArray_del_2cycle[i][j-iprc])
                    Phaseshift = min(abs(Phasepost - Phasepre), abs(Phasepre - Phasepost))
                    whole_phaseL.append(Phasepre)
                    whole_phaseR.append(Phasepost)
                    if Phaseshift < 0.7:
                        PhAve[j%length] = PhAve[j%length] + Phaseshift
                        E[j%length] = E[j%length] + 1       ## jの値を 5～length+5 で動かしてるので j%length にする
                        whole_phaseshift.append(Phaseshift)
                        continue

        ####ここまでPhaseShift


        for i in range (length):
            if E[i] > (turnnum / 4):
                PhAve[i] = abs (PhAve[i] / E[i])
            else:
                PhAve[i] = 0
        PhAve24 = PhAve
        for i in range (len (delete)):
             PhAve24 = np.insert (PhAve24, delete[i], 0)  ##24行配列に直すために解析しなかったペグに0を入れた

########################################################################################################################
        # ここからinterval shift

        IntAve = np.zeros (length)
        F = np.zeros (length)
        for i in range (1, turnnum-1):
            for j in range (5, length+5):
                if WaterOnArray[int (PegTouch2DArray_del_2cycle[i][j])] == 0:
                    continue
                if (low < (PegTouch2DArray_del_2cycle[i,j] - PegTouch2DArray_del_2cycle[i,j-2]) < Inthigh) & (
                        low < (PegTouch2DArray_del_2cycle[i,j+2] - PegTouch2DArray_del_2cycle[i,j]) < Inthigh):
                    Intervalpost = (PegTouch2DArray_del_2cycle[i,j+2] - PegTouch2DArray_del_2cycle[i,j])
                    whole_interval.append(Intervalpost)
                    Intervalpre = (PegTouch2DArray_del_2cycle[i,j] - PegTouch2DArray_del_2cycle[i,j-2])
                    Intervalshift = abs (Intervalpost - Intervalpre)
                    IntAve[j%length] = IntAve[j%length] + Intervalshift
                    F[j%length] = F[j%length] + 1
                    whole_intervalshift.append(Intervalshift)
                    

        for i in range (length):
            if F[i] > (turnnum / 4):
                IntAve[i] = abs (IntAve[i] / F[i])
            else:
                IntAve[i] = 0
        IntAve24 = IntAve
        for i in range (len (delete)):
            IntAve24 = np.insert (IntAve24, delete[i], 0)
            

        Block_start = []
        Block_end = []
        for i in range (1, length - 1):
            for j in range (1, length - 1):
                if (Var[i, j] == 0)or(Var[i, j] != Var[i, j]):
                    Var[i, j] = (Var[i - 1, j] + Var[i + 1, j] + Var[
                        i, j - 1] + Var[i, j + 1]) / 4
        avevar = np.mean(Var[np.where(Var>1)])
        stdvar = np.std(Var[np.where(Var>1)])
        numchunk = np.zeros(length)
        counter = 0
        for i in range (length):
            for j in range (length):
                # if(Var[i, j]  < avevar-stdvar*0.2)and(Var[i, j] > 1):
                if(Var[i, j]  < avevar)and(Var[i, j] > 1):
                    counter = counter + 1
                else:
                    numchunk[i] = counter
                    counter = 0
                    break
        numchunk_expand = np.zeros(length*3)
        for i in range(length):
            numchunk_expand[i] = numchunk[i]
            numchunk_expand[length + i] = numchunk[i]
            numchunk_expand[length*2 + i] = numchunk[i]

        counter = np.zeros(length)  ##counter: Chunkに属したとみなされるペグに1を入れていく
        for i in range(length):
            for j in range(length):
                if numchunk[j] < 3:   ##3以下はそもそも考えない
                    continue
                if (j in Block_start) or (j in Block_end):  ##既にカウントされたペグは考えない
                    continue
                if (numchunk[j] >= length - i):  ##chunkの大きさをlength(最大)から順に下げていき、大きなChunkから採用できるようにする。
                    if len(Block_end) == 0:               ##配列が0なら取り敢えず入れる。
                        Block_start = np.append(Block_start, int(j))
                        Block_end = np.append(Block_end, int((j + numchunk[j])%length))  ##startからchunknumだけ後ろに行ったペグmod(length)の値で処理
                        for k in range(int(numchunk[j])):
                            counter[(k+j)%length] = 1     ## j番目のペグからChunkの中のペグに値１を入れる
                        continue
                    else:
                        counter1 = 0
                        for k in range(int(numchunk[j])):
                            if counter[(j+k)%length] == 1:
                                counter1 = counter1 + 1  ## これまでのchunkと今回のchunkが被るとcounter1には１が足されていく
                        if counter1 > 0:                 ## counter1が０以上なら被り　→　棄却
                            continue
                        else:
                            Block_start = np.append(Block_start, int(j))
                            Block_end = np.append(Block_end, int(j))
                            for k in range(int(numchunk[j])):
                                counter[(k+j)%length] = 1  ## 採用したChunkに１を入れていく
        print(Block_start)
        print(Block_end)
        
        for i in range (1, length - 1):
            if Covar[i] == 0:
                Covar[i] = (Covar[i - 1] + Covar[i + 1]) / 2
        if Covar[0] == 0:
            Covar[0] = (Covar[length - 1] + Covar[1]) / 2
        if Covar[length - 1] == 0:
            Covar[length - 1] = (Covar[length - 2] + Covar[0]) / 2
        if len (Covar[np.nonzero (Covar)]) > 2:
            for i in range (11):
                for j in Block_start:
                    j = int(j)
                    if j + i - 5 > length - 1:
                        k = j + i - 4 - length
                        if Covar[k] != 0:
                            Cost[i].append (Covar[k])
                        continue
                    elif j + i - 5 < 0:
                        k = j + i - 6 + length
                        if Covar[k] != 0:
                            Cost[i].append (Covar[k])
                        continue
                    else:
                        k = j + i - 5
                        if Covar[k] != 0:
                            Cost[i].append (Covar[k])
        if len (Covar[np.nonzero (Covar)]) > 5:
            for i in range (11):
                for j in Block_end:
                    j = int(j)
                    if j + i - 5 > length - 1:
                        k = j - 1 + i - 4 - length
                        if Covar[k] != 0:
                            Coend[i].append (Covar[k])
                        continue
                    elif j + i - 5 < 0:
                        k = j - 1 + i - 6 + length
                        if Covar[k] != 0:
                            Coend[i].append (Covar[k])
                        continue
                    else:
                        k = j + i - 5
                        if Covar[k] != 0:
                            Coend[i].append (Covar[k])

        if len (PhAve[np.nonzero(PhAve)]) > 2:
            for i in range (11):
                for j in Block_start:
                    j = int(j)
                    if j + i - 5 > length - 1:
                        k = j + i - 4 - length
                        if PhAve[k] != 0:
                            Phst[i].append (PhAve[k])
                        continue
                    elif j + i - 5 < 0:
                        k = j + i - 6 + length
                        if PhAve[k] != 0:
                            Phst[i].append (PhAve[k])
                        continue
                    else:
                        k = j + i - 5
                        if PhAve[k] != 0:
                            Phst[i].append (PhAve[k])
        if len (PhAve[np.nonzero (PhAve)]) > 5:
            for i in range (11):
                for j in Block_end:
                    j = int(j)
                    if j + i - 5 > length - 1:
                        k = j - 1 + i - 4 - length
                        if PhAve[k] != 0:
                            Phend[i].append (PhAve[k])
                        continue
                    elif j + i - 5 < 0:
                        k = j - 1 + i - 6 + length
                        if PhAve[k] != 0:
                            Phend[i].append (PhAve[k])
                        continue
                    else:
                        k = j + i - 5
                        if PhAve[k] != 0:
                            Phend[i].append (PhAve[k])

        for i in range (1, length - 1):
            if IntAve[i] == 0:
                IntAve[i] = (IntAve[i - 1] + IntAve[i + 1]) / 2
        if IntAve[0] == 0:
            IntAve[0] = (IntAve[length - 1] + IntAve[1]) / 2
        if IntAve[length - 1] == 0:
            IntAve[length - 1] = (IntAve[length - 2] + IntAve[0]) / 2
        for i in range (15):
            for j in Block_start:
                j = int(j)
                if j + i - 7 > length - 1:
                    k = j + i - 6 - length
                    if IntAve[k] != 0:
                        Intst[i].append (IntAve[k])
                    continue
                elif j + i - 7 < 0:
                    k = j + i - 8 + length
                    if IntAve[k] != 0:
                        Intst[i].append (IntAve[k])
                    continue
                else:
                    k = j + i - 7
                    if IntAve[k] != 0:
                        Intst[i].append (IntAve[k])

        for i in range (15):
            for j in Block_end:
                j = int(j)
                if j + i - 7 > length - 1:
                    k = j + i - 6 - length
                    if IntAve[k] != 0:
                        Intend[i].append (IntAve[k])
                    continue
                elif j + i - 7 < 0:
                    k = j + i - 8 + length
                    if IntAve[k] != 0:
                        Intend[i].append (IntAve[k])
                    continue
                else:
                    k = j + i - 7
                    if IntAve[k] != 0:
                        Intend[i].append (IntAve[k])

#         Intend_each_R = [[] for i in range (15)]
#         Intend_each_L = [[] for i in range (15)]
#         for i in range (15):
#             for j in Block_end:
#                 j = int(j)
#                 if j in Rpeg_number:
#                     if j + i - 7 > length - 1:
#                         k = j + i - 6 - length
#                         if IntAve[k] != 0:
#                             Intend_R[i].append (IntAve[k])
#                             Intend_each_R[i].append (IntAve[k])
#                         continue
#                     elif j + i - 7 < 0:
#                         k = j + i - 8 + length
#                         if IntAve[k] != 0:
#                             Intend_R[i].append (IntAve[k])
#                             Intend_each_R[i].append (IntAve[k])
#                         continue
#                     else:
#                         k = j + i - 7
#                         if IntAve[k] != 0:
#                             Intend_R[i].append (IntAve[k])
#                             Intend_each_R[i].append (IntAve[k])
#                 elif :
#                     if j + i - 7 > length - 1:
#                         k = j + i - 6 - length
#                         if IntAve[k] != 0:
#                             Intend_L[i].append (IntAve[k])
#                             Intend_each_L[i].append (IntAve[k])
#                         continue
#                     elif j + i - 7 < 0:
#                         k = j + i - 8 + length
#                         if IntAve[k] != 0:
#                             Intend_L[i].append (IntAve[k])
#                             Intend_each_L[i].append (IntAve[k])
#                         continue
#                     else:
#                         k = j + i - 7
#                         if IntAve[k] != 0:
#                             Intend_L[i].append (IntAve[k])
#                             Intend_each_L[i].append (IntAve[k])

        Intend_each_R = [[] for i in range (15)]
        Intend_each_L = [[] for i in range (15)]
        for i in range (15):
            for j in Block_end:
                j = int(j)
                if j in Rpeg_number:
                    if j + i - 7 > length - 1:
                        k = j + i - 6 - length
                        if IntAve[k] != 0:
                            Intend_R[i].append(IntAve[k])
                            Intend_each_R[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intend_L[i].append(IntAve[(k+1)%length])
                            Intend_each_L[i].append (IntAve[(k+1)%length])
                            
                        continue
                    elif j + i - 7 < 0:
                        k = j + i - 8 + length
                        if IntAve[k] != 0:
                            Intend_R[i].append (IntAve[k])
                            Intend_each_R[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intend_L[i].append (IntAve[(k+1)%length])
                            Intend_each_L[i].append (IntAve[(k+1)%length])
                        continue
                    else:
                        k = j + i - 7
                        if IntAve[k] != 0:
                            Intend_R[i].append (IntAve[k])
                            Intend_each_R[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intend_L[i].append (IntAve[(k+1)%length])
                            Intend_each_L[i].append (IntAve[(k+1)%length])
                elif j in Lpeg_number:
                    if j + i - 7 > length - 1:
                        k = j + i - 6 - length
                        if IntAve[k] != 0:
                            Intend_L[i].append (IntAve[k])
                            Intend_each_L[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intend_R[i].append (IntAve[(k+1)%length])
                            Intend_each_R[i].append (IntAve[(k+1)%length])
                        continue
                    elif j + i - 7 < 0:
                        k = j + i - 8 + length
                        if IntAve[k] != 0:
                            Intend_L[i].append (IntAve[k])
                            Intend_each_L[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intend_R[i].append (IntAve[(k+1)%length])
                            Intend_each_R[i].append (IntAve[(k+1)%length])
                        continue
                    else:
                        k = j + i - 7
                        if IntAve[k] != 0:
                            Intend_L[i].append (IntAve[k])
                            Intend_each_L[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intend_R[i].append (IntAve[(k+1)%length])
                            Intend_each_R[i].append (IntAve[(k+1)%length])

#         Intst_each_R = [[] for i in range (15)]
#         Intst_each_L = [[] for i in range (15)]
#         for i in range (15):
#             for j in Block_start:
#                 j = int(j)
#                 if j in Rpeg_number:
#                     if j + i - 7 > length - 1:
#                         k = j + i - 6 - length
#                         if IntAve[k] != 0:
#                             Intst_R[i].append (IntAve[k])
#                             Intst_each_R[i].append (IntAve[k])
#                         continue
#                     elif j + i - 7 < 0:
#                         k = j + i - 8 + length
#                         if IntAve[k] != 0:
#                             Intst_R[i].append (IntAve[k])
#                             Intst_each_R[i].append (IntAve[k])
#                         continue
#                     else:
#                         k = j + i - 7
#                         if IntAve[k] != 0:
#                             Intst_R[i].append (IntAve[k])
#                             Intst_each_R[i].append (IntAve[k])
#                 else:
#                     if j + i - 7 > length - 1:
#                         k = j + i - 6 - length
#                         if IntAve[k] != 0:
#                             Intst_L[i].append (IntAve[k])
#                             Intst_each_L[i].append (IntAve[k])
#                         continue
#                     elif j + i - 7 < 0:
#                         k = j + i - 8 + length
#                         if IntAve[k] != 0:
#                             Intst_L[i].append (IntAve[k])
#                             Intst_each_L[i].append (IntAve[k])
#                         continue
#                     else:
#                         k = j + i - 7
#                         if IntAve[k] != 0:
#                             Intst_L[i].append (IntAve[k])
#                             Intst_each_L[i].append (IntAve[k])

        Intst_each_R = [[] for i in range (15)]
        Intst_each_L = [[] for i in range (15)]
        for i in range (15):
            for j in Block_start:
                j = int(j)
                if j in Rpeg_number:
                    if j + i - 7 > length - 1:
                        k = j + i - 6 - length
                        if IntAve[k] != 0:
                            Intst_R[i].append (IntAve[k])
                            Intst_each_R[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intst_L[i].append (IntAve[(k+1)%length])
                            Intst_each_L[i].append (IntAve[(k+1)%length])
                        continue
                    elif j + i - 7 < 0:
                        k = j + i - 8 + length
                        if IntAve[k] != 0:
                            Intst_R[i].append (IntAve[k])
                            Intst_each_R[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intst_L[i].append (IntAve[(k+1)%length])
                            Intst_each_L[i].append (IntAve[(k+1)%length])
                        continue
                    else:
                        k = j + i - 7
                        if IntAve[k] != 0:
                            Intst_R[i].append (IntAve[k])
                            Intst_each_R[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intst_L[i].append (IntAve[(k+1)%length])
                            Intst_each_L[i].append (IntAve[(k+1)%length])
                elif j in Lpeg_number:
                    if j + i - 7 > length - 1:
                        k = j + i - 6 - length
                        if IntAve[k] != 0:
                            Intst_L[i].append (IntAve[k])
                            Intst_each_L[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intst_R[i].append (IntAve[(k+1)%length])
                            Intst_each_R[i].append (IntAve[(k+1)%length])
                        continue
                    elif j + i - 7 < 0:
                        k = j + i - 8 + length
                        if IntAve[k] != 0:
                            Intst_L[i].append (IntAve[k])
                            Intst_each_L[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intst_R[i].append (IntAve[(k+1)%length])
                            Intst_each_R[i].append (IntAve[(k+1)%length])
                        continue
                    else:
                        k = j + i - 7
                        if IntAve[k] != 0:
                            Intst_L[i].append (IntAve[k])
                            Intst_each_L[i].append (IntAve[k])
                        if IntAve[(k+1)%length] != 0:
                            Intst_R[i].append (IntAve[(k+1)%length])
                            Intst_each_R[i].append (IntAve[(k+1)%length])

        for i in range (11):
            for j in Block_end:
                j = int(j)
                if j in Rpeg_number:
                    if j + i - 5 > length - 1:
                        k = j + i - 4 - length
                        if PhAve[k] != 0:
                            Phend_R[i].append (PhAve[k])
                        continue
                    elif j + i - 5 < 0:
                        k = j + i - 6 + length
                        if PhAve[k] != 0:
                            Phend_R[i].append (PhAve[k])
                        continue
                    else:
                        k = j + i - 5
                        if PhAve[k] != 0:
                            Phend_R[i].append (PhAve[k])
                else:
                    if j + i - 5 > length - 1:
                        k = j + i - 4 - length
                        if PhAve[k] != 0:
                            Phend_L[i].append (PhAve[k])
                        continue
                    elif j + i - 5 < 0:
                        k = j + i - 6 + length
                        if PhAve[k] != 0:
                            Phend_L[i].append (PhAve[k])
                        continue
                    else:
                        k = j + i - 5
                        if PhAve[k] != 0:
                            Phend_L[i].append (PhAve[k])

        Intend_ave_eachR = np.zeros (5)
        Intend_ave_eachL = np.zeros (5)
        for i in range (5):
            Intend_ave_eachR[i] = np.nanmean (Intend_each_R[i * 2 + 3])
            Intend_ave_eachL[i] = np.nanmean (Intend_each_L[i * 2 + 3])

        if (max (Intend_ave_eachR[1:4]) - min (Intend_ave_eachR)) < (
                max (Intend_ave_eachL[1:4]) - min (Intend_ave_eachL)):
            for i in range (15):
                Intend_ideal[i] = np.append (Intend_ideal[i], Intend_each_L[i])
        else:
            for i in range (15):
                Intend_ideal[i] = np.append (Intend_ideal[i], Intend_each_R[i])
        
        Intend_R_modulation = max(Intend_ave_eachR[1:4])/min(Intend_ave_eachR[0], Intend_ave_eachR[4])
        Intend_L_modulation = max(Intend_ave_eachL[1:4])/min(Intend_ave_eachL[0], Intend_ave_eachL[4])
        modulation[number].append(Intend_R_modulation)
        modulation[number].append(Intend_L_modulation)

        Intst_ave_eachR = np.zeros (5)
        Intst_ave_eachL = np.zeros (5)
        for i in range (5):
            Intst_ave_eachR[i] = np.nanmean (Intst_each_R[i * 2 + 3])
            Intst_ave_eachL[i] = np.nanmean (Intst_each_L[i * 2 + 3])

        if (max (Intst_ave_eachR[1:4]) - min (Intst_ave_eachR)) < (
                max (Intst_ave_eachL[1:4]) - min (Intst_ave_eachL)):
            for i in range (15):
                Intst_ideal[i] = np.append (Intst_ideal[i], Intst_each_L[i])
        else:
            for i in range (15):
                Intst_ideal[i] = np.append (Intst_ideal[i], Intst_each_R[i])

        Intst_R_modulation = max(Intst_ave_eachR[1:4])/min(Intst_ave_eachR[0], Intst_ave_eachR[4])
        Intst_L_modulation = max(Intst_ave_eachL[1:4])/min(Intst_ave_eachL[0], Intst_ave_eachL[4])
        modulation[number].append(Intst_R_modulation)
        modulation[number].append(Intst_L_modulation)
        plt.figure()
        x = np.arange(1, length+1)
        y = PhAve
        plt.bar(x, y, color="k")
        plt.title("Phase Shift each trial")
        PATH = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/Figure'
        os.chdir (PATH)
        plt.savefig (str (number) + '_' + str (num) + '_Phase.png')

        Chunk = np.zeros(length) ## Chunkのinsideなら1, outsideなら0
        if len(Block_start) > 0:
            for i in range(len(Block_start)):
                if Block_start[i]<Block_end[i]:
                    for j in range(int(Block_end[i]-Block_start[i])):
                        Chunk[int(Block_start[i]+j)]=1
                else:
                    for j in range(int(Block_end[i]+length-Block_end[i])):
                        Chunk[int(Block_start[i]+j)%length]=1
            
            for i in range(length):
                if Chunk[i] == 1:
                    Intervalshift_in.append(IntAve[i])
                    Phaseshift_in.append(PhAve[i])
                else:
                    Intervalshift_out.append(IntAve[i])
                    Phaseshift_out.append(PhAve[i])

            for i in range(len(Block_start)):
                for j in range(3):
                    Intervalshift_out_boundary.append(IntAve[int(Block_start[i])-1-j])
                    Phaseshift_out_boundary.append(PhAve[int(Block_start[i])-1-j])
            for i in range(len(Block_end)):
                for j in range(3):
                    Intervalshift_out_boundary.append(IntAve[(int(Block_end[i])+1+j)%length])
                    Phaseshift_out_boundary.append(PhAve[(int(Block_start[i])+1+j)%length])
            
                

        
        
        ###それぞれのトライアルごとにInterval shiftを右足、左足に分けてStart/End Pegでアラインした。
        # PATH = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/Figure'
        # os.chdir (PATH)

        # plt.figure()
        # y = IntAve
        # plt.bar(x, y, color="k")
        # plt.title ("Interval Shift each trial")
        # plt.savefig (str (number) + '_' + str (num) + '_Interval.png')
        
        # plt.figure()
        # y = Intst_ave_eachR
        # x = np.arange(-2,3)
        # plt.bar(x,y)
        # plt.title("start interval Right each trial")
        # plt.savefig (str (number) + '_' + str (num) + '_Intst_R.png')
        
        # plt.figure()
        # y = Intst_ave_eachL
        # x = np.arange(-2,3)
        # plt.bar(x,y)
        # plt.title("start interval Left each trial")
        # plt.savefig (str (number) + '_' + str (num) + '_Intst_L.png')
        
        # plt.figure()
        # y = Intend_ave_eachR
        # x = np.arange(-2,3)
        # plt.bar(x,y)
        # plt.title("end interval Right each trial")
        # plt.savefig (str (number) + '_' + str (num) + '_Intend_R.png')
        
        # plt.figure()
        # y = Intend_ave_eachL
        # x = np.arange(-2,3)
        # plt.bar(x,y)
        # plt.title("end interval Left each trial")
        # plt.savefig (str (number) + '_' + str (num) + '_Intend_L.png')

        for i in range(5):
            if not np.isnan(Intst_ave_eachR[i]):
                Intst_R_mouse[i].append(Intst_ave_eachR[i])
            if not np.isnan(Intst_ave_eachL[i]):
                Intst_L_mouse[i].append(Intst_ave_eachL[i])
            if not np.isnan(Intend_ave_eachR[i]):
                Intend_R_mouse[i].append(Intend_ave_eachR[i])
            if not np.isnan(Intend_ave_eachL[i]):
                Intend_L_mouse[i].append(Intend_ave_eachL[i])


        HMwin = np.arange (0.5, 1.5, 0.02)
        HMhist = np.zeros ((length, length, len (HMwin)))
        for i in range (turnnum - 2):
            for j in range (length-1):
                for k in range (j + 1, length):
                    if (WaterOnArray[int (PegTouch2DArray_del[i, k])] == 0) or (
                            WaterOnArray[int (PegTouch2DArray_del[i, j])] == 0):
                        continue
                    elif (0 < (PegTouch2DArray_del[i, k] - PegTouch2DArray_del[i, j])) & (
                            (PegTouch2DArray_del[i, k] - PegTouch2DArray_del[i, j]) < (
                            cut * 2 + (Peg2DArray_del[1, k] - Peg2DArray_del[1, j]))):
                        for l in range (len (HMwin) - 1):
                            if (Ave[j, k] == 0):
                                continue
                            else:
                                if (HMwin[l] <= (PegTouch2DArray_del[i, k] - PegTouch2DArray_del[i, j]) / Ave[j, k]
                                        < HMwin[l + 1]):
                                    HMhist[j, k, l] = HMhist[j, k, l] + 1
        plt.clf()
        plt.close()
    
    plt.figure()
    plt.title("start interval Right each mouse")
    x = np.arange(-2, 3)
    y = np.zeros(5)
    for i in range(5):
        y[i] = np.nanmean(Intst_R_mouse[i])
    sem_y = np.zeros(5)
    for i in range(5):
        sem_y[i] = np.nanstd(Intst_R_mouse[i])/[np.sqrt(len(Intst_R_mouse[i]))]
    plt.bar(x, y, yerr=sem_y)
    PATH = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/Figure'
    os.chdir (PATH)
    plt.savefig (str (number) + '_Intst_R_mouse.png')
    tukey_hsd(list("ABCDE"),Intst_R_mouse[0] ,Intst_R_mouse[1], Intst_R_mouse[2], Intst_R_mouse[3],Intst_R_mouse[4])
    

    
    plt.figure()
    plt.title("start interval Left each mouse")
    x = np.arange(-2, 3)
    y = np.zeros(5)
    for i in range(5):
        y[i] = np.nanmean(Intst_L_mouse[i])
    sem_y = np.zeros(5)
    for i in range(5):
        sem_y[i] = np.nanstd(Intst_L_mouse[i])/[np.sqrt(len(Intst_L_mouse[i]))]
    plt.bar(x, y, yerr=sem_y)
    PATH = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/Figure'
    os.chdir (PATH)
    plt.savefig (str (number) + '_Intst_L_mouse.png')
    tukey_hsd(list("ABCDE"),Intst_L_mouse[0] ,Intst_L_mouse[1], Intst_L_mouse[2], Intst_L_mouse[3],Intst_L_mouse[4])
    
    
                       
    plt.figure()
    plt.title("end interval Right each mouse")
    x = np.arange(-2, 3)
    y = np.zeros(5)
    for i in range(5):
        y[i] = np.nanmean(Intend_R_mouse[i])
    sem_y = np.zeros(5)
    for i in range(5):
        sem_y[i] = np.nanstd(Intend_R_mouse[i])/[np.sqrt(len(Intend_R_mouse[i]))]
    plt.bar(x, y, yerr=sem_y)
    PATH = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/Figure'
    os.chdir (PATH)
    plt.savefig (str (number) + '_Intend_R_mouse.png')
    tukey_hsd(list("ABCDE"),Intend_R_mouse[0] ,Intend_R_mouse[1], Intend_R_mouse[2], Intend_R_mouse[3],Intend_R_mouse[4])
    
    
                       
    plt.figure()
    plt.title("end interval Left each mouse")
    x = np.arange(-2, 3)
    y = np.zeros(5)
    for i in range(5):
        y[i] = np.nanmean(Intend_L_mouse[i])
    sem_y = np.zeros(5)
    for i in range(5):
        sem_y[i] = np.nanstd(Intend_L_mouse[i])/[np.sqrt(len(Intend_L_mouse[i]))]
    plt.bar(x, y, yerr=sem_y)
    PATH = '/Users/kouhiro1006/Desktop/解析用/wheel_analyze/Figure'
    os.chdir (PATH)
    plt.savefig (str (number) + '_Intend_L_mouse.png')
    tukey_hsd(list("ABCDE"),Intend_L_mouse[0] ,Intend_L_mouse[1], Intend_L_mouse[2], Intend_L_mouse[3],Intend_L_mouse[4])

    plt.clf()
    plt.close()
deleteline = [0,1,3,5,7,9,11,13]
Intst_delete = np.delete(Intst, deleteline, axis=0)
Intend_delete = np.delete(Intend, deleteline, axis=0)

os.chdir (PATH)

plt.figure()
plt.title("Start Peg Aligned Phase Shift")
x = np.arange(-4, 5)
Phst_ave = np.zeros(9)
Phst.insert(0, Phst[10])
for i in range(9):
    Phst_ave[i] = np.mean(Phst[i+1])
y = Phst_ave
sem_y = np.zeros(9)
for i in range(9):
    sem_y[i] = np.std(Phst[i+1])/[np.sqrt(len(Phst[i+1]))]
plt.bar(x, y,yerr=sem_y)
result = scipy.stats.f_oneway(Phst[1], Phst[2], Phst[3], Phst[4], Phst[5], Phst[6], Phst[7], Phst[8], Phst[9])
print("Phase shift startpeg", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDEFGHI"),Phst[1] ,Phst[2], Phst[3], Phst[4], Phst[5], Phst[6], Phst[7], Phst[8], Phst[9])
plt.savefig ('PhaseShift_Start.png')


plt.figure()
plt.title("Start Peg Aligned Interval Shift")
x = np.arange(-7, 8)
Intst_ave = np.zeros(15)
for i in range(15):
    Intst_ave[i] = np.mean(Intst[i])
y = Intst_ave
sem_y = np.zeros(15)
for i in range(15):
    sem_y[i] = np.std(Intst[i])/[np.sqrt(len(Intst[i]))]
plt.bar(x, y, yerr=sem_y)
plt.savefig ('IntervalShift_Start.png')


plt.figure()
plt.title("Start Peg Aligned Interval Shift(2Step)")
x = np.arange(-2, 3)
Intst_ave = np.zeros(5)
for i in range(5):
    Intst_ave[i] = np.mean(Intst[i*2+3])
y = Intst_ave
sem_y = np.zeros(5)
for i in range(5):
    sem_y[i] = np.std(Intst[i*2+3])/[np.sqrt(len(Intst[i*2+3]))]
plt.bar(x, y, yerr=sem_y)
result = scipy.stats.f_oneway(Intst[3], Intst[5], Intst[7], Intst[9], Intst[11])
print("Interval shift startpeg", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDE"), Intst[3], Intst[5], Intst[7], Intst[9], Intst[11])
plt.savefig ('IntervalShift_Start_2Step.png')

plt.figure()
plt.title("End Peg Aligned Phase Shift")
x = np.arange(-4, 5)
Phend_ave = np.zeros(9)
for i in range(9):
    Phend_ave[i] = np.mean(Phend[i+1])
y = Phend_ave
sem_y = np.zeros(9)
for i in range(9):
    sem_y[i] = np.std(Phend[i+1])/[np.sqrt(len(Phend[i+1]))]
plt.bar(x, y, yerr=sem_y)
result = scipy.stats.f_oneway(Phend[1], Phend[2], Phend[3], Phend[4], Phend[5], Phend[6], Phend[7], Phend[8], Phend[9])
#result = scipy.stats.f_oneway(Phend[1], Phend[2], Phend[3], Phend[4], Phend[5], Phend[6], Phend[7])
#result = scipy.stats.f_oneway(Phend[2], Phend[3], Phend[4], Phend[5], Phend[6])
print("Phase shift endpeg", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDEFGHI"),Phend[1], Phend[2], Phend[3], Phend[4], Phend[5], Phend[6], Phend[7], Phend[8], Phend[9])
    #tukey_hsd (list ("ABCDEFG"), Phend[1], Phend[2], Phend[3], Phend[4], Phend[5], Phend[6], Phend[7])
    #tukey_hsd (list ("ABCDE"), Phend[2], Phend[3], Phend[4], Phend[5], Phend[6])
plt.savefig ('PhaseShift_End.png')



plt.figure()
plt.title("End Peg Aligned Interval Shift")
x = np.arange(-7, 8)
Intend_ave = np.zeros(15)
for i in range(15):
    Intend_ave[i] = np.mean(Intend[i])
y = Intend_ave
sem_y = np.zeros(15)
for i in range(15):
    sem_y[i] = np.std(Intend[i])/[np.sqrt(len(Intend[i]))]
plt.bar(x, y, yerr=sem_y)
plt.savefig ('IntervalShift_End.png')


plt.figure()
plt.title("End Peg Aligned Interval Shift(2Step)")
x = np.arange(-2, 3)
Intend_ave = np.zeros(5)
for i in range(5):
    Intend_ave[i] = np.mean(Intend[i*2+3])
y = Intend_ave
sem_y = np.zeros(5)
for i in range(5):
    sem_y[i] = np.std(Intend[i*2+3])/[np.sqrt(len(Intend[i*2+3]))]
plt.bar(x, y, yerr=sem_y)
result = scipy.stats.f_oneway(Intend[3], Intend[5], Intend[7], Intend[9], Intend[11])
print("Interval shift endpeg", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDE"), Intend[3], Intend[5], Intend[7], Intend[9], Intend[11])
plt.savefig ('IntervalShift_End_2Step.png')


plt.figure()
plt.title("End Peg Aligned Interval Shift Right(2Step)")
x = np.arange(-2, 3)
Intend_ave_R = np.zeros(5)
for i in range(5):
    Intend_ave_R[i] = np.mean(Intend_R[i*2+3])
y = Intend_ave_R
sem_y = np.zeros(5)
for i in range(5):
    sem_y[i] = np.std(Intend_R[i*2+3])/[np.sqrt(len(Intend_R[i*2+3]))]
plt.bar(x, y, yerr=sem_y)
result = scipy.stats.f_oneway(Intend_R[3], Intend_R[5], Intend_R[7], Intend_R[9], Intend_R[11])
print("Interval shift endpeg Right", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDE"), Intend_R[3], Intend_R[5], Intend_R[7], Intend_R[9], Intend_R[11])
plt.savefig ('IntervalShift_Start_2Step.png')


plt.figure()
plt.title("End Peg Aligned Interval Shift Left(2Step)")
x = np.arange(-2, 3)
Intend_ave_L = np.zeros(5)
for i in range(5):
    Intend_ave_L[i] = np.mean(Intend_L[i*2+3])
y = Intend_ave_L
sem_y = np.zeros(5)
for i in range(5):
    sem_y[i] = np.std(Intend_L[i*2+3])/[np.sqrt(len(Intend_L[i*2+3]))]
plt.bar(x, y, yerr=sem_y)
result = scipy.stats.f_oneway(Intend_L[3], Intend_L[5], Intend_L[7], Intend_L[9], Intend_L[11])
print("Interval shift endpeg Left", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDE"), Intend_L[3], Intend_L[5], Intend_L[7], Intend_L[9], Intend_L[11])
plt.savefig ('IntervalShift_End_2Step.png')


plt.figure()
plt.title("End Peg Aligned Phase Shift Right")
x = np.arange(-4, 5)
Phend_ave = np.zeros(9)
for i in range(9):
    Phend_ave[i] = np.mean(Phend_R[i+1])
y = Phend_ave
sem_y = np.zeros(9)
for i in range(9):
    sem_y[i] = np.std(Phend_R[i+1])/[np.sqrt(len(Phend_R[i+1]))]
plt.bar(x, y, yerr=sem_y)
result = scipy.stats.f_oneway(Phend_R[1], Phend_R[2], Phend_R[3], Phend_R[4], Phend_R[5], Phend_R[6], Phend_R[7], Phend_R[8], Phend_R[9])
print("Phase shift endpeg right", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDEFGHI"), Phend_R[1], Phend_R[2], Phend_R[3], Phend_R[4], Phend_R[5], Phend_R[6], Phend_R[7], Phend_R[8], Phend_R[9])
plt.savefig ('PhaseShift_End_R.png')


plt.figure()
plt.title("End Peg Aligned Phase Shift Left")
x = np.arange(-4, 5)
Phend_ave = np.zeros(9)
for i in range(9):
    Phend_ave[i] = np.mean(Phend_L[i+1])
y = Phend_ave
sem_y = np.zeros(9)
for i in range(9):
    sem_y[i] = np.std(Phend_L[i+1])/[np.sqrt(len(Phend_L[i+1]))]
plt.bar(x, y, yerr=sem_y)
result = scipy.stats.f_oneway(Phend_L[1], Phend_L[2], Phend_L[3], Phend_L[4], Phend_L[5], Phend_L[6], Phend_L[7], Phend_L[8], Phend_L[9])
print("Phase shift endpeg Left", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDEFGHI"), Phend_L[1], Phend_L[2], Phend_L[3], Phend_L[4], Phend_L[5], Phend_L[6], Phend_L[7], Phend_L[8], Phend_L[9])
plt.savefig ('PhaseShift_End_L.png')


plt.figure()
plt.title("End Peg Aligned Interval Shift ideal")
x = np.arange(-2, 3)
Intend_ave_ideal = np.zeros(5)
for i in range(5):
    Intend_ave_ideal[i] = np.mean(Intend_ideal[i*2+3])
y = Intend_ave_ideal
sem_y = np.zeros(5)
for i in range(5):
    sem_y[i] = np.std(Intend_ideal[i*2+3])/[np.sqrt(len(Intend_ideal[i*2+3]))]
plt.bar(x, y, yerr=sem_y)
result = scipy.stats.f_oneway(Intend_ideal[3], Intend_ideal[5], Intend_ideal[7], Intend_ideal[9], Intend_ideal[11])
print("Interval shift endpeg ideal", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDE"), Intend_ideal[3], Intend_ideal[5], Intend_ideal[7], Intend_ideal[9], Intend_ideal[11])
plt.savefig ('IntervalShift_End_ideal.png')


plt.figure()
plt.title("Start Peg Aligned Interval Shift ideal")
x = np.arange(-2, 3)
Intst_ave_ideal = np.zeros(5)
for i in range(5):
    Intst_ave_ideal[i] = np.mean(Intst_ideal[i*2+3])
y = Intst_ave_ideal
sem_y = np.zeros(5)
for i in range(5):
    sem_y[i] = np.std(Intst_ideal[i*2+3])/[np.sqrt(len(Intst_ideal[i*2+3]))]
plt.bar(x, y, yerr=sem_y)
result = scipy.stats.f_oneway(Intst_ideal[3], Intst_ideal[5], Intst_ideal[7], Intst_ideal[9], Intst_ideal[11])
print("Interval shift startpeg ideal", result.pvalue)
if result.pvalue < 0.05:
    tukey_hsd(list("ABCDE"), Intst_ideal[3], Intst_ideal[5], Intst_ideal[7], Intst_ideal[9], Intst_ideal[11])
plt.savefig ('IntervalShift_Start_ideal.png')


plt.figure()
plt.title("Interval Histogram")
x = np.arange(0, 800, 8)
int_hist = np.zeros(100)
for i in range(len(whole_interval)):
    for j in range(100):
        if j*8 < whole_interval[i] < (j+1)*8:
            int_hist[j] = int_hist[j] + 1
for i in range(3):
    int_hist[i] = 0
for i in range(94):
    int_hist[i+3] = np.average(int_hist[i: i+6])
y = int_hist
plt.plot(x,y)
plt.savefig('Interval_histogram.png')


plt.figure()
x = np.arange(0, 200, 2)
intshift_hist = np.zeros(100)
for i in range(len(whole_intervalshift)):
    for j in range(100):
        if j*2 < whole_intervalshift[i] < (j+1)*2:
            intshift_hist[j] = intshift_hist[j] + 1
for i in range(3):
    intshift_hist[i] = 0
for i in range(94):
    intshift_hist[i+3] = np.average(intshift_hist[i: i+6])
y = intshift_hist
plt.plot(x,y)
plt.savefig ('IntervalShift_hisotgram.png')



plt.figure()
plt.title("R Phase Histogram")
x = np.arange(0, 1, 0.01)
ph_hist = np.zeros(100)
for i in range(len(whole_phaseR)):
    for j in range(100):
        if j*0.01 < whole_phaseR[i] < (j+1)*0.01:
            ph_hist[j] = ph_hist[j] + 1
for i in range(3):
    ph_hist[i] = 0
for i in range(94):
    ph_hist[i+3] = np.average(ph_hist[i: i+6])
y = ph_hist
plt.plot(x,y)
plt.savefig ('Phase_HIstogram_R.png')


plt.figure()
plt.title("L Phase Histogram")
x = np.arange(0, 1, 0.01)
ph_hist = np.zeros(100)
for i in range(len(whole_phaseL)):
    for j in range(100):
        if j*0.01 < whole_phaseL[i] < (j+1)*0.01:
            ph_hist[j] = ph_hist[j] + 1
for i in range(3):
    ph_hist[i] = 0
for i in range(94):
    ph_hist[i+3] = np.average(ph_hist[i: i+6])
y = ph_hist
plt.plot(x,y)
plt.savefig ('Phase_Histogram_L.png')


whole_phaseR = np.array(whole_phaseR)
whole_phaseL = np.array(whole_phaseL)
whole_phase = np.append(whole_phaseR, whole_phaseL)
plt.figure()
plt.title("Phase Histogram")
x = np.arange(0, 1, 0.01)
ph_hist = np.zeros(100)
for i in range(len(whole_phase)):
    for j in range(100):
        if j*0.01 < whole_phase[i] < (j+1)*0.01:
            ph_hist[j] = ph_hist[j] + 1
for i in range(3):
    ph_hist[i] = 0
for i in range(94):
    ph_hist[i+3] = np.average(ph_hist[i: i+6])
y = ph_hist
plt.plot(x,y)
plt.savefig ('Phase_Histogram.png')


plt.figure()
x = np.arange(0, 1, 0.01)
phshift_hist = np.zeros(100)
for i in range(len(whole_phaseshift)):
    for j in range(100):
        if j*0.01 < whole_phaseshift[i] < (j+1)*0.01:
            phshift_hist[j] = phshift_hist[j] + 1
for i in range(3):
    phshift_hist[i] = 0
for i in range(94):
    phshift_hist[i+3] = np.average(phshift_hist[i: i+6])
y = phshift_hist
plt.plot(x,y)
plt.savefig ('PhaseShift_Histogram.png')


print("modulation of R/L Aligned Interval shift")
print(modulation)


plt.figure()
x = ["Chunk_in","Chunk_out"]
intervalshift_in_average = np.mean(Intervalshift_in)
intervalshift_out_average = np.mean(Intervalshift_out)
in_error = np.std(Intervalshift_in)/np.sqrt(len(Intervalshift_in))
out_error = np.std(Intervalshift_out)/np.sqrt(len(Intervalshift_out))
y = [intervalshift_in_average, intervalshift_out_average]
sem_y = [in_error,out_error]
plt.title("interval Shift Chunkin/out")
plt.ylim(120,160)
plt.bar(x,y,yerr=sem_y)
_, pvalue = stats.ttest_ind(Intervalshift_in, Intervalshift_out)
print(pvalue)


plt.figure()
x = ["Chunk_in","Chunk_out"]
phaseshift_in_average = np.mean(Phaseshift_in)
phaseshift_out_average = np.mean(Phaseshift_out)
in_error = np.std(Phaseshift_in)/np.sqrt(len(Phaseshift_in))
out_error = np.std(Phaseshift_out)/np.sqrt(len(Phaseshift_out))
y = [phaseshift_out_average, phaseshift_in_average]
sem_y = [in_error, out_error]
plt.title("pahse shift Chunkin/out")
plt.bar(x,y,yerr=sem_y)
_, pvalue = stats.ttest_ind(Phaseshift_in, Phaseshift_out)
print(pvalue)


plt.figure()
x = ["Boundary_in","Boundary_out"]
intervalshift_in_average = np.mean(Intervalshift_in)
intervalshift_out_average = np.mean(Intervalshift_out_boundary)
in_error = np.std(Intervalshift_in)/np.sqrt(len(Intervalshift_in))
out_error = np.std(Intervalshift_out_boundary)/np.sqrt(len(Intervalshift_out_boundary))
y = [intervalshift_in_average, intervalshift_out_average]
sem_y = [in_error,out_error]
plt.title("interval Shift Boundaryin/out")
plt.ylim(120,160)
plt.bar(x,y,yerr=sem_y)
_, pvalue = stats.ttest_ind(Intervalshift_in, Intervalshift_out_boundary)
print(pvalue)


plt.figure()
x = ["Boundary_in","Boundary_out"]
phaseshift_in_average = np.mean(Phaseshift_in)
phaseshift_out_average = np.mean(Phaseshift_out_boundary)
in_error = np.std(Phaseshift_in)/np.sqrt(len(Phaseshift_in))
out_error = np.std(Phaseshift_out_boundary)/np.sqrt(len(Phaseshift_out_boundary))
y = [phaseshift_out_average, phaseshift_in_average]
sem_y = [in_error, out_error]
plt.title("pahse shift Boundaryin/out")
plt.bar(x,y,yerr=sem_y)
_, pvalue = stats.ttest_ind(Phaseshift_in, Phaseshift_out_boundary)
print(pvalue)

plt.show()