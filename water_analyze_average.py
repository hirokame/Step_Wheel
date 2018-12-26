import scipy.io as sci
import matplotlib.pyplot as plt
import os

win = []
z = 16
for num in range(z):
    win.append(300 + num * 20)

A = list(range(len(win)-1))
#B = list(range(len(win)-1))
C = list(range(len(win)-1))
#D = list(range(len(win)-1))
for i in range(0, z-1):
    A[i] = 0
    #B[i] = 0
    C[i] = 0
    #D[i] = 0


for num in range(15):
    print(num)              ## print processing number
    PATH = 'C:\\Users\\hirokane\\PycharmProjects\\Wheel_Analyze\\wateranalysis_Agra'
    os.chdir(PATH)
    matdata = sci.loadmat(str(num+1))    #create file name
    #print(matdata.keys())
    rpeg = matdata["RpegTimeArray"]
    lpeg = matdata["LpegTimeArray"]                   #read .mat data
    drink = matdata["DrinkOnArray"]
    #print(rpeg)
    #print(lpeg)
    #print(drink)


    n = (len(rpeg)) - 1
    m = (len(lpeg)) - 1
    dr = (len(drink)) - 1

    rs = list(range(n + 1))
    ls = list(range(m + 1))

    for i in range(1, n):  # counting putative rpeg touch number
        rs[i] = rpeg[i] - rpeg[i - 1]
        for k in range(0, z-1):
            if win[k] < rs[i] < win[k + 1]:  # filtering interval
                A[k] = A[k] + 1
                for j in range(dr):
                    if (rpeg[i-1]) < drink[j] < (rpeg[i]):
                        C[k] = C[k] + 1

    for i in range(1, m):  # counting putative lpeg touch number
        ls[i] = lpeg[i] - lpeg[i - 1]
        for k in range(0, z-1):
            if win[k] < ls[i] < win[k + 1]:  # filtering interval
                A[k] = A[k] + 1
                for j in range(dr):
                    if (lpeg[i-1]) < drink[j] < (lpeg[i]):
                        C[k] = C[k] + 1

    print(A)
#print(B)
    print(C)
#print(D)

for i in range(0, z-1):
    if not(A[i] == 0):
        C[i] = (C[i]/A[i])
        #D[i] = (D[i]/B[i])
    else:
        C[i] = 0
        #D[i] = 0

win.pop(-1)

plt.figure()
x = win
y = C
#y_L = D
plt.plot(x, y, marker="x")
#plt.plot(x, y_L, marker="x")
plt.title("Drink_Histogram")
plt.ylim(0, 5)
plt.show()


###### 5s Inter Peg Interval = about400ms with avarage 2 drink on.  So,10s IPI = about800ms with 4 drink on
###### therefore, need 5 samples each of 10s 8s 6s 5s 4s trials, for prediction drink on nember 4,3,2,2,2 of each trial
###### 水飲みが整数倍になっているか、全てのマウス、全ての日にちで平均したもの
###### マウスごとのtransitionの位置が違うのか、母数が多ければ多いほど線形になってしまう