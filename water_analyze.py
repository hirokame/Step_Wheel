import scipy.io as sci
import matplotlib.pyplot as plt
import os

mn = 20
mx = 200

for num in range(15):

    win = []
    z = 17
    for number in range(z):
        win.append(20 + number * 30)

    print(num)              ## print processing number
    PATH = 'C:\\Users\\hirokane\\PycharmProjects\\Wheel_Analyze\\wateranalysis_Agra'
    os.chdir(PATH)
    matdata = sci.loadmat(str(num+1))    #create file name
    peg = matdata["RLPegTouchAll"]
    drink = matdata["DrinkOnArray"]

    n = (len (peg)) - 1
    dr = (len (drink)) - 1
    s = list (range (n + 1))

    A = list (range (z))
    # B = list(range(len(win)-1))
    C = list (range (z))
    # D = list(range(len(win)-1))

    for i in range(z):
        A[i] = 0
        # B[i] = 0
        C[i] = 0
        # D[i] = 0

    for i in range(1, n):  # counting putative rpeg touch number
        s[i] = peg[i] - peg[i-1]
        if s[i] < 600:
            for k in range(0, z-1):
                counter = 0
                if win[k] < s[i] < win[k + 1]:  # filtering interval
                    for j in range(2, dr):
                        if ((peg[i-1]) < drink[j] < (peg[i]))and(mn < drink[j+2]-drink[j+1] < mx)and(mn < drink[j+1]-drink[j] < mx)and(mn < drink[j]-drink[j-1] < mx)and(mn < drink[j-1]-drink[j-2] < mx):
                            C[k] = C[k] + 1
                            print(drink[j])
                            counter = counter +1
                if 0.5 < counter:
                    A[k] = A[k] + 1
    print(A)
    print(C)
    print(win)

    for i in range(z):
        if not(A[i] == 0):
            C[i] = (C[i]/A[i])
            #D[i] = (D[i]/B[i])
        else:
            C[i] = 0
            #D[i] = 0

    win.pop(-1)
    C.pop(-1)
    plt.figure()
    x = win
    y = C
    plt.plot(x, y)
    plt.title("Drink_Histogram")
    plt.ylim(0, 6)
    plt.xlabel("interval(ms)")
    plt.ylabel("drink frequency (time/interval)")
    plt.show()


###### 5s Inter Peg Interval = about400ms with avarage 2 drink on.  So,10s IPI = about800ms with 4 drink on
###### therefore, need 5 samples each of 10s 8s 6s 5s 4s trials, for prediction drink on nember 4,3,2,2,2 of each trial.
###### マウス1匹づつ階段グラフを出したもの。