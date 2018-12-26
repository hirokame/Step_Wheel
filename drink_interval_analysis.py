import scipy.io as sci
import matplotlib.pyplot as plt
import os

win = []
mn = 40
mx  =200
z = 15
for num in range(z):
    win.append(100 + num * 30)

drinkhist = list(range(z))
C = list(range(len(win)-1))

for num in range(15):
    print(num)              ## print processing number
    PATH = 'C:\\Users\\hirokane\\PycharmProjects\\Wheel_Analyze\\wateranalysis_Agra'
    os.chdir(PATH)
    matdata = sci.loadmat(str(num+1))    #create file name
    peg = matdata["RLPegTouchAll"]
    drink = matdata["DrinkOnArray"]

    n = (len(peg)) - 1
    dr = (len(drink)) - 1
    s = list(range(n + 1))

    for i in range (0, z - 1):
        C[i] = 0
        drinkhist[i] = 0

    for i in range(1, n):  # counting putative rpeg touch number
        s[i] = peg[i] - peg[i - 1]
        for k in range(z-1):
            if win[k] < s[i] < win[k + 1]:  # filtering interval
                for j in range(dr-1):
                    if ((peg[i-1]) < drink[j] < (peg[i])) and((peg[i]) < drink[j+1] < (peg[i+1])):
                        if (mn < drink[j+2]-drink[j+1] < mx)and(mn < drink[j+1]-drink[j] < mx)and(mn < drink[j]-drink[j-1] < mx)and(mn < drink[j-1]-drink[j-2] < mx):
                            C[k] = C[k] + 1
                            drinkhist[k] = drinkhist[k] + (drink[j+1]-drink[j])

    for i in range(z-1):
        if not(C[i] == 0):
            drinkhist[i] = (drinkhist[i]/C[i])
            #D[i] = (D[i]/B[i])
        else:
            drinkhist[i] = 0
            #D[i] = 0



    plt.figure()
    x = win
    y = drinkhist
    plt.plot(x, y)
    plt.title("Drink interval histogram")
    plt.xlabel("Running interval (ms)")
    plt.ylabel("Averaged drink interval (ms)")
    plt.ylim(0, 400)
    plt.show()

### 各歩行インターバルのビンごとのIDI(Inter Drink Interval)の平均
### 歩行のスピードが変わっても、階段の平行部分では安定したIDIで水を飲んでいることが示したい
###多分プログラムがミスっている