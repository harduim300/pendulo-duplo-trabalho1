from init import initValue
import numpy as np
import matplotlib.pyplot as grafico
# from prettytable import PrettyTable
o1 = [initValue.o1i]
o2 = [initValue.o2i]
w1 = [initValue.w1i]
w2 = [initValue.w2i]
delt = initValue.deltT
def systemPenduloEXE():
    t = [initValue.initT]
    # tab = PrettyTable()
    totalEnerg = [initValue.energy(o1[0],o2[0],w1[0],w2[0])]
    # tab.field_names = ["T","O1", "O2", "W1", "W2","E"]
    for i in range(initValue.itere):
        theta = contRungeKutta(t[i],o1[i],o2[i],w1[i],w2[i])
        o1.append(o1[i] + theta[0]*delt)
        o2.append(o2[i] + theta[1]*delt)
        w1.append(w1[i] + theta[2]*delt)
        w2.append(w2[i] + theta[3]*delt)
        t.append(t[i] + delt)
        totalEnerg.append(initValue.energy(o1[i+1],o2[i+1],w1[i+1],w2[i+1]))
        # tab.add_row([round(t[i],2),o1[i],o2[i],w1[i],w2[i],totalEnerg[i]])
    # print(tab)
    x1 = initValue.l1*np.sin(o1)
    y1 = -initValue.l1*np.cos(o1)
    x2 = x1 + initValue.l2*np.sin(o2)
    y2 = y1 - initValue.l2*np.cos(o2)
    
    plotarGrafXY(x1,y1,x2,y2)
    plotfuncaoAngulo(t,o1)
    plotfuncaoAngulo(t,o2)
    plotFase(o1,w1)
    plotFase(o2,w2)
    plotEnergy(t,totalEnerg)

def contRungeKutta(t,o1,o2,w1,w2):
    # ----------------------------------------------
    k1 = np.array(k(t,o1,o2,w1,w2))
    # ----------------------------------------------
    k2 =  np.array(k(t +(1/3)*delt,
    o1 +(1.0/3.0)*k1[0]*delt,
    o2 +(1.0/3.0)*k1[1]*delt,
    w1 +(1.0/3.0)*k1[2]*delt,
    w2 +(1.0/3.0)*k1[3]*delt))
    # ----------------------------------------------
    k3 =  np.array(k(t +(2/3)*delt,
    o1 -(1.0/3.0)*k1[0]*delt + k2[0]*delt,
    o2 -(1.0/3.0)*k1[1]*delt + k2[1]*delt,
    w1 -(1.0/3.0)*k1[2]*delt + k2[2]*delt,
    w2 -(1.0/3.0)*k1[3]*delt + k2[3]*delt))
    # ----------------------------------------------
    k4 =  np.array(k(t + delt,
    o1 +k1[0]*delt -k2[0]*delt +k3[0]*delt,
    o2 +k1[1]*delt -k2[1]*delt +k3[1]*delt,
    w1 +k1[2]*delt -k2[2]*delt +k3[2]*delt,
    w2 +k1[3]*delt -k2[3]*delt +k3[3]*delt))
    
    # ----------------------------------------------
    return (1.0/8.0)*( k1 +  3*k2 +  3*k3 + k4 )
    # ----------------------------------------------

def k(t,o1,o2,w1,w2):
    # ----------------------------------------------
        return initValue.fO1(w1
        ),initValue.fO2(w2
        ),initValue.fW1(o1,o2,w1,w2
        ),initValue.fW2(o1,o2,w1,w2) 

def plotarGrafXY(x1,y1,x2,y2):
    grafico.ylabel('Y')
    grafico.xlabel('X')
    grafico.axis(ymin=-2,ymax=-0.8)
    grafico.axis(xmin=-1,xmax=1)
    grafico.title("Movimento do Pêndulo Duplo")
    grafico.plot(x1,y1,label='1')
    grafico.plot(x2,y2,label='2')
    grafico.legend()
    grafico.show()

def plotfuncaoAngulo(x,y):
    grafico.ylabel(' Oi / radianos')
    grafico.xlabel('Tempo')
    grafico.title("Oi em função do tempo")
    grafico.plot(x,y)
    grafico.show()

def plotFase(x,y):
    grafico.ylabel('Wi')
    grafico.xlabel('Oi')
    grafico.title("Espaço de Fase")
    grafico.plot(x,y)
    grafico.show()

def plotEnergy(t,totalEnerg):
    grafico.ylabel('E')
    grafico.xlabel('Tempo')
    grafico.title("Energia Total")
    grafico.axis(ymin=-60,ymax=60)
    grafico.plot(t,totalEnerg)
    grafico.show()