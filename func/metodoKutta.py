from init import initValue
import numpy as np
import matplotlib.pyplot as grafico
#------------------------------------------------------------------------
# from prettytable import PrettyTable
#------------------------------------------------------------------------
# Biblioteca da Tabela
#------------------------------------------------------------------------
o1 = [initValue.o1i]
o2 = [initValue.o2i]
w1 = [initValue.w1i]
w2 = [initValue.w2i]
delt = initValue.deltT
def systemPenduloEXE():
    #--------------------------------------------------
    # inicializa os valores T (tempo) e E(energia)
    #--------------------------------------------------
    t = [initValue.initT]
    totalEnerg = [initValue.energy(o1[0],o2[0],w1[0],w2[0])]
    #--------------------------------------------------
    # Inicializa a tabela com todos os dados do problema
    #--------------------------------------------------
    # tab = PrettyTable()
    #--------------------------------------------------
    # Constrói as colunadas da tabela
    #--------------------------------------------------
    # tab.field_names = ["T","\u03B81", "\u03B82", "\u03C91", "\u03C92","\u03B5"]
    #-----------------------------------------------------------------------
    # Metodo de Runge-Kutta clássico 3/8
    #-----------------------------------------------------------------------
    for i in range(initValue.itere):
    #--------------------------------------------------
    # Primeiro calulamos o O(t,O1,O2,W1,W2)
    #--------------------------------------------------
        theta = contRungeKutta(t[i],o1[i],o2[i],w1[i],w2[i])
    #--------------------------------------------------
    # O que o metodo retorna é um array O's que serão usados
    # Pra dar valores a cada termo do sitema, sendo calcualdos
    # toda a espressão : 
    # O(t,O1,O2,W1,W2) = [1/8*(k1+ 3k2 + 3k3 +k4)]
    # dentro do método
    #--------------------------------------------------
        o1.append(o1[i] + theta[0]*delt)
        o2.append(o2[i] + theta[1]*delt)
        w1.append(w1[i] + theta[2]*delt)
        w2.append(w2[i] + theta[3]*delt)
        t.append(t[i] + delt)
        totalEnerg.append(initValue.energy(o1[i+1],o2[i+1],w1[i+1],w2[i+1]))
    #------------------------------------------------------------------------
    # Adiciona os valores na tabela(Faz parte do loop for())
    #------------------------------------------------------------------------
        # tab.add_row([round(t[i],2),o1[i],o2[i],w1[i],w2[i],totalEnerg[i]])
    #------------------------------------------------------------------------
    # Impressão da Tabela
    #------------------------------------------------------------------------
    # print(tab)
    #------------------------------------------------------------------------
    # Calculo dos Valores Xi e Yi para fazer o grafico xy das posiçoes
    #------------------------------------------------------------------------
    x1 = initValue.l1*np.sin(o1)
    y1 = -initValue.l1*np.cos(o1)
    x2 = x1 + initValue.l2*np.sin(o2)
    y2 = y1 - initValue.l2*np.cos(o2)
    #--------------------------------------------------
    # Plotagem dos graficos do Problema
    #--------------------------------------------------
    #--------------------------------------------------
    # Movimento do Pêndulo Duplo
    #--------------------------------------------------
    plotarGrafXY(x1,y1,x2,y2)
    #--------------------------------------------------
    # theta i em função do tempo
    #--------------------------------------------------
    plotfuncaoAngulo(t,o1,1)
    plotfuncaoAngulo(t,o2,2)
    #--------------------------------------------------
    # Espaço de Fase
    #--------------------------------------------------
    plotFase(o1,w1,1)
    plotFase(o2,w2,2)
    #--------------------------------------------------
    # Energia
    #--------------------------------------------------
    plotEnergy(t,totalEnerg)
    #--------------------------------------------------
    # Fim do programa
    #--------------------------------------------------


#------------------------------------------------------------------------
#   O(t,O1,O2,W1,W2) = [1/8*(k1+ 3k2 + 3k3 +k4)]
#------------------------------------------------------------------------
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
#------------------------------------------------------------------------

def k(t,o1,o2,w1,w2):
    # ----------------------------------------------
        return initValue.fO1(w1
        ),initValue.fO2(w2
        ),initValue.fW1(o1,o2,w1,w2
        ),initValue.fW2(o1,o2,w1,w2) 
#------------------------------------------------------------------------
def plotarGrafXY(x1,y1,x2,y2):
    grafico.ylabel('Y')
    grafico.xlabel('X')
    grafico.title("Movimento do Pêndulo Duplo")
    grafico.plot(x1,y1,label='1')
    grafico.plot(x2,y2,label='2')
    grafico.legend()
    grafico.show()
#------------------------------------------------------------------------
def plotfuncaoAngulo(x,y,i):
    grafico.ylabel(' \u03B8'+ str(i))
    grafico.xlabel('Tempo')
    grafico.title("\u03B8"+ str(i) +" em função do tempo")
    grafico.plot(x,y)
    grafico.show()
#------------------------------------------------------------------------
def plotFase(x,y,i):
    grafico.ylabel('\u03C9'+ str(i))
    grafico.xlabel('\u03B8'+ str(i))
    grafico.title("Espaço de Fase")
    grafico.plot(x,y)
    grafico.show()
#------------------------------------------------------------------------
def plotEnergy(t,totalEnerg):
    grafico.ylabel('\u03B5')
    grafico.xlabel('Tempo')
    grafico.title("Energia Total")
    grafico.axis(ymin=-60,ymax=60)
    grafico.plot(t,totalEnerg)
    grafico.show()
#------------------------------------------------------------------------