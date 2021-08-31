# Bibliotecas necessárias
# pip install prettytable
# pip install matplotlib
# pip install sympy
# pip install numpy
import numpy as np
from sympy import *
# -------------------------------
# Inicialização das variaveis
# -------------------------------
# Comprimento das barra
# -------------------------------
l1 = 1.0
l2 = 1.0
# -------------------------------
# Massa das particulas
# -------------------------------
m1 = 1.0
m2 = 1.0
# -------------------------------
# Angulo 𝜃1 e 𝜃2 iniciais(recebendo
# em graus e convertendo em radianos)
# -------------------------------
o1i = np.radians(30)
o2i = np.radians(0)
# -------------------------------
#  velocidades angulares W1 e W2 
#  iniciais
# -------------------------------
w1i = 0.0
w2i = 0.0
# -------------------------------
# delta T e tempo inicial e final
# -------------------------------
deltT = 0.01
initT = 0.0
finalT = 50.0
itere = int(round(finalT/deltT,1))
# --------------------------------
a = (m1 + m2)*l1
d = m2*l2


def fO1(w1):
    return w1
def fO2(w2):
    return w2
def fW1(o1,o2,w1,w2):
    b = m2*l2*np.cos(o1 - o2)
    c = m2*l1*np.cos(o1 - o2)
    e = -m2*l2*(w2**2)*np.sin(o1 - o2) - 9.8*(m1 + m2)*np.sin(o1)
    f = m2*l1*(w1**2)*np.sin(o1 - o2) - m2*9.8*np.sin(o2)
    return (e*d -b*f)/(a*d - c*b)
def fW2(o1,o2,w1,w2):
    b = m2*l2*np.cos(o1 - o2)
    c = m2*l1*np.cos(o1 - o2)
    e = -m2*l2*(w2**2)*np.sin(o1 - o2) - 9.8*(m1 + m2)*np.sin(o1)
    f = m2*l1*(w1**2)*np.sin(o1 - o2) - m2*9.8*np.sin(o2)
    return (a*f - c*e)/(a*d - c*b)
def energy(o1,o2,w1,w2):
    P = -(m1 + m2)*9.8*l1*np.cos(o1) - m2*9.8*l2*np.cos(o2)
    C = (1.0/2.0)*m1*(l1**2)*(w1**2) + (1.0/2.0)*m2*((l1**2)*(w1**2) + (l2**2)*(w2**2) + 2*l1*l2*w1*w2*np.cos(o1 -o2))
    return P + C