import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import math 
from scipy.integrate import odeint  

#Faisons rentrer les constantes
q=-1.6*(10**(-19))
#Masse de l'electron
m=9.1*10**(-31)
#Champ magnetique initial
B0=1.0 
#Pulsation cyclotron
Pc=abs(q)*(B0/m)
#Periode cyclotron
Tc=2*(math.pi)/Pc
#Vitesse initial de l'electron
v0=2.0*(10**8)
#distance initiale
d=0.02
#Nombre de pas
N=500
#Rayon de la lentille
a=0.001

#Partie 1 : Champ magnetique constant

#Question 1 
#Definissons une liste qui contient N+1 elements uniformement repartis entre 0 et 5TC

listeT=[i*(5*(Tc/N)) for i in range(N+1)]

#Question 2 : Methode EXACTE
#a
listeX_exacte=[]
listeY_exacte=[]
listeZ_exacte=[]

for t in listeT:
    x = (v0/math.sqrt(2))*(m/(q*B0))*math.sin((q*B0)*t/m)
    listeX_exacte.append(x)
    
    y = (v0/math.sqrt(2))*(m/(q*B0))*math.cos((q*B0)*t/m) - (v0*m)/(math.sqrt(2)*q*B0)
    listeY_exacte.append(y)
    
    z = -d + (t*v0)/math.sqrt(2)
    listeZ_exacte.append(z)

#b Trajectoire
fig = plt.figure(1)
fig.suptitle('Trace des differentes trajectoires', fontsize=20, fontweight ='bold') #anticipation qst 5
plt.gca(projection= '3d')
plt.plot(listeX_exacte,listeY_exacte,listeZ_exacte, 'y', label = "Methode EXACTE")
plt.legend()    

#Question 3 : Methode d EULER 
#a
listeX_euler=[0]
listeY_euler=[0]
listeZ_euler=[-d]

listeVX_euler=[v0/math.sqrt(2)]
listeVY_euler=[0]
listeVZ_euler=[v0/math.sqrt(2)]


for i in range (N):
    listeX_euler.append(listeX_euler[i] + (listeT[i+1] - listeT[i])*(listeVX_euler[i]))
    listeY_euler.append(listeY_euler[i] + (listeT[i+1] - listeT[i])*(listeVY_euler[i]))
    listeZ_euler.append(listeZ_euler[i] + (listeT[i+1] - listeT[i])*(listeVZ_euler[i]))
    
    listeVX_euler.append(listeVX_euler[i] + (listeT[i+1] - listeT[i])*((q*B0/m)*listeVY_euler[i]))
    listeVY_euler.append(listeVY_euler[i] + (listeT[i+1] - listeT[i])*((-q*B0/m)*listeVX_euler[i]))
    listeVZ_euler.append(listeVZ_euler[i])

#b Trajectoire 
plt.gca(projection='3d')
plt.plot(listeX_euler,listeY_euler,listeZ_euler, 'b', label ="Methode EULER")
plt.legend()    

#Question 4 : Methode ODEINT 
#a definition de la fonction F
def F(X,t):
    vx = X[3]
    vy = X[4]
    vz = X[5]
    fvx = (q*B0/m)*vy
    fvy=(-q*B0/m)*vx
    fvz = 0
    return [vx, vy, vz, fvx, fvy, fvz]  
      
#b Resolution de l equation avec odeint
    #Conditions initiales 

X0=[0, 0 , -d, v0/math.sqrt(2), 0, v0/math.sqrt(2) ]

    #Resolution 
X= odeint(F, X0, listeT)

listeX= [liste[0] for liste in X]
listeY= [liste[1] for liste in X]
listeZ= [liste[2] for liste in X]

#c Trace de la trajectoire
ax = plt.gca(projection='3d')
plt.plot(listeX,listeY,listeZ, 'r-.', label ="Methode ODEINT") 
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.legend()    
plt.show()

#Question 5 
# La solution calculee par ODEINT est celle qui se rapproche le plus de la solution exacte

#Partie 2: Une lentille magnetique 
#Question 1

def Bx(x,y,z):
    return 3*B0* (x*z/(2*(a**2))) *(1+z**2/a**2)**(-5/2)
def By(x,y,z):
    return 3*B0* (y*z/(2*(a**2))) *(1+z**2/a**2)**(-5/2)
def Bz(x,y,z):
    return B0*(1+z**2/a**2)**(-3/2)

#Question 2

B=[Bx, By, Bz]  

#Definition de la fonction F 
def F(X,t):
    x = X[0]
    y = X[1]
    z = X[2]
    vx = X[3]
    vy = X[4]
    vz = X[5]
    V= [vx, vy, vz]
    vB = [Bx(x,y,z), By(x,y,z), Bz(x,y,z)]
    p = produit_vect(V,vB)
    fv = [(q/m)*pelement for pelement in p]
    [fvx, fvy, fvz]=fv
    
    return [vx, vy, vz, fvx, fvy, fvz]
    
#Definition du produit vectoriel 
def produit_vect(V,B):
    x = V[1]*B[2] - V[2]*B[1]
    y = V[2]*B[0] - V[0]*B[2]
    z = V[0]*B[1] - V[1]*B[0]
    return [x, y, z]
    
#Solution en utilisant la methode ODEINT 
alphaSolutions = []
for k in range(-10,11):
    alpha = k*math.pi/500
    X0=[0, 0, -d, v0*math.sin(alpha), 0, v0*math.cos(alpha)]
    X=odeint(F, X0, listeT)
    listeX= [liste[0] for liste in X]
    listeY= [liste[1] for liste in X]
    listeZ= [liste[2] for liste in X]
    
    alphaSolutions.append([listeX,listeY,listeZ])
    fig = plt.figure(2)
    fig.suptitle('Ensemble des trajectoires', fontsize=20, fontweight='bold')
    ax = plt.gca(projection='3d')
    plt.plot(listeX,listeY,listeZ) 
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    plt.show()

#Question 3
#a Calcul des z_alpha


def minimum(x): #fct qui retourne le minimum d une liste 
    mini = x[0]
    for i in x:
        if i < mini:
            mini = i
    return mini
    

def minimumZPositive(listeR, listeZ): #fct qui retourne z_alpha
    #
    newLZ = []
    newLR = []
    for i in range(len(listeZ)):
        if listeZ[i] >= 0:
            newLZ.append(listeZ[i])
            newLR.append(listeR[i])
            
    mini= minimum(newLR)
    indice = newLR.index(mini)
    
    return newLZ[indice]

#On cherche la liste des valeurs z alpha
listeZ_alpha=[]    
for solution in alphaSolutions:
    listeR=[]
    listeX = solution[0]
    listeY = solution[1]    
    listeZ = solution[2]
    for i in range(N+1):
        listeR.append((listeX[i])**2+(listeY[i])**2)
        
    z_alpha = minimumZPositive(listeR, listeZ)
    listeZ_alpha.append(z_alpha)
    
    
print ' Les valeurs z_alpha les plus proches de laxe sont:'
print listeZ_alpha   

#b Valeur moyennne des altitudes 
val_moy=np.mean(listeZ_alpha)
print  'la valeur moyenne des altitudes est: '"{}".format(val_moy)

#c : calcul de la distance focale 
v = 1/d +1/val_moy
f = 1/v
print f