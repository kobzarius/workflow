# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 23:59:02 2019

@author: olegk
"""
from scipy import special
import os
import scipy as sp

import matplotlib.pyplot as plt
import numpy as np
import pylab


from bokeh.plotting import figure, output_file, show
from bokeh.io import output_file, show
from bokeh.layouts import widgetbox
from bokeh.models.widgets import Slider


def out(string1, result,string2=''): #выпод вычислений с комментариями
    space='  '
    print(string1 + space + str(result) + space + string2 + '\n')





#Инегральная экспонента и логарифмическая аппроксимация

Ei=sp.special.expi  #f=-Ei(-x), используется в уравнении решения фильтрации
E1=sp.special.exp1 #f=E1(x);   E1(x)=-Ei(-x)=-ln(x)-gammd
ln=np.log

def LogApprox(x): #логарифмическая аппроксимация, только для числа, не списка
    gamma=0.57
    #gamma=0.809
    #result=(-1)*math.log(x)-gamma
    result=(-1)*ln(x)-gamma
    
    return result


def genvaluesx(k,x0=0.001): #генерация значений x c множителем
    values=[float(x0)] #начальное значение
    for i in range(20): #число элементов 
        values.append(values[i]*k)
    return values


def save(name='', fmt='png'): #в текущей директории должна быть папка pictures
    pwd = os.getcwd()
    iPath = './pictures/{}'.format(fmt)
    if not os.path.exists(iPath):
        os.mkdir(iPath)
    os.chdir(iPath)
    plt.savefig('{}.{}'.format(name, fmt), fmt='png')
    os.chdir(pwd)
    #plt.close()

k_md=10
mu_sp=1
ct_1atm=5*10**(-5)
rw_m=0.1
q_m3day=10
B_m3m3=1.2
Pi_bar=250
m=0.2
h_m=1

def td(t,k=k_md,m=m,mu=mu_sp,ct=ct_1atm,rw=rw_m):
    result=0.00036*k*t/(m*mu*ct*rw*rw)
    #print(result)
    return result

def pd_E1(rd,td):
    result=1/2*E1(rd**2/4/td)
    #print(result)
    return result

def rd(r,rw=rw_m):
    result=r/rw
    #print(result)
    return result

def pwf(pd,k=k_md,h=h_m,q=q_m3day,b=B_m3m3,mu=mu_sp,pi=Pi_bar):
    result=pi-18.41*pd*q*b*mu/k/h
    #print(result)
    return result

def pd_ln(rd,td):
    #result=1/2*LogApprox(rd**2/4/td)
    result=1/2*(ln(td/rd**2)+0.80907)
    #print(result)
    return result
    
def p(r,t,mod=1):
    td1=td(t)
    rd1=rd(r)
    if mod==1: #решение линейного стока с помощью интегральной экспоненты
        pd1=pd_E1(rd1,td1)
        result=pwf(pd1)
        return result
    if mod==2:
        pd1=pd_ln(rd1,td1)
        result=pwf(pd1)
        return result       


    
td1=td(83.336111111)
rd1=rd(0.1)
x=0.001
k=1.4 #множитель по оси x
'''pd1=pd_E1(rd1,td1)
pd2=pd_ln(rd1,td1)
pwf(pd1)
pwf(pd2)'''

valuesx=genvaluesx(k)
valuesx=np.asanyarray(valuesx)
print(p(0.1,0.00001))
print(p(0.1,0.00001,2))
solvE1=p(0.1,valuesx)
solvLog=p(0.1,valuesx,2)

out('-Ei(-x)=',-Ei(-x), 'Пример значения интегральной экспоненты при x={}'.format(str(x)))
out('E1(x)=',E1(x), 'E1(x)=-Ei(-x)=-ln(x)-gamma при x={}'.format(str(x)))
out('-ln(x)-gamma=',LogApprox(x),'логарифмическая аппроксимация интегральной экспоненты при x={}'.format(str(x)))


#figsize = (4,3)
#fig=plt.figure(figsize=figsize, dpi=150, facecolor='green', frameon=True)
fig=plt.figure()
#plt.xscale('log') 
plt.plot(valuesx,solvE1, label='E1')
plt.plot(valuesx,solvLog, label='Ln')

plt.title('Интегральная экспонента и ее аппроксимация')   # заголовок
plt.xlabel('x')   # подпись оси OX
plt.ylabel('E1(x); -ln(x)-0.57')   # подпись оси OY
plt.legend()
plt.grid(True)

plt.show()

