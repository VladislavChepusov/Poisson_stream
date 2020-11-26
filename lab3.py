import scipy 
import random
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import re
import copy
from scipy.stats import chi2

T1 = 21
T2 = 121
L1 = 29/45
L2 = 30/46

del_t=(T2-T1)/25
#более красивый вывод
def mprint(che):
    for x in che: print(x)

 #Функция для отрисовки
def draw(mass1,mass2,mass3,mass4):
    dpi = 100
    fig = plt.figure(dpi = dpi, figsize = (1024 / dpi, 640 / dpi) )
    mpl.rcParams.update({'font.size': 10})
    plt.axis([20, 122, -0.1,0.1])#Координатная плоскость
    plt.title(f'График')
    plt.xlabel('x')
    plt.ylabel('Y')
    xs = []
    cos_vals = []
    x = 0
    while x < 122.0:
        cos_vals += [0]
        xs += [x]
        x += 1
    plt.plot(xs, cos_vals, color = 'black', linestyle = 'solid',label = '!!!!!!!!')
    #Точки для 1 массива пока что
    size1=len(mass1)
    size2=len(mass2)
    size3=len(mass3)
    size4=len(mass4)
    #Построение точек
    for i in range(size1):
        plt.scatter(mass1[i], 0.05, color ='red', s = 15, marker = '*')
        plt.vlines(mass1[i], -1, 1,color = 'b',linewidth = 0.3,linestyle = '--')#Вертикальные линии
    for i in range(size2):
        plt.scatter(mass2[i], 0.025 , color ='red', s = 15, marker = '*')
        plt.vlines(mass2[i], -1, 1,color = 'b',linewidth = 0.3,linestyle = '--')
    for i in range(size3):
        plt.scatter(mass3[i], 0, color ='red', s = 15, marker = '*')
    for i in range(size4):
        plt.scatter(mass4[i], -0.025, color ='red', s = 15, marker = '*')
    #Горизонтальные линии
    plt.hlines(0.05, -1, 200,color = 'black', linewidth = 1, linestyle = '-')
    plt.hlines(0.025, -1, 200, color = 'black',linewidth = 1,linestyle = '-')
    plt.hlines(-0.025, 1, 200, color = 'black',linewidth = 1,linestyle = '-')
    #Надписи 
    plt.text(60, 0.055, "Поток X1(t) с интенсивностью λ1")
    plt.text(60,  0.03, "Поток X2(t) с интенсивностью λ2")
    plt.text(60, 0.005, "Поток X(t) с интенсивностью  λ1+λ2")
    plt.text(60, -0.02, "Поток Xпр(t) как сумма потоков X1(t) и X2(t) ")

    plt.show()                                               
    #fig.savefig(f'График.png')

#Реализайия алгоритма 4
#Алгоритм построения паусоновского потока (согласно теореме Смирнова)
def alg_four(coef):
    T_loc = T1
   # u_mass=[]
    t_mass = []
    i = 0
    while (T_loc<=T2):
        e  =random.random()
        u = -((math.log(e))/coef)
        t = T_loc+u
        T_loc = t
       # u_mass.append(u)
        t_mass.append(t)
        i += 1
    if t_mass[-1] > T2:
        t_mass.pop()
    return  t_mass
    #return u_mass,t_mass

#Массив на 25  интервалов
def interval(mass):
    mass.sort()
    new_nass = []
    result = []
    T1 = 21
    shag = 25
    T2 = 121
    col_vo=[0]*25
    j = 0
    chek = True
    while (T2 > shag):
        for i in mass:
            if (T1 <= i < shag):
                new_nass.append(i) 
                check = True
            else:
                col_vo[j]+=len(new_nass)
                result.append(new_nass)
                while not(T1 <= i < shag):
                    T1 = shag
                    shag += 4
                    j += 1
                new_nass = []
                new_nass.append(i)
                check = False
        shag = 123
        if chek:
            result.append(new_nass)
            col_vo[j] += len(new_nass)
    return result,col_vo



#ВЫВОДИТ МАССИВА ИНТЕРВАЛО И КОЛ-ВО ЗНАЧЕНИЙ В ИНТЕРВАЛЕ
def variation(mass):
    z=0
    col=[0]*len(mass)
    while z != len(mass):
        mass[z],col[z] = interval(mass[z])
        z += 1
    return mass,col

#Функция получает массив вариантов и массив частот.
def chastotik(mass):
    #chastot-список частот
    #lst=уникальные значнеия(варианты)
    #масс просто из 2 в 1 массив
    lst = []
    for line in mass:
        lst += line
    mass = lst
    lst = list(set(lst))
    chastot = [0]*len(lst)
    i = 0
    while i != len(lst):
        for elem in mass:
            if elem == lst[i]:
                chastot[i] += 1
        i+=1
    return lst,chastot


#Функия подсчета характеристик (Обьем,матожидание,дисперсия)
def charact1(variat,chastot):
    #N=Обьем выбоки
    #n_nl=матрица перемножения значений частот и вариантов 
    #Mn-математическое ожидание(выборочная интенсивность наступления)
    #Dn-дисперсия
    N = sum(chastot)
   # print(f"Обьм выборки={N}")
    n_nl = np.array(variat)*np.array(chastot)
    #print(f"Умножение 2х массивов={n_nl}")
    Mn = sum(n_nl)/N
    #print(f"Мат ожидание ={Mn}")
    Dn = sum(((np.array(variat)-Mn)**2)*np.array(chastot))/N
    print(f"Дисперсия ={Dn}")
    #Дисперсия пока не возвращается 
    return N,Mn


#Функия подсчета характеристик (оценка теоретической вероятности  по формуле Паусона и  Оценка теоретический частот)
def charact2(variation,N,Mn):
    #variation-массива вариантов
    #Pl-оценка теоретической вероятности  по формуле Паусона
    Pl = [0]*len(variation)
    i=0
    while i !=len(variation):
        Pl[i]=((Mn**variation[i])/math.factorial(variation[i]))*math.exp(-Mn)
        i += 1
    NlTeor=np.array(Pl)*N
    return NlTeor

#Функция вичисления  Хи квадрата практ и квантиля хи квадрата
def xi(chastot,NlTeor):
    xi_square = 0
    xi_square = sum(((NlTeor-np.array(chastot))**2)/NlTeor)
    print(f"\nХи-квадрат(практ)={xi_square}")
    xi_krit = chi2.ppf(0.95, 23)
    print(f"\nКвантиль Хи-Квадрата(крит)={xi_krit}")
    if xi_square < xi_krit:
        print("\nУсловие выполнено ")
    else:
        print("\nУсловие не выполнено")
 
        
########
#Генерация 4х выборок(ралзичных интенсивностей)  по 50 реализаций
########
X1_L1 = np.zeros((50), 'object')
X1_L2 = np.zeros((50), 'object')
X1_L1L2 = np.zeros((50), 'object')
X1_L1NaL2 = np.zeros((50), 'object')
for i in range(50):
    one = alg_four(L1)
    two = alg_four(L2)
    X1_L1[i] = one
    X1_L2[i] = two
    X1_L1L2[i] = alg_four(L1+L2)
    z = (one+two)
    z.sort()
    X1_L1NaL2[i] = z
    
#Отрисовака графичекой интерпритации потока событий
draw(X1_L1[0],X1_L2[0],X1_L1L2[0],X1_L1NaL2[0])   
########
#######  Разбитие матрица на интервалы и подсчет попададаний значений в тот или иной интервал
########
peremen1,col_vo_X1_L1 = variation(X1_L1)
peremen2,col_vo_X1_L2 = variation(X1_L2)
peremen3,col_vo_X1_L1L2  = variation(X1_L1L2)
peremen4,col_vo_X1_L1NaL2 = variation(X1_L1NaL2)

########
######## #Вычисление вариантов и частот вариантов
######### 
variation1,chastot1=chastotik(col_vo_X1_L1)
variation2,chastot2=chastotik(col_vo_X1_L2)
variation3,chastot3=chastotik(col_vo_X1_L1L2)
variation4,chastot4=chastotik(col_vo_X1_L1NaL2)
print(f"\nТаблица вариантов для Потока X1(t) с интенсивностью λ1\n {variation1}   ")  
print(f"\nТаблица частот для Потока X1(t) с интенсивностью λ1\n {chastot1} ")  
print(f"\nТаблица вариантов для Потока X2(t) с интенсивностью λ2\n {variation2}   ")  
print(f"\nТаблица частот для Потока X2(t) с интенсивностью λ2\n {chastot2} ")  
print(f"\nТаблица вариантов для Потока X(t) с интенсивностью  λ1+λ2\n {variation3}   ")  
print(f"\nТаблица частот для Поток X(t) с интенсивностью  λ1+λ2\n {chastot3} ")  
print(f"\nТаблица вариантов для  Потока Xпр(t) как сумма потоков X1(t) и X2(t)\n {variation4}   ")  
print(f"\nТаблица частот для Потока Xпр(t) как сумма потоков X1(t) и X2(t)\n {chastot4} ")  
########
######## Подсчет обьема и матожидания(там еще дисперсия)
######### 
N1,Mn1 = charact1(variation1,chastot1)
N2,Mn2 = charact1(variation2,chastot2)
N3,Mn3 = charact1(variation3,chastot3)
N4,Mn4 = charact1(variation4,chastot4)
print(f"\nОбьем выборки Потока X1(t) с интенсивностью λ1 = {N1}  ")  
print(f"\nМат ожидание Потока X1(t) с интенсивностью λ1 = {Mn1} ") 
print(f"\nОбьем выборки Потока X2(t) с интенсивностью λ2 = {N2}  ")  
print(f"\nМат ожидание Потока X2(t) с интенсивностью λ2 = {Mn2} ") 
print(f"\nОбьем выборки Потока X(t) с интенсивностью  λ1+λ2 = {N3}  ")  
print(f"\nМат ожидание  Потока X(t) с интенсивностью  λ1+λ2 = {Mn3} ") 
print(f"\nОбьем выборки Потока Xпр(t) как сумма потоков X1(t) и X2(t) = {N4}  ")  
print(f"\nМат ожидание Потока Xпр(t) как сумма потоков X1(t) и X2(t) = {Mn4} ") 

########
######## Подсчет оценки теоретических частот
######### 
NlTeor1 = charact2(variation1,N1,Mn1)
NlTeor2 = charact2(variation2,N2,Mn2)
NlTeor3 = charact2(variation3,N3,Mn3)
NlTeor4 = charact2(variation4,N4,Mn4)
print(f"\nОценка теоретический частот Потока X1(t) с интенсивностью λ1 = {NlTeor1} ")
print(f"\nОценка теоретический частот Потока X2(t) с интенсивностью λ2 = {NlTeor2} ")
print(f"\nОценка теоретический частот Потока X(t) с интенсивностью  λ1+λ2 = {NlTeor3} ")
print(f"\nОценка теоретический частот Потока Xпр(t) как сумма потоков X1(t) и X2(t) = {NlTeor4} ")
########
######## вычисление Хи-квадрата практ и крит
######### 
xi(chastot4,NlTeor4)
########
######## Сравнение интенсивностей
#########
print(f"Сравнение интенсивности выборочных и теоретических потоков")
print(f"λ1 = {L1}")
print(f"ᴧ1 = {Mn1/del_t}")
print(f"\nλ2 = {L2}")
print(f"ᴧ2 = {Mn2/del_t}")
########
########  Сумма интенсивностей примерно равна интенсивности суммарного потока
#########
print(f"\n Интенсивность суммарного потока")
print(f"λ4 (сумма интенсивностей) = {L1+L2}")
print(f"ᴧ4 для Xпр(t) = {Mn4/del_t}")