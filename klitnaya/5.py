from math import sin, cos, exp, pi
import math
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.figsize'] = (12,12)
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.wedgets import Button

'''def yavnaya(u,a,k,h,t,j):
    return (u[k+1][j] - u[k][j])/t +
'''
#аналитические решения
def analitic(x,t):
    return exp(-4*(pi**2)*a*t)*sin(2*pi*x)

#начальные условия
def u_x_0(x):
    return sin(2*pi*x)

#граничные условия, в одном уравнении тк их возможно обобщить
def u_gr():
    return 0
#Аналитическое решение на сетке
def analitic_on_net(t_, lc, h, a):
    x_c = int(pi/h) + 1
    group = np.zeros((lc, x_c))

    for i in range (lc):
        group[i] = np.array([analitic(j * h, i*t_) for j in range (x_c)])
    return group

#Явная схема
def explicit_solution(t_, lc, h, a):
    x_count = int(pi/ h) + 1
    group = np.zeros((lc, x_count))
    alpha = np.zeros(x_count - 1)
    beta = np.zeros(x_count - 1)
    
    aa = -a * t_/ (h**2)
    bb = 1 + 2 * a * t_/ (h**2)
    cc = aa #-a * t_/ (h**2)
    
    for i in range(lc):
        if i == 0:
            group_ = [u_x_0(j * h) for j in range(x_count)]
        else:
            group_ = np.zeros(x_count)
            beta[0] = u_gr()
            for j in range(1, x_count - 1):
                alpha[j] = -aa/ (bb + cc * alpha[j-1])
                beta[j] = ((group[i-1][j] - aa * (group[i-1][j+1] - 2 * group[i-1][j] + group[i-1][j-1]))) / (bb + cc * alpha[j-1])
            group_[x_count - 1] = u_gr()
            for j in range(x_count - 2, -1, -1):
                group_[j] = group_[j+1] * alpha[j] + beta[j]
        group[i] = np.array(group_)
    return group

'''
    x_c = int(pi/h) + 1
    group = np.zeros((lc, x_c))

    for i in range (lc):
        if i == 0:
            group_ = [u_x_0(j*h) for j in range (x_c)]
        else:
            group_ = [u_gr()] + [group[i-1][j] + a * t_ * (group[i-1][j+1] - 2 * group[i-1][j] + group[i-1][j-1]) / (h ** 2) for j in range(1, x_c - 1)] + [u_gr()]
            group[i] = np.array(group_)
        return group            
'''


#HeЯвная схема
def implicit_solve(t_, lc, h, a):
    x_count = int(pi/ h) + 1
    layers = np.zeros((lc, x_count))
    alpha = np.zeros(x_count - 1)
    beta = np.zeros(x_count - 1)
    
    aa = -a * t_/ (h**2)
    bb = 1 + 2 * a * t_/ (h**2)
    cc = -a * t_/ (h**2)
    
    for i in range(lc):
        if i == 0:
            layer = [u_x_0(j * h) for j in range(x_count)]
        else:
            layer = np.zeros(x_count)
            beta[0] = u_gr()
            for j in range(1, x_count - 1):
                alpha[j] = -aa/ (bb + cc * alpha[j-1])
                beta[j] = (layers[i-1][j] - cc * beta[j-1]) / (bb + cc * alpha[j-1])
            layer[x_count - 1] = u_gr()
            for j in range(x_count - 2, -1, -1):
                layer[j] = layer[j+1] * alpha[j] + beta[j]
        layers[i] = np.array(layer)
    return layers

#схема Кранка-Николсона
def KN_solve(t_, lc, h, a):
    x_count = int(pi/ h) + 1
    group = np.zeros((lc, x_count))
    alpha = np.zeros(x_count - 1)
    beta = np.zeros(x_count - 1)
    
    aa = -a * t_/ (2 * h**2)
    bb = 1 + a * t_/ (h**2)
    cc = aa
    
    for i in range(lc):
        if i == 0:
            group_ = [u_x_0(j * h) for j in range(x_count)]
        else:
            group_ = np.zeros(x_count)
            beta[0] = u_gr()
            for j in range(1, x_count - 1):
                alpha[j] = -aa/ (bb + cc * alpha[j-1])
                beta[j] = ((group[i-1][j] - aa * (group[i-1][j+1] - 2 * group[i-1][j] + group[i-1][j-1])) - cc * beta[j-1]) / (bb + cc * alpha[j-1])
            group_[x_count - 1] = u_gr()
            for j in range(x_count - 2, -1, -1):
                group_[j] = group_[j+1] * alpha[j] + beta[j]
        group[i] = np.array(group_)
    return group


#погрешность
def error(explicit, analitic):
    er = sum((explicit - analitic)**2)/ len(explicit)
    return er

#для кнопок
def button_callback(event):
    plt.clf()
    if label == 'Явная схема.':
        name = 'Явная схема.'
        print ("Запишем решение в явном виде ")
        res = explicit_solution(t_, lc, h, a)
        plot_function( name, analitic, res, t_, lc, h)
    elif label == 'HeЯвная схема.':
        name = 'HeЯвная схема.'
        res = implicit_solve(t_, lc, h, a)
        plot_function( name, analitic, res, t_, lc, h)
    elif label == 'схема Кранка-Николсона':
        name = 'Схема Кранка-Николсона.'
        res = KN_solve(t_, lc, h, a)
        plot_function( name, analitic, res, t_, lc, h)
   
def update_plot(analitic, t_, lc, h):
    fst.clear()
    scd.clear()
    thd.clear()
    frth.clear()
    fvth.clear()
    sxth.clear()

    
    cur_display_mode = display_modes[display_index]
    if cur_display_mode == 1:
        name = 'Явная схема.'
        print ("Запишем решение в явном виде ")
        res = explicit_solution(t_, lc, h, a)
        plot_function( name, analitic, res, t_, lc, h)
    elif cur_display_mode == 2:
        name = 'HeЯвная схема.'
        res = implicit_solve(t_, lc, h, a)
        plot_function( name, analitic, res, t_, lc, h)
    elif cur_display_mode == 3:
        name = 'Схема Кранка-Николсона.'
        res = KN_solve(t_, lc, h, a)
        plot_function( name, analitic, res, t_, lc, h)
    fig.canvas.draw()


def on_key(event):
    global display_index

    if event.key == 'right':
        display_index = (display_index + 1) % len(display_modes)
        update_plot(analitic, t_, lc, h)
    elif event.key == 'left':
        display_index = (display_index - 1) % len(display_modes)
        update_plot(analitic, t_, lc, h)




#функция визуализации
def plot_function( name, analiyic, res, t_, lc, h):
    x_c = int(pi/h)+1
    X = np.array([i * h for i in range(x_c)])

    for layer in analitic:
        fst.scatter(X, layer, marker='x')
    fst.legend(['t = {}'.format(i * t_) for i in range(lc)], fontsize=5, loc='upper right')
    
    scd.set_title(name + 'Результаты сетке', fontsize=10)
    scd.set_xlabel('x', fontsize=5)
    scd.set_ylabel('U(x, t)', fontsize=5)

    for layer in res:
        scd.scatter(X, layer, marker='x')
    scd.legend(['t = {}'.format(i * t_) for i in range(lc)], fontsize=5, loc='upper right')
    T = np.array([i * t_ for i in range(lc)])
    MSE_error = np.array([error(i, j) for i, j in zip(res, analitic)])
    thd.set_title('MSE: h = {}, tau = {}'.format(h, t_), fontsize=10)
    thd.set_xlabel('t', fontsize=5)
    thd.set_ylabel('mse_error', fontsize=5)
    thd.plot(T, MSE_error)
    thd.scatter(T, MSE_error, marker='o', c='r', s=50)
    thd.legend(['Значение ошибки в разные моменты времени'], fontsize=5, loc='upper right')

    X, T = np.meshgrid(X, T)
    U_real = analitic.ravel().reshape(X.shape)
    frth.set_title('Аналитическое решение. Результаты сетке, 3D', fontsize=10)
    frth.set_xlabel('x', fontsize=5)
    frth.set_ylabel('t', fontsize=5)
    frth.set_zlabel('u', fontsize=5)
    surf_real = frth.plot_surface(X, T, U_real, rstride=1, cstride=1, cmap=cm.coolwarm,
                 linewidth=0, antialiased=False)
    #fig.colorbar(surf_real, ax=frth, shrink=0.5, aspect=10)

    U_kn = res.ravel().reshape(X.shape)
    fvth.set_title(name +'Результаты на сетке, 3D', fontsize=10)
    fvth.set_xlabel('x', fontsize=10)
    fvth.set_ylabel('t', fontsize=10)
    fvth.set_zlabel('u', fontsize=10)
    surf_kn = fvth.plot_surface(X, T, U_kn, rstride=1, cstride=1, cmap=cm.twilight,
                linewidth=0, antialiased=False)
    #fig.colorbar(surf_kn, ax=fvth, shrink=0.5, aspect=10)

    U_real = analitic.ravel().reshape(X.shape)
    U_kn = res.ravel().reshape(X.shape)
    U_err = (U_real - U_kn)**2

    sxth.set_title('Ошибка. Результаты на сетке, 3D', fontsize=10)
    sxth.set_xlabel('x', fontsize=10)
    sxth.set_ylabel('t', fontsize=10)
    sxth.set_zlabel('u', fontsize=10)
    surf_err= sxth.plot_surface(X, T, U_err, rstride=1, cstride=1, cmap=cm.hsv,
                linewidth=0, antialiased=False)
    #fig.colorbar(surf_err, ax=sxth, shrink=0.5, aspect=10)


    
    


a = 0.5#float(input()) #0.5
t_= 0.01#float(input()) #0.1
h = 0.1#float(input()) #0.5
lc = 5#float(input()) #3


print ("Наша система:")
print ('\tdu/dt = a* d^2u/dx^2, a>0',
       'u(0,t) = 0,',
       'u(1,t)=0,',
       'u(x,0)=sin(2pix)', sep = '\n\t')
print ("")
o = False
analitic = analitic_on_net(t_, lc, h, a)
fig = plt.figure(figsize=(50, 35))
fst = fig.add_subplot(2, 3, 1)
scd = fig.add_subplot(2, 3, 2)
thd = fig.add_subplot(2, 3, 3)  
frth = fig.add_subplot(2, 3, 4, projection='3d')
fvth = fig.add_subplot(2, 3, 5, projection='3d')
sxth = fig.add_subplot(2, 3, 6, projection='3d')
    
fst.set_title('Аналитическое решение. Результаты сетке', fontsize=10)
fst.set_xlabel('x', fontsize=5)
fst.set_ylabel('U(x, t)', fontsize=5)
        
'''        
while (o != True):
    bar = int(input())
    if bar == 0 :
        while ((t_/h**2) > (1/ (2*a))):
            print("неверно введено h, !!условие устойчивости нарушено!!")
            print(t_/h**2, ' ', 1/ (2*a))
            h = float(input())
        name = 'Явная схема.'
        print ("Запишем решение в явном виде ")
        res = explicit_solution(t_, lc, h, a)
        plot_function( name, analitic, res, t_, lc, h)
        bar = int(input())
    if bar == 1:
        while ((2* t_ * a /h**2) >= (1+ 2* a * t_/ (h**2))):
            print("неверно введено h, !!условие устойчивости нарушено!!")
            print(t_/h**2, ' ', 1/ (2*a))
            h = float(input())
        print ("")
        name = 'HeЯвная схема.'
        res = implicit_solve(t_, lc, h, a)
        plot_function( name, analitic, res, t_, lc, h)
    if bar == 2:
        while ((2* t_ * a /h**2) >= (1+ 2* a * t_/ (h**2))):
            print("неверно введено h, !!условие устойчивости нарушено!!")
            print(t_/h**2, ' ', 1/ (2*a))
            h = float(input())
        print ("")
        name = 'Схема Кранка-Николсона.'
        res = KN_solve(t_, lc, h, a)
        plot_function( name, analitic, res, t_, lc, h)
        '''
display_modes  = [1, 2, 3]
display_index = 0
# Привязка обработчика нажатий клавиш к схеме
fig.canvas.mpl_connect('key_press_event', on_key)
# Первоначальное отображение графиков
update_plot(analitic, t_, lc, h)

# Отображение схемы
plt.show()
