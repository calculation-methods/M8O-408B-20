from math import sin, cos, exp, pi, log
import math
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.figsize'] = (12,12)
from mpl_toolkits.mplot3d import Axes3D

#аналитические решения
def analitic(a,x,t,y, mu1, mu2):
    return cos(x) * cos(y)* exp((-(mu1**2 + mu2**2))*a*t)
#Граничные условия.
def u_x0(y,t,a,mu1,mu2):
    return cos(mu2*y) * exp(-(mu1**2 + mu2**2) * a*t)

def u_xpi(y,t,a,mu1,mu2):
    return (-1)**(mu1) * cos(mu2*y) * exp(-(mu1**2 + mu2**2)*a*t)

def u_y0(x,t,a,mu1,mu2):
    return cos(mu1*x)*exp(-(mu1**2 + mu2**2)*a*t)

def u_ypi(x,t,a,mu1,mu2):
    return (-1)**(mu2) * cos(mu1*x)*exp(-(mu1**2 + mu2**2)*a*t)


#Начальное условие.
def u_xy(x,y,a,mu1,u2):
    return cos(mu1*x)*cos(mu2*y)

#Аналитическое решение на сетке

def analitic_on_net(a, dx, dy, t_, lc, mu1, mu2):
    x_count = int(pi/ (4*dx)) + 1
    y_count = int(log(2)/ dy) + 1
    layers = np.zeros((lc, y_count, x_count))
    
    for i in range(lc):
        for j in range(y_count):
            layers[i][j] = np.array([analitic(k * dx, j * dy, i * t_, a, mu1, mu2) for k in range(x_count)])
    return layers

#погрешность
def error(explicit, analitic):
    er = sum((explicit - analitic)**2)/ len(explicit)
    return er

#Метод переменных направлений
def variable_direction_method(a, dx, dy, t_, lc, mu1,mu2):
    
    x_count = int(pi/ (4*dx)) + 1
    y_count = int(log(2)/dy) + 1
    
    layers = np.zeros((lc, y_count, x_count))
    layers[0] = np.array([[ u_xy(i * dx, j * dy, a,mu1,mu2) for i in range(x_count)] for j in range(y_count)])
    
    for i in range(1, lc):
        
        aa = -a * t_/ (2 * dx**2)
        bb = 1 + a * t_/ (dx**2)
        cc = aa
        
        fract_layer = np.zeros((y_count, x_count))
        
        for k in range(1, y_count - 1):
           
            alpha = np.zeros(x_count - 1)
            beta = np.zeros(x_count - 1)
            beta[0] = u_x0(k * dy, t_ * (i - 0.5),a,mu1,mu2)
            
            for j in range(1, x_count - 1):
                alpha[j] = -aa/ (bb + cc * alpha[j-1])
                xi_jk = layers[i-1][k][j] + a * t_/ 2 * (layers[i-1][k+1][j] - 2 * layers[i-1][k][j] + layers[i-1][k-1][j])/ (dy**2)
                beta[j] = (xi_jk - cc * beta[j-1]) / (bb + cc * alpha[j-1])
            
            fract_layer[k][x_count - 1] = u_xpi(k * dy, t_ * (i - 0.5), a,mu1,mu2)
            
            for j in range(x_count - 2, -1, -1):
                fract_layer[k][j] = fract_layer[k][j+1] * alpha[j] + beta[j]
        fract_layer[0] = np.array([u_y0(k * dx, t_ * (i - 0.5),a,mu1,mu2) for k in range(x_count)])
        fract_layer[y_count - 1] = np.array([u_ypi(k * dx, t_ * (i - 0.5),a, mu1,mu2) for k in range(x_count)])
        
        aa = -a * t_/ (2 * dy**2)
        bb = 1 + a * t_/ (dy**2)
        cc = aa
        
        layer = np.zeros((y_count, x_count))
        
        for k in range(1, x_count - 1):
           
            alpha = np.zeros(y_count - 1)
            beta = np.zeros(y_count - 1)
            beta[0] = u_y0(k * dx, t_ * i,a, mu1,mu2)
            
            for j in range(1, y_count - 1):
                alpha[j] = -aa/ (bb + cc * alpha[j-1])
                xi_jk = fract_layer[j][k] + a * t_/ 2 * (fract_layer[j][k+1] - 2 * fract_layer[j][k] + fract_layer[j][k-1])/ (dx**2)
                beta[j] = (xi_jk - cc * beta[j-1]) / (bb + cc * alpha[j-1])
            
            layer[y_count - 1][k] = u_ypi(k * dx, t_ * i,a,mu1,mu2)
            
            for j in range(y_count - 2, -1, -1):
                layer[j][k] = layer[j+1][k] * alpha[j] + beta[j]
                
        layer[:, 0] = np.array([u_x0(dy * j, t_ * i,a, mu1,mu2) for j in range(y_count)])
        layer[:, x_count - 1] = np.array([u_xpi(dy * j, t_ * i,a, mu1,mu2) for j in range(y_count)])
        layers[i] = layer
        
    return layers


#Метод дробных шагов
def double_steps(a, dx, dy, t_, lc,mu1,mu2 ):
    x_count = int(pi/ (4*dx)) + 1
    y_count = int(log(2)/dy) + 1
    
    layers = np.zeros((lc, y_count, x_count))
    layers[0] = np.array([[ u_xy(i * dx, j * dy,a, mu1,mu2) for i in range(x_count)] for j in range(y_count)])
    
    for i in range(1, lc):
        aa = -a * t_/ (dx**2)
        bb = 1 + 2 * a * t_/ (dx**2)
        cc = aa
        
        fract_layer = np.zeros((y_count, x_count))
        
        for k in range(1, y_count - 1):
            alpha = np.zeros(x_count - 1)
            beta = np.zeros(x_count - 1)
            beta[0] = u_x0(k * dy, t_ * (i - 0.5),a,mu1,mu2)
            
            for j in range(1, x_count - 1):
                alpha[j] = -aa/ (bb + cc * alpha[j-1])
                beta[j] = (layers[i-1][k][j] - cc * beta[j-1]) / (bb + cc * alpha[j-1])
            
            fract_layer[k][x_count - 1] = u_xpi(k * dy, t_ * (i - 0.5),a,mu1,mu2)
            
            for j in range(x_count - 2, -1, -1):
                fract_layer[k][j] = fract_layer[k][j+1] * alpha[j] + beta[j]
                
        fract_layer[0] = np.array([u_y0(k * dx, t_ * (i - 0.5),a,mu1,mu2) for k in range(x_count)])
        fract_layer[y_count - 1] = np.array([u_ypi(k * dx, t_ * (i - 0.5),a,mu1,mu2) for k in range(x_count)])
        
        aa = -a * t_/ (dy**2)
        bb = 1 + 2 * a * t_/ (dy**2)
        cc = aa
        
        layer = np.zeros((y_count, x_count))
        
        for k in range(1, x_count - 1):
            alpha = np.zeros(y_count - 1)
            beta = np.zeros(y_count - 1)
            beta[0] = u_y0(k * dx, t_ * i,a,mu1,mu2)
            
            for j in range(1, y_count - 1):
                alpha[j] = -aa/ (bb + cc * alpha[j-1])
                beta[j] = (fract_layer[j][k] - cc * beta[j-1]) / (bb + cc * alpha[j-1])
            
            layer[y_count - 1][k] = u_ypi(k * dx, t_ * i,a,mu1,mu2)
            
            for j in range(y_count - 2, -1, -1):
                layer[j][k] = layer[j+1][k] * alpha[j] + beta[j]
                
        layer[:, 0] = np.array([u_x0(dy * j, t_ * i,a,mu1,mu2) for j in range(y_count)])
        layer[:, x_count - 1] = np.array([u_xpi(dy * j, t_ * i,a,mu1,mu2) for j in range(y_count)])
        layers[i] = layer
    return layers
        
'''    
#для кнопок
def button_callback(event):
    plt.clf()
    if label == 'Метод дробных шагов':
        name = 'Метод дробных шагов'
        res = double_steps(a, dx, dy, t_, lc,mu1,mu2 )
        plot_function( name, analitic, res, t_, lc, h)
    elif label == 'Схема переменных направлений':
        name = 'Схема переменных направлений'
        res = variable_direction_method(a, dx,dy, t_, lc,mu1,mu2)
        plot_function( name, analitic, res, a,dx,dy, t_, lc,mu1,mu2)
    
def update_plot(analitic, dx,dy, t_, lc,a,mu1,mu2):
    fst.clear()
    scd.clear()
    thd.clear()
    frth.clear()
    fvth.clear()
    sxth.clear()

    
    cur_display_mode = display_modes[display_index]
    if cur_display_mode == 1:
        name = 'Метод дробных шагов'
        res = double_steps(a, dx, dy, t_, lc,mu1,mu2 )
        plot_function( name, analitic, res, a,dx,dy, t_, lc,mu1,mu2)
    elif cur_display_mode == 2:
        name = 'Схема переменных направлений'
        res = variable_direction_method(a, dx,dy, t_, lc,mu1,mu2)
        plot_function( name, analitic, res, a,dx,dy, t_, lc,mu1,mu2)
    fig.canvas.draw()


def on_key(event):
    global display_index

    if event.key == 'right':
        display_index = (display_index + 1) % len(display_modes)
        update_plot( analitic, dx, dy, t_, lc,a,mu1,mu2)
    elif event.key == 'left':
        display_index = (display_index - 1) % len(display_modes)
        update_plot( analitic, dx, dy, t_, lc, a, mu1,mu2)




#функция визуализации
def plot_function( name, analitic, res, a,dx,dy, t_, lc,mu1,mu2):

    x_c = int(pi/(4*dx))+1
    X = np.array([i * dx for i in range(x_c)])

    y_c = int(log(2)/dy)+1
    Y = np.array([i * dy for i in range(y_c)])


    y_c = int(y_c - 1)
    t_c = lc -1

    fst.set_title('Аналитическое решение. Результаты сетке', fontsize=10)
    fst.set_xlabel('x', fontsize=5)
    fst.set_ylabel('U(x,{}, t)'.format(y_c*dy), fontsize=5)

    
    for layer in analitic:
        fst.scatter(X, layer[y_c], marker='x')
    fst.legend(['t = {}'.format(i * t_) for i in range(lc)], fontsize=5, loc='upper right')
    
    scd.set_title(name + 'Результаты сетке', fontsize=10)
    scd.set_xlabel('x', fontsize=5)
    scd.set_ylabel('U(x,{}, t)'.format(y_c * dy), fontsize = 5)

    for layer in res:
        scd.scatter(X, layer[y_c], marker='x')
    scd.legend(['t = {}'.format(i * t_) for i in range(lc)], fontsize=5, loc='upper right')

    T = np.array([i * t_ for i in range(lc)])
    
    MSE_error = np.array([error(i[y_c], j[y_c]) for i, j in zip(res, analitic)])
    thd.set_title('MSE: dx = {},dy ={}, tau = {}'.format(dx,dy, t_), fontsize=10)
    thd.set_xlabel('t', fontsize=5)
    thd.set_ylabel('mse_error', fontsize=5)
    thd.plot(T, MSE_error)
    thd.scatter(T, MSE_error, marker='o', c='r', s=50)
    thd.legend(['Значение ошибки в разные моменты времени'.format(y_c*dy)], fontsize=5, loc='upper right')
    

    X, Y = np.meshgrid(X, Y)
    
   # U_real = analitic.ravel().reshape(X.shape)
    frth.set_title('Аналитическое решение. Результаты сетке, 3D', fontsize=10)
    frth.set_xlabel('x', fontsize=5)
    frth.set_ylabel('t', fontsize=5)
    frth.set_zlabel('u', fontsize=5)

    for t in range(lc):
        U_real= analitic[t].ravel().reshape(X.shape)
        sur_real = frth.plot_surface(X, Y, U_real, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    #fig.colorbar(sur_real, ax=frth, shrink=0.5, aspect=10)
    #surf_real = frth.plot_surface(X, T, U_real, rstride=1, cstride=1, cmap=cm.coolwarm,
     #            linewidth=0, antialiased=False)
    #fig.colorbar(surf_real, ax=frth, shrink=0.5, aspect=10)

    #U_kn = res.ravel().reshape(X.shape)
     
    fvth.set_title(name +'Результаты на сетке, 3D', fontsize=10)
    fvth.set_xlabel('x', fontsize=10)
    fvth.set_ylabel('t', fontsize=10)
    fvth.set_zlabel('u', fontsize=10)
    for t in range(lc):
        U_kn = res[t].ravel().reshape(X.shape)
        surf_kn = fvth.plot_surface(X, Y, U_kn, rstride=1, cstride=1, cmap=cm.twilight,
                linewidth=0, antialiased=False)
    #fig.colorbar(surf_kn, ax=fvth, shrink=0.5,aspect=10)
    #fig.colorbar(surf_kn, ax=fvth, shrink=0.5, aspect=10)

    
    

    sxth.set_title('Ошибка. Результаты на сетке, 3D', fontsize=10)
    sxth.set_xlabel('x', fontsize=10)
    sxth.set_ylabel('t', fontsize=10)
    sxth.set_zlabel('u', fontsize=10)
    
    for t in range(lc):
        U_real = analitic.ravel().reshape(X.shape)
        U_kn = res.ravel().reshape(X.shape)
        U_err = (U_real - U_kn)**2
        surf_err=sxth.plot_surface(X,Y,U_err, rstride=1,cstride=1,cmap=cm.hsv,
                           linewidth=0, antialiased=False)
    #fig.colorbar(surf_err, ax=sxth, shrink=0.5, aspect=10)
    
'''
print("введите параметры а, t_, dx, dy, lc ")
a = float(input()) #0.5
t_= float(input()) #0.1
dx = float(input())#0.1
dy=float(input())#0.1
lc = int(input())#5


print("введите номер 1, 2 или 3 для определения параметров мю1 и мю2")

order = int(input())
#order = 1
if (order == 1):
    mu1 = 1
    mu2 = 1
elif(order == 2):
    mu1 = 2
    mu2 = 1
elif(order == 3):
    mu1 = 1
    mu2 = 2


    
real_res = analitic_on_net(a, dx, dy, t_, lc, mu1, mu2)
var_res = double_steps(a, dx, dy, t_, lc,mu1,mu2 )#variable_direction_method(a, dx, dy, t_, lc, mu1,mu2)

layers_count = lc

'''display_modes  = [1, 2]
display_index = 0
# Привязка обработчика нажатий клавиш к схеме
fig.canvas.mpl_connect('key_press_event', on_key)
# Первоначальное отображение графиков
update_plot( analitic, dx, dy, t_, lc, a, mu1,mu2)

# Отображение схемы
plt.show()
        
    '''

x_count = int(pi/ (4*dx)) + 1
y_count = int(log(2)/ dy) + 1

X = np.array([i * dx for i in range(x_count)])
Y = np.array([i * dy for i in range(y_count)])

fig = plt.figure(figsize=(50, 35))
fst = fig.add_subplot(2, 3, 1)
scd = fig.add_subplot(2, 3, 2)
thd = fig.add_subplot(2, 3, 3)
frth = fig.add_subplot(2, 3, 4, projection='3d')
fvth = fig.add_subplot(2, 3, 5, projection='3d')
sxth = fig.add_subplot(2, 3, 6, projection='3d')

y_const = int(y_count - 1)
t_const = layers_count - 1

fst.set_title('Аналитическое решение.y ={}'.format(y_const * dy), fontsize=10)
fst.set_xlabel('x', fontsize=10)
fst.set_ylabel('U(x, {}, t)'.format(y_const * dy), fontsize=10)

for layer in real_res:
    fst.scatter(X, layer[y_const], marker='x')
fst.legend(['t = {}'.format(i * t_) for i in range(layers_count)], fontsize=5, loc='upper right')

scd.set_title('Схема переменных направлений. y = {}'.format(y_const * dy), fontsize=7)
scd.set_xlabel('x', fontsize=10)
scd.set_ylabel('U(x, {}, t)'.format(y_const * dy), fontsize=10)

for layer in var_res:
    scd.scatter(X, layer[y_const], marker='x')
scd.legend(['t = {}'.format(i * t_) for i in range(layers_count)], fontsize=5, loc='upper right')

T = np.array([i * t_ for i in range(layers_count)])

MSE_error = np.array([error(i[y_const], j[y_const]) for i, j in zip(var_res, real_res)])
thd.set_title('MSE: dx = {}, dy = {}, tau = {}'.format(dx, dy, t_), fontsize=10)
thd.set_xlabel('t', fontsize=10)
thd.set_ylabel('mse_error', fontsize=10)
thd.plot(T, MSE_error)
thd.scatter(T, MSE_error, marker='o', c='r', s=50)
thd.legend(['Значение ошибки при y = {}'.format(y_const * dy)], fontsize=5, loc='upper right')

X, Y = np.meshgrid(X, Y)
                      
frth.set_title('Аналитическое решение. Результаты на всей сетке, 3D', fontsize=10)
frth.set_xlabel('x', fontsize=10)
frth.set_ylabel('y', fontsize=10)
frth.set_zlabel('u', fontsize=10)
                      
for t in range(layers_count):
    U_real = real_res[t].ravel().reshape(X.shape)
    surf_real = frth.plot_surface(X, Y, U_real, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
fig.colorbar(surf_real, ax=frth, shrink=0.5, aspect=10)
#frth.legend(['t = {}'.format(i * tau) for i in range(layers_count)], fontsize=20, loc='upper right')
                      
#U_num = var_res.ravel().reshape(X.shape)
fvth.set_title('Схема переменных направлений. Результаты на сетке 3D', fontsize=7)
fvth.set_xlabel('x', fontsize=5)
fvth.set_ylabel('y', fontsize=5)
fvth.set_zlabel('u', fontsize=5)

for t in range(layers_count):
    U_num = var_res[t].ravel().reshape(X.shape)
    surf_num = fvth.plot_surface(X, Y, U_num, rstride=1, cstride=1, cmap=cm.twilight,
                           linewidth=0, antialiased=False)
fig.colorbar(surf_num, ax=fvth, shrink=0.5, aspect=10)
#fvth.legend(['t = {}'.format(i * tau) for i in range(layers_count)], fontsize=20, loc='upper right')

sxth.set_title('Ошибка. Результаты на 3D', fontsize=10)
sxth.set_xlabel('x', fontsize=5)
sxth.set_ylabel('y', fontsize=5)
sxth.set_zlabel('u', fontsize=5)                      

for t in range(layers_count):
    U_real = real_res[t].ravel().reshape(X.shape)
    U_num = var_res[t].ravel().reshape(X.shape)
    U_err = (U_real - U_num)**2
    surf_err= sxth.plot_surface(X, Y, U_err, rstride=1, cstride=1, cmap=cm.hsv,
                           linewidth=0, antialiased=False)
fig.colorbar(surf_err, ax=sxth, shrink=0.5, aspect=10)


#plt.show()

real_res = analitic_on_net(a, dx, dy, t_, lc, mu1, mu2)
var_res = variable_direction_method(a, dx, dy, t_, lc, mu1,mu2)

layers_count = lc
   
x_count = int(pi/ (4*dx)) + 1
y_count = int(log(2)/ dy) + 1

X = np.array([i * dx for i in range(x_count)])
Y = np.array([i * dy for i in range(y_count)])

fig1 = plt.figure(figsize=(50, 35))
fst = fig1.add_subplot(2, 3, 1)
scd = fig1.add_subplot(2, 3, 2)
thd = fig1.add_subplot(2, 3, 3)
frth = fig1.add_subplot(2, 3, 4, projection='3d')
fvth = fig1.add_subplot(2, 3, 5, projection='3d')
sxth = fig1.add_subplot(2, 3, 6, projection='3d')

y_const = int(y_count - 1)
t_const = layers_count - 1

fst.set_title('Аналитическое решение. y = {}'.format(y_const * dy), fontsize=10)
fst.set_xlabel('x', fontsize=10)
fst.set_ylabel('U(x, {}, t)'.format(y_const * dy), fontsize=10)

for layer in real_res:
    fst.scatter(X, layer[y_const], marker='x')
fst.legend(['t = {}'.format(i * t_) for i in range(layers_count)], fontsize=5, loc='upper right')

scd.set_title('Схема дробных шагов.y = {}'.format(y_const * dy), fontsize=10)
scd.set_xlabel('x', fontsize=10)
scd.set_ylabel('U(x, {}, t)'.format(y_const * dy), fontsize=10)

for layer in var_res:
    scd.scatter(X, layer[y_const], marker='x')
scd.legend(['t = {}'.format(i * t_) for i in range(layers_count)], fontsize=10, loc='upper right')

T = np.array([i * t_ for i in range(layers_count)])

MSE_error = np.array([error(i[y_const], j[y_const]) for i, j in zip(var_res, real_res)])
thd.set_title('MSE: dx = {}, dy = {}, tau = {}'.format(dx, dy, t_), fontsize=10)
thd.set_xlabel('t', fontsize=10)
thd.set_ylabel('mse_error', fontsize=10)
thd.plot(T, MSE_error)
thd.scatter(T, MSE_error, marker='o', c='r', s=50)
thd.legend(['Значение ошибки при y = {}'.format(y_const * dy)], fontsize=5, loc='upper right')

X, Y = np.meshgrid(X, Y)
                      
frth.set_title('Аналитическое решение. Результаты на 3D', fontsize=10)
frth.set_xlabel('x', fontsize=5)
frth.set_ylabel('y', fontsize=5)
frth.set_zlabel('u', fontsize=5)
                      
for t in range(layers_count):
    U_real = real_res[t].ravel().reshape(X.shape)
    surf_real = frth.plot_surface(X, Y, U_real, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
fig1.colorbar(surf_real, ax=frth, shrink=0.5, aspect=10)
#frth.legend(['t = {}'.format(i * tau) for i in range(layers_count)], fontsize=20, loc='upper right')
                      
#U_num = var_res.ravel().reshape(X.shape)
fvth.set_title('Схема дробных шагов.3D', fontsize=10)
fvth.set_xlabel('x', fontsize=10)
fvth.set_ylabel('y', fontsize=10)
fvth.set_zlabel('u', fontsize=10)

for t in range(layers_count):
    U_num = var_res[t].ravel().reshape(X.shape)
    surf_num = fvth.plot_surface(X, Y, U_num, rstride=1, cstride=1, cmap=cm.twilight,
                           linewidth=0, antialiased=False)
fig1.colorbar(surf_num, ax=fvth, shrink=0.5, aspect=10)
#fvth.legend(['t = {}'.format(i * tau) for i in range(layers_count)], fontsize=20, loc='upper right')

sxth.set_title('Ошибка. 3D', fontsize=10)
sxth.set_xlabel('x', fontsize=10)
sxth.set_ylabel('y', fontsize=10)
sxth.set_zlabel('u', fontsize=10)                      

for t in range(layers_count):
    U_real = real_res[t].ravel().reshape(X.shape)
    U_num = var_res[t].ravel().reshape(X.shape)
    U_err = (U_real - U_num)**2
    surf_err= sxth.plot_surface(X, Y, U_err, rstride=1, cstride=1, cmap=cm.hsv,
                           linewidth=0, antialiased=False)
fig1.colorbar(surf_err, ax=sxth, shrink=0.5, aspect=10)


plt.show()

