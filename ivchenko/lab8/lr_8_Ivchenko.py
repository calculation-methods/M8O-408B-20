import math as m
import matplotlib.pyplot as plt
import numpy as np
import copy

n = 10
x_end = m.pi/2
y_end = m.pi
t_end = 3
acoef = 2
bcoef = 1
mu = 1
x = np.linspace(0,x_end,n)
y = np.linspace(0,y_end,n)
t = np.linspace(0,t_end,n)
x_plt, y_plt = np.meshgrid(x, y)
h1 = x_end/n
h2 = y_end/n
tau = t_end/n

# Аналитическое решение
def isxF(x,y,t):
    u = np.zeros((n,n,n))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                u[i][j][k] = np.sin(x[i])*np.sin(y[j])*np.sin(mu*t[k])

    return u

U = isxF(x,y,t)
fig = plt.figure()
ax = plt.axes(projection ='3d')
ax.plot_surface(x_plt,y_plt,np.array(U[:,:,3]))
plt.show()

# МДШ
def mdsh(n):
    u = np.zeros((n,n,n))
    for i in range(n):
        for j in range(n):
            u[i][j][0] = 0

    for k in range(1,n):
        u1 = np.zeros((n,n))
        u2 = np.zeros((n,n))
        tau2 = t[k-1] + tau/2
        for i in range(n):
            u1[i][0] = 0
            u1[0][i] = 0
            u1[-1][i] = np.sin(y[i])*np.sin(mu*tau2)
            u2[i][0] = 0

        for j in range(n-1):
            a = np.zeros(n)
            b = np.zeros(n)
            c = np.zeros(n)
            d = np.zeros(n)
            a[0] = 0
            b[0] = 1
            c[0] = 0
            d[0] = 0
            for i in range(1,n-1):
                a[i] = acoef/h1/h1
                b[i] = -(1/tau + 2*acoef/h1/h1)
                c[i] = acoef/h1/h1
                d[i] = -(f(x[i],y[j],t[k-1])/2+u[i][j][k-1]/tau)
            a[-1] = 0
            b[-1] = 1
            c[-1] = 0
            d[-1] = np.sin(y[j])*np.sin(mu*tau2)

            result = slau([a,b,c,d])
            #print(result)
            for i in range(n):
                u1[i][j] = result[i]
                u1[i][-1] = -np.sin(x[i])*np.sin(mu*tau2)*h1 + u1[i][-2]


        for i in range(n-1):
            a = np.zeros(n)
            b = np.zeros(n)
            c = np.zeros(n)
            d = np.zeros(n)
            a[0] = 0
            b[0] = 1
            c[0] = 0
            d[0] = 0
            for j in range(1,n-1):
                a[j] = bcoef/h2/h2
                b[j] = -(1/tau + 2*bcoef/h2/h2)
                c[j] = bcoef/h2/h2
                d[j] = -(f(x[i],y[j],t[k])/2+u1[i][j]/tau)
            a[-1] = -1/h2
            b[-1] = 1/h2
            c[-1] = 0
            d[-1] = -np.sin(x[i])*np.sin(mu*t[k])

            result2 = slau([a,b,c,d])
            for j in range(n):
                u2[i][j] = result2[j]
                u2[0][j] = 0
                u2[-1][j] = np.sin(y[j])*np.sin(t[k])

        for i in range(n):
            u2[i][-1] = -np.sin(x[i])*np.sin(mu*t[k])*h2 + u2[i][-2]

        for i in range(n):
            for j in range(n):
                u[i][j][k] = u2[i][j]
    return u

U2 = mdsh(n)
fig = plt.figure()
ax = plt.axes(projection ='3d')
ax.plot_surface(x_plt,y_plt,np.array(U2[:,:,3]))
plt.show()

# МПН
def mpn(n):
    u = np.zeros((n,n,n))
    for i in range(n):
        for j in range(n):
            u[i][j][0] = 0

    for k in range(1,n):
        u1 = np.zeros((n,n))
        u2 = np.zeros((n,n))
        tau2 = t[k-1] + tau/2
        for i in range(n):
            u1[i][0] = 0
            u1[0][i] = 0
            u1[-1][i] = np.sin(y[i])*np.sin(mu*tau2)
            u2[i][0] = 0

        for j in range(n-1):
            a = np.zeros(n)
            b = np.zeros(n)
            c = np.zeros(n)
            d = np.zeros(n)
            a[0] = 0
            b[0] = 1
            c[0] = 0
            d[0] = 0
            for i in range(1,n-1):
                a[i] = acoef/h1/h1*tau/2
                b[i] = -(2/tau + 2*acoef/h1/h1)*tau/2
                c[i] = acoef/h1/h1*tau/2
                d[i] = -(2*u[i][j][k-1]/tau + bcoef/h2/h2*(u[i][j+1][k-1]-2*u[i][j][k-1]+u[i][j-1][k-1])-f(x[i],y[j],tau2))*tau/2
            a[-1] = 0
            b[-1] = 1
            c[-1] = 0
            d[-1] = np.sin(y[j])*np.sin(mu*tau2)

            result = slau([a,b,c,d])
            for i in range(n):
                u1[i][j] = result[i]
                u1[i][-1] = -np.sin(x[i])*np.sin(mu*tau2)*h1 + u1[i][-2]

        for i in range(n-1):
            a = np.zeros(n)
            b = np.zeros(n)
            c = np.zeros(n)
            d = np.zeros(n)
            a[0] = 0
            b[0] = 1
            c[0] = 0
            d[0] = 0
            for j in range(1,n-1):
                a[j] = bcoef/h2/h2*tau/2
                b[j] = -(2/tau + 2*bcoef/h2/h2)*tau/2
                c[j] = bcoef/h2/h2*tau/2
                d[j] = -(u1[i][j]/tau/2 + acoef/h1/h1*(u1[i+1][j]-2*u1[i][j]+u1[i-1][j]) + f(x[i],y[j],tau2))*tau/2
            a[-1] = -1/h2
            b[-1] = 1/h2
            c[-1] = 0
            d[-1] = -np.sin(x[i])*np.sin(mu*t[k])

            result2 = slau([a,b,c,d])
            #print(result2)
            for j in range(n):
                u2[i][j] = result2[j]
                u2[0][j] = 0
                u2[-1][j] = np.sin(y[j])*np.sin(t[k])

        for i in range(n):
            u2[i][-1] = -np.sin(x[i])*np.sin(mu*t[k])*h2 + u2[i][-2]

        for i in range(n):
            for j in range(n):
                u[i][j][k] = u2[i][j]
    return u

U1 = mpn(n)
fig = plt.figure()
ax = plt.axes(projection ='3d')
ax.plot_surface(x_plt,y_plt,np.array(U1[:,:,3]))
plt.show()


def f(x,y,t):
    return np.sin(x)*np.sin(y)*(mu*np.cos(mu*t)+(a+b)*np.sin(mu*t))

def progonka(a,b,c,d):
    P = [0]*n
    Q = [0]*n

    P[0] = -c[0]/b[0]
    for i in range(1,n):
        P[i] = -c[i] / (b[i] + a[i] * P[i - 1])

    Q[0] = d[0] / b[0]
    for i in range(1,n):
        Q[i] = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1])

    result = [0]*n
    result[n-1] = Q[n-1]
    for i in range(n-2,0):
        result[i] = P[i] * result[i + 1] + Q[i]

    return result

def slau(a_):

    a = a_[0]
    b = a_[1]
    c = a_[2]
    d = a_[3]

    for i in range(len(a)):
        if m.fabs(b[i]) < m.fabs(a[i]) + m.fabs(c[i]):
            raise Exception("LOL"+str(i))



    P,Q = [-c[0]/b[0]],[d[0]/b[0]]
    for i in range(1,len(a)):
        P.append(-c[i]/(b[i]+a[i]*P[i-1]))
        Q.append((d[i]-a[i]*Q[i-1])/(b[i]+a[i]*P[i-1]))

    x = [Q[-1]]
    for i in range(1,len(a)):
        x.append(P[len(a)-i-1]*x[i-1]+Q[len(a)-i-1])


    x = list(reversed(x))

    return x

