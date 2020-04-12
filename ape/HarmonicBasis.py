# -*- coding: utf-8 -*-
import sys
import numpy as np
from numpy import exp
from scipy.special import erf
from math import factorial as fact
from math import pi, sqrt

def Hermite(m,x):
    if m < 0:
        value = 0
    elif m == 0:
        value = 1
    elif m == 1:
        value = 2*x
    else:
        i = 1
        P1 = 2*x
        P2 = 1
        while True:
            i += 1
            value = 2*x*P1-2*(i-1)*P2
            P2 = P1
            P1 = value
            if i == m:
                break
    return value

def IndefInt(m,n,x):
    x_max = np.sqrt(np.log(sys.float_info[0]))
    value = 0
    if x > x_max:
        if m == n:
            value += pow(2.,n)*fact(n)*sqrt(pi)/2
            
    elif x < -x_max:
        if m == n:
            value -= pow(2.,n)*fact(n)*sqrt(pi)/2

    else:
        if m != n:
            value = exp(-pow(x,2))*(n*Hermite(m,x)*Hermite(n-1,x)-m*Hermite(n,x)*Hermite(m-1,x))/(m-n)
        else:
            value += pow(2.,n)*sqrt(pi)/2*erf(x)
            for i in range(n):
                value -= exp(-pow(x,2))*Hermite(i+1,x)*Hermite(i,x)*pow(2.,n-1-i)/fact(i+1)
            value *= fact(n)
    return value

def IntHmHnexp(m,n,x1,x2):
    return (IndefInt(m,n,x2)- IndefInt(m,n,x1))

def IntXHmHnexp(m,n,x1,x2,a):
    y=0
    if (a[0]!=0): y += a[0]*IntHmHnexp(m,n,x1,x2)
    if (a[1]!=0): y += a[1]*(1/2.0*IntHmHnexp(m,n+1,x1,x2)+n*IntHmHnexp(m,n-1,x1,x2))
    if (a[2]!=0): y += a[2]*(1/4.0*IntHmHnexp(m,n+2,x1,x2)+(n+1/2.0)*IntHmHnexp(m,n,x1,x2)\
            +n*(n-1)*IntHmHnexp(m,n-2,x1,x2))
    if (a[3]!=0): y += a[3]*(1/8.0*IntHmHnexp(m,n+3,x1,x2)+3/4.0*(n+1)*IntHmHnexp(m,n+1,x1,x2)\
            +3/2.0*n*n*IntHmHnexp(m,n-1,x1,x2) + n*(n-1)*(n-2)*IntHmHnexp(m,n-3,x1,x2))
    return y
