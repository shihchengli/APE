# -*- coding: utf-8 -*-
import sys
import numpy as np
from numpy import exp
from scipy.special import erf
from math import factorial as fact
from math import pi, sqrt

# fix overflow problem
from decimal import Decimal as D
from decimal import getcontext
getcontext().prec = 15

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
    return D(value)

def IndefInt(m,n,x,f):
    f1 = f[0]
    f2 = f[1]
    f3 = f[2]
    x_max = np.sqrt(np.log(sys.float_info[0]))
    value = D(0)
    if x > x_max:
        if m == n:
            value += D(pow(2.,n))/f1*D(fact(n))*D(sqrt(pi)/2)/f2/f3
            
    elif x < -x_max:
        if m == n:
            value -= D(pow(2.,n))/f1*D(fact(n))*D(sqrt(pi)/2)/f2/f3

    else:
        if m != n:
            value = D(exp(-pow(x,2)))*(D(n)*Hermite(m,x)/f1*Hermite(n-1,x)/f2-D(m)*Hermite(n,x)/f1*Hermite(m-1,x)/f2)/D(m-n)/f3
        else:
            value += D(pow(2.,n)*sqrt(pi)/2)/f1*D(erf(x))/f2/f3
            for i in range(n):
                value -= D(exp(-pow(x,2)))*Hermite(i+1,x)/f1*Hermite(i,x)/f2*D(pow(2.,n-1-i))/D(fact(i+1))/f3
            value *= D(fact(n))
    return value

def IntHmHnexp(m,n,x1,x2,f):
    return (IndefInt(m,n,x2,f)- IndefInt(m,n,x1,f))

def IntXHmHnexp(m,n,x1,x2,a,f):
    y=0
    if (a[0]!=0): y += D(a[0])*IntHmHnexp(m,n,x1,x2,f)
    if (a[1]!=0): y += D(a[1])*(D(1/2.0)*IntHmHnexp(m,n+1,x1,x2,f)+D(n)*IntHmHnexp(m,n-1,x1,x2,f))
    if (a[2]!=0): y += D(a[2])*(D(1/4.0)*IntHmHnexp(m,n+2,x1,x2,f)+D(n+1/2.0)*IntHmHnexp(m,n,x1,x2,f)\
            +D(n*(n-1))*IntHmHnexp(m,n-2,x1,x2,f))
    if (a[3]!=0): y += D(a[3])*(D(1/8.0)*IntHmHnexp(m,n+3,x1,x2,f)+D(3/4.0*(n+1))*IntHmHnexp(m,n+1,x1,x2,f)\
            +D(3/2.0*n*n)*IntHmHnexp(m,n-1,x1,x2,f) + D(n*(n-1)*(n-2))*IntHmHnexp(m,n-3,x1,x2,f))
    return y