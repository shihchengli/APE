# -*- coding: utf-8 -*-
import numpy as np
from math import sqrt
from math import pi, sin, cos

def IndefIntxPhimPhin(m,n,x,L,x_power):
    if m == 0 or n == 0:
        if m == n: #const*const
            y = pow(x,x_power+1)/(x_power+1)
            y /= (2*L)
        else:
            if n == 0:
                n = m
                m = 0
            if n%2 == 1: #const*cos
                n = (n+1)/2
                k = n*pi/L
                if x_power == 0:
                    y = sin(k*x)/k
                elif x_power == 1:
                    y = (k*x*sin(k*x)+cos(k*x))/pow(k,2)
                elif x_power == 2:
                    y = ((pow(k,2)*pow(x,2)-2)*sin(k*x)+2*k*x*cos(k*x))/pow(k,3)
                elif x_power == 3:
                    y = (k*x*(pow(k,2)*pow(x,2)-6)*sin(k*x)+3*(pow(k,2)*pow(x,2)-2)*cos(k*x))/pow(k,4)
                y /= L*sqrt(2)
            else: #const*sin
                n = n/2
                k = n*pi/L
                if x_power == 0:
                    y = -cos(k*x)/k
                elif x_power == 1:
                    y = (-k*x*cos(k*x)+sin(k*x))/pow(k,2)
                elif x_power == 2:
                    y = ((2-pow(k,2)*pow(x,2))*cos(k*x)+2*k*x*sin(k*x))/pow(k,3)
                elif x_power == 3:
                    y = (-k*x*(pow(k,2)*pow(x,2)-6)*cos(k*x)+3*(pow(k,2)*pow(x,2)-2)*sin(k*x))/pow(k,4)
                y /= L*sqrt(2)
    
    elif m%2 == 1 and n%2 == 1: #cos*cos
        m = (m+1)/2
        n = (n+1)/2
        k1 = m*pi/L
        k2 = n*pi/L
        if m == n:
            if x_power == 0:
                y = (2*k1*x+sin(2*k1*x))/(4*k1)
            elif x_power == 1:
                y = (2*k1*x*(k1*x+sin(2*k1*x))+cos(2*k1*x))/(8*pow(k1,2))
            elif x_power == 2:
                y = (4*pow(k1,3)*pow(x,3)+(6*pow(k1,2)*pow(x,2)-3)*sin(2*k1*x)+6*k1*x*cos(2*k1*x))/(24*pow(k1,3))
            elif x_power == 3:
                y = (2*pow(k1,4)*pow(x,4)+2*k1*x*(2*pow(k1,2)*pow(x,2)-3)*sin(2*k1*x)+(6*pow(k1,2)*pow(x,2)-3)*cos(2*k1*x))/(16*pow(k1,4))
        else:
            if x_power == 0:
                y = (k1*sin(k1*x)*cos(k2*x)-k2*cos(k1*x)*sin(k2*x))/(pow(k1,2)-pow(k2,2))
            elif x_power == 1:
                y = 0.5*(x*sin(x*(k1-k2))/(k1-k2)+x*sin(x*(k1+k2))/(k1+k2)\
                +cos(x*(k1-k2))/pow(k1-k2,2)+cos(x*(k1+k2))/pow(k1+k2,2))
            elif x_power == 2:
                y = 0.5*((pow(k1,2)*pow(x,2)-2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*sin(x*(k1-k2))/pow(k1-k2,3)\
                + (pow(k1,2)*pow(x,2)+2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*sin(x*(k1+k2))/pow(k1+k2,3)\
                + 2*x*cos(x*(k1-k2))/pow(k1-k2,2) + 2*x*cos(x*(k1+k2))/pow(k1+k2,2))
            elif x_power == 3:
                y = 0.5*(x*(pow(k1,2)*pow(x,2)-2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-6)*sin(x*(k1-k2))/pow(k1-k2,3)\
                + x*(pow(k1,2)*pow(x,2)+2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-6)*sin(x*(k1+k2))/pow(k1+k2,3)\
                + 3*(pow(k1,2)*pow(x,2)-2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*cos(x*(k1-k2))/pow(k1-k2,4)\
                + 3*(pow(k1,2)*pow(x,2)+2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*cos(x*(k1+k2))/pow(k1+k2,4))
        y /= L
    elif m%2 == 0 and n%2 == 0: #sin*sin
        m = m/2
        n = n/2
        k1 = m*pi/L
        k2 = n*pi/L
        if m == n:
            if x_power == 0:
                y = (2*k1*x-sin(2*k1*x))/(4*k1)
            elif x_power == 1:
                y = -(2*k1*x*(-k1*x+sin(2*k1*x))+cos(2*k1*x))/(8*pow(k1,2))
            elif x_power == 2:
                y = (4*pow(k1,3)*pow(x,3)+(-6*pow(k1,2)*pow(x,2)+3)*sin(2*k1*x)-6*k1*x*cos(2*k1*x))/(24*pow(k1,3));
            elif x_power == 3:
                y = (2*pow(k1,4)*pow(x,4)+2*k1*x*(-2*pow(k1,2)*pow(x,2)+3)*sin(2*k1*x)+(-6*pow(k1,2)*pow(x,2)+3)*cos(2*k1*x))/(16*pow(k1,4))
        else:
            if x_power == 0:
                y = (k2*sin(k1*x)*cos(k2*x)-k1*cos(k1*x)*sin(k2*x))/(pow(k1,2)-pow(k2,2))
            elif x_power == 1:
                y = 0.5*(x*sin(x*(k1-k2))/(k1-k2)-x*sin(x*(k1+k2))/(k1+k2)\
                +cos(x*(k1-k2))/pow(k1-k2,2)-cos(x*(k1+k2))/pow(k1+k2,2))
            elif x_power == 2:
                y = 0.5*((pow(k1,2)*pow(x,2)-2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*sin(x*(k1-k2))/pow(k1-k2,3)\
                - (pow(k1,2)*pow(x,2)+2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*sin(x*(k1+k2))/pow(k1+k2,3)\
                + 2*x*cos(x*(k1-k2))/pow(k1-k2,2) - 2*x*cos(x*(k1+k2))/pow(k1+k2,2))
            elif x_power == 3:
                y = 0.5*(x*(pow(k1,2)*pow(x,2)-2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-6)*sin(x*(k1-k2))/pow(k1-k2,3)\
                - x*(pow(k1,2)*pow(x,2)+2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-6)*sin(x*(k1+k2))/pow(k1+k2,3)\
                + 3*(pow(k1,2)*pow(x,2)-2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*cos(x*(k1-k2))/pow(k1-k2,4)\
                - 3*(pow(k1,2)*pow(x,2)+2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*cos(x*(k1+k2))/pow(k1+k2,4))
        y /= L
    else: #cos*sin
        if m%2 == 0:
            temp = m
            m = n
            n = temp
        m = (m+1)/2
        n = n/2
        k1 = m*pi/L
        k2 = n*pi/L
        if m == n:
            if x_power == 0:
                y = -pow(cos(k1*x),2)/(2*k1)
            elif x_power == 1:
                y = (sin(2*k1*x)-2*k1*x*cos(2*k1*x))/(8*pow(k1,2))
            elif x_power == 2:
                y = ((1-2*pow(k1,2)*pow(x,2))*cos(2*k1*x)+2*k1*x*sin(2*k1*x))/(8*pow(k1,3))
            elif x_power == 3:
                y = ((6*k1*x-4*pow(k1,3)*pow(x,3))*cos(2*k1*x)+3*(2*pow(k1,2)*pow(x,2)-1)*sin(2*k1*x))/(16*pow(k1,4))
        else:
            if x_power == 0:
                y = (k1*sin(k1*x)*sin(k2*x)+k2*cos(k1*x)*cos(k2*x))/(pow(k1,2)-pow(k2,2))
            elif x_power == 1:
                y = 0.5*(x*cos(x*(k1-k2))/(k1-k2)-x*cos(x*(k1+k2))/(k1+k2)\
                -sin(x*(k1-k2))/pow(k1-k2,2)+sin(x*(k1+k2))/pow(k1+k2,2))
            elif x_power == 2:
                y = 0.5*((pow(k1,2)*pow(x,2)-2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*cos(x*(k1-k2))/pow(k1-k2,3)\
                - (pow(k1,2)*pow(x,2)+2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*cos(x*(k1+k2))/pow(k1+k2,3)\
                - 2*x*sin(x*(k1-k2))/pow(k1-k2,2) + 2*x*sin(x*(k1+k2))/pow(k1+k2,2))
            elif x_power == 3:
                y = 0.5*(x*(pow(k1,2)*pow(x,2)-2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-6)*cos(x*(k1-k2))/pow(k1-k2,3)\
                - x*(pow(k1,2)*pow(x,2)+2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-6)*cos(x*(k1+k2))/pow(k1+k2,3)\
                - 3*(pow(k1,2)*pow(x,2)-2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*sin(x*(k1-k2))/pow(k1-k2,4)\
                + 3*(pow(k1,2)*pow(x,2)+2*k1*k2*pow(x,2)+pow(k2,2)*pow(x,2)-2)*sin(x*(k1+k2))/pow(k1+k2,4))
        y /= L
    return y

def IntxPhimPhin(m,n,x1,x2,L,x_power):
    return (IndefIntxPhimPhin(m,n,x2,L,x_power) - IndefIntxPhimPhin(m,n,x1,L,x_power))

def IntXPhimPhin(m,n,x1,x2,L,a):
    y = 0
    if (a[0]!=0): y += a[0]*IntxPhimPhin(m,n,x1,x2,L,0)
    if (a[1]!=0): y += a[1]*IntxPhimPhin(m,n,x1,x2,L,1)
    if (a[2]!=0): y += a[2]*IntxPhimPhin(m,n,x1,x2,L,2)
    if (a[3]!=0): y += a[3]*IntxPhimPhin(m,n,x1,x2,L,3)
    return y

# Fourier series