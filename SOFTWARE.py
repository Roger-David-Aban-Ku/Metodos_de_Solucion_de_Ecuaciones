#!/usr/bin/env python
# coding: utf-8

# # Materia: Metodos Numericos
# ### Ingenieria en Sistemas Computacionales
# ### Docente: MM. Jorger Pool Cel
# ### UNIDAD 2
# ### 4°A
# ## Integrantes del Equipo:
#  **Roger David Aban Ku**
#  
#  **Axel Salvador Aguilar Nuñez**
#  
#  **Jesus Armando Cabrera Piña**

# # METODO DE BISECCIÓN
# $$f(x)=\cos(x) - x$$
# **Intervalos** = $[0,\pi/2]$
# 
# **Tolerancia** = $0.0001$

# In[1]:


import math
def fx(x):
    return math.cos(x) - x

def xr(xi,xs):
    return (xi+xs)/2

def error(xanterior, xactual):
    return abs((xanterior - xactual)/xanterior)


# In[6]:


def biseccion(xi, xs, error):
    ecalculado=1
    xra=xr(xi, xs)
    print("\033[1;32m","El intervalo es [ ","\033[1;35m","     | Xi |","\t\033[1;31m","                 | Xs |","\t\033[1;34m","            | xra |","\t\033[1;32m","       | Error |]")
    print("\033[1;32m","El intervalo es [","\033[1;35m",xi,"\t\033[1;31m",xs,"\t\033[1;34m",xra,"]")
    k=1;
    while(ecalculado>error):
        #xra=xr(xi, xs)
        if((fx(xi)*fx(xra))<0):
            xs=xra
        if((fx(xi)*fx(xra))>0):
            xi=xra
        if ((fx(xi)*fx(xra))==0):
            ecalculado=error
        
        xold=xra
        xra=xr(xi, xs)
        ecalculado=abs((xold - xra)/xold)
        if ecalculado>=error:
            print("\033[1;32m","El intervalo es [","\033[1;35m",xi,"\t\033[1;31m",xs,"\t\033[1;34m",xra,"\t\033[1;32m",ecalculado,"]")
            k=k+1;
    
    print("")
    print("\t\033[1;30m",k," ITERACIONES para aproximar a: ",error,"; ","\033[1;32m",xra,"\033[1;30m","la aproximacion")


# In[7]:


xs=math.pi/2
biseccion(0,xs,0.0001)


# In[ ]:





# # METODO DE PUNTO FIJO
# $$f(x)=e^{-x} - x$$
# 
# $$e^{-x} - x = 0$$
# **Despejar x**
# $$x=e^{-x}$$
# **Aproximación Inicial** = $[0]$
# 
# **Tolerancia** = $0.001$

# In[8]:


import math

def fxdespejada(x):
    return math.e**(-x)

def errorabsoluto(xactual,xanterior):
    return abs(((xactual-xanterior)/xactual))


# In[9]:


def punto_fijo(x,cota):
    error=1
    k=1;
    print("\033[1;35m","           | xi | ","\t\033[1;31m","                   | Error | ")
    print("\033[1;35m","Raiz: ",x)
    while error>cota:
        xactual=fxdespejada(x)
        error=errorabsoluto(xactual,x)
        x=xactual
        if error>=cota:
            print("\033[1;35m","Raiz: ",xactual,"\t\033[1;31m", "   Error: ",error)
            k=k+1;
    
    print("")
    print("\t\033[1;30m",k," ITERACIONES para aproximar a: ",cota,"; ","\033[1;32m",xactual,"\033[1;30m","la aproximacion")


# In[10]:


punto_fijo(0,0.001)


# # METODO DE NEWTON RAPHSON
# $$f(x)=3(x+1)^{-1} - 1.8$$
# 
# $$f'(x)=-3(x+1)^{-2}$$
# **Aproximación Inicial** = $[0]$
# 
# **Tolerancia** = $0.001$

# In[12]:


def fx(x):
    return 3*(x+1)**(-1) - 1.8

def derivadafx(x):
    return -3*(x+1)**(-2)

def Error(xactual, xanterior):
    return abs(xactual-xanterior)


# In[13]:


def newton(x, error):
    ecalculado=1
    k=1;
    print("\033[1;35m","           | xi | ","\t\033[1;31m","                   | Error | ")
    print("\033[1;35m","Raiz: ",x)
    while(ecalculado>error):
        
        xactual=x-(fx(x)/derivadafx(x))
        
        ecalculado=Error(xactual,x)
        
        x=xactual
        if ecalculado>=error:
            print("\033[1;35m","Raiz: ",xactual,"\t\033[1;31m", "   Error: ",ecalculado)
            k=k+1;
        
    print("")
    print("\t\033[1;30m",k," ITERACIONES para aproximar a: ",error,"; ","\033[1;32m",xactual,"\033[1;30m","la aproximacion")


# In[14]:


newton(0,0.001)


# In[ ]:





# # METODO DE LA SECANTE
# $$f(x)=\cos(x+1) - \sin(x+1) + 0.8$$
# **Aproximación Inicial** = $[0]$
# 
# **Tolerancia** = $0.001$

# In[16]:


import math

def funcionx(x):
    return math.cos(x+1) - math.sin(x+1) + 0.8
def fun(xia, xi):
    return (xi-((funcionx(xi) * (xia - xi)) / (funcionx(xia) - funcionx(xi))))

def error(xiactual, xianterior):
    return abs((xiactual-xianterior)/xiactual)


# In[17]:


def secante(xia, xi, tol):
    xinew=fun(xia, xi)
    ecal=error(xinew,xi)
    print("\033[1;35m","(xi-1): ",xia,"\t\033[1;31m","  (xi): ",xi,"\t\033[1;32m","   (Error): ",ecal)
    
    k=1;
    
    while(ecal>tol):
        
        xinuevo=fun(xia, xi)
        
        xia=xi
        xi=xinuevo
        ecal=error(xi,xia)
        
        if ecal>=tol:
            print("\033[1;35m","(xi-1): ",xia,"\t\033[1;31m","  (xi): ",xi,"\t\033[1;32m","   (Error): ",ecal)
            k=k+1;
        
    print("")
    print("\t\033[1;30m",k," ITERACIONES para aproximar a: ",tol,"; ","\033[1;32m",xia,"\033[1;30m","la aproximacion")


# In[18]:


secante(0,1,0.001)

