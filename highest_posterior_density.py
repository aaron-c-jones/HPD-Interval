import scipy as sp
import numpy as np
import math as m
import pandas as pd

def normal(x,mu=0.,sig=1.):
    return (1/m.sqrt(2*m.pi*m.pow(sig,2)))*m.exp((-1/(2*m.pow(sig,2)))*m.pow(x-mu,2))
    
def inv_normal(y,mu=0.,sig=1.):
    l=mu-m.sqrt((-2*m.pow(sig,2))*m.log(y*m.sqrt(2*m.pi*m.pow(sig,2))))
    u=mu+m.sqrt((-2*m.pow(sig,2))*m.log(y*m.sqrt(2*m.pi*m.pow(sig,2))))
    return l,u
    
def highest_posterior_density(func1,func2,desired=0.95):
    max_x=sp.optimize.fmin(lambda x: -func1(x),0)
    max_func=func1(max_x)
    
    search_grid=np.arange(0,max_func,0.000001)[::-1][:-1]

    bounds=list(map(func2,search_grid))
    
    rows_list=[]
    for i in range(len(search_grid)):
        l,u=bounds[i]
        area,error=sp.integrate.quad(func1,l,u)
        diff=abs(desired-area)
        rows_list.append([area,l,u,diff])
    df=pd.DataFrame(rows_list,columns=['Area','Lower Bound','Upper Bound','Difference'])
    
    result=df.iloc[df['Difference'].argmin()]
    return result
highest_posterior_density(func1=normal,func2=inv_normal)