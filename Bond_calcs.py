# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 13:25:37 2023

@author: Gilberto
"""

""" Get yield-to-maturity of a bond """
import scipy.optimize as optimize


def bond_ytm(price, par, T, coup, freq=2, guess=0.05):
    freq = float(freq)
    periods = T * freq
    coupon = coup / 100. * par / freq
    dt = [(i + 1) / freq for i in range(int(periods))]
    ytm_func = lambda y: \
        sum([coupon / (1 + y / freq) ** (freq * t) for t in dt]) + \
        par / (1 + y / freq) ** (freq * T) - price

    return optimize.newton(ytm_func, guess)
 
    

ytm = bond_ytm(95.0428, 100, 1.5, 5.75, 2)
print (ytm)

""" Get bond price from YTM """
def bond_price(par, T, ytm, coup, freq=2):
     freq = float(freq)
     periods = T*freq
     coupon = coup/100.*par/freq
     dt = [(i+1)/freq for i in range(int(periods))]
     price = sum([coupon/(1+ytm/freq)**(freq*t) for t in dt]) + \
         par/(1+ytm/freq)**(freq*T)
     return price
 
    
##Bond Pricing##
bond_price(100, 1.5, ytm, 5.75, 2)


def bond_mod_duration(price, par, T, coup, freq, dy=0.01):
     ytm = bond_ytm(price, par, T, coup, freq)
    
     ytm_minus = ytm - dy 
     price_minus = bond_price(par, T, ytm_minus, coup, freq)
    
     ytm_plus = ytm + dy
     price_plus = bond_price(par, T, ytm_plus, coup, freq)
    
     mduration = (price_minus-price_plus)/(2*price*dy)
     return mduration 
 
    

duration = bond_mod_duration(95.04, 100, 1.5, 5.75, 2, 0.01)   
print(duration) #Modified bond duration#

""" Calculate convexity of a bond """

def bond_convexity(price, par, T, coup, freq, dy=0.01):
    ytm = bond_ytm(price, par, T, coup, freq)
    ytm_minus = ytm - dy
    price_minus = bond_price(par, T, ytm_minus, coup, freq)
    
    ytm_plus = ytm + dy
    price_plus = bond_price(par, T, ytm_plus, coup, freq)
    
    convexity = (price_minus+price_plus-2*price)/(price*dy**2)
    return convexity
    
convexity = bond_convexity(95.0428, 100, 1.5, 5.75, 2)
print(convexity)    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 