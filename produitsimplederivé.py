#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 23:33:54 2024

@author: thomasroujou
"""

import numpy as np
from scipy.stats import norm
from math import *

import math
from scipy.stats import norm

# Black-Scholes pour call
def black_scholes_call(S, K, T, r, sigma):
    d1 = (math.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)
    return S * norm.cdf(d1) - K * math.exp(-r * T) * norm.cdf(d2)

# Black-Scholes pour put
def black_scholes_put(S, K, T, r, sigma):
    d1 = (math.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)
    return K * math.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)

# Bull Call Spread
def bull_call_spread(S, K1, K2, T, r, sigma):
    """
    Stratégie : Acheter un call à K1, vendre un call à K2
    """
    call_long = black_scholes_call(S, K1, T, r, sigma)
    call_short = black_scholes_call(S, K2, T, r, sigma)
    return call_long - call_short

# Bear Put Spread
def bear_put_spread(S, K1, K2, T, r, sigma):
    """
    Stratégie : Acheter un put à K2, vendre un put à K1
    """
    put_long = black_scholes_put(S, K2, T, r, sigma)
    put_short = black_scholes_put(S, K1, T, r, sigma)
    return put_long - put_short

# Box Spread
def box_spread(S, K1, K2, T, r, sigma):
    bull_call = bull_call_spread(S, K1, K2, T, r, sigma)
    bear_put = bear_put_spread(S, K1, K2, T, r, sigma)
    return bull_call + bear_put

# Calendar Spread
def calendar_spread(S, K, T1, T2, r, sigma):
    """
    Stratégie : Acheter une option à long terme (T2), vendre une option à court terme (T1)
    """
    call_long = black_scholes_call(S, K, T2, r, sigma)
    call_short = black_scholes_call(S, K, T1, r, sigma)
    return call_long - call_short

# Exemple d'utilisation
if __name__ == "__main__":
    # Paramètres
    S = 100       # Prix actuel du sous-jacent
    K1 = 95       # Strike 1 (bas)
    K2 = 105      # Strike 2 (haut)
    T = 1         # Maturité (en années)
    T1 = 0.5      # Maturité courte (calendar spread)
    T2 = 1        # Maturité longue (calendar spread)
    r = 0.05      # Taux sans risque (5%)
    sigma = 0.2   # Volatilité (20%)
    
    # Calculs
    bull_call_price = bull_call_spread(S, K1, K2, T, r, sigma)
    bear_put_price = bear_put_spread(S, K1, K2, T, r, sigma)
    box_spread_price = box_spread(S, K1, K2, T, r, sigma)
    calendar_price = calendar_spread(S, K1, T1, T2, r, sigma)
    
    # Résultats
    print(f"Bull Call Spread Price: {bull_call_price:.2f}")
    print(f"Bear Put Spread Price: {bear_put_price:.2f}")
    print(f"Box Spread Price: {box_spread_price:.2f}")
    print(f"Calendar Spread Price: {calendar_price:.2f}")


# Paramètres communs
S = 100      # Prix actuel du sous-jacent
T = 1        # Temps jusqu'à l'échéance (1 an)
r = 0.05     # Taux sans risque
sigma = 0.2  # Volatilité

# ** Vertical Spread (Bull Call Spread) **
K1 = 95  # Strike du call acheté
K2 = 105 # Strike du call vendu

# Calcul des prix des options
call_K1 = black_scholes_call(S, K1, T, r, sigma)
call_K2 = black_scholes_call(S, K2, T, r, sigma)

# Prix du vertical spread
vertical_spread_price = call_K1 - call_K2
print(f"Prix du Bull Call Spread : {vertical_spread_price:.2f}")

# ** Butterfly Spread **
K3 = 80 # Strike du premier call acheté
K4 = 110 # Strike du second call acheté

# Calcul des prix des options pour le butterfly spread
call_K3 = black_scholes_call(S, K3, T, r, sigma)
call_K4 = black_scholes_call(S, K4, T, r, sigma)

# Prix du butterfly spread
butterfly_spread_price = call_K3 - 2 * call_K1 + call_K4
print(f"Prix du Butterfly Spread : {butterfly_spread_price:.2f}")

import numpy as np
import matplotlib.pyplot as plt

# Bull Call Spread Payoff
def bull_call_spread_payoff(S, K1, K2, premium1, premium2):
    payoff = np.maximum(S - K1, 0) - premium1 - (np.maximum(S - K2, 0) - premium2)
    return payoff

# Bear Put Spread Payoff
def bear_put_spread_payoff(S, K1, K2, premium1, premium2):
    payoff = np.maximum(K2 - S, 0) - premium1 - (np.maximum(K1 - S, 0) - premium2)
    return payoff

# Box Spread Payoff
def box_spread_payoff(S, K1, K2, premium):
    payoff = np.where(S < K1, 0, np.where(S < K2, S - K1, K2 - K1)) - premium
    return payoff

# Calendar Spread Payoff (simplifié pour visualisation)
def calendar_spread_payoff(S, K, premium1, premium2):
    payoff = (np.maximum(S - K, 0) - premium2) - (np.maximum(S - K, 0) - premium1)
    return payoff

# Visualisation des payoffs
def plot_payoff():
    # Paramètres
    S = np.linspace(50, 150, 500)  # Prix du sous-jacent
    K1 = 80
    K2 = 120
    premium_call1 = 4  # Prime du call long
    premium_call2 = 2 # Prime du call short
    premium_put1 = 4  # Prime du put long
    premium_put2 = 2  # Prime du put short
    
    # Calcul des payoffs
    payoff_bull_call = bull_call_spread_payoff(S, K1, K2, premium_call1, premium_call2)
    payoff_bear_put = bear_put_spread_payoff(S, K1, K2, premium_put1, premium_put2)
    payoff_box = payoff_bull_call + payoff_bear_put
    payoff_calendar = calendar_spread_payoff(S, K1, premium_call1, premium_call2)
    
    # Tracés des graphiques
    plt.figure(figsize=(12, 8))
    plt.plot(S, payoff_bull_call, label='Bull Call Spread', color='blue')
    plt.plot(S, payoff_bear_put, label='Bear Put Spread', color='red')
    plt.plot(S, payoff_box, label='Box Spread', color='green')
    plt.plot(S, payoff_calendar, label='Calendar Spread', color='purple')
    
    # Légendes et personnalisation
    plt.axhline(0, color='black', linestyle='--', linewidth=0.7)
    plt.title("Évolution des Payoffs - Différentes Stratégies d'Options")
    plt.xlabel("Prix du Sous-jacent (S)")
    plt.ylabel("Payoff (€)")
    plt.legend()
    plt.grid()
    plt.show()

# Afficher les graphiques
plot_payoff()
