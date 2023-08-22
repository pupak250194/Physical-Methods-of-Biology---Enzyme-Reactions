import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint



# COSTANTI

kf = 10
kb = 9
kcat = 1

kM = (kb + kcat) / kf

St = kM
Et = 0.1*(kM + St)

# SISTEMA DI EDO RIDOTTO - C NON TRASCURABILE

def Full_Model(reactancts, t):

    [C, P] = reactancts

    dC = kf * (St - C - P) * (Et - C) - kb * C - kcat * C
    dP = kcat * C

    return [dC, dP]

# QSSA - DA VERIFICARE LE CONDIZIONI DENTRO AL SISTEMA IN PARTICLARE LA TRASCURABILITA' DI C

def QSSA(P, t):

    S =  St - P

    dP_dt = (kcat * Et * S) / (S + kM)

    return dP_dt

# tQSSA - DA VEDERE SE SI PUO' SCRIVERE MEGLIO

def tQSSA (P, t):

    S_hat = St - P
    square_Root_arg = (Et + S_hat + kM)**2 - (4 * Et * S_hat)

    dP_dt = (2 * kcat * Et * S_hat) / (Et + S_hat + kM + np.sqrt(square_Root_arg))
    
    return dP_dt

# CONDIZIONI INIZIALI

C_initial = 0
P_initial = 0

initial_conditions = [C_initial, P_initial]

# INTERVALLO TEMPORALE

time_span = np.linspace(0, 50, 1000)

# RISOLUZIONE ODE

Full_Model_Solution = odeint (Full_Model, initial_conditions, time_span)
QSSAsolution = odeint (QSSA, P_initial, time_span)
tQSSAsolution = odeint (tQSSA, P_initial, time_span)

Full_Model_P_values = Full_Model_Solution[:, 1]
QSSA_P_values = QSSAsolution[:, 0]
tQSSA_P_values = tQSSAsolution[:, 0]


# Calcolo del rapporto P/S con i tre metodi

Full_Model_ratio_values = Full_Model_P_values / St
QSSA_ratio_values = QSSA_P_values / St
tQSSA_ratio_values = tQSSA_P_values / St

# GRAFICI

plt.plot (time_span, Full_Model_ratio_values, label ='Full Model - P/Stot', color = 'yellow', linewidth = 3.5)
plt.plot (time_span, QSSA_ratio_values, label = 'QSSA - P/Stot', linewidth = 2)
plt.plot (time_span, tQSSA_ratio_values, label = 'tQSSA - P/Stot', linestyle = '--', color = 'black')

plt.xlabel('Time')
plt.ylabel('P/Stot')
plt.legend()
plt.show()
