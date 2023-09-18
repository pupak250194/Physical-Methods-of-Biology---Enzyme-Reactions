import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# COSTANTI

def Constants (kf, kb, kcat):

#kf = 10
#kb = 9
#kcat = 1

    kM = (kb + kcat) / kf
    St = kM
    Et = 0.1*(kM + St)

    return kf, kb, kcat, kM, St, Et

# SISTEMA DI EDO RIDOTTO - C NON TRASCURABILE

def Full_Model(reactancts, t, kf, kb, kcat, St, Et):

    [C, P] = reactancts
    
    # Call the Constants function and pass the parameters

    dC_dt = kf * (St - C - P) * (Et - C) - kb * C - kcat * C
    dP_dt = kcat * C

    return [dC_dt, dP_dt]

# QSSA - DA VERIFICARE LE CONDIZIONI DENTRO AL SISTEMA IN PARTICLARE LA TRASCURABILITA' DI C

def QSSA(P, t, kcat, kM, St, Et):

    S =  St - P

    dP_dt = (kcat * Et * S) / (S + kM)

    return dP_dt

# tQSSA - DA VEDERE SE SI PUO' SCRIVERE MEGLIO

def tQSSA (P, t,  kcat, kM, St, Et):

    S_hat = St - P
    square_Root_arg = (Et + S_hat + kM)**2 - (4 * Et * S_hat)

    dP_dt = (2 * kcat * Et * S_hat) / (Et + S_hat + kM + np.sqrt(square_Root_arg))

    return dP_dt

# CONDIZIONI INIZIALI

# Collect user input for kf
inputkf = float(input('Enter a value for kf: '))

# Collect user input for kb
inputkb = float(input('Enter a value for kb: '))

# Collect user input for kcat
inputkcat = float(input('Enter a value for kcat: '))

# Collect user input for C
inputC = float(input('Enter a value for the initial concentration of C: '))

# Collect user input for P
inputP = float(input('Enter a value for the initial concentration of P: '))
 

Globkf, Globkb, Globkcat, GlobkM, GlobSt, GlobEt = Constants(inputkf, inputkb, inputkcat)

C_initial = inputC
P_initial = inputP

initial_conditions = [C_initial, P_initial]

# INTERVALLO TEMPORALE

time_span = np.linspace(0, 50, 1000)

# RISOLUZIONE ODE

Full_Model_Solution = odeint(Full_Model, initial_conditions, time_span, args=(Globkf, Globkb, Globkcat, GlobSt, GlobEt))
QSSAsolution = odeint(QSSA, initial_conditions, time_span, args=(Globkcat, GlobkM, GlobSt, GlobEt))
tQSSAsolution = odeint(tQSSA, initial_conditions, time_span, args=(Globkcat, GlobkM, GlobSt, GlobEt))


Full_Model_P_values = Full_Model_Solution[:, 1]
QSSA_P_values = QSSAsolution[:, 1]
tQSSA_P_values = tQSSAsolution[:, 1]


# Calcolo del rapporto P/S con i tre metodi

Full_Model_ratio_values = Full_Model_P_values / GlobSt
QSSA_ratio_values = QSSA_P_values / GlobSt
tQSSA_ratio_values = tQSSA_P_values / GlobSt

# GRAFICI

plt.plot (time_span, Full_Model_ratio_values, label ='Full Model', color = 'yellow', linewidth = 3.5)
plt.plot (time_span, QSSA_ratio_values, label = 'QSSA - P/Stot', linewidth = 2)
plt.plot (time_span, tQSSA_ratio_values, label = 'tQSSA - P/S-hat', linestyle = '--', color = 'black')

plt.xlabel('Time')
plt.ylabel('P/S')
plt.legend()
plt.show()
