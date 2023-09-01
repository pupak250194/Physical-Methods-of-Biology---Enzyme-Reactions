# simulation_functions.py

import numpy as np

def Constants(kf, kb, kcat):
    kM = (kb + kcat) / kf
    return kf, kb, kcat, kM

def Full_Model(reactancts, t, kf, kb, kcat, St, Et):

    [C, P] = reactancts
    dC_dt = kf * (St - C - P) * (Et - C) - kb * C - kcat * C
    dP_dt = kcat * C
    return [dC_dt, dP_dt]

def QSSA(P, t, kcat, kM, St, Et):
    S = St - P
    dP_dt = (kcat * Et * S) / (S + kM)
    return dP_dt

def tQSSA(P, t, kcat, kM, St, Et):
    S_hat = St - P
    square_Root_arg = (Et + S_hat + kM) ** 2 - (4 * Et * S_hat)
    dP_dt = (2 * kcat * Et * S_hat) / (Et + S_hat + kM + np.sqrt(square_Root_arg))
    return dP_dt