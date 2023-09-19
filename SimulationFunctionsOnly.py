# simulation_functions.py

import numpy as np

def MMConstant(kf, kb, kcat):
    kM = (kb + kcat) / kf
    return kM 

def Full_Model(x, t, kf, kb, kcat, St, Et):
    C, P = x
    dC_dt = kf * (St - C - P) * (Et - C) - kb * C - kcat * C
    dP_dt = kcat * C
    return dC_dt, dP_dt

def QSSA(P, t, kcat, kM, St, Et):
    S = St - P
    dP_dt = (kcat * Et * S) / (S + kM)
    return dP_dt

def tQSSA(P, t, kcat, kM, St, Et):
    S_hat = St - P
    square_Root_arg = (Et + S_hat + kM) ** 2 - (4 * Et * S_hat)
    dP_dt = (2 * kcat * Et * S_hat) / (Et + S_hat + kM + np.sqrt(square_Root_arg))
    return dP_dt

# SWITCH ----------------------------------------------------------------
# -----------------------------------------------------------------------

def Constants_Switch(kfe, kbe, ke, kfd, kbd, kd):
     kME = (kbe + ke) / kfe
     kMD = (kbd + kd) / kfd
     return kME, kMD

def Full_Model_Switch(x, t, St, Et, Dt, kfe, kbe, ke, kfd, kbd, kd):
    Sp, C, Cp = x
    dSp_dt = -kfd * (Dt - Cp) * Sp + kbd * Cp + ke * C
    dC_dt = kfe * (Et - C) * (St - Sp - C - Cp) - kbe * C - ke * C
    dCp_dt = kfd * (Dt - Cp) * Sp - kbd * Cp - kd * Cp
    return dSp_dt, dC_dt, dCp_dt

def QSSA_Switch(Sp, t, St, Et, Dt, ke, kd, kME, kMD):
    S = St - Sp
    dSp_dt = (ke * Et * S) / (S + kME) - (kd * Dt * Sp) / (Sp + kMD)
    return dSp_dt

def deltaFunctions(Sp_hat, St, Et, Dt, kME, kMD):
    S_hat = St - Sp_hat
    delta = (Et + S_hat + kME)**2 - (4 * Et * S_hat)
    deltap = (Dt + Sp_hat + kMD)**2 - (4 * Dt * Sp_hat)
    return delta, deltap

def tQSSA_Switch(Sp_hat, t, St, Et, Dt, ke, kd, kME, kMD):
    S_hat = St - Sp_hat
    delta, deltap = deltaFunctions(Sp_hat, St, Et, Dt, kME, kMD)
    dSp_dt = (2 * ke * Et * S_hat) / (Et + S_hat + kME + np.sqrt(delta))
    dSp_dt = dSp_dt - (2 * kd * Dt * Sp_hat) / (Dt + Sp_hat + kMD + np.sqrt(deltap))
    return dSp_dt
    

  


