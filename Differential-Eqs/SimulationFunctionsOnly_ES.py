# simulation_functions.py

import numpy as np


def Constants (kf, kb, kcat):
    kM = (kb + kcat) / kf
    return kf, kb, kcat, kM 

def Full_Model (reactancts, t, kf, kb, kcat, St, Et):
    [C, P] = reactancts
    dC_dt = kf * (St - C - P) * (Et - C) - kb * C - kcat * C
    dP_dt = kcat * C
    return [dC_dt, dP_dt]

def QSSA (P, t, kcat, kM, St, Et):
    S = St - P
    dP_dt = (kcat * Et * S) / (S + kM)
    return dP_dt

def tQSSA (P, t, kcat, kM, St, Et):
    S_hat = St - P
    square_Root_arg = (Et + S_hat + kM) ** 2 - (4 * Et * S_hat)
    dP_dt = (2 * kcat * Et * S_hat) / (Et + S_hat + kM + np.sqrt(square_Root_arg))
    return dP_dt

# SWITCH ----------------------------------------------------------------
# -----------------------------------------------------------------------

#def compute_Sphat_over_St(E, D, ke, kd, Sh, kME, kMD, Delta, Delta_p):
    numerator = 2 * ke * E * Sh * D + 2 * ke * E * Sh**2 + 2 * ke * E * Sh * kMD + 2 * ke * E * np.sqrt(Delta_p) - 2 * kd * D * Sh * E - 2 * kd * D * Sh**2 - 2 * kd * D * Sh * kME - 2 * kd * D * np.sqrt(Delta)
    denominator = E + Sh + kME + np.sqrt(Delta)
    Sphat_over_St = (numerator / denominator) * (ke * E / (kd * D))
    return Sphat_over_St

def Constants_Switch (ke, kbe, kfe, kbd, kfd, kd):
     kME = (kbe + ke) / kfe
     kMD = (kbd + kd) / kfd
     #kD = (kMD * kfd) - kbd
     return kME, kMD

def Full_Model_Switch (reactantcts, t, St, Et, Dt, kfd, kbd, ke, kfe, kbe, kd):
    [Sp, C, Cp] = reactantcts
    dSp_dt = -kfd * (Dt - Cp) * Sp + kbd * Cp + ke * C
    dC_dt = kfe * (Et - C) * (St - Sp - C - Cp) - kbe * C - ke * C
    dCp_dt = kfd * (Dt - Cp) * Sp - kbd * Cp - kd * Cp
    return [dSp_dt, dC_dt, dCp_dt]

def QSSA_Switch (Sp, t, St, Et, Dt, ke, kd, kME, kMD):
    S = St - Sp
    dSp_dt = (ke * Et * S) / (S + kME) - (kd * Dt * Sp) / (Sp + kMD)
    return dSp_dt

def deltaFunction (Sp_hat, St, Et, Dt, kME, kMD):
    S_hat = St - Sp_hat
    delta = (Et + S_hat + kME)**2 - (4 * Et * S_hat)
    deltap = (Dt + Sp_hat + kMD)**2 - (4 * Dt * Sp_hat)
    return [delta, deltap]


def tQSSA_Switch (Sp_hat, t, St, Et, Dt, ke, kd, kME, kMD):
    S_hat = St - Sp_hat
    delta = deltaFunction(Sp_hat, St, Et, Dt, kME, kMD)[0]
    deltap = deltaFunction(Sp_hat, St, Et, Dt, kME, kMD)[1]
    dSp_dt = (2 * ke * Et * S_hat) / (Et + S_hat + kME + np.sqrt(delta))
    dSp_dt = dSp_dt - (2 * kd * Dt * Sp_hat) / (Dt + Sp_hat + kMD + np.sqrt(deltap))
    return dSp_dt


# PROPENSITY FUNCTIONS - Switch 

# QSSA ------------------------

def propensity_function_AFE_QSSA (Sp, St, Et, C, Cp, kfe):
    afe = kfe * (Et - C) * (St - Sp - C - Cp)
    return afe

def propensity_function_ABE_QSSA (C, kbe):
    abe = kbe * C 
    return abe

def propensity_function_AE_QSSA (C, ke):
    ae = ke * C 
    return ae

def propensity_function_AFD_QSSA (Sp, Dt, Cp, kfd):
    afd = kfd * (Dt - Cp) * Sp
    return afd

def propensity_function_ABD_QSSA (C, kbd):
    abd = kbd * C 
    return abd

def propensity_function_AD_QSSA (C, kd):
    ad = kd * C 
    return ad

# tQSSA -----------------------

def propensity_function_AE_tQSSA (Sp, St, Et, ke, kME):
    S_hat = St - Sp
    delta = (Et + S_hat + kME)**2 - (4 * Et * S_hat)
    ae = (2 * ke * Et * S_hat) / (Et + S_hat + kME + np.sqrt(delta))
    return ae

def propensity_function_AD_tQSSA (Sp, St, Et, Dt, kd, kME, kMD):
    S_hat = St - Sp
    Sp_hat = St + S_hat
    delta = (Et + S_hat + kME)**2 - (4 * Et * S_hat)
    ad = (2 * kd * Dt * Sp_hat) / (Dt + Sp_hat + kMD + np.sqrt(delta))
    return ad
    
    
# TEST Functions --------------------------------------------------------
# -----------------------------------------------------------------------

# CME Simulations

def weighted_ave_sd(arr, weights):
	average = np.average(arr, weights=weights)
	variance = np.average((arr-average)**2, weights=weights)
	return average, np.sqrt(variance)
		



    

  


