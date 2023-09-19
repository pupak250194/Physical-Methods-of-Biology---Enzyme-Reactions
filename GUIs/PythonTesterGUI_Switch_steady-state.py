import tkinter as tk
import matplotlib.pyplot as plt
import sys
import numpy as np

from tkinter import ttk
from scipy.optimize import fsolve
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import SimulationFunctionsOnly # assuming Function file is in the GUI folder, otherwise use sys.path()...


def start_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    DT = float(entry_DT.get())
    ST = float(entry_ST.get())

    initial_SP = 0
    initial_C  = 0
    initial_CP = 0

    initial_conditions = [initial_SP, initial_C, initial_CP]


    # Getting the user input from the entry widgets and setting constants

    kfe = float(entry_kfe.get())  # kfe value
    kbe = float(entry_kbe.get())  # kbe value
    ke  = float(entry_ke.get())   # ke value
    kfd = float(entry_kfd.get())  # kfd value
    kbd = float(entry_kbd.get())  # kbd value
    kd  = float(entry_kd.get())   # kd value

    kME, kMD = SimulationFunctionsOnly.Constants_Switch(kfe, kbe, ke, kfd, kbd, kd)
    
    # SETTING TIME SPAN

    # Getting the user input for maximum ke ET / kd DT ratio
    max_ratio = float(entry_max_ratio.get())

    # Creating the ke ET / kd DT ratio span
    ratio_span = np.linspace(0, max_ratio, 1000)

    # Creating the ET span
    ET_span = ratio_span * kd * DT / ke

    # solving odes

    Full_Model_Solution = []
    QSSAsolution = []
    tQSSAsolution = []
    for ET in ET_span:
        Full_Model_Solution.append(fsolve(SimulationFunctionsOnly.Full_Model_Switch, initial_conditions, args=(0, ST, ET, DT, kfe, kbe, ke, kfd, kbd, kd)))
        QSSAsolution.append(fsolve(SimulationFunctionsOnly.QSSA_Switch, initial_SP, args=(0, ST, ET, DT, ke, kd, kME, kMD)))
        tQSSAsolution.append(fsolve(SimulationFunctionsOnly.tQSSA_Switch, initial_SP+initial_CP, args=(0, ST, ET, DT, ke, kd, kME, kMD)))

    # convert lists to numpy arrays

    Full_Model_Solution = np.array(Full_Model_Solution)
    QSSAsolution = np.array(QSSAsolution)
    tQSSAsolution = np.array(tQSSAsolution)

    Full_Model_SP_hat = Full_Model_Solution[:, 0] + Full_Model_Solution[:, 2]
    QSSA_SP_hat = QSSAsolution[:, 0]
    tQSSA_SP_hat = tQSSAsolution[:, 0]

    # Ratio SP_hat/ST according to the 3 models

    Full_Model_SP_hat_ST = Full_Model_SP_hat / ST
    QSSA_SP_hat_ST = QSSA_SP_hat / ST
    tQSSA_SP_hat_ST = tQSSA_SP_hat / ST


    # Plots

    # Create a new figure for the plot
    fig, ax = plt.subplots()

    # Plot the data
    ax.plot(ratio_span, Full_Model_SP_hat_ST, label='Full Model', color='yellow', linewidth=3.5)
    ax.plot(ratio_span, QSSA_SP_hat_ST, label='QSSA', linewidth=2)
    ax.plot(ratio_span, tQSSA_SP_hat_ST, label='tQSSA', linestyle='--', color='black')

    ax.set_xlabel('$k_e E_T / k_d D_T$')
    ax.set_ylabel('steady-state $\hat{S}_P/S_T$')
    ax.legend()

    # Embed the plot in the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=rootSwitch)
    canvas.get_tk_widget().grid(row=5, columnspan=6)


rootSwitch = tk.Tk()
rootSwitch.title("Kinetic Reaction Simulator - GK Switch steady state")

# ET
label_max_ratio = ttk.Label(rootSwitch, text="Max ke ET / kd DT ratio:")
label_max_ratio.grid(row=0, column=0)

entry_max_ratio = ttk.Entry(rootSwitch)
entry_max_ratio.insert(0, "2")
entry_max_ratio.grid(row=0, column=1)

# DT
label_DT = ttk.Label(rootSwitch, text="DT:")
label_DT.grid(row=0, column=2)

entry_DT = ttk.Entry(rootSwitch)
entry_DT.insert(0, "100")
entry_DT.grid(row=0, column=3)

# ST
label_ST = ttk.Label(rootSwitch, text="ST:")
label_ST.grid(row=0, column=4)

entry_ST = ttk.Entry(rootSwitch)
entry_ST.insert(0, "100")
entry_ST.grid(row=0, column=5)

# kfe value
label_kfe = ttk.Label(rootSwitch, text="kfe value:")
label_kfe.grid(row=1, column=0)

entry_kfe = ttk.Entry(rootSwitch)
entry_kfe.insert(0, "10")
entry_kfe.grid(row=1, column=1)

# kbe value
label_kbe = ttk.Label(rootSwitch, text="kbe value:")
label_kbe.grid(row=1, column=2)

entry_kbe = ttk.Entry(rootSwitch)
entry_kbe.insert(0, "8.3")
entry_kbe.grid(row=1, column=3)

# ke value
label_ke = ttk.Label(rootSwitch, text="ke value:")
label_ke.grid(row=1, column=4)

entry_ke = ttk.Entry(rootSwitch)
entry_ke.insert(0, "1.7")
entry_ke.grid(row=1, column=5)

# kfd value
label_kfd = ttk.Label(rootSwitch, text="kfd value:")
label_kfd.grid(row=2, column=0)

entry_kfd = ttk.Entry(rootSwitch)
entry_kfd.insert(0, "10")
entry_kfd.grid(row=2, column=1)

# kbd value
label_kbd = ttk.Label(rootSwitch, text="kbd value:")
label_kbd.grid(row=2, column=2)

entry_kbd = ttk.Entry(rootSwitch)
entry_kbd.insert(0, "8.3")
entry_kbd.grid(row=2, column=3)

# kd value
label_kd = ttk.Label(rootSwitch, text="kd value:")
label_kd.grid(row=2, column=4)

entry_kd = ttk.Entry(rootSwitch)
entry_kd.insert(0, "1.7")
entry_kd.grid(row=2, column=5)


start_button = ttk.Button(rootSwitch, text="Start Simulation", command=start_simulation)
start_button.grid(row=3, columnspan=6)

label_result = ttk.Label(rootSwitch, text="CONCENTRATION PLOTS")
label_result.grid(row=4, columnspan=6)

rootSwitch.mainloop()
