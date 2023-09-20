import tkinter as tk
import matplotlib.pyplot as plt
import numpy as np

from tkinter import ttk
from scipy.integrate import odeint
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import sys
sys.path.append('..')
import SimulationFunctionsOnly


def start_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    ET = float(entry_ET.get())
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

    # Getting the user input for simulation time
    simulation_time = float(entry_time.get())  # Simulation time

    # Creating the time span based on user-defined simulation time
    time_span = np.linspace(0, simulation_time, 1000)

    # solving odes

    Full_Model_Solution = odeint(SimulationFunctionsOnly.Full_Model_Switch, initial_conditions, time_span, args=(ST, ET, DT, kfe, kbe, ke, kfd, kbd, kd))
    QSSAsolution = odeint(SimulationFunctionsOnly.QSSA_Switch, initial_SP, time_span, args=(ST, ET, DT, ke, kd, kME, kMD))
    tQSSAsolution = odeint(SimulationFunctionsOnly.tQSSA_Switch, initial_SP+initial_CP, time_span, args=(ST, ET, DT, ke, kd, kME, kMD))

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
    ax.plot(time_span, Full_Model_SP_hat_ST, label='Full Model', color='yellow', linewidth=3.5)
    ax.plot(time_span, QSSA_SP_hat_ST, label='QSSA', linewidth=2)
    ax.plot(time_span, tQSSA_SP_hat_ST, label='tQSSA', linestyle='--', color='black')

    ax.set_xlabel('Time')
    ax.set_ylabel('$\hat{S_P}/S_T$')
    ax.legend()

    # Embed the plot in the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=rootSwitch)
    canvas.get_tk_widget().grid(row=5, columnspan=6)


rootSwitch = tk.Tk()
rootSwitch.title("Kinetic Reaction Simulator - GK Switch")

# ET
label_ET = ttk.Label(rootSwitch, text="ET:")
label_ET.grid(row=0, column=0)

entry_ET = ttk.Entry(rootSwitch)
entry_ET.insert(0, "100")
entry_ET.grid(row=0, column=1)

# DT
label_DT = ttk.Label(rootSwitch, text="DT:")
label_DT.grid(row=0, column=2)

entry_DT = ttk.Entry(rootSwitch)
entry_DT.insert(0, "95")
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


# Label for Simulation Time
label_time = ttk.Label(rootSwitch, text="Simulation Time:")
label_time.grid(row=3, column=0)

entry_time = ttk.Entry(rootSwitch)
entry_time.insert(0, "10")
entry_time.grid(row=3, column=1)

start_button = ttk.Button(rootSwitch, text="Start Simulation", command=start_simulation)
start_button.grid(row=3, columnspan=6)

label_result = ttk.Label(rootSwitch, text="CONCENTRATION PLOTS")
label_result.grid(row=4, columnspan=6)

rootSwitch.mainloop()
