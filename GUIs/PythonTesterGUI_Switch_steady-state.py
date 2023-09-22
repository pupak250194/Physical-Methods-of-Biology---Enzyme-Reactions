
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk

from scipy.optimize import fsolve
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk

import sys
sys.path.append('..')

import SimulationFunctionsOnly 

def run_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    DT = float(input_entries[1].get())
    ST = float(input_entries[2].get())

    initial_SP = 0
    initial_C  = 0
    initial_CP = 0

    initial_conditions = [initial_SP, initial_C, initial_CP]


    # Getting the user input from the entry widgets and setting constants

    kfe = float(input_entries[3].get())  # kfe value
    kbe = float(input_entries[4].get())  # kbe value
    ke  = float(input_entries[5].get())   # ke value
    kfd = float(input_entries[6].get())  # kfd value
    kbd = float(input_entries[7].get())  # kbd value
    kd  = float(input_entries[8].get())   # kd value

    kME, kMD = SimulationFunctionsOnly.Constants_Switch(kfe, kbe, ke, kfd, kbd, kd)
    
    # SETTING TIME SPAN

    # Getting the user input for maximum ke ET / kd DT ratio
    max_ratio = float(input_entries[0].get())

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

    ax.clear()

    # Plot the data
    ax.plot(ratio_span, Full_Model_SP_hat_ST, label='Full Model', color='yellow', linewidth=3.5)
    ax.plot(ratio_span, QSSA_SP_hat_ST, label='QSSA', linewidth=2)
    ax.plot(ratio_span, tQSSA_SP_hat_ST, label='tQSSA', linestyle='--', color='black')

    ax.set_xlabel('$k_e E_T / k_d D_T$')
    ax.set_ylabel('steady-state $\hat{S}_P/S_T$')
    ax.legend()

    # Update the Tkinter canvas

    canvas.draw()


root = tk.Tk()
root.title('Kinetic Reaction Simulator - GK Switch steady state')

fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=8, column=0, columnspan=6)

# Create and layout input widgets
input_labels = ['Max ke ET / kd DT ratio:', 'DT value:', 'ST value:', 'kfe value:', 'kbe value:', 'ke value:', 'kfd value:', 'kbd value:', 'kd value:']
input_entries = []
default_entries = ['2', '100', '100', '10', '8.3', '1.7', '10', '8.3', '1.7']

for i, label_text in enumerate(input_labels):
    ttk.Label(root, text=label_text).grid(row=(i//3)*2, column=i%3)
    entry = ttk.Entry(root)
    entry.insert(0, default_entries[i])
    entry.grid(row=(i//3)*2+1, column=i%3)
    input_entries.append(entry)


start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=6, columnspan=6)

label_result = ttk.Label(root, text='CONCENTRATION PLOTS')
label_result.grid(row=7, columnspan=6)

root.mainloop()
