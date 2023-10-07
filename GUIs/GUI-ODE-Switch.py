import tkinter as tk
import matplotlib.pyplot as plt
import numpy as np

from tkinter import ttk
from scipy.integrate import odeint
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import sys
sys.path.append('..')
import SimulationFunctionsOnly


def run_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    ET = float(input_entries[0].get())
    DT = float(input_entries[1].get())
    ST = float(input_entries[2].get())

    initial_SP = 0
    initial_C  = 0
    initial_CP = 0

    initial_conditions = [initial_SP, initial_C, initial_CP]


    # Getting the user input from the entry widgets and setting constants

    kfe = float(input_entries[3].get())  # kfe value
    kbe = float(input_entries[4].get())  # kbe value
    ke  = float(input_entries[5].get())  # ke value
    kfd = float(input_entries[6].get())  # kfd value
    kbd = float(input_entries[7].get())  # kbd value
    kd  = float(input_entries[8].get())  # kd value

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

    ax.clear()

    # Plot the data
    ax.plot(time_span, Full_Model_SP_hat_ST, label='Full Model', color='yellow', linewidth=3.5)
    ax.plot(time_span, QSSA_SP_hat_ST, label='QSSA', linewidth=2)
    ax.plot(time_span, tQSSA_SP_hat_ST, label='tQSSA', linestyle='--', color='black')

    ax.set_xlabel('Time')
    ax.set_ylabel('$\hat{S_P}/S_T$')
    ax.legend()

    # Update the Tkinter canvas

    canvas.draw()
    
    # Save plot
    
    plt.savefig('Switch.png')


root = tk.Tk()
root.title('Kinetic Reaction Simulator - GK Switch')

fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=9, column=0, columnspan=6)

# Create and layout input widgets
input_labels = ['ET value:', 'DT value:', 'ST value:', 'kfe value:', 'kbe value:', 'ke value:', 'kfd value:', 'kbd value:', 'kd value:']
input_entries = []
default_entries = ['100', '95', '100', '10', '8.3', '1.7', '10', '8.3', '1.7']

for i, label_text in enumerate(input_labels):
    ttk.Label(root, text=label_text).grid(row=(i//3)*2, column=i%3)
    entry = ttk.Entry(root)
    entry.insert(0, default_entries[i])
    entry.grid(row=(i//3)*2+1, column=i%3)
    input_entries.append(entry)


# Label for Simulation Time
label_time = ttk.Label(root, text='Simulation Time:')
label_time.grid(row=6, column=0)

entry_time = ttk.Entry(root)
entry_time.insert(0, '10')
entry_time.grid(row=6, column=1)

start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=7, columnspan=6)

label_result = ttk.Label(root, text='CONCENTRATION PLOTS')
label_result.grid(row=8, columnspan=6)

root.mainloop()
