import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk

from scipy.integrate import odeint
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk
 
import sys
sys.path.append('..')
import SimulationFunctionsOnly 

def run_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    ET = int(input_entries[0].get())
    ST = int(input_entries[1].get())

    initial_C = input_entries[2].get()
    initial_P = input_entries[3].get()
    initial_conditions = [initial_C, initial_P]
    
    # Getting the user input from the entry widgets and setting constants

    kf = float(input_entries[4].get())    # kf value
    kb = float(input_entries[5].get())    # kb value
    kcat = float(input_entries[6].get())  # kcat value

    kM = SimulationFunctionsOnly.MMConstant(kf, kb, kcat)
    
    # SETTING TIME SPAN

    # Getting the user input for simulation time
    simulation_time = float(entry_time.get())  # Simulation time
    
    # Creating the time span based on user-defined simulation time
    time_span = np.linspace(0, simulation_time, 1000)

    # solving odes
    
    Full_Model_Solution = odeint(SimulationFunctionsOnly.Full_Model, initial_conditions, time_span, args=(kf, kb, kcat, ST, ET))
    QSSAsolution = odeint(SimulationFunctionsOnly.QSSA, initial_P, time_span, args=(kcat, kM, ST, ET))
    tQSSAsolution = odeint(SimulationFunctionsOnly.tQSSA, initial_P, time_span, args=(kcat, kM, ST, ET))

    Full_Model_P_values = Full_Model_Solution[:, 1]
    QSSA_P_values = QSSAsolution[:, 0]
    tQSSA_P_values = tQSSAsolution[:, 0]

    # Ratio P/S according to the 3 models

    Full_Model_ratio_values = Full_Model_P_values / ST
    QSSA_ratio_values = QSSA_P_values / ST
    tQSSA_ratio_values = tQSSA_P_values / ST

    # Plots

    ax.clear()

    # Plot the data
    ax.plot(time_span, Full_Model_ratio_values, label='Full Model', color='yellow', linewidth=3.5)
    ax.plot(time_span, QSSA_ratio_values, label='QSSA', linewidth=2)
    ax.plot(time_span, tQSSA_ratio_values, label='tQSSA', linestyle='--', color='black')

    ax.set_xlabel('Time')
    ax.set_ylabel('$P/S_T$')
    ax.legend()

    # Update the Tkinter canvas

    canvas.draw()
    
    # Save plot
    
    plt.savefig('Full-Model.png')

    # Replace the print statement with your simulation code
    print(f"Running simulation with parameters: ET={ET}, ST={ST}, C0={initial_C}, P0={initial_P}, kf={kf}, kb={kb}, kcat={kcat}")


root = tk.Tk()
root.title('Kinetic Reaction Simulator')

fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=7, column=0, columnspan=6)

# Create and layout input widgets
input_labels = ['ET value:', 'ST value:', 'C0 value:', 'P0 value:', 'kf value:', 'kb value:', 'kcat value:']
input_entries = []
default_entries = ['10', '9', '0', '0', '10', '9', '1']

for i, label_text in enumerate(input_labels):
    ttk.Label(root, text=label_text).grid(row=(i//4)*2, column=i%4)
    entry = ttk.Entry(root)
    entry.insert(0, default_entries[i])
    entry.grid(row=(i//4)*2+1, column=i%4)
    input_entries.append(entry)

# Label for Simulation Time
label_time = ttk.Label(root, text='Simulation Time:')
label_time.grid(row=4, column=0)

entry_time = ttk.Entry(root)
entry_time.insert(0, '10')
entry_time.grid(row=4, column=1)

start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=5, columnspan=6)

label_result = ttk.Label(root, text='CONCENTRATION PLOTS')
label_result.grid(row=6, columnspan=6)

root.mainloop()
