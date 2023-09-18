import tkinter as tk
import matplotlib.pyplot as plt
import sys
import numpy as np

from tkinter import ttk
from scipy.integrate import odeint
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import SimulationFunctionsOnly


def start_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    ST = float(entry_ST.get())
    ET = float(entry_ET.get())

    initial_C = entry_C0.get()
    initial_P = entry_P0.get()
    initial_conditions = [initial_C, initial_P]
    
    # Getting the user input from the entry widgets and setting constants

    kf = float(entry_kf.get())      # kf value
    kb = float(entry_kb.get())      # kb value
    kcat = float(entry_kcat.get())  # kcat value  

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

    # Create a new figure for the plot
    fig, ax = plt.subplots()
    
    # Plot the data
    ax.plot(time_span, Full_Model_ratio_values, label='Full Model', color='yellow', linewidth=3.5)
    ax.plot(time_span, QSSA_ratio_values, label='QSSA', linewidth=2)
    ax.plot(time_span, tQSSA_ratio_values, label='tQSSA', linestyle='--', color='black')

    ax.set_xlabel('Time')
    ax.set_ylabel('P/ST')
    ax.legend()

    # Embed the plot in the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().grid(row=5, columnspan=6)

    # Replace the print statement with your simulation code
    print(f"Running simulation with parameters: ET={ET}, ST={ST}, C0={initial_C}, P0={initial_P}, kf={kf}, kb={kb}, kcat={kcat}")


root = tk.Tk()
root.title("Kinetic Reaction Simulator")

# ET value
label_ET = ttk.Label(root, text="ET:")
label_ET.grid(row=0, column=0)

entry_ET = ttk.Entry(root)
entry_ET.insert(0, "10")
entry_ET.grid(row=0, column=1)

# ST value
label_ST = ttk.Label(root, text="ST:")
label_ST.grid(row=0, column=2)

entry_ST = ttk.Entry(root)
entry_ST.insert(0, "9")
entry_ST.grid(row=0, column=3)

# Initial C
label_C0 = ttk.Label(root, text="C0:")
label_C0.grid(row=1, column=0)

entry_C0 = ttk.Entry(root)
entry_C0.insert(0, "0")
entry_C0.grid(row=1, column=1)

# Initial P
label_P0 = ttk.Label(root, text="P0:")
label_P0.grid(row=1, column=2)

entry_P0 = ttk.Entry(root)
entry_P0.insert(0, "0")
entry_P0.grid(row=1, column=3)

# kf value
label_kf = ttk.Label(root, text="kf value:")
label_kf.grid(row=2, column=0)

entry_kf = ttk.Entry(root)
entry_kf.insert(0, "10")
entry_kf.grid(row=2, column=1)

# kb value
label_kb = ttk.Label(root, text="kb value:")
label_kb.grid(row=2, column=2)

entry_kb = ttk.Entry(root)
entry_kb.insert(0, "9")
entry_kb.grid(row=2, column=3)

# kcat value
label_kcat = ttk.Label(root, text="kcat value:")
label_kcat.grid(row=2, column=4)

entry_kcat = ttk.Entry(root)
entry_kcat.insert(0, "1")
entry_kcat.grid(row=2, column=5)

# Label for Simulation Time
label_time = ttk.Label(root, text="Simulation Time:")
label_time.grid(row=3, column=0)

entry_time = ttk.Entry(root)
entry_time.insert(0, "10")
entry_time.grid(row=3, column=1)

start_button = ttk.Button(root, text="Start Simulation", command=start_simulation)
start_button.grid(row=4, columnspan=6)

label_result = ttk.Label(root, text="CONCENTRATION PLOTS")
label_result.grid(row=5, columnspan=6)

root.mainloop()
