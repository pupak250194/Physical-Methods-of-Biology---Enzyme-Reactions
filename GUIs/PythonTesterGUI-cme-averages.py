import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk

import sys
sys.path += ['path/to/SimulationFunctionsOnly','path/to/cme']

import SimulationFunctionsOnly
import cme

def start_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    initial_S = int(entry1.get()) 
    initial_E = int(entry2.get()) 
    
    # Getting the user input from the entry widgets and setting constants

    kf_value = float(entry3.get())    # kf value
    kb_value = float(entry4.get())    # kb value
    kcat_value = float(entry5.get())  # kcat value  
    
    simkM = SimulationFunctionsOnly.MMConstants(kf_value, kb_value, kcat_value)
    
    # SETTING TIME SPAN

    # Getting the user input for simulation time
    simulation_time = float(entry6.get())  # Simulation time
    
    
    # Running simulation
    
    sim = ['Exact', 'tQSSA', 'sQSSA']
    c = {}

    c['Exact'] = cme.single_substrate(kf=kf_value, kb=kb_value, kcat=kcat_value, ET=initial_E, ST=initial_S)
    c['tQSSA'] = cme.single_substrate_tqssa(kM=simkM, kcat=kcat_value, ET=initial_E, ST=initial_S)
    c['sQSSA'] = cme.single_substrate_sqssa(kM=simkM, kcat=kcat_value, ET=initial_E, ST=initial_S)

    avePs = {}
    msqPs = {}
    ts = {}
    marginal_axes = {'Exact': (int(c['Exact'].species.C)+1), 'tQSSA': (), 'sQSSA': ()}
    max_t = simulation_time
    possible_Ps = np.arange(0, initial_S+1)

    for s in sim:
        prob, ts[s] = c[s].simulate(dt=1e-4, t_final=max_t, n_sampling=100)
        marginal_prob = np.sum(prob, axis=marginal_axes[s])
        avePs[s] = np.tensordot(marginal_prob, possible_Ps, axes=1)
        msqPs[s] = np.tensordot(marginal_prob, possible_Ps**2, axes=1)
	    
    # PLOTTING SIMULATION RESULTS
    
    # Create a figure
    fig, ax = plt.subplots()
    
    # Manually set the y-axis limit based on the maximum value you want to display
    ax.set_xlim(0, max_t)  # Adjust as needed
    ax.set_ylim(0, 10)  # Adjust as needed (replace 10 with your desired maximum value)

    ax.set_xticks(np.arange(0, max_t+1, step=1))  # Adjust as needed
    ax.set_yticks(np.arange(0, 11, step=1))  # Adjust as needed

    ax.set_xlabel('Time')
    ax.set_ylabel('Products count (average and st. dev., $P$)')

    # Add plots to the figure
    for s in sim:
        error = np.sqrt(np.maximum(msqPs[s] - avePs[s]**2, 0))
        ax.plot(ts[s], avePs[s], label=s)
        ax.fill_between(ts[s], 0, error, alpha=.3)

    # Embed the figure in the tkinter window
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().grid(row=6, columnspan=6)

    print(f"Running simulation with parameters: Initial S={initial_S}, Initial E={initial_E}, kf={kf_value}, kb={kb_value}, kcat={kcat_value}")    
    

root = tk.Tk()
root.title("Kinetic Reaction Simulator")

# Initial C
label1 = ttk.Label(root, text="Initial S:")
label1.grid(row=0, column=0)

entry1 = ttk.Entry(root)
entry1.grid(row=0, column=1)

# Initial P
label2 = ttk.Label(root, text="Initial E:")
label2.grid(row=0, column=2)

entry2 = ttk.Entry(root)
entry2.grid(row=0, column=3)

# kf value
label3 = ttk.Label(root, text="kf value:")
label3.grid(row=1, column=0)

entry3 = ttk.Entry(root)
entry3.grid(row=1, column=1)

# kb value
label4 = ttk.Label(root, text="kb value:")
label4.grid(row=1, column=2)

entry4 = ttk.Entry(root)
entry4.grid(row=1, column=3)

# kcat value
label5 = ttk.Label(root, text="kcat value:")
label5.grid(row=1, column=4)

entry5 = ttk.Entry(root)
entry5.grid(row=1, column=5)

# Label for Simulation Time
label6 = ttk.Label(root, text="Simulation Time:")
label6.grid(row=2, column=0)

entry6 = ttk.Entry(root)
entry6.grid(row=2, column=1)

start_button = ttk.Button(root, text="Start Simulation", command=start_simulation)
start_button.grid(row=3, columnspan=6)

result_label = ttk.Label(root, text="SIMULATION RESULTS")
result_label.grid(row=4, columnspan=6)

root.mainloop()
