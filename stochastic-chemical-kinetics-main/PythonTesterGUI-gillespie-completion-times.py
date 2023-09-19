import tkinter as tk
import matplotlib.pyplot as plt
import sys
import numpy as np

from tkinter import ttk
from scipy.integrate import odeint
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


sys.path += ['C:/Users/39333/Desktop/Physical-Methods-of-Biology---Enzyme-Reactions-main/Differential-Eqs','C:/Users/39333/Desktop/Physical-Methods-of-Biology---Enzyme-Reactions-main/stochastic-chemical-kinetics-main/pybind']

import SimulationFunctionsOnly
import gillespie

def start_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    initial_S = int(entry1.get())
    initial_E = int(entry2.get())

    # Getting the user input from the entry widgets and setting constants

    kf_value = float(entry3.get())    # kf value
    kb_value = float(entry4.get())    # kb value
    kcat_value = float(entry5.get())  # kcat value

    simkf, simkb, simkcat, simkM = SimulationFunctionsOnly.Constants(kf_value, kb_value, kcat_value)

    simulationtime = int(entry6.get())

    # Running simulation

    sim = ['Exact', 'tQSSA', 'sQSSA']
    g = {}

    g['Exact'] = gillespie.single_substrate(kf=float(simkf), kb=float(simkb), kcat=float(simkcat), ET=initial_E, ST=initial_S)
    g['tQSSA'] = gillespie.single_substrate_tqssa(kM=float(simkM), kcat=float(simkcat), ET=initial_E, ST=initial_S)
    g['sQSSA'] = gillespie.single_substrate_sqssa(kM=float(simkM), kcat=float(simkcat), ET=initial_E, ST=initial_S)

    n_simulations = 20000
    max_t = simulationtime

    completion_times = {
	'Exact': np.empty((n_simulations,)),
	'tQSSA': np.empty((n_simulations,)),
	'sQSSA': np.empty((n_simulations,))
     }

    init_conditions = {'Exact': [0, 0], 'tQSSA': [0], 'sQSSA': [0]}

    for s in sim:
        for i in range(n_simulations):
            g[s].x = init_conditions[s]
            g[s].t = 0
            g[s].simulate(t_final=max_t, noreturn=True)
            completion_times[s][i] = g[s].t

    for s in sim:
        print(s, np.mean(completion_times[s]), '+/-', np.std(completion_times[s]))

    bins = np.linspace(0, 9, 19)
    
    fig, ax = plt.subplots(figsize=(5, 4), dpi=100)  # Create a matplotlib figure
    ax.hist(completion_times['Exact'], bins, label='Exact', density=True, color='blue', histtype='step', fill=False)
    ax.hist(completion_times['tQSSA'], bins, label='tQSSA', density=True, color='red', alpha=.6)
    ax.hist(completion_times['sQSSA'], bins, label='sQSSA', density=True, color='gray', alpha=.3)

    ax.set_ylim(0, .4)
    ax.set_xlabel('Completion time (τ)')
    ax.set_ylabel('Probability density')
    ax.legend()

    canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
    canvas.draw()

    canvas.get_tk_widget().grid(row=5, column=0, columnspan=6)


root = tk.Tk()
root.title("Kinetic Reaction Simulator")

# Initial S
label1 = ttk.Label(root, text="Initial S:")
label1.grid(row=0, column=0)

entry1 = ttk.Entry(root)
entry1.grid(row=0, column=1)

# Initial E
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

#time entry
label6 = ttk.Label(root, text="Simulation time:")
label6.grid(row=2, column=0)

entry6 = ttk.Entry(root)
entry6.grid(row=2, column=1)


start_button = ttk.Button(root, text="Start Simulation", command=start_simulation)
start_button.grid(row=3, columnspan=6)

result_label = ttk.Label(root, text="SIMULATION RESULTS")
result_label.grid(row=4, columnspan=6)

#output_text = tk.Text(root, height=5, width=40)
#output_text.grid(row=4, columnspan=6)

root.mainloop()
