import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk

import sys
sys.path += ['..', '../stochastic-chemical-kinetics/pybind']

import SimulationFunctionsOnly
import gillespie

def run_simulation():
    """Run the simulation and plot the results."""
    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    ET = int(input_entries[0].get()) 
    ST = int(input_entries[1].get()) 
    
    # Getting the user input from the entry widgets and setting constants

    kf = float(input_entries[2].get())    # kf value
    kb = float(input_entries[3].get())    # kb value
    kcat = float(input_entries[4].get())  # kcat value

    kM = SimulationFunctionsOnly.MMConstant(kf, kb, kcat)

    max_t = int(entry_time.get())

    # Running simulation

    sim = ['Exact', 'tQSSA', 'sQSSA']
    g = {}

    g['Exact'] = gillespie.single_substrate(kf=kf, kb=kb, kcat=kcat, ET=ET, ST=ST)
    g['tQSSA'] = gillespie.single_substrate_tqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)
    g['sQSSA'] = gillespie.single_substrate_sqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)

    n_simulations = 20000

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
            g[s].simulate(t_final=max_t+1, noreturn=True)
            completion_times[s][i] = g[s].t

    for s in sim:
        print(s, np.mean(completion_times[s]), '+/-', np.std(completion_times[s]))

    bins = np.linspace(0, max_t, 19)

    ax.clear()

    ax.hist(completion_times['Exact'], bins, label='Exact', density=True, color='blue', histtype='step', fill=False)
    ax.hist(completion_times['tQSSA'], bins, label='tQSSA', density=True, color='red', alpha=.6)
    ax.hist(completion_times['sQSSA'], bins, label='sQSSA', density=True, color='gray', alpha=.3)

    ax.set_ylim(0, .4)
    ax.set_xlabel('Completion time (τ)')
    ax.set_ylabel('Probability density')
    ax.legend()

    canvas.draw()
    
    # Save plot
    
    plt.savefig('.\Simulator-Outputs\Gillespie-completion-times.png')


root = tk.Tk()
root.title("Kinetic Reaction Simulator - Completion times for Gillespie")

fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=5, column=0, columnspan=6)

# Create and layout input widgets
input_labels = ['ET value:', 'ST value:', 'kf value:', 'kb value:', 'kcat value:']
input_entries = []
default_entries = ['10', '9', '10', '9', '1']

for i, label_text in enumerate(input_labels):
    ttk.Label(root, text=label_text).grid(row=0, column=i)
    entry = ttk.Entry(root)
    entry.insert(0, default_entries[i])
    entry.grid(row=1, column=i)
    input_entries.append(entry)

# Label for Max Time
ttk.Label(root, text='Max Time:').grid(row=2, column=0)
entry_time = ttk.Entry(root)
entry_time.insert(0, '9')
entry_time.grid(row=2, column=1)

start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=3, columnspan=6)

result_label = ttk.Label(root, text="SIMULATION RESULTS")
result_label.grid(row=4, columnspan=6)


root.mainloop()
