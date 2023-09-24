import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk

from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

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

    avePs = {}
    msqPs = {}
    
    hists = {'Exact': [], 'tQSSA': [], 'sQSSA': []}
    n_simulations = 5000
    init_conditions = {'Exact': [0, 0], 'tQSSA': [0], 'sQSSA': [0]}

    ndiv = int(max_t*50 + 1)
    bins = np.linspace(0, max_t, ndiv)

    hists = {s: [] for s in sim}  # Initialize hists as an empty list for each simulation

    for s in sim:
        for _ in range(n_simulations):
            g[s].x = init_conditions[s]
            g[s].t = 0
            x, t = g[s].simulate(t_final=max_t)
            dprods = np.diff(x[:, g[s].species.P])
            dhist, _ = np.histogram(t[1:], bins, weights=dprods)
            hists[s].append(np.cumsum(dhist))

        hists[s] = np.array(hists[s])


    for s in sim:
        avePs[s] = np.sum(hists[s], axis=0) / hists[s].shape[0]
        msqPs[s] = np.sum(hists[s]**2, axis=0) / hists[s].shape[0]
    
    # PLOTTING SIMULATION RESULTS
    
    t = (bins[1:] + bins[:-1]) / 2

    ax.clear()

    for s in sim:
        error = np.sqrt(np.maximum(msqPs[s] - avePs[s]**2, 0))
        ax.plot(t, avePs[s], label=s)
        ax.fill_between(t, 0, error, alpha=.3)

    ax.set_ylim(0)
    ax.set_xlabel('Time')
    ax.set_ylabel('Products count (average and st. dev., $P$)')
    ax.legend()

    # Update the Tkinter canvas

    canvas.draw()
    
    print(f"Running simulation with parameters: ST={ST}, ET={ET}, kf={kf}, kb={kb}, kcat={kcat}")    


root = tk.Tk()
root.title("Kinetic Reaction Simulator")

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
