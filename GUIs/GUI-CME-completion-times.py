import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from scipy.integrate import odeint
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk
import sys

sys.path += ['..', '../stochastic-chemical-kinetics/pybind']
import SimulationFunctionsOnly
import cme

def weighted_average_and_std(arr, weights):
    """Calculate weighted average and standard deviation."""
    average = np.average(arr, weights=weights)
    variance = np.average((arr - average) ** 2, weights=weights)
    return average, np.sqrt(variance)

def run_simulation():
    """Run the simulation and plot the results."""
    global canvas
    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS
    ET = int(input_entries[0].get())
    ST = int(input_entries[1].get())
    kf = float(input_entries[2].get())
    kb = float(input_entries[3].get())
    kcat = float(input_entries[4].get())

    kM = SimulationFunctionsOnly.MMConstant(kf, kb, kcat)

    max_t = float(entry_time.get())  # Simulation time
    
    # Running simulation
    sim = ['Exact', 'tQSSA', 'sQSSA']
    c = {}

    c['Exact'] = cme.single_substrate(kf=kf, kb=kb, kcat=kcat, ET=ET, ST=ST)
    c['tQSSA'] = cme.single_substrate_tqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)
    c['sQSSA'] = cme.single_substrate_sqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)

    bins = np.linspace(0, max_t, 19)
    
    completion_time_weights = {
        'Exact': np.empty_like(bins),
        'tQSSA': np.empty_like(bins),
        'sQSSA': np.empty_like(bins)
    }
    
    completion_states = {'Exact': (0, ST), 'tQSSA': ST, 'sQSSA': ST}

    for s in sim:
        for i, bin in enumerate(bins):
            c[s].simulate(dt=1e-4, t_final=bin, noreturn=True)
            completion_time_weights[s][i] = c[s].p[completion_states[s]]

    for s in sim:
        completion_time_weights[s] = np.diff(completion_time_weights[s])

    t = (bins[1:] + bins[:-1]) / 2

    for s in sim:
        average, sd = weighted_average_and_std(t, weights=completion_time_weights[s])
        print(s, average, '+/-', sd)

    ax.clear()
    ax.hist(t, bins, weights=completion_time_weights['Exact'], label='Exact', density=True, color='blue', histtype='step', fill=False)
    ax.hist(t, bins, weights=completion_time_weights['tQSSA'], label='tQSSA', density=True, color='red', alpha=.6)
    ax.hist(t, bins, weights=completion_time_weights['sQSSA'], label='sQSSA', density=True, color='gray', alpha=.3)

    ax.set_ylim(0, 0.4)
    ax.set_xlabel('Completion time (τ)')
    ax.set_ylabel('Probability density')
    ax.legend()

    canvas.draw()
    
    # Save plot
    
    plt.savefig('.\Simulator-Outputs\CME-completion-times.png')

root = tk.Tk()
root.title('Kinetic Reaction Simulator')

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

# Start button
start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=2, columnspan=6)

result_label = ttk.Label(root, text='SIMULATION RESULTS')
result_label.grid(row=3, columnspan=6)

root.mainloop()