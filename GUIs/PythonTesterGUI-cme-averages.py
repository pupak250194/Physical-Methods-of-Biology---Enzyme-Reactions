import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk

import sys
sys.path += ['..', '../stochastic-chemical-kinetics/pybind']

import SimulationFunctionsOnly
import cme

def run_simulation():
    """Run the simulation and plot the results."""
    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS
    ET = int(input_entries[0].get())
    ST = int(input_entries[1].get())
    kf = float(input_entries[2].get())
    kb = float(input_entries[3].get())
    kcat = float(input_entries[4].get())

    kM = SimulationFunctionsOnly.MMConstant(kf, kb, kcat)
    
    # SETTING TIME SPAN

    # Getting the user input for simulation time
    simulation_time = float(entry_time.get())  # Simulation time

    # Running simulation
    
    sim = ['Exact', 'tQSSA', 'sQSSA']
    c = {}

    c['Exact'] = cme.single_substrate(kf=kf, kb=kb, kcat=kcat, ET=ET, ST=ST)
    c['tQSSA'] = cme.single_substrate_tqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)
    c['sQSSA'] = cme.single_substrate_sqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)

    avePs = {}
    msqPs = {}
    ts = {}
    marginal_axes = {'Exact': (int(c['Exact'].species.C)+1), 'tQSSA': (), 'sQSSA': ()}
    max_t = simulation_time
    possible_Ps = np.arange(0, ST+1)

    for s in sim:
        prob, ts[s] = c[s].simulate(dt=1e-4, t_final=max_t, n_sampling=100)
        marginal_prob = np.sum(prob, axis=marginal_axes[s])
        avePs[s] = np.tensordot(marginal_prob, possible_Ps, axes=1)
        msqPs[s] = np.tensordot(marginal_prob, possible_Ps**2, axes=1)
	    
    # PLOTTING SIMULATION RESULTS

    ax.clear()

    ax.set_xlim(0, max_t)
    ax.set_ylim(0, ST+1)

    ax.set_xlabel('Time')
    ax.set_ylabel('Products count (average and st. dev., $P$)')
    

    # Add plots to the figure
    for s in sim:
        error = np.sqrt(np.maximum(msqPs[s] - avePs[s]**2, 0))
        ax.plot(ts[s], avePs[s], label=f'{s}')
        ax.fill_between(ts[s], 0, error, alpha=.3)
        
    ax.legend()  # Add this line to display the legend

    # Update the Tkinter canvas

    canvas.draw()
    
    # Save plot
    
    plt.savefig('CME-averages.png')

    print(f'Running simulation with parameters: ST={ST}, ET={ET}, kf={kf}, kb={kb}, kcat={kcat}')   
    
    
     


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

# Label for Simulation Time
ttk.Label(root, text='Simulation Time:').grid(row=2, column=0)
entry_time = ttk.Entry(root)
entry_time.insert(0, '10')
entry_time.grid(row=2, column=1)

start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=3, columnspan=6)

result_label = ttk.Label(root, text='SIMULATION RESULTS')
result_label.grid(row=4, columnspan=6)

root.mainloop()
