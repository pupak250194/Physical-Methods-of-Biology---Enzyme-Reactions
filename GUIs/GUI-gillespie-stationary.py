
'''
    Stochastic enzyme kinetics: CME python example
    Copyright (C) 2023 Alessandro Lo Cuoco (alessandro.locuoco@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

# have mmpeg installed to run this script. Using conda:

#   conda install -c conda-forge ffmpeg


import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk

from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import sys
sys.path += ['..', '../stochastic-chemical-kinetics/pybind']

import SimulationFunctionsOnly
import gillespie

def weighted_ave_sd(arr, weights):
	average = np.average(arr, weights=weights)
	variance = np.average((arr-average)**2, weights=weights)
	return average, np.sqrt(variance)

def run_simulation():
    """Run the simulation and plot the results."""
    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    ET = int(input_entries[0].get())
    DT = int(input_entries[1].get())
    ST = int(input_entries[2].get())

    # Getting the user input from the entry widgets and setting constants

    kfe = float(input_entries[3].get())  # kfe value
    kbe = float(input_entries[4].get())  # kbe value
    ke  = float(input_entries[5].get())   # ke value
    kfd = float(input_entries[6].get())  # kfd value
    kbd = float(input_entries[7].get())  # kbd value
    kd  = float(input_entries[8].get())   # kd value

    kME, kMD = SimulationFunctionsOnly.Constants_Switch(kfe, kbe, ke, kfd, kbd, kd)

    # SETTING TIME SPAN

    # Getting the user input for simulation time
    simulation_time = float(entry_time.get())  # Simulation time

    # not using exact formulation due to very long computational time
    sim = ['Exact', 'tQSSA', 'sQSSA']
    g = {}

    g['Exact'] = gillespie.goldbeter_koshland(kfe=kfe, kbe=kbe, ke=ke, kfd=kfd, kbd=kbd, kd=kd, ET=ET, DT=DT, ST=ST)
    g['tQSSA'] = gillespie.goldbeter_koshland_tqssa(kME=kME, ke=ke, kMD=kMD, kd=kd, ET=ET, DT=DT, ST=ST)
    g['sQSSA'] = gillespie.goldbeter_koshland_sqssa(kME=kME, ke=ke, kMD=kMD, kd=kd, ET=ET, DT=DT, ST=ST)

    SP_hats = {}
    weights = {}
    init_conditions = {'Exact': [ST//2, 0, 0], 'tQSSA': [ST//2], 'sQSSA': [ST//2]}
    max_t = simulation_time

    tracked_species = {
        'Exact': np.array([g['Exact'].species.SP, g['Exact'].species.CP], dtype=int),
        'tQSSA': np.array([g['tQSSA'].species.SP_hat], dtype=int),
        'sQSSA': np.array([g['sQSSA'].species.SP], dtype=int)
    }

    for s in sim:
        g[s].x = init_conditions[s]
        g[s].t = 0
        x, t = g[s].simulate(t_final=max_t, max_steps=int(5e7))
        SP_hats[s] = np.sum(x[:-1,tracked_species[s]], axis=1)
        weights[s] = np.diff(t)

    for s in sim:
        average, sd = weighted_ave_sd(SP_hats[s], weights=weights[s])
        print(s, average, '+/-', sd)

    if always_div.get():
        bins = np.linspace(0, ST, 21)
    else:
        bins = np.linspace(0, ST+1, ST+2) - 0.5

    # Plots

    ax.clear()

    # Add your histograms

    ax.hist(SP_hats['Exact'], bins, weights=weights['Exact'], label='Exact', density=True, color='blue', histtype='step', fill=False)
    ax.hist(SP_hats['tQSSA'], bins, weights=weights['tQSSA'], label='tQSSA', density=True, color='red', alpha=.6)
    ax.hist(SP_hats['sQSSA'], bins, weights=weights['sQSSA'], label='sQSSA', density=True, color='gray', alpha=.3)

    ax.set_xlabel('Steady-state phosphorylated substrate count ($\hat{S}_P$)')
    ax.set_ylabel('Probability density')
    ax.legend()

    canvas.draw()
    
    # Save plot
    
    plt.savefig('.\Simulator-Outputs\Gillespie-stationary.png')


root = tk.Tk()
root.title('Kinetic Reaction Simulator - Stationary for Gillespie')

fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=9, column=0, columnspan=6)

# Create and layout input widgets
input_labels = ['ET value:', 'DT value:', 'ST value:', 'kfe value:', 'kbe value:', 'ke value:', 'kfd value:', 'kbd value:', 'kd value:']
input_entries = []
default_entries = ['100', '100', '100', '10', '8.3', '1.7', '10', '8.3', '1.7']

for i, label_text in enumerate(input_labels):
    ttk.Label(root, text=label_text).grid(row=(i//3)*2, column=i%3)
    entry = ttk.Entry(root)
    entry.insert(0, default_entries[i])
    entry.grid(row=(i//3)*2+1, column=i%3)
    input_entries.append(entry)

# Label for Simulation Time
ttk.Label(root, text='Simulation Time:').grid(row=6, column=0)
entry_time = ttk.Entry(root)
entry_time.insert(0, '10000')
entry_time.grid(row=6, column=1)

# Checkbox: always divide into 20 bins
always_div = tk.BooleanVar()
checkbox = ttk.Checkbutton(root, text="Always divide into 20 bins", variable=always_div)
checkbox.grid(row=6, column=2)

start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=7, columnspan=6)

result_labelcmestat= ttk.Label(root, text='SIMULATION RESULTS')
result_labelcmestat.grid(row=8, columnspan=6)

root.mainloop()