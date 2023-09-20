
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
        x, _ = g[s].simulate(t_final=max_t)
        SP_hats[s] = np.sum(x[:,tracked_species[s]], axis=1)

    for s in sim:
        print(s, np.mean(SP_hats[s]), '+/-', np.std(SP_hats[s]))

    bins = np.linspace(0, ST, 21)

    # Plots

    ax.clear()

    # Add your histograms

    ax.hist(SP_hats['Exact'], bins, label='Exact', density=True, color='blue', histtype='step', fill=False)
    ax.hist(SP_hats['tQSSA'], bins, label='tQSSA', density=True, color='red', alpha=.6)
    ax.hist(SP_hats['sQSSA'], bins, label='sQSSA', density=True, color='gray', alpha=.3)

    ax.set_xlabel('Steady-state phosphorylated substrate count ($\hat{S}_P$)')
    ax.set_ylabel('Probability density')
    ax.legend()

    canvas.draw()


root = tk.Tk()
root.title('Kinetic Reaction Simulator - cme stat')

fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=9, column=0, columnspan=6)

# Create and layout input widgets
input_labels = ['ET value:', 'DT value:', 'ST value:', 'kfe value:', 'kbe value:', 'ke value:', 'kfd value:', 'kbd value:', 'kd value:']
input_entries = []
default_entries = ['100', '100', '100', '8.7', '1.3', '1', '8.7', '1.3', '1']

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

start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=7, columnspan=6)

result_labelcmestat= ttk.Label(root, text='SIMULATION RESULTS')
result_labelcmestat.grid(row=8, columnspan=6)

root.mainloop()