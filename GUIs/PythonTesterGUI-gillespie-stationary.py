
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
sys.path += ['path/to/SimulationFunctionsOnly','path/to/gillespie']


import gillespie
import SimulationFunctionsOnly


def start_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    initial_S = int(entry1cmestat.get()) 
    initial_E = int(entry2cmestat.get()) 
    initial_D = int(entry3cmestat.get()) 

    # Getting the user input from the entry widgets and setting constants

    ke_value = float(entry4cmestat.get())    # ke value
    kbe_value = float(entry5cmestat.get())    # kbe value
    kfe_value = float(entry6cmestat.get())  # kfe value 
    kbd_value = float(entry7cmestat.get())  # kbd value  
    kfd_value = float(entry8cmestat.get())  # kfd value  
    kd_value = float(entry9cmestat.get())  # kfd value 

    simkME, simkMD = SimulationFunctionsOnly.Constants_Switch (ke_value, kbe_value, kfe_value, kbd_value, kfd_value, kd_value)
    
    # SETTING TIME SPAN

    # Getting the user input for simulation time
    simulation_time = float(entry10cmestat.get())  # Simulation time
    
    # not using exact formulation due to very long computational time
    sim = ['Exact', 'tQSSA', 'sQSSA']
    g = {}

    g['Exact'] = gillespie.goldbeter_koshland(kfe=kfe_value, kbe=kbe_value, ke=ke_value, kfd=kfd_value, kbd=kbd_value, kd=kd_value, ET=initial_E, DT=initial_D, ST=initial_S)
    g['tQSSA'] = gillespie.goldbeter_koshland_tqssa(kME=simkME, ke=ke_value, kMD=simkMD, kd=kd_value, ET=initial_E, DT=initial_D, ST=initial_S)
    g['sQSSA'] = gillespie.goldbeter_koshland_sqssa(kME=simkME, ke=ke_value, kMD=simkMD, kd=kd_value, ET=initial_E, DT=initial_D, ST=initial_S)

    SP_hats = {}
    init_conditions = {'Exact': [initial_S//2, 0, 0], 'tQSSA': [initial_S//2], 'sQSSA': [initial_S//2]}
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

    bins = np.linspace(0, initial_S, 21)
    
    # Plots

    # Create a new figure for the plot
    
    fig = plt.figure()
    
    # Add your histograms
    
    plt.hist(SP_hats['Exact'], bins, label='Exact', density=True, color='blue', histtype='step', fill=False)
    plt.hist(SP_hats['tQSSA'], bins, label='tQSSA', density=True, color='red', alpha=.6)
    plt.hist(SP_hats['sQSSA'], bins, label='sQSSA', density=True, color='gray', alpha=.3)

    plt.xlabel('Steady-state phosphorylated substrate count ($\hat{S}_P$)')
    plt.ylabel('Probability density')
    plt.legend()
    
    # Create the FigureCanvasTkAgg widget
    canvas = FigureCanvasTkAgg(fig, master=rootcmestat)
    canvas.draw()

    # Add the widget to your GUI
    canvas.get_tk_widget().grid(row=5, column=0, columnspan=6)    

    # Replace the print statement with your simulation code
    #print(f"Running simulation with parameters: Initial S={initial_S}, Initial E={initial_E}, kf={kf_value}, kb={kb_value}, kcat={kcat_value}")


rootcmestat = tk.Tk()
rootcmestat.title("Kinetic Reaction Simulator - cme stat")

# Initial S
label1cmestat = ttk.Label(rootcmestat, text="Initial S:")
label1cmestat.grid(row=0, column=0)

entry1cmestat = ttk.Entry(rootcmestat)
entry1cmestat.grid(row=0, column=1)

# Initial E
label2cmestat = ttk.Label(rootcmestat, text="Initial E:")
label2cmestat.grid(row=0, column=2)

entry2cmestat = ttk.Entry(rootcmestat)
entry2cmestat.grid(row=0, column=3)

# Initial D
label3cmestat = ttk.Label(rootcmestat, text="Initial D:")
label3cmestat.grid(row=0, column=4)

entry3cmestat = ttk.Entry(rootcmestat)
entry3cmestat.grid(row=0, column=5)

# ke value
label4cmestat = ttk.Label(rootcmestat, text="ke value:")
label4cmestat.grid(row=1, column=0)

entry4cmestat = ttk.Entry(rootcmestat)
entry4cmestat.grid(row=1, column=1)

# kbe value
label5cmestat = ttk.Label(rootcmestat, text="kbe value:")
label5cmestat.grid(row=1, column=2)

entry5cmestat = ttk.Entry(rootcmestat)
entry5cmestat.grid(row=1, column=3)

# kfe value
label6cmestat = ttk.Label(rootcmestat, text="kfe value:")
label6cmestat.grid(row=1, column=4)

entry6cmestat = ttk.Entry(rootcmestat)
entry6cmestat.grid(row=1, column=5)

# kbd value
label7cmestat = ttk.Label(rootcmestat, text="kbd value:")
label7cmestat.grid(row=2, column=0)

entry7cmestat = ttk.Entry(rootcmestat)
entry7cmestat.grid(row=2, column=1)

# kfd value
label8cmestat = ttk.Label(rootcmestat, text="kfd value:")
label8cmestat.grid(row=2, column=2)

entry8cmestat = ttk.Entry(rootcmestat)
entry8cmestat.grid(row=2, column=3)

# kfd value
label9cmestat = ttk.Label(rootcmestat, text="kd value:")
label9cmestat.grid(row=2, column=4)

entry9cmestat = ttk.Entry(rootcmestat)
entry9cmestat.grid(row=2, column=5)


# Label for Simulation Time
label10cmestat = ttk.Label(rootcmestat, text="Simulation Time:")
label10cmestat.grid(row=3, column=0)

entry10cmestat = ttk.Entry(rootcmestat)
entry10cmestat.grid(row=3, column=1)

start_button = ttk.Button(rootcmestat, text="Start Simulation", command=start_simulation)
start_button.grid(row=3, columnspan=6)

result_labelcmestat= ttk.Label(rootcmestat, text="SIMULATION RESULTS")
result_labelcmestat.grid(row=4, columnspan=6)

rootcmestat.mainloop()