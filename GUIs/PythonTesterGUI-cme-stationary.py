
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

import sys

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk

sys.path += ['path/to/SimulationFunctionsOnly','path/to/cme']


import cme
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
    sim = ['tQSSA', 'sQSSA']
    c = {}

    c['tQSSA'] = cme.goldbeter_koshland_tqssa(kME=simkME, ke=ke_value, kMD=simkMD, kd=kd_value, ET=initial_E, DT=initial_D, ST=initial_S)
    c['sQSSA'] = cme.goldbeter_koshland_sqssa(kME=simkME, ke=ke_value, kMD=simkMD, kd=kd_value, ET=initial_E, DT=initial_D, ST=initial_S)

    SP_hat_dists = {}

    init_conditions = {}
    
    for s in sim:
       init_conditions[s] = np.zeros_like(c[s].p)
       
    init_conditions['tQSSA'][initial_S//2] = 1
    init_conditions['sQSSA'][initial_S//2] = 1

    max_t = simulation_time
    ave = {}
    msq = {}
    possible_SP_hats = np.arange(0, initial_S+1)

    for s in sim:
       c[s].p = init_conditions[s]
       c[s].t = 0
       c[s].simulate(dt=1e-3, t_final=max_t, noreturn=True)
       SP_hat_dists[s] = c[s].p
       ave[s] = np.dot(c[s].p, possible_SP_hats)
       msq[s] = np.dot(c[s].p, possible_SP_hats**2)

    for s in sim: 
        print(s, ave[s], '+/-', np.sqrt(np.maximum(msq[s] - ave[s]**2, 0)))

    bins = np.linspace(0, initial_S, 21)
    bins_all = np.linspace(0, initial_S+1, initial_S+2)
    x = (bins_all[1:] + bins_all[:-1])/2

    # Plots

    # Create a new figure for the plot
    fig = plt.figure()

    # Add your histograms
    plt.hist(x, bins, weights=SP_hat_dists['tQSSA'], label='tQSSA', density=True, color='red', alpha=.6)
    plt.hist(x, bins, weights=SP_hat_dists['sQSSA'], label='sQSSA', density=True, color='gray', alpha=.3)

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


