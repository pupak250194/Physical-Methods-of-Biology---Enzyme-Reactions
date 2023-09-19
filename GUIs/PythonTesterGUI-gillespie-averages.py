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
    
    # SETTING TIME SPAN

    # Getting the user input for simulation time
    simulation_time = float(entry6.get())  # Simulation time
    
    
    # Running simulation
    
    sim = ['Exact', 'tQSSA', 'sQSSA']
    g = {}

    g['Exact'] = gillespie.single_substrate(kf=float(simkf), kb=float(simkb), kcat=float(simkcat), ET=initial_E, ST=initial_S)
    g['tQSSA'] = gillespie.single_substrate_tqssa(kM=float(simkM), kcat=float(simkcat), ET=initial_E, ST=initial_S)
    g['sQSSA'] = gillespie.single_substrate_sqssa(kM=float(simkM), kcat=float(simkcat), ET=initial_E, ST=initial_S)

    avePs = {}
    msqPs = {}
    max_t = simulation_time
    
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

    # Creating a new figure for the plots
    fig, ax = plt.subplots(figsize=(8, 6))

    for s in sim:
        error = np.sqrt(np.maximum(msqPs[s] - avePs[s]**2, 0))
        t = np.linspace(0, max_t, len(avePs[s]))  # Ensure t has the correct dimension
        ax.plot(t, avePs[s], label=s)
        ax.fill_between(t, 0, error, alpha=.3)

    ax.set_ylim(0)
    ax.set_xlabel('Time')
    ax.set_ylabel('Products count (average and st. dev., $P$)')
    ax.legend()

    # Embedding the Matplotlib figure in  Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().grid(row=5, column=0, columnspan=6)
    
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

#output_text = tk.Text(root, height=5, width=40)
#output_text.grid(row=5, columnspan=6)

root.mainloop()
