import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk
import sys

sys.path += ['..', '../stochastic-chemical-kinetics/pybind']
import SimulationFunctionsOnly
import cme
import gillespie

def weighted_average_and_std(arr, weights):
    """Calculate weighted average and standard deviation."""
    average = np.average(arr, weights=weights)
    variance = np.average((arr - average) ** 2, weights=weights)
    return average, np.sqrt(variance)

def save_figure():
    plt.savefig('./Simulator-Outputs/tQSSA-completion-times-overlay.png')

def run_cme_and_gillespie_simulation():
    """Run both CME and Gillespie simulations and overlay the results for tQSSA."""
    global canvas
    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS
    ET = int(input_entries[0].get())
    ST = int(input_entries[1].get())
    kf = float(input_entries[2].get())
    kb = float(input_entries[3].get())
    kcat = float(input_entries[4].get())

    kM = SimulationFunctionsOnly.MMConstant(kf, kb, kcat)

    max_t = float(entry_time.get())  # Simulation time
    
    # Running tQSSA simulation for CME
    cme_sim = cme.single_substrate_tqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)
    bins_cme = np.linspace(0, max_t, 19)
    
    completion_time_weights_cme = np.empty_like(bins_cme)
    completion_states_cme = ST

    for i, bin in enumerate(bins_cme):
        cme_sim.simulate(dt=1e-4, t_final=bin, noreturn=True)
        completion_time_weights_cme[i] = cme_sim.p[completion_states_cme]

    completion_time_weights_cme = np.diff(completion_time_weights_cme)
    t_cme = (bins_cme[1:] + bins_cme[:-1]) / 2

    average_cme, sd_cme = weighted_average_and_std(t_cme, weights=completion_time_weights_cme)
    print('CME tQSSA', average_cme, '+/-', sd_cme)

    # Running tQSSA simulation for Gillespie
    gillespie_sim = gillespie.single_substrate_tqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)
    n_simulations = 20000

    completion_times_gillespie = np.empty((n_simulations,))

    init_conditions_gillespie = [0]

    for i in range(n_simulations):
        gillespie_sim.x = init_conditions_gillespie
        gillespie_sim.t = 0
        gillespie_sim.simulate(t_final=max_t, noreturn=True)
        completion_times_gillespie[i] = gillespie_sim.t

    average_gillespie, sd_gillespie = weighted_average_and_std(completion_times_gillespie, weights=np.ones_like(completion_times_gillespie))
    print('Gillespie tQSSA', average_gillespie, '+/-', sd_gillespie)

    # Overlaying both plots
    ax.clear()
    ax.hist(t_cme, bins_cme, weights=completion_time_weights_cme, label='CME tQSSA', density=True, color='blue', histtype='step', fill=False)
    ax.hist(completion_times_gillespie, bins_cme, label='Gillespie tQSSA', density=True, color='red', alpha=.6)

    ax.set_ylim(0, 0.4)
    ax.set_xlabel('Completion time (τ)')
    ax.set_ylabel('Probability density')
    ax.legend()

    canvas.draw()
    
    # Save plot
    plt.savefig('./Simulator-Outputs/tQSSA-completion-times-CMEvsGillespie.png')
    print('Figure saved as \'./Simulator-Outputs/tQSSA-completion-times-CMEvsGillespie.png\'')


root = tk.Tk()
root.title('Combined Kinetic Reaction Simulator - Completion Times for Gillespie and CME - tQSSA')

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

# Buttons for running simulations
start_button = ttk.Button(root, text='Run simulation', command=run_cme_and_gillespie_simulation)
start_button.grid(row=2, column=0, columnspan=6, pady=10)

result_label = ttk.Label(root, text='SIMULATION RESULTS')
result_label.grid(row=4, columnspan=6)

# Create a button for saving the figure
save_button = ttk.Button(root, text='Save Figure', command=save_figure)
save_button.grid(row=6, column=0, columnspan=6)

root.mainloop()
