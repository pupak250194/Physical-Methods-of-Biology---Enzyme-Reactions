import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import sys

sys.path += ['..', '../stochastic-chemical-kinetics/pybind']

import cme
import SimulationFunctionsOnly
import gillespie

def save_figure():
    # Save the figure in the same directory as the Python code
    fig.savefig('.\Simulator-Outputs\Stationary-CMEvsGillespie-tQSSA.png')
    print('Figure saved as \'Stationary-CMEvsGillespie-tQSSA.png\'')

def weighted_ave_sd(arr, weights):
	average = np.average(arr, weights=weights)
	variance = np.average((arr-average)**2, weights=weights)
	return average, np.sqrt(variance)


def run_simulation():
    ET = int(input_entries[0].get())
    DT = int(input_entries[1].get())
    ST = int(input_entries[2].get())

    kfe = float(input_entries[3].get())  # kfe value
    kbe = float(input_entries[4].get())  # kbe value
    ke = float(input_entries[5].get())   # ke value
    kfd = float(input_entries[6].get())  # kfd value
    kbd = float(input_entries[7].get())  # kbd value
    kd = float(input_entries[8].get())   # kd value

    kME, kMD = SimulationFunctionsOnly.Constants_Switch(kfe, kbe, ke, kfd, kbd, kd)

    g_sim = ['Exact', 'tQSSA', 'sQSSA']
    c_sim = ['tQSSA', 'sQSSA']

    g = {}
    c = {}

    g['Exact'] = gillespie.goldbeter_koshland(kfe=kfe, kbe=kbe, ke=ke, kfd=kfd, kbd=kbd, kd=kd, ET=ET, DT=DT, ST=ST)
    g['tQSSA'] = gillespie.goldbeter_koshland_tqssa(kME=kME, ke=ke, kMD=kMD, kd=kd, ET=ET, DT=DT, ST=ST)
    g['sQSSA'] = gillespie.goldbeter_koshland_sqssa(kME=kME, ke=ke, kMD=kMD, kd=kd, ET=ET, DT=DT, ST=ST)

    c['tQSSA'] = cme.goldbeter_koshland_tqssa(kME=kME, ke=ke, kMD=kMD, kd=kd, ET=ET, DT=DT, ST=ST)
    c['sQSSA'] = cme.goldbeter_koshland_sqssa(kME=kME, ke=ke, kMD=kMD, kd=kd, ET=ET, DT=DT, ST=ST)

    SP_hats = {}
    SP_hat_dists = {}
    weights = {}

    init_conditions_g = {'Exact': [ST // 2, 0, 0], 'tQSSA': [ST // 2], 'sQSSA': [ST // 2]}
    init_conditions_c = {'tQSSA': np.zeros_like(c['tQSSA'].p), 'sQSSA': np.zeros_like(c['sQSSA'].p)}
    init_conditions_c['tQSSA'][ST // 2] = 1
    init_conditions_c['sQSSA'][ST // 2] = 1

    max_t_g = float(entry_time_gillespie.get())
    max_t_c = float(entry_time_cme.get())

    tracked_species_g = {
        'Exact': np.array([g['Exact'].species.SP, g['Exact'].species.CP], dtype=int),
        'tQSSA': np.array([g['tQSSA'].species.SP_hat], dtype=int),
        'sQSSA': np.array([g['sQSSA'].species.SP], dtype=int)
    }

    possible_SP_hats = np.arange(0, ST+1)
    
    for s in g_sim:
        g[s].x = init_conditions_g[s]
        g[s].t = 0
        x, t = g[s].simulate(t_final=max_t_g, max_steps=int(5e7))
        SP_hats[s] = np.sum(x[:-1, tracked_species_g[s]], axis=1)
        weights[s] = np.diff(t)

    for s in c_sim:
        c[s].p = init_conditions_c[s]
        c[s].t = 0
        c[s].simulate(dt=1e-3, t_final=max_t_c, noreturn=True)
        SP_hat_dists[s] = c[s].p
    
    print('\nGillespie stats\n')
    for s in g_sim:
        average, sd = weighted_ave_sd(SP_hats[s], weights=weights[s])
        print(s, average, '+/-', sd)

    print('\nCME stats\n')
    for s in c_sim:
        ave = np.dot(c[s].p, possible_SP_hats)
        msq = np.dot(c[s].p, possible_SP_hats**2)
        print(s, ave, '+/-', np.sqrt(np.maximum(msq - ave**2, 0)))

    bins = np.linspace(0, ST, 21)

    ax.clear()

    plt.hist(SP_hats['tQSSA'], bins, weights=weights['tQSSA'], label='tQSSA (Gillespie)', density=True, color='blue', alpha=.6, histtype='step')
    plt.hist(possible_SP_hats, bins, weights=SP_hat_dists['tQSSA'], label='tQSSA (CME)', density=True, color='red', alpha=.6)

    ax.set_xlabel('Steady-state phosphorylated substrate count ($\hat{S}_P$)')
    ax.set_ylabel('Probability density')
    ax.legend()

    canvas.draw()

root = tk.Tk()
root.title('Kinetic Reaction Simulator - Stationary for Gillespie and CME - tQSSA')

fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.2)  # Increase space at the bottom

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=10, column=0, columnspan=6, pady=(10, 0))

# Create and layout input widgets
input_labels = ['ET value:', 'DT value:', 'ST value:', 'kfe value:', 'kbe value:', 'ke value:',
                'kfd value:', 'kbd value:', 'kd value:']
input_entries = []
default_entries = ['100', '100', '100', '10', '8.3', '1.7', '10', '8.3', '1.7']

for i, label_text in enumerate(input_labels):
    ttk.Label(root, text=label_text).grid(row=(i // 3) * 2, column=i % 3)
    entry = ttk.Entry(root)
    entry.insert(0, default_entries[i])
    entry.grid(row=(i // 3) * 2 + 1, column=i % 3)
    input_entries.append(entry)

# Label for Simulation Time
ttk.Label(root, text='Simulation Time (Gillespie):').grid(row=6, column=0)
entry_time_gillespie = ttk.Entry(root)
entry_time_gillespie.insert(0, '10000')
entry_time_gillespie.grid(row=7, column=0)

ttk.Label(root, text='Simulation Time (CME):').grid(row=6, column=1)
entry_time_cme = ttk.Entry(root)
entry_time_cme.insert(0, '100')
entry_time_cme.grid(row=7, column=1)

start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=8, column=0, columnspan=3, pady=(10, 0))  # Increased vertical space

result_label_cmestat = ttk.Label(root, text="SIMULATION RESULTS")
result_label_cmestat.grid(row=9, column=0, columnspan=3, pady=(10, 0)) 

save_button = ttk.Button(root, text='Save Figure', command=save_figure)
save_button.grid(row=11, columnspan=6)

root.mainloop()