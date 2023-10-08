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

def save_figure():
    # Save the figure in the same directory as the Python code
    fig.savefig('.\Simulator-Outputs\Averages-CMEvsGillespie-tQSSA.png')
    print('Figure saved as \'Averages-CMEvsGillespie-tQSSA.png\'')
    
    
def run_simulation():
    # Extracting CME parameters
    ET_cme = int(input_entries[0].get())
    ST_cme = int(input_entries[1].get())
    kf_cme = float(input_entries[2].get())
    kb_cme = float(input_entries[3].get())
    kcat_cme = float(input_entries[4].get())
    max_t_cme = float(entry_time_cme.get())

    kM_cme = SimulationFunctionsOnly.MMConstant(kf_cme, kb_cme, kcat_cme)

    sim_cme = ['Exact', 'tQSSA', 'sQSSA']
    c_cme = {}

    c_cme['Exact'] = cme.single_substrate(kf=kf_cme, kb=kb_cme, kcat=kcat_cme, ET=ET_cme, ST=ST_cme)
    c_cme['tQSSA'] = cme.single_substrate_tqssa(kM=kM_cme, kcat=kcat_cme, ET=ET_cme, ST=ST_cme)
    c_cme['sQSSA'] = cme.single_substrate_sqssa(kM=kM_cme, kcat=kcat_cme, ET=ET_cme, ST=ST_cme)

    avePs_cme = {}
    msqPs_cme = {}
    ts_cme = {}
    marginal_axes_cme = {'Exact': (int(c_cme['Exact'].species.C)+1), 'tQSSA': (), 'sQSSA': ()}
    possible_Ps_cme = np.arange(0, ST_cme+1)

    for s in sim_cme:
        prob, ts_cme[s] = c_cme[s].simulate(dt=1e-4, t_final=max_t_cme, n_sampling=100)
        marginal_prob = np.sum(prob, axis=marginal_axes_cme[s])
        avePs_cme[s] = np.tensordot(marginal_prob, possible_Ps_cme, axes=1)
        msqPs_cme[s] = np.tensordot(marginal_prob, possible_Ps_cme**2, axes=1)

    # Extracting Gillespie parameters
    ET_gillespie = int(input_entries_gillespie[0].get())
    ST_gillespie = int(input_entries_gillespie[1].get())
    kf_gillespie = float(input_entries_gillespie[2].get())
    kb_gillespie = float(input_entries_gillespie[3].get())
    kcat_gillespie = float(input_entries_gillespie[4].get())
    max_t_gillespie = float(entry_time_gillespie.get())

    kM_gillespie = SimulationFunctionsOnly.MMConstant(kf_gillespie, kb_gillespie, kcat_gillespie)
    sim_gillespie = ['Exact', 'tQSSA', 'sQSSA']
    g = {}

    g['Exact'] = gillespie.single_substrate(kf=kf_gillespie, kb=kb_gillespie, kcat=kcat_gillespie, ET=ET_gillespie, ST=ST_gillespie)
    g['tQSSA'] = gillespie.single_substrate_tqssa(kM=kM_gillespie, kcat=kcat_gillespie, ET=ET_gillespie, ST=ST_gillespie)
    g['sQSSA'] = gillespie.single_substrate_sqssa(kM=kM_gillespie, kcat=kcat_gillespie, ET=ET_gillespie, ST=ST_gillespie)

    avePs_gillespie = {}
    msqPs_gillespie = {}
    hists_gillespie = {'Exact': [], 'tQSSA': [], 'sQSSA': []}
    n_simulations = 5000
    init_conditions = {'Exact': [0, 0], 'tQSSA': [0], 'sQSSA': [0]}

    ndiv = int(max_t_gillespie*50 + 1)
    bins_gillespie = np.linspace(0, max_t_gillespie, ndiv)

    hists_gillespie = {s: [] for s in sim_gillespie}  # Initialize hists as an empty list for each simulation

    for s in sim_gillespie:
        for _ in range(n_simulations):
            g[s].x = init_conditions[s]
            g[s].t = 0
            x, t = g[s].simulate(t_final=max_t_gillespie)
            dprods = np.diff(x[:, g[s].species.P])
            dhist, _ = np.histogram(t[1:], bins_gillespie, weights=dprods)
            hists_gillespie[s].append(np.cumsum(dhist))

        hists_gillespie[s] = np.array(hists_gillespie[s])

    for s in sim_gillespie:
        avePs_gillespie[s] = np.sum(hists_gillespie[s], axis=0) / hists_gillespie[s].shape[0]
        msqPs_gillespie[s] = np.sum(hists_gillespie[s]**2, axis=0) / hists_gillespie[s].shape[0]

    # PLOTTING SIMULATION RESULTS
    ax.clear()
    ax.set_xlim(0, max_t_cme)
    ax.set_ylim(0, max(ST_cme, ST_gillespie)+1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Products count (average and st. dev., $P$)')

 # Add plots for tQSSA model only
    s = 'tQSSA'
    error = np.sqrt(np.maximum(msqPs_cme[s] - avePs_cme[s]**2, 0))
    ax.plot(ts_cme[s], avePs_cme[s], label=f'CME {s}')
    ax.fill_between(ts_cme[s], 0, error, alpha=.3)

    t_gillespie = (bins_gillespie[1:] + bins_gillespie[:-1]) / 2
    error = np.sqrt(np.maximum(msqPs_gillespie[s] - avePs_gillespie[s]**2, 0))
    ax.plot(t_gillespie, avePs_gillespie[s], label=f'Gillespie {s}')
    ax.fill_between(t_gillespie, 0, error, alpha=.3)

    ax.legend()
    canvas.draw()
    
    # Enable the "Save Figure" button after simulation
    save_button['state'] = 'normal'

root = tk.Tk()
root.title('Combined Kinetic Reaction Simulator - Averages for Gillespie and CME - tQSSA model')

# Create and layout input widgets for CME
input_labels = ['ET value (CME):', 'ST value (CME):', 'kf value (CME):', 'kb value (CME):', 'kcat value (CME):', 'Max Time (CME):']
input_entries = []
default_entries = ['10', '9', '10', '9', '1', '10']

for i, label_text in enumerate(input_labels):
    ttk.Label(root, text=label_text).grid(row=i, column=0)
    entry = ttk.Entry(root)
    entry.insert(0, default_entries[i])
    entry.grid(row=i, column=1)
    input_entries.append(entry)

entry_time_cme = ttk.Entry(root)
entry_time_cme.insert(0, '10')
entry_time_cme.grid(row=len(input_labels), column=1)

# Create and layout input widgets for Gillespie
input_labels_gillespie = ['ET value (Gillespie):', 'ST value (Gillespie):', 'kf value (Gillespie):', 'kb value (Gillespie):', 'kcat value (Gillespie):', 'Max Time (Gillespie):']
input_entries_gillespie = []
default_entries_gillespie = ['10', '9', '10', '9', '1', '9']

for i, label_text in enumerate(input_labels_gillespie):
    ttk.Label(root, text=label_text).grid(row=i, column=2)
    entry = ttk.Entry(root)
    entry.insert(0, default_entries_gillespie[i])
    entry.grid(row=i, column=3)
    input_entries_gillespie.append(entry)

entry_time_gillespie = ttk.Entry(root)
entry_time_gillespie.insert(0, '9')
entry_time_gillespie.grid(row=len(input_labels_gillespie), column=3)

# Button for running simulation and updating the plot
start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=len(input_labels_gillespie)+1, columnspan=4)

result_label = ttk.Label(root, text='SIMULATION RESULTS')
result_label.grid(row=len(input_labels_gillespie)+2, columnspan=4)

# Create a figure for plotting
fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=0, column=4, rowspan=len(input_labels_gillespie)+3)

# Button for saving figure 

save_button = ttk.Button(root, text='Save Figure', command=save_figure, state='disabled')
save_button.grid(row=len(input_labels_gillespie)+3, columnspan=4)


root.mainloop()
