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
    ET = int(input_entries[0].get())
    ST = int(input_entries[1].get())
    kf = float(input_entries[2].get())
    kb = float(input_entries[3].get())
    kcat = float(input_entries[4].get())
    max_t = float(entry_time.get())

    kM = SimulationFunctionsOnly.MMConstant(kf, kb, kcat)

    sim = ['Exact', 'tQSSA', 'sQSSA']
    c = {}

    c['Exact'] = cme.single_substrate(kf=kf, kb=kb, kcat=kcat, ET=ET, ST=ST)
    c['tQSSA'] = cme.single_substrate_tqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)
    c['sQSSA'] = cme.single_substrate_sqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)

    avePs_c = {}
    msqPs_c = {}
    ts_c = {}
    marginal_axes_c = {'Exact': (int(c['Exact'].species.C)+1), 'tQSSA': (), 'sQSSA': ()}
    possible_Ps_c = np.arange(0, ST+1)

    for s in sim:
        prob, ts_c[s] = c[s].simulate(dt=1e-4, t_final=max_t, n_sampling=100)
        marginal_prob = np.sum(prob, axis=marginal_axes_c[s])
        avePs_c[s] = np.tensordot(marginal_prob, possible_Ps_c, axes=1)
        msqPs_c[s] = np.tensordot(marginal_prob, possible_Ps_c**2, axes=1)

    g = {}

    g['Exact'] = gillespie.single_substrate(kf=kf, kb=kb, kcat=kcat, ET=ET, ST=ST)
    g['tQSSA'] = gillespie.single_substrate_tqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)
    g['sQSSA'] = gillespie.single_substrate_sqssa(kM=kM, kcat=kcat, ET=ET, ST=ST)

    avePs_g = {}
    msqPs_g = {}
    hists_g = {'Exact': [], 'tQSSA': [], 'sQSSA': []}
    n_simulations = 5000
    init_conditions = {'Exact': [0, 0], 'tQSSA': [0], 'sQSSA': [0]}

    ndiv = int(max_t*50 + 1)
    bins_g = np.linspace(0, max_t, ndiv)

    hists_g = {s: [] for s in sim}  # Initialize hists as an empty list for each simulation

    for s in sim:
        for _ in range(n_simulations):
            g[s].x = init_conditions[s]
            g[s].t = 0
            x, t = g[s].simulate(t_final=max_t)
            dprods = np.diff(x[:, g[s].species.P])
            dhist, _ = np.histogram(t[1:], bins_g, weights=dprods)
            hists_g[s].append(np.cumsum(dhist))

        hists_g[s] = np.array(hists_g[s])

    for s in sim:
        avePs_g[s] = np.sum(hists_g[s], axis=0) / hists_g[s].shape[0]
        msqPs_g[s] = np.sum(hists_g[s]**2, axis=0) / hists_g[s].shape[0]

    # PLOTTING SIMULATION RESULTS
    ax.clear()
    ax.set_xlim(0, max_t)
    ax.set_ylim(0, ST+1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Products count (average and st. dev., $P$)')

 # Add plots for tQSSA model only
    s = 'tQSSA'
    error = np.sqrt(np.maximum(msqPs_c[s] - avePs_c[s]**2, 0))
    ax.plot(ts_c[s], avePs_c[s], label=f'CME {s}')
    ax.fill_between(ts_c[s], 0, error, alpha=.3)

    t_g = (bins_g[1:] + bins_g[:-1]) / 2
    error = np.sqrt(np.maximum(msqPs_g[s] - avePs_g[s]**2, 0))
    ax.plot(t_g, avePs_g[s], label=f'Gillespie {s}')
    ax.fill_between(t_g, 0, error, alpha=.3)

    ax.legend()
    canvas.draw()
    
    # Enable the "Save Figure" button after simulation
    save_button['state'] = 'normal'

root = tk.Tk()
root.title('Combined Kinetic Reaction Simulator - Averages for Gillespie and CME - tQSSA model')

# Create and layout input widgets
input_labels = ['ET value:', 'ST value:', 'kf value:', 'kb value:', 'kcat value:']
input_entries = []
default_entries = ['10', '9', '10', '9', '1']

for i, label_text in enumerate(input_labels):
    ttk.Label(root, text=label_text).grid(row=i, column=0)
    entry = ttk.Entry(root)
    entry.insert(0, default_entries[i])
    entry.grid(row=i, column=1)
    input_entries.append(entry)

entry_time = ttk.Entry(root)
entry_time.insert(0, '10')
entry_time.grid(row=len(input_labels), column=1)

# Button for running simulation and updating the plot
start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=len(input_labels)+1, columnspan=4)

result_label = ttk.Label(root, text='SIMULATION RESULTS')
result_label.grid(row=len(input_labels)+2, columnspan=4)

# Create a figure for plotting
fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=0, column=4, rowspan=len(input_labels)+3)

# Button for saving figure 

save_button = ttk.Button(root, text='Save Figure', command=save_figure, state='disabled')
save_button.grid(row=len(input_labels)+3, columnspan=4)


root.mainloop()
