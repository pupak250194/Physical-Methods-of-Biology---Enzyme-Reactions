
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

#   conda install ffmpeg


import matplotlib.pyplot as plt
import tkinter as tk

from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import sys
sys.path.append('../stochastic-chemical-kinetics/pybind')

import gillespie

def run_simulation():
    """Run the simulation and plot the results."""
    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    ET = int(input_entries[0].get()) 
    ST = int(input_entries[1].get()) 
    
    # Getting the user input from the entry widgets and setting constants

    kf = float(input_entries[2].get())    # kf value
    kb = float(input_entries[3].get())    # kb value
    kcat = float(input_entries[4].get())  # kcat value  
     
    # Running simulation
    
    g = gillespie.single_substrate(kf=kf, kb=kb, kcat=kcat, ET=ET , ST=ST)
    x, t = g.simulate()

    ax.clear()

    ax.plot(t, x[:,g.species.C], drawstyle='steps-post', label='C')
    ax.plot(t, x[:,g.species.P], drawstyle='steps-post', label='P')

    ax.set_xlabel('Time')
    ax.set_ylabel('Population')
    ax.legend()
  
    canvas.draw()
    
    # Save plot
    
    plt.savefig('Gillespie-concentrations.png')
    
          
root = tk.Tk()
root.title("Kinetic Reaction Simulator")

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

start_button = ttk.Button(root, text='Run simulation', command=run_simulation)
start_button.grid(row=3, columnspan=6)

result_label = ttk.Label(root, text="SIMULATION RESULTS")
result_label.grid(row=4, columnspan=6)

root.mainloop()




