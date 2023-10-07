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

import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import animation, colors

from tkinter import ttk
import sys

sys.path.append('../stochastic-chemical-kinetics/pybind')
import cme

# have mmpeg installed to run this script. Using conda:

#   conda install ffmpeg

# Create the initial tick labels (assuming a 10x10 grid)
init_x_max = 10
init_y_max = 10
x_ticks = np.arange(0, init_x_max, 1)
y_ticks = np.arange(0, init_y_max, 1)
x_tick_labels = [str(int(label)) for label in x_ticks]
y_tick_labels = [str(int(label)) for label in y_ticks]

def initialize_simulation():
    ET = int(input_entries[0].get())
    ST = int(input_entries[1].get())
    kf = float(input_entries[2].get())
    kb = float(input_entries[3].get())
    kcat = float(input_entries[4].get())

    c = cme.single_substrate(kf=kf, kb=kb, kcat=kcat, ET=ET, ST=ST)
    
    return c

def init_animation():
    c = initialize_simulation()
    im.set_data(c.p)
    annot.set_text('Time: 0.000')
    return im, annot

def animate_frame(i, c):
    current_t = i * 0.002
    c.simulate(dt=1e-4, t_final=current_t, noreturn=True)
    im.set_data(c.p)
    annot.set_text(f'Time: {current_t:.3f}')

    return im, annot

def start_simulation():
    global anim  # Declare anim as a global variable
    if anim:
        anim.event_source.stop()  # Stop the previous animation if it exists
    
    c = initialize_simulation()
    
    # Get the current data shape
    x_max, y_max = c.p.shape[1], c.p.shape[0]
    
    # Calculate tick positions for the x-axis and set the labels
    x_tick_positions = np.arange(0, x_max, 1)
    x_tick_labels = [str(int(label)) for label in x_tick_positions]
    x_tick_positions = (x_tick_positions + .5) * init_x_max / x_max - .5
    ax.set_xticks(x_tick_positions)
    ax.set_xticklabels(x_tick_labels)
    
    # Calculate tick positions for the y-axis and set the labels
    y_tick_positions = np.arange(0, y_max, 1)
    y_tick_labels = [str(int(label)) for label in y_tick_positions]
    y_tick_positions = (y_tick_positions + .5) * init_y_max / y_max - .5
    ax.set_yticks(y_tick_positions)
    ax.set_yticklabels(y_tick_labels)

    fps = 30
    anim = animation.FuncAnimation(fig, animate_frame, fargs=(c,), init_func=init_animation, interval=1000/fps, blit=False, cache_frame_data=False)
    canvas.draw()

def stop_simulation():
    if anim:
        anim.event_source.stop()
        print('Simulation stopped manually.')

root = tk.Tk()
root.title('Kinetic Reaction Simulator')

fig, ax = plt.subplots()
im = ax.imshow(np.zeros(shape=(10, 10)), norm=colors.SymLogNorm(1e-20), origin='lower')
plt.colorbar(im, ax=ax, label='Probability')
im.set_clim(1e-10, 1)
annot = ax.annotate('Time: 0.000', (0.09, 0.92), xycoords='figure fraction')
plt.xlabel('Product count ($P$)')
plt.ylabel('Enzyme-substrate complex count ($C$)')

anim = None

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=5, column=0, columnspan=6)

# Set initial tick labels
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.set_xticklabels(x_tick_labels)
ax.set_yticklabels(y_tick_labels)

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

# Start and stop buttons
start_button = ttk.Button(root, text='Start Simulation', command=start_simulation)
start_button.grid(row=3, columnspan=6)

stop_button = ttk.Button(root, text='Stop Simulation', command=stop_simulation)
stop_button.grid(row=6, column=0, columnspan=6)

root.mainloop()
