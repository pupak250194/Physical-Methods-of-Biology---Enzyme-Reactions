
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

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib import animation, colors
import tkinter as tk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk

import sys
sys.path.append(r'path/to/cme')
import cme

def start_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    initial_S = int(entry1.get()) 
    initial_E = int(entry2.get()) 
    
    # Getting the user input from the entry widgets and setting constants

    kf_value = float(entry3.get())    # kf value
    kb_value = float(entry4.get())    # kb value
    kcat_value = float(entry5.get())  # kcat value  
     
    # Running simulation
    
    c = cme.single_substrate(kf=kf_value, kb=kb_value, kcat= kcat_value, ET=initial_E , ST=initial_S)

    def animate(i):
      current_t = i*0.002
      c.simulate(dt=1e-4, t_final=current_t, noreturn=True)
      im.set_data(c.p)
      annot.set_text('Time: {:.3f}'.format(current_t))
      return im, annot


    fig = plt.figure()
    ax = plt.axes()
    im = ax.imshow(c.p, norm=colors.SymLogNorm(1e-20), origin='lower')
    plt.colorbar(im, ax=ax, label='Probability')
    im.set_clim(1e-10, 1)

    annot = ax.annotate('Time: 0.000', (0.09, 0.92), xycoords='figure fraction')

    fps = 30
    anim = animation.FuncAnimation(fig, animate, init_func=None, frames=5000, interval=1000/fps, blit=False)
   
    # Embedding the animation in the Tkinter window...
  
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().grid(row=5, column=0, columnspan=6)
    
 
    plt.xlabel('Product count ($P$)')
    plt.ylabel('Enzyme-substrate complex count ($C$)')

    def update_animation():
        try:
            anim.__next__()
            root.after(1000 // fps, update_animation)  # Update at the desired FPS
        except StopIteration:
            pass

    update_animation()

    
    if len(sys.argv) > 1:
	      if sys.argv[1] == 's' or sys.argv[1] == 'save':
		        anim.save('cme_example.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])
          
    

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

start_button = ttk.Button(root, text="Start Simulation", command=start_simulation)
start_button.grid(row=3, columnspan=6)

result_label = ttk.Label(root, text="SIMULATION RESULTS")
result_label.grid(row=4, columnspan=6)

#output_text = tk.Text(root, height=5, width=40)
#output_text.grid(row=5, columnspan=6)

root.mainloop()




