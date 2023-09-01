import tkinter as tk
import matplotlib.pyplot as plt
import sys
import numpy as np

from tkinter import ttk
from scipy.integrate import odeint
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


sys.path.append(r"C:/Users/lacri/Desktop/UNI/MAGISTRALE/PhysicalMethodsForBiology/EPaci/Progetto/BozzeCodes")
import SimulationFunctionsOnlyES


def start_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    initial_S = float(entry1.get()) 
    initial_E = float(entry2.get()) 

    initial_C = 0
    initial_P = 0 
    initial_conditions = [initial_C, initial_P]
    
    # Getting the user input from the entry widgets and setting constants

    kf_value = float(entry3.get())    # kf value
    kb_value = float(entry4.get())    # kb value
    kcat_value = float(entry5.get())  # kcat value  

    simkf, simkb, simkcat, simkM = SimulationFunctionsOnlyES.Constants(kf_value, kb_value, kcat_value)
    
    # SETTING TIME SPAN

    # Getting the user input for simulation time
    simulation_time = float(entry6.get())  # Simulation time
    
    # Creating the time span based on user-defined simulation time
    time_span = np.linspace(0, simulation_time, 1000)

    # solving odes
    
    Full_Model_Solution = odeint (SimulationFunctionsOnlyES.Full_Model, initial_conditions, time_span, args=(simkf, simkb, simkcat, initial_S, initial_E))
    QSSAsolution = odeint (SimulationFunctionsOnlyES.QSSA, initial_conditions, time_span, args=(simkcat, simkM, initial_S, initial_E))
    tQSSAsolution = odeint (SimulationFunctionsOnlyES.tQSSA, initial_conditions, time_span, args=(simkcat, simkM, initial_S, initial_E))

    Full_Model_P_values = Full_Model_Solution[:, 1]
    QSSA_P_values = QSSAsolution[:, 1]
    tQSSA_P_values = tQSSAsolution[:, 1]

    # Ratio P/S according to the 3 models

    Full_Model_ratio_values = Full_Model_P_values / initial_S
    QSSA_ratio_values = QSSA_P_values / initial_S
    tQSSA_ratio_values = tQSSA_P_values / initial_S

    # Plots

    # Create a new figure for the plot
    fig, ax = plt.subplots()
    
    # Plot the data
    ax.plot(time_span, Full_Model_ratio_values, label='Full Model', color='yellow', linewidth=3.5)
    ax.plot(time_span, QSSA_ratio_values, label='QSSA - P/Stot', linewidth=2)
    ax.plot(time_span, tQSSA_ratio_values, label='tQSSA - P/S-hat', linestyle='--', color='black')

    ax.set_xlabel('Time')
    ax.set_ylabel('P/S')
    ax.legend()

    # Embed the plot in the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().grid(row=5, columnspan=6) 

    #plt.show()


    # Replace the print statement with your simulation code
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

result_label = ttk.Label(root, text="CONCENTRATION PLOTS")
result_label.grid(row=4, columnspan=6)

#output_text = tk.Text(root, height=5, width=40)
#output_text.grid(row=4, columnspan=6)

root.mainloop()
