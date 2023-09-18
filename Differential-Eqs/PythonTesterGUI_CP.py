import tkinter as tk
import matplotlib.pyplot as plt
import sys
import numpy as np

from tkinter import ttk
from scipy.integrate import odeint
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


sys.path.append(r"C:/Users/lacri/Desktop/UNI/MAGISTRALE/PhysicalMethodsForBiology/EPaci/Progetto/BozzeCodes")
import SimulationFunctionsOnly


def start_simulation():

    # Getting the user input from the entry widgets and setting initial conditions

    initial_C = float(entry1.get())  # Initial C
    initial_P = float(entry2.get())  # Initial P
    initial_conditions = [initial_C, initial_P]
    
    # Getting the user input from the entry widgets and setting constants

    kf_value = float(entry3.get())    # kf value
    kb_value = float(entry4.get())    # kb value
    kcat_value = float(entry5.get())  # kcat value   
    simkf, simkb, simkcat, simkM, simSt, simEt = SimulationFunctionsOnly.Constants(kf_value, kb_value, kcat_value)
    
    # Setting time span

    time_span = np.linspace(0, 50, 1000)

    # solving odes
    
    Full_Model_Solution = odeint (SimulationFunctionsOnly.Full_Model, initial_conditions, time_span, args=(simkf, simkb, simkcat, simSt, simEt))
    QSSAsolution = odeint (SimulationFunctionsOnly.QSSA, initial_conditions, time_span, args=(simkcat, simkM, simSt, simEt))
    tQSSAsolution = odeint (SimulationFunctionsOnly.tQSSA, initial_conditions, time_span, args=(simkcat, simkM, simSt, simEt))

    Full_Model_P_values = Full_Model_Solution[:, 1]
    QSSA_P_values = QSSAsolution[:, 1]
    tQSSA_P_values = tQSSAsolution[:, 1]

    # Ratio P/S according to the 3 models

    Full_Model_ratio_values = Full_Model_P_values / simSt
    QSSA_ratio_values = QSSA_P_values / simSt
    tQSSA_ratio_values = tQSSA_P_values / simSt

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
    print(f"Running simulation with parameters: Initial C={initial_C}, Initial P={initial_P}, kf={kf_value}, kb={kb_value}, kcat={kcat_value}")


root = tk.Tk()
root.title("Kinetic Reaction Simulator")

# Initial C
label1 = ttk.Label(root, text="Initial C:")
label1.grid(row=0, column=0)

entry1 = ttk.Entry(root)
entry1.grid(row=0, column=1)

# Initial P
label2 = ttk.Label(root, text="Initial P:")
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
start_button.grid(row=2, columnspan=6)

result_label = ttk.Label(root, text="Simulation Results:")
result_label.grid(row=3, columnspan=6)

#output_text = tk.Text(root, height=5, width=40)
#output_text.grid(row=4, columnspan=6)

root.mainloop()
