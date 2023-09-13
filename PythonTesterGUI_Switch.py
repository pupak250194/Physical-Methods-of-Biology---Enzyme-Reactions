import tkinter as tk
import matplotlib.pyplot as plt
import sys
import numpy as np

from tkinter import ttk
from scipy.integrate import odeint
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


sys.path.append(r"C:/Users/lacri/Desktop/UNI/MAGISTRALE/PhysicalMethodsForBiology/EPaci/Progetto/BozzeCodes")
import SimulationFunctionsOnly_ES


def start_simulation():

    # ESTABLISHING INITIAL CONDITIONS AND PARAMETERS

    initial_S = float(entry1Switch.get()) 
    initial_E = float(entry2Switch.get()) 
    initial_D = float(entry3Switch.get()) 

    initial_Sp = 0
    initial_C = 0
    initial_Cp = 0
    
    initial_input = [initial_S, initial_E, initial_D]
    initial_conditions = [initial_Sp, initial_C, initial_Cp]

    
    # Getting the user input from the entry widgets and setting constants

    ke_value = float(entry4Switch.get())    # ke value
    kbe_value = float(entry5Switch.get())    # kbe value
    kfe_value = float(entry6Switch.get())  # kfe value 
    kbd_value = float(entry7Switch.get())  # kbd value  
    kfd_value = float(entry8Switch.get())  # kfd value  
    kd_value = float(entry9Switch.get())  # kfd value 

    simkME, simkMD = SimulationFunctionsOnly_ES.Constants_Switch (ke_value, kbe_value, kfe_value, kbd_value, kfd_value, kd_value)
    
    # SETTING TIME SPAN

    # Getting the user input for simulation time
    simulation_time_Switch = float(entry10Switch.get())  # Simulation time
    
    # Creating the time span based on user-defined simulation time
    time_span_Switch = np.linspace(0, simulation_time_Switch, 1000)

    # solving odes


    Full_Model_Solution_Switch = odeint (SimulationFunctionsOnly_ES.Full_Model_Switch, initial_conditions, time_span_Switch, args=(*initial_input, kfd_value, kbd_value, ke_value, kfe_value, kbe_value, kd_value))
    QSSAsolution_Switch = odeint (SimulationFunctionsOnly_ES.QSSA_Switch, initial_conditions[0], time_span_Switch, args=(*initial_input, ke_value, kd_value, simkME, simkMD))
    tQSSAsolution_Switch = odeint (SimulationFunctionsOnly_ES.tQSSA_Switch, initial_conditions[0], time_span_Switch, args=(*initial_input, ke_value, kd_value, simkME, simkMD))

    Full_Model_Switch_values = Full_Model_Solution_Switch[:, 0] + Full_Model_Solution_Switch[:, 2]
    QSSA_Switch_values = QSSAsolution_Switch[:, 0]
    tQSSA_Switch_values = tQSSAsolution_Switch[:, 0]

    # Ratio Sp/S according to the 3 models

    Full_Model_Switch_ratio_values = Full_Model_Switch_values / initial_S
    QSSA_Switch_ratio_values = QSSA_Switch_values / initial_S
    tQSSA_Switch_ratio_values = tQSSA_Switch_values / initial_S

    # Plot as a function of ke, kd, E and D ...


    # Plots

    # Create a new figure for the plot
    fig_Switch, ax_Switch = plt.subplots()
    
    # Plot the data
    ax_Switch.plot(time_span_Switch, Full_Model_Switch_ratio_values, label='Full Model', color='yellow', linewidth=3.5)
    ax_Switch.plot(time_span_Switch, QSSA_Switch_ratio_values , label='QSSA - Sp/Stot', linewidth=2)
    ax_Switch.plot(time_span_Switch, tQSSA_Switch_ratio_values , label='tQSSA - Sp/S-hat', linestyle='--', color='black')

    ax_Switch.set_xlabel('Time')
    ax_Switch.set_ylabel('Sp/S')
    ax_Switch.legend()

    # Embed the plot in the Tkinter window
    canvasSwitch = FigureCanvasTkAgg(fig_Switch, master=rootSwitch)
    canvasSwitch.get_tk_widget().grid(row=5, columnspan=6) 

    #plt.show()


    # Replace the print statement with your simulation code
    #print(f"Running simulation with parameters: Initial S={initial_S}, Initial E={initial_E}, kf={kf_value}, kb={kb_value}, kcat={kcat_value}")


rootSwitch = tk.Tk()
rootSwitch.title("Kinetic Reaction Simulator - Switch")

# Initial S
label1Switch = ttk.Label(rootSwitch, text="Initial S:")
label1Switch.grid(row=0, column=0)

entry1Switch = ttk.Entry(rootSwitch)
entry1Switch.grid(row=0, column=1)

# Initial E
label2Switch = ttk.Label(rootSwitch, text="Initial E:")
label2Switch.grid(row=0, column=2)

entry2Switch = ttk.Entry(rootSwitch)
entry2Switch.grid(row=0, column=3)

# Initial D
label3Switch = ttk.Label(rootSwitch, text="Initial D:")
label3Switch.grid(row=0, column=4)

entry3Switch = ttk.Entry(rootSwitch)
entry3Switch.grid(row=0, column=5)

# ke value
label4Switch = ttk.Label(rootSwitch, text="ke value:")
label4Switch.grid(row=1, column=0)

entry4Switch = ttk.Entry(rootSwitch)
entry4Switch.grid(row=1, column=1)

# kbe value
label5Switch = ttk.Label(rootSwitch, text="kbe value:")
label5Switch.grid(row=1, column=2)

entry5Switch = ttk.Entry(rootSwitch)
entry5Switch.grid(row=1, column=3)

# kfe value
label6Switch = ttk.Label(rootSwitch, text="kfe value:")
label6Switch.grid(row=1, column=4)

entry6Switch = ttk.Entry(rootSwitch)
entry6Switch.grid(row=1, column=5)

# kbd value
label7Switch = ttk.Label(rootSwitch, text="kbd value:")
label7Switch.grid(row=2, column=0)

entry7Switch = ttk.Entry(rootSwitch)
entry7Switch.grid(row=2, column=1)

# kfd value
label8Switch = ttk.Label(rootSwitch, text="kfd value:")
label8Switch.grid(row=2, column=2)

entry8Switch = ttk.Entry(rootSwitch)
entry8Switch.grid(row=2, column=3)

# kfd value
label9Switch = ttk.Label(rootSwitch, text="kd value:")
label9Switch.grid(row=2, column=4)

entry9Switch = ttk.Entry(rootSwitch)
entry9Switch.grid(row=2, column=5)


# Label for Simulation Time
label10Switch = ttk.Label(rootSwitch, text="Simulation Time:")
label10Switch.grid(row=3, column=0)

entry10Switch = ttk.Entry(rootSwitch)
entry10Switch.grid(row=3, column=1)

start_button = ttk.Button(rootSwitch, text="Start Simulation", command=start_simulation)
start_button.grid(row=3, columnspan=6)

result_labelSwitch= ttk.Label(rootSwitch, text="CONCENTRATION PLOTS")
result_labelSwitch.grid(row=4, columnspan=6)

#output_text = tk.Text(root, height=5, width=40)
#output_text.grid(row=4, columnspan=6)

rootSwitch.mainloop()
