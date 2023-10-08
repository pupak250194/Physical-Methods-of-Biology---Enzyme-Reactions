
import tkinter as tk
from tkinter import ttk
import subprocess

def show_gui(selected_gui):
    try:
        subprocess.Popen(["python", f"{selected_gui}.py"])
    except FileNotFoundError:
        print(f"GUI {selected_gui} not found.")

def on_select():
    selected_gui = combo.get()
    show_gui(selected_gui)

# Create the main window
root = tk.Tk()
root.title("GUI Selector")
root.geometry("350x200")

# Create a label
label = ttk.Label(root, text="Select an experiment to start simulation:")
label.pack(pady=10)

# Create a combobox (dropdown menu) to select GUI
options = [
    "GUI-ODE-single-substrate",
    "GUI-ODE-gk-Switch",
    "GUI-ODE-gk-Switch-steady-state",
    "GUI-CME-averages",
    "GUI-CME-completion-times",
    "GUI-CME-concentrations",
    "GUI-CME-stationary",
    "GUI-gillespie-averages",
    "GUI-gillespie-completion-times",
    "GUI-gillespie-concentrations",
    "GUI-gillespie-stationary",
    "GUI-Averages-CMEvsGillespie"
]  

combo = ttk.Combobox(root, values=options)
combo.pack(pady=10)
combo.current(0)  # Set the default selection

# Set the width of the dropdown menu
combo.config(width=max(len(name) for name in options) + 2)

start_button = ttk.Button(root, text="Start simulation", command=on_select)
start_button.pack(pady=10)

# Start the main event loop
root.mainloop()