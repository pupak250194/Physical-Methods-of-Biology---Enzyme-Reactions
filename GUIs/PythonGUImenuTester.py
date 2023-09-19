
import tkinter as tk
from tkinter import ttk
import subprocess

def show_gui(selected_gui):
    gui_file = f"{selected_gui}.py"
    gui_path = f"path\\to\\GUIs\\{gui_file}"

    try:
        subprocess.Popen(["python", gui_path])
    except FileNotFoundError:
        print(f"GUI {selected_gui} not found.")

def on_select(event):
    selected_gui = combo.get()
    show_gui(selected_gui)
    
# Main window
root = tk.Tk()
root.title("GUI Selector for Simulation")
root.geometry("350x200")

# Main label
label = ttk.Label(root, text="Select a GUI to start simulation:")
label.pack(pady=10)

# Dropdown Menu

options = [
    "PythonTesterGUI",
    "PythonTesterGUI_Switch_steady-state",
    "PythonTesterGUI-cme-averages",
    "PythonTesterGUI-cme-completion-times",
    "PythonTesterGUI-cme-concentrations",
    "PythonTesterGUI-cme-stationary",
    "PythonTesterGUI-gillespie-averages",
    "PythonTesterGUI-gillespie-completion-times",
    "PythonTesterGUI-gillespie-concentrations",
    "PythonTesterGUI-gillespie-stationary"
]  
combo = ttk.Combobox(root, values=options)
combo.pack(pady=10)

combo.current(0) 
combo.bind("<<ComboboxSelected>>", on_select)

# Dropdown menu width
combo.config(width=max(len(name) for name in options) + 2)

# Start button
start_button = ttk.Button(root, text="Start Simulation", command=lambda: on_select(None))
start_button.pack(pady=10)

# Start the main event loop
root.mainloop()
