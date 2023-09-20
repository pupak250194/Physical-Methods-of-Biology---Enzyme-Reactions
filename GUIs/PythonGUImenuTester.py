
import tkinter as tk
from tkinter import ttk
import subprocess

def show_gui(selected_gui):
    gui_file = f"{selected_gui}.py"
    gui_path = f"C:\\Users\\39333\\Desktop\\Physical-Methods-of-Biology---Enzyme-Reactions-main\\GUIs\\{gui_file}"

    try:
        subprocess.Popen(["python", gui_path])
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
label = ttk.Label(root, text="Select a GUI to start simulation:")
label.pack(pady=10)

# Create a combobox (dropdown menu) to select GUI
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
combo.current(0)  # Set the default selection

# Set the width of the dropdown menu
combo.config(width=max(len(name) for name in options) + 2)

# Create a "Start Simulation" button
start_button = ttk.Button(root, text="Start Simulation", command=on_select)
start_button.pack(pady=10)

# Start the main event loop
root.mainloop()