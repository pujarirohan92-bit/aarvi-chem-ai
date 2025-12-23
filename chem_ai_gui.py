import tkinter as tk
from tkinter import messagebox

def predict_reaction():
    reactants = entry.get().lower().replace(" ", "")

    if "hcl" in reactants and "naoh" in reactants:
        result = "Neutralization Reaction\nHCl + NaOH → NaCl + H₂O"
    elif "ch4" in reactants or "methane" in reactants:
        result = "Combustion Reaction\nCH₄ + O₂ → CO₂ + H₂O"
    elif "zn" in reactants and "hcl" in reactants:
        result = "Single Displacement\nZn + HCl → ZnCl₂ + H₂"
    elif "caco3" in reactants:
        result = "Decomposition\nCaCO₃ → CaO + CO₂"
    else:
        result = "Reaction not found ❌\nTry: HCl+NaOH, CH4, Zn+HCl"

    messagebox.showinfo("Result", result)

root = tk.Tk()
root.title("Chem AI – Reaction Predictor")
root.geometry("420x260")

tk.Label(root, text="Chemical Reaction Predictor",
         font=("Arial", 14, "bold")).pack(pady=10)

tk.Label(root, text="Enter Reactants (example: HCl + NaOH)").pack()
entry = tk.Entry(root, width=40)
entry.pack(pady=6)

tk.Button(root, text="Predict Reaction",
          command=predict_reaction,
          bg="#1e8e3e", fg="white",
          font=("Arial", 10, "bold")).pack(pady=15)

root.mainloop()
