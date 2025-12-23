import tkinter as tk
from engine import predict_reactions

def run_retro():
    smiles = entry.get().strip()
    output.delete("1.0", "end")

    if not smiles:
        output.insert("end", "Please enter a SMILES string\n")
        return

    reactions, fgs = predict_reactions(smiles)

    output.insert("end", "Detected Functional Groups:\n")
    output.insert("end", ", ".join(fgs) + "\n\n")

    if not reactions:
        output.insert("end", "No matching reactions found\n")
        return

    output.insert("end", "Possible Reactions (LEVEL-2):\n\n")

    for i, r in enumerate(reactions, start=1):
        output.insert("end", f"{i}. {r['name']}\n")
        output.insert("end", f"   Retro logic : {r['retro']}\n")
        output.insert("end", f"   Conditions : {r['conditions']}\n\n")

# ---------------- GUI LAYOUT ----------------

root = tk.Tk()
root.title("Aarvi Chem AI – LEVEL 2 Reaction Predictor")
root.geometry("650x450")

tk.Label(
    root,
    text="Enter Product SMILES",
    font=("Arial", 12, "bold")
).pack(pady=5)

entry = tk.Entry(root, width=60, font=("Arial", 11))
entry.pack(pady=5)

tk.Button(
    root,
    text="Predict Reactions",
    command=run_retro,
    bg="#2E86C1",
    fg="white",
    font=("Arial", 11, "bold")
).pack(pady=10)

output = tk.Text(root, height=18, width=75, font=("Consolas", 10))
output.pack(pady=5)

root.mainloop()
