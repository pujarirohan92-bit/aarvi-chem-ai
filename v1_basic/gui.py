import tkinter as tk
from tkinter import messagebox
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import ImageTk
from engine import predict_reaction


def run_reaction():
    r1 = entry_r1.get().strip()
    r2 = entry_r2.get().strip()

    if not r1 or not r2:
        messagebox.showerror("Error", "Please enter both reactants")
        return

    product = predict_reaction(r1, r2)
    product_label.config(text=f"Product SMILES:\n{product}")

    try:
        mol = Chem.MolFromSmiles(product)
        img = Draw.MolToImage(mol, size=(300, 300))
        img_tk = ImageTk.PhotoImage(img)
        canvas.image = img_tk
        canvas.create_image(150, 150, image=img_tk)
    except:
        messagebox.showinfo("Info", "Structure not drawable")


# ---------- GUI ----------
root = tk.Tk()
root.title("Aarvi chem Ai")
root.geometry("420x600")

tk.Label(root, text="Reactant 1 (SMILES)").pack()
entry_r1 = tk.Entry(root, width=45)
entry_r1.pack()

tk.Label(root, text="Reactant 2 (SMILES)").pack()
entry_r2 = tk.Entry(root, width=45)
entry_r2.pack()

tk.Button(root, text="Predict Product", command=run_reaction,
          bg="green", fg="white").pack(pady=10)

product_label = tk.Label(root, text="Product SMILES:", wraplength=380)
product_label.pack(pady=10)

canvas = tk.Canvas(root, width=300, height=300, bg="white")
canvas.pack()

root.mainloop()
