import tkinter as tk
from tkinter import ttk, messagebox
import itertools
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

#Diccionario IUPAC para las letras de codones ambiguos
IUPAC_DNA = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
    'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
}

TABLA_CODONES = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

PESOS_AMINOACIDOS = {
    'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 
    'Q': 146, 'E': 147, 'G': 75,  'H': 155, 'I': 131, 
    'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115, 
    'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117, 
    '_': 0, '?': 0
}

def validar_secuencia(secuencia):
    secuencia = secuencia.upper().strip().replace(" ", "")
    if not secuencia:
        return None, "Secuencia vacía."

    invalidos = [b for b in secuencia if b not in IUPAC_DNA]
    if invalidos:
        return None, f"Caracteres inválidos: {', '.join(sorted(set(invalidos)))}"

    resto = len(secuencia) % 3
    if resto != 0:
        secuencia = secuencia[:-resto]

    return secuencia, None


def algoritmo_voraz_pesado(secuencia):
    dna_final = ""
    proteina_final = ""
    peso_total = 0
    historial = []
    for i in range(0, len(secuencia), 3):
        codon_ambiguo = secuencia[i:i+3]

        opciones = [list(IUPAC_DNA[b]) for b in codon_ambiguo]
        combinaciones = itertools.product(*opciones)

        mejor_codon = ""
        mejor_aa = ""
        max_peso = -1

        for comb in combinaciones:
            codon = "".join(comb)
            aa = TABLA_CODONES.get(codon, '?')
            peso = PESOS_AMINOACIDOS.get(aa, 0)

            if peso > max_peso:
                max_peso = peso
                mejor_codon = codon
                mejor_aa = aa

        dna_final += mejor_codon
        proteina_final += mejor_aa
        peso_total += max_peso

        historial.append((codon_ambiguo, mejor_codon, mejor_aa, max_peso))

    return dna_final, proteina_final, peso_total, historial


#GUI
class DashboardVoraz:
    def __init__(self, root):
        self.root = root
        self.root.title("Decode ways genetic con voraces")
        self.root.geometry("1100x650")
        self.panel_izq = tk.Frame(root)
        self.panel_izq.pack(side="left", fill="both", expand=True, padx=10, pady=5)
        self.panel_der = tk.Frame(root)
        self.panel_der.pack(side="right", fill="y", padx=10, pady=5)
        lf_entrada = tk.LabelFrame(self.panel_izq, text="Buscador de Proteína Más Pesada")
        lf_entrada.pack(fill="x", pady=5)
        tk.Label(lf_entrada, text="Secuencia (IUPAC):").pack(side="left", padx=5)
        self.entrada_secuencia = tk.Entry(lf_entrada, font=("Arial", 11))
        self.entrada_secuencia.pack(side="left", fill="x", expand=True, padx=5, pady=10)
        self.entrada_secuencia.bind('<Return>', lambda e: self.ejecutar())
        self.btn_ejecutar = tk.Button(
            lf_entrada,
            text="EJECUTAR GREEDY",
            bg="#F57C00",
            fg="white",
            font=("Arial", 9, "bold"),
            padx=15,
            command=self.ejecutar
        )
        self.btn_ejecutar.pack(side="right", padx=10, pady=10)
        lf_result = tk.LabelFrame(self.panel_izq, text="Ruta Voraz (Decisión Local Paso a Paso)")
        lf_result.pack(fill="both", expand=True, pady=5)
        columnas = ("input", "eleccion", "aa", "peso")
        self.tree = ttk.Treeview(lf_result, columns=columnas, show="headings")
        self.tree.heading("input", text="Codón Original")
        self.tree.heading("eleccion", text="Decisión Voraz (ADN)")
        self.tree.heading("aa", text="Aminoácido")
        self.tree.heading("peso", text="Peso (Daltons)")

        self.tree.column("input", width=100, anchor="center")
        self.tree.column("eleccion", width=150, anchor="center")
        self.tree.column("aa", width=150, anchor="center")
        self.tree.column("peso", width=100, anchor="center")

        sb = ttk.Scrollbar(lf_result, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=sb.set)
        sb.pack(side="right", fill="y")
        self.tree.pack(side="left", fill="both", expand=True)
        frame_status = tk.Frame(self.panel_izq)
        frame_status.pack(fill="x", pady=5)
        self.lbl_resumen = tk.Label(
            frame_status,
            text="Resultado Final: ...",
            font=("Arial", 10, "bold"),
            fg="#E65100"
        )
        self.lbl_resumen.pack(side="left")
        lf_glosario = tk.LabelFrame(self.panel_der, text="Referencia de Pesos", width=350, height=250)
        lf_glosario.pack(fill="x", pady=5)
        lf_glosario.pack_propagate(False)
        tree_ref = ttk.Treeview(lf_glosario, columns=("aa", "w"), show="headings")
        tree_ref.heading("aa", text="Aminoácido")
        tree_ref.heading("w", text="Peso Mol.")
        tree_ref.column("aa", width=150)
        tree_ref.column("w", width=80, anchor="center")

        sb_ref = ttk.Scrollbar(lf_glosario, orient="vertical", command=tree_ref.yview)
        tree_ref.configure(yscroll=sb_ref.set)
        sb_ref.pack(side="right", fill="y")
        tree_ref.pack(side="left", fill="both", expand=True)
        for aa, peso in sorted(PESOS_AMINOACIDOS.items(), key=lambda x: x[1], reverse=True):
            tree_ref.insert("", "end", values=(aa, peso))

    def ejecutar(self):
        sec = self.entrada_secuencia.get()
        sec, err = validar_secuencia(sec)

        if err:
            messagebox.showerror("Error", err)
            return

        dna, prot, peso, historial = algoritmo_voraz_pesado(sec)

        for item in self.tree.get_children():
            self.tree.delete(item)

        for codon_in, codon_out, aa, w in historial:
            self.tree.insert("", "end", values=(codon_in, codon_out, aa, w))

        self.lbl_resumen.config(
            text=f"Resultado Final → ADN: {dna}  |  Proteína: {prot}  |  Peso Total: {peso}"
        )

root = tk.Tk()
app = DashboardVoraz(root)
root.mainloop()
