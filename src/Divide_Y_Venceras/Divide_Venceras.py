import tkinter as tk
from tkinter import ttk, messagebox
import itertools
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

#Diccionario IUPAC para las letras ambiguas
IUPAC_ADN = {
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

def validar_secuencia(secuencia):
    secuencia = secuencia.upper().strip().replace(" ", "")
    if not secuencia:
        return None, "Secuencia vacía."

    invalidos = [b for b in secuencia if b not in IUPAC_ADN]
    if invalidos:
        return None, f"Caracteres inválidos: {', '.join(sorted(set(invalidos)))}"

    resto = len(secuencia) % 3
    if resto != 0:
        secuencia = secuencia[:-resto]

    return secuencia, None


def divide_y_venceras(secuencia):
    if len(secuencia) <= 3:
        opciones = [list(IUPAC_ADN[b]) for b in secuencia]
        resultados = []
        for comb in itertools.product(*opciones):
            codon = "".join(comb)
            aa = TABLA_CODONES.get(codon, '?')
            if aa == '_':
                aa = ""
            resultados.append((codon, aa))
        return resultados

    mitad = len(secuencia) // 2
    rem = mitad % 3
    if rem != 0:
        mitad -= rem
    if mitad == 0:
        mitad = 3

    izq = divide_y_venceras(secuencia[:mitad])
    der = divide_y_venceras(secuencia[mitad:])
    combinados = []
    limite = 20000
    for (c1, a1) in izq:
        for (c2, a2) in der:
            combinados.append((c1 + c2, a1 + a2))
            if len(combinados) >= limite:
                return combinados

    return combinados

#GUI
class DashboardADN:
    def __init__(self, root):
        self.root = root
        self.root.title("Decode ways con Divide y venceras")
        self.root.geometry("1100x650")
        self.panel_izq = tk.Frame(root)
        self.panel_izq.pack(side="left", fill="both", expand=True, padx=10, pady=5)
        self.panel_der = tk.Frame(root)
        self.panel_der.pack(side="right", fill="y", padx=10, pady=5)

        #Entrada
        lf_entrada = tk.LabelFrame(self.panel_izq, text="Entrada de Secuencia")
        lf_entrada.pack(fill="x", pady=5)
        tk.Label(lf_entrada, text="Secuencia (IUPAC):").pack(side="left", padx=5)

        self.entry_secuencia = tk.Entry(lf_entrada, font=("Arial", 11))
        self.entry_secuencia.pack(side="left", fill="x", expand=True, padx=5, pady=10)
        self.entry_secuencia.bind('<Return>', lambda e: self.analizar())

        self.btn_analizar = tk.Button(lf_entrada, text="ANALIZAR", bg="#1976D2", fg="white", 
                                      font=("Arial", 9, "bold"), padx=15, command=self.analizar)
        self.btn_analizar.pack(side="right", padx=10, pady=10)
        lf_resultados = tk.LabelFrame(self.panel_izq, text="Resultados Detallados")
        lf_resultados.pack(fill="both", expand=True, pady=5)

        columnas = ("adn", "prot")
        self.tree = ttk.Treeview(lf_resultados, columns=columnas, show="headings")
        self.tree.heading("adn", text="Variante de ADN")
        self.tree.heading("prot", text="Proteína Resultante")

        self.tree.column("adn", width=300)
        self.tree.column("prot", width=200)

        barra = ttk.Scrollbar(lf_resultados, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=barra.set)
        barra.pack(side="right", fill="y")
        self.tree.pack(side="left", fill="both", expand=True)
        frame_status = tk.Frame(self.panel_izq)
        frame_status.pack(fill="x", pady=5)

        self.lbl_var = tk.Label(frame_status, text="Variantes ADN: 0", font=("Arial", 9))
        self.lbl_var.pack(side="left")

        self.lbl_prot = tk.Label(frame_status, text="Proteínas Únicas: 0", font=("Arial", 9, "bold"), fg="#2E7D32")
        self.lbl_prot.pack(side="right")

        # Glosario
        lf_glosario = tk.LabelFrame(self.panel_der, text="Glosario Aminoácidos", width=350, height=250)
        lf_glosario.pack(fill="x", pady=5)
        lf_glosario.pack_propagate(False)

        tree_ref = ttk.Treeview(lf_glosario, columns=("l", "n"), show="headings")
        tree_ref.heading("l", text="Letra")
        tree_ref.heading("n", text="Nombre")
        tree_ref.column("l", width=50, anchor="center")
        tree_ref.column("n", width=180)
        sb_ref = ttk.Scrollbar(lf_glosario, orient="vertical", command=tree_ref.yview)
        tree_ref.configure(yscroll=sb_ref.set)
        sb_ref.pack(side="right", fill="y")
        tree_ref.pack(side="left", fill="both", expand=True)

        ref_data = [
            ('A','Alanina'),('C','Cisteína'),('D','Ác. Aspártico'),('E','Ác. Glutámico'),
            ('F','Fenilalanina'),('G','Glicina'),('H','Histidina'),('I','Isoleucina'),
            ('K','Lisina'),('L','Leucina'),('M','Metionina'),('N','Asparagina'),
            ('P','Prolina'),('Q','Glutamina'),('R','Arginina'),('S','Serina'),
            ('T','Treonina'),('V','Valina'),('W','Triptófano'),('Y','Tirosina'),('_','STOP')
        ]
        for item in ref_data:
            tree_ref.insert("", "end", values=item)

        
        self.lf_grafica = tk.LabelFrame(self.panel_der, text="Rendimiento Algorítmico")
        self.lf_grafica.pack(fill="both", expand=True, pady=5)

        self.generar_grafica()

    
    #Graficación
    def generar_grafica(self):
        x = [1, 5, 10, 15, 20]
        y_fb = [4**i for i in x]
        y_dyv = [i * 50 for i in x]

        self.fig, self.ax = plt.subplots(figsize=(4, 3.5), dpi=100)
        self.ax.plot(x, y_fb, label="Fuerza Bruta (4ⁿ)")
        self.ax.plot(x, y_dyv, label="Divide y Vencerás (n)")

        self.ax.set_title("Complejidad Comparada")
        self.ax.set_xlabel("Cantidad de 'N'")
        self.ax.set_ylabel("Operaciones")
        self.ax.legend()

        canvas = FigureCanvasTkAgg(self.fig, master=self.lf_grafica)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

    def analizar(self):
        sec = self.entry_secuencia.get()
        sec, err = validar_secuencia(sec)

        if err:
            messagebox.showerror("Error", err)
            return

        inicio = time.time()
        resultados = divide_y_venceras(sec)
        duracion = time.time() - inicio

        self.tree.delete(*self.tree.get_children())

        prot_unicas = set()
        for adn, prot in resultados:
            prot_unicas.add(prot)
            self.tree.insert("", "end", values=(adn, prot))

        self.lbl_var.config(text=f"Variantes ADN: {len(resultados)}")
        self.lbl_prot.config(text=f"Proteínas Únicas: {len(prot_unicas)}")

        messagebox.showinfo("Completado", f"Análisis terminado en {duracion:.3f} segundos.")

root = tk.Tk()
app = DashboardADN(root)
root.mainloop()
