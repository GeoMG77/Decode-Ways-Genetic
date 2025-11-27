import tkinter as tk 
from tkinter import ttk, messagebox
import itertools
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


#Lógica para fuerza bruta
def analizar_variantes_adn_proteina(secuencia_ambigua):
    #Diccionario IUPAC para las letras ambiguas
    iupac_adn = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
        'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
        'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
    }
    #Diccionario de codones a aminoacidos
    tabla_codones = {
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

    secuencia_ambigua = secuencia_ambigua.upper().strip()
    if not secuencia_ambigua:
        return set(), [], 0 

    caracteres_invalidos = [b for b in secuencia_ambigua if b not in iupac_adn]
    if caracteres_invalidos:
        raise ValueError(f"Caracteres no válidos: {', '.join(sorted(set(caracteres_invalidos)))}")

    residuo = len(secuencia_ambigua) % 3
    if residuo != 0:
        secuencia_ambigua = secuencia_ambigua[:-residuo]

    opciones_posibles = [list(iupac_adn[base]) for base in secuencia_ambigua]
    todas_las_secuencias = itertools.product(*opciones_posibles)
    proteinas_unicas = set()
    lista_limitada = []
    total = 0
    limite_visual = 10000
    for nucleotidos in todas_las_secuencias:
        total += 1
        adn_generado = "".join(nucleotidos)
        proteina_actual = ""
        for i in range(0, len(adn_generado), 3):
            codon = adn_generado[i:i+3]
            aminoacido = tabla_codones.get(codon, '?')
            if aminoacido == '_':
                break
            proteina_actual += aminoacido

        proteinas_unicas.add(proteina_actual)
        if len(lista_limitada) < limite_visual:
            lista_limitada.append((adn_generado, proteina_actual))

    return proteinas_unicas, lista_limitada, total

#Interfaz de Usuario GUI
class AplicacionDNA:
    def __init__(self, root):
        self.root = root
        self.root.title("Decode ways con Fuerza Bruta")
        self.root.geometry("1100x650")
        self.panel_izq = tk.Frame(root)
        self.panel_izq.pack(side="left", fill="both", expand=True, padx=10, pady=5)
        self.panel_der = tk.Frame(root)
        self.panel_der.pack(side="right", fill="y", padx=10, pady=5)
        lf_entrada = tk.LabelFrame(self.panel_izq, text="Entrada de Secuencia")
        lf_entrada.pack(fill="x", pady=5)
        tk.Label(lf_entrada, text="Secuencia (IUPAC):").pack(side="left", padx=5)
        self.entrada_seq = tk.Entry(lf_entrada, font=("Arial", 11))
        self.entrada_seq.pack(side="left", fill="x", expand=True, padx=5, pady=10)
        self.entrada_seq.bind('<Return>', lambda e: self.ejecutar_analisis())
        btn_analizar = tk.Button(lf_entrada, text="ANALIZAR", bg="#1976D2", fg="white", 
                                 font=("Arial", 9, "bold"), padx=15, command=self.ejecutar_analisis)
        btn_analizar.pack(side="right", padx=10, pady=10)
        lf_resultados = tk.LabelFrame(self.panel_izq, text="Generación Exhaustiva")
        lf_resultados.pack(fill="both", expand=True, pady=5)
        columns = ("adn", "proteina")
        self.tree = ttk.Treeview(lf_resultados, columns=columns, show="headings")
        self.tree.heading("adn", text="Variante de ADN")
        self.tree.heading("proteina", text="Proteína Resultante")
        self.tree.column("adn", width=350)
        self.tree.column("proteina", width=250)
        scrollbar = ttk.Scrollbar(lf_resultados, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        self.tree.pack(side="left", fill="both", expand=True)

        frame_stats = tk.Frame(self.panel_izq)
        frame_stats.pack(fill="x", pady=5)
        self.lbl_total_adn = tk.Label(frame_stats, text="Total Variantes: 0", font=("Arial", 9))
        self.lbl_total_adn.pack(side="left")
        self.lbl_total_prot = tk.Label(frame_stats, text="Proteínas Únicas: 0", font=("Arial", 9, "bold"), fg="#2E7D32")
        self.lbl_total_prot.pack(side="right")

        lf_glosario = tk.LabelFrame(self.panel_der, text="Glosario Aminoácidos", width=350, height=300)
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
            ('A', 'Alanina'), ('C', 'Cisteína'), ('D', 'Ác. Aspártico'),
            ('E', 'Ác. Glutámico'), ('F', 'Fenilalanina'), ('G', 'Glicina'),
            ('H', 'Histidina'), ('I', 'Isoleucina'), ('K', 'Lisina'),
            ('L', 'Leucina'), ('M', 'Metionina'), ('N', 'Asparagina'),
            ('P', 'Prolina'), ('Q', 'Glutamina'), ('R', 'Arginina'),
            ('S', 'Serina'), ('T', 'Treonina'), ('V', 'Valina'),
            ('W', 'Triptófano'), ('Y', 'Tirosina'), ('_', 'STOP')
        ]
        for item in ref_data:
            tree_ref.insert("", "end", values=item)
        #Grafica del algoritmo
        self.lf_grafica = tk.LabelFrame(self.panel_der, text="Complejidad Fuerza Bruta")
        self.lf_grafica.pack(fill="both", expand=True, pady=5)
        self.generar_grafica_teorica()

    def generar_grafica_teorica(self):
        fig, ax = plt.subplots(figsize=(4, 3.5), dpi=100)
        x = list(range(1, 11))
        y = [4**n for n in x]
        ax.plot(x, y)
        ax.set_title("Crecimiento Exponencial de Fuerza Bruta")
        ax.set_xlabel("Longitud De secuencia(n)")
        ax.set_ylabel("Variantes posibles")

        canvas = FigureCanvasTkAgg(fig, master=self.lf_grafica)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)


    def ejecutar_analisis(self):
        self.tree.delete(*self.tree.get_children())
        secuencia = self.entrada_seq.get().strip()
        try:
            proteinas_unicas, lista_adn, total = analizar_variantes_adn_proteina(secuencia)
        except ValueError as e:
            messagebox.showerror("Error", str(e))
            return

        for adn, prot in lista_adn:
            self.tree.insert("", "end", values=(adn, prot))

        self.lbl_total_adn.config(text=f"Total Variantes: {total}")
        self.lbl_total_prot.config(text=f"Proteínas Únicas: {len(proteinas_unicas)}")

root = tk.Tk()
app = AplicacionDNA(root)
root.mainloop()
