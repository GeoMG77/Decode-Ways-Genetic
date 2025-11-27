import tkinter as tk
from tkinter import ttk, messagebox
import itertools
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

#Diccionario para ambiguedades
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


#Logica de la validacion y generacion dinamica
def validar_secuencia(secuencia):
    #Normaliza y valida la secuencia (remueve espacios, mayúsculas).
    #Devuelve (secuencia_valida, None) o (None, mensaje_error).
    sec = secuencia.upper().strip().replace(" ", "")
    if not sec:
        return None, "Secuencia vacía."
    invalidos = [b for b in sec if b not in IUPAC_DNA]
    if invalidos:
        return None, f"Caracteres inválidos: {', '.join(sorted(set(invalidos)))}"
    resto = len(sec) % 3
    if resto != 0:
        sec = sec[:-resto]
    return sec, None


def generar_variantes_dinamica(secuencia, limite_seguridad=500000):
    #Construcción iterativa (enfoque dinámico) que mantiene una lista de variantes acumuladas (adn_acumulado, prot_acumulado)
    #y la va expandiendo codón por codón, si la lista intermedia supera  el limite_seguridad, devuelve lo acumulado y True.
    variantes_actuales = [("", "")]
    longitud = len(secuencia)

    for i in range(0, longitud, 3):
        triplete = secuencia[i:i+3]
        # Generamos todas las combinaciones de ese codón según IUPAC
        opciones_bases = [list(IUPAC_DNA[b]) for b in triplete]
        codones_generados = []
        for comb in itertools.product(*opciones_bases):
            codon = "".join(comb)
            aa = TABLA_CODONES.get(codon, "?")
            if aa == "_":
                aa = ""  
            codones_generados.append((codon, aa))

        nuevas_variantes = []
        for adn_prev, prot_prev in variantes_actuales:
            for adn_nuevo, prot_nuevo in codones_generados:
                if len(nuevas_variantes) >= limite_seguridad:
                    return nuevas_variantes, True
                nuevas_variantes.append((adn_prev + adn_nuevo, prot_prev + prot_nuevo))

        variantes_actuales = nuevas_variantes

    return variantes_actuales, False


#Interfaz Gráfica
class DashboardDPList:
    def __init__(self, root):
        self.root = root
        self.root.title("Decode ways con Programación Dinámica")
        self.root.geometry("1100x650")
        self.TOPE_CALCULO = 500000
        self.TOPE_VISUAL = 5000

        self.panel_izq = tk.Frame(root)
        self.panel_izq.pack(side="left", fill="both", expand=True, padx=10, pady=5)
        self.panel_der = tk.Frame(root)
        self.panel_der.pack(side="right", fill="y", padx=10, pady=5)

        lf_entrada = tk.LabelFrame(self.panel_izq, text=f"Generador Seguro (Tope: {self.TOPE_CALCULO:,})")
        lf_entrada.pack(fill="x", pady=5)
        tk.Label(lf_entrada, text="Secuencia (IUPAC):").pack(side="left", padx=5)

        self.entry_seq = tk.Entry(lf_entrada, font=("Arial", 11))
        self.entry_seq.pack(side="left", fill="x", expand=True, padx=5, pady=10)
        self.entry_seq.bind('<Return>', lambda e: self.analizar())

        self.btn_analizar = tk.Button(
            lf_entrada, text="GENERAR LISTA", bg="#5E35B1", fg="white",
            font=("Arial", 9, "bold"), padx=15, command=self.analizar
        )
        self.btn_analizar.pack(side="right", padx=10, pady=10)
        # Resultados 
        lf_resultados = tk.LabelFrame(self.panel_izq, text="Resultados Generados")
        lf_resultados.pack(fill="both", expand=True, pady=5)
        columnas = ("adn", "prot")
        self.tree = ttk.Treeview(lf_resultados, columns=columnas, show="headings")
        self.tree.heading("adn", text="Variante de ADN")
        self.tree.heading("prot", text="Proteína Resultante")
        self.tree.column("adn", width=300)
        self.tree.column("prot", width=200)

        sb = ttk.Scrollbar(lf_resultados, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=sb.set)
        sb.pack(side="right", fill="y")
        self.tree.pack(side="left", fill="both", expand=True)

        frame_status = tk.Frame(self.panel_izq)
        frame_status.pack(fill="x", pady=5)
        self.lbl_adn = tk.Label(frame_status, text="Variantes: 0", font=("Arial", 9))
        self.lbl_adn.pack(side="left")
        self.lbl_prot = tk.Label(frame_status, text="Proteínas Únicas: 0", font=("Arial", 9, "bold"), fg="#2E7D32")
        self.lbl_prot.pack(side="right")

        #glosario y la gráfica
        lf_glosario = tk.LabelFrame(self.panel_der, text="Glosario", width=350, height=250)
        lf_glosario.pack(fill="x", pady=5)
        lf_glosario.pack_propagate(False)

        self.tree_ref = ttk.Treeview(lf_glosario, columns=("l", "n"), show="headings")
        self.tree_ref.heading("l", text="L")
        self.tree_ref.heading("n", text="Aminoácido")
        self.tree_ref.column("l", width=40, anchor="center")
        self.tree_ref.column("n", width=150)

        sb_ref = ttk.Scrollbar(lf_glosario, orient="vertical", command=self.tree_ref.yview)
        self.tree_ref.configure(yscroll=sb_ref.set)
        sb_ref.pack(side="right", fill="y")
        self.tree_ref.pack(side="left", fill="both", expand=True)

        glosario = [
            ('A','Alanina'),('C','Cisteína'),('D','Ác. Aspártico'),
            ('E','Ác. Glutámico'),('F','Fenilalanina'),('G','Glicina'),
            ('H','Histidina'),('I','Isoleucina'),('K','Lisina'),
            ('L','Leucina'),('M','Metionina'),('N','Asparagina'),
            ('P','Prolina'),('Q','Glutamina'),('R','Arginina'),
            ('S','Serina'),('T','Treonina'),('V','Valina'),
            ('W','Triptófano'),('Y','Tirosina'),('_','STOP')
        ]
        for item in glosario:
            self.tree_ref.insert("", "end", values=item)
        self.lf_grafica = tk.LabelFrame(self.panel_der, text="Curva de Complejidad")
        self.lf_grafica.pack(fill="both", expand=True, pady=5)
        self.generar_grafica_teorica()

    def generar_grafica_teorica(self):

        fig, ax = plt.subplots(figsize=(4, 3.5), dpi=100)
        x = [1, 5, 8, 9, 10]
        y_teorico = [4**i for i in x]
        y_limite = [self.TOPE_CALCULO for _ in x]
        ax.plot(x, y_teorico, color='#D32F2F', linestyle='--', marker='o', label='Variantes Reales')
        ax.plot(x, y_limite, color='#1976D2', linestyle='-', label=f'Tope Seguridad ({self.TOPE_CALCULO:,})', linewidth=2)
        ax.set_yscale('log')
        ax.set_title("Crecimiento vs Límite", fontsize=10)
        ax.set_xlabel("Letras 'N'", fontsize=8)
        ax.set_ylabel("Variantes (Log)", fontsize=8)
        ax.grid(True, which="both", ls="--", alpha=0.3)
        ax.legend(fontsize=8, loc="upper left")
        ax.tick_params(labelsize=8)

        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=self.lf_grafica)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
    # Método principal 
    def analizar(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.root.update()
        raw = self.entry_seq.get()
        seq, err = validar_secuencia(raw)
        if err:
            messagebox.showerror("Error", err)
            return
        lista_resultados, hubo_corte = generar_variantes_dinamica(seq, self.TOPE_CALCULO)
        mostrar = lista_resultados[:self.TOPE_VISUAL]
        for adn, prot in mostrar:
            self.tree.insert("", "end", values=(adn, prot))
        total_calc = len(lista_resultados)
        msg_adn = f"Calculados: {total_calc:,}"
        if hubo_corte:
            msg_adn += " (TRUNCADO POR SEGURIDAD)"
            self.lbl_adn.config(fg="#D32F2F")
        else:
            self.lbl_adn.config(fg="black")

        if total_calc > self.TOPE_VISUAL:
            msg_adn += f" | Vista: Primeros {self.TOPE_VISUAL}"

        self.lbl_adn.config(text=msg_adn)
        proteinas_unicas = len(set([p for (_, p) in lista_resultados]))
        self.lbl_prot.config(text=f"Proteínas Únicas: {proteinas_unicas:,}")

        if hubo_corte:
            messagebox.showwarning(
                "Límite de Seguridad Activado",
                f"La secuencia genera más de {self.TOPE_CALCULO:,} variantes.\n"
                "El programa detuvo la generación para proteger la memoria RAM.\n\n"
                "Se muestran los resultados hasta el punto de corte."
            )

root = tk.Tk()
app = DashboardDPList(root)
root.mainloop()
