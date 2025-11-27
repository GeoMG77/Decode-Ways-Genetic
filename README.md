# Decode-Ways-Genetic
Proyecto final de análisis de algoritmos 25b
Decode Ways Genetic: Análisis de Variantes de ADN
Proyecto Universitario: Comparativa de Algoritmos (Fuerza Bruta vs Divide y Vencerás vs Voraz vs Dinámica).
Este proyecto explora cómo aplicar conceptos clásicos de Algoritmia para resolver un problema real de Bioinformática: la traducción de secuencias de ADN con ambigüedad (código IUPAC).
El Problema
En biología, las secuenciaciones de ADN no siempre son perfectas. A veces obtenemos una N (que puede ser A, C, G o T) u otras letras ambiguas.
El problema es similar al algoritmo "Decode Ways":
•	Dada una secuencia ATG NNN, ¿cuántas proteínas distintas se pueden formar?
•	¿Cuál es la más probable o la más pesada molecularmente?
•	¿Cómo evitamos que la computadora colapse al calcular millones de combinaciones?
Módulos del Proyecto
El software está dividido en 4 módulos independientes, cada uno utilizando una técnica algorítmica diferente para resolver el mismo problema:
1. Fuerza Bruta (Exhaustivo)
•	Archivo: fuerza_bruta_dashboard.py
•	Descripción: Genera todas las combinaciones posibles mediante iteración pura.
•	Aprendizaje: Demuestra la explosión combinatoria ($O(4^n)$) y cómo la memoria RAM se satura rápidamente con secuencias largas.
2. Divide y Vencerás
•	Archivo: divide_y_venceras_dashboard.py
•	Descripción: Divide la cadena de ADN en mitades recursivamente, resuelve los fragmentos y combina los resultados.
•	Aprendizaje: Muestra el poder de la recursividad y cómo se visualiza la complejidad logarítmica vs exponencial.
3. Algoritmo Voraz (Greedy)
•	Archivo: algoritmo_voraz.py
•	Descripción: En lugar de buscar todas las proteínas, busca solo una: la de mayor peso molecular.
•	Aprendizaje: Demuestra eficiencia extrema ($O(n)$). Toma la decisión local más "pesada" en cada codón sin mirar atrás. Es ideal para cadenas infinitas.
4. Programación Dinámica (Generativa con Tope)
•	Archivo: dinamica_generativa.py
•	Descripción: Construye la lista de variantes iterativamente, reutilizando cálculos previos. Incluye un sistema de seguridad que detiene el proceso si se detectan más de 500,000 variantes.
•	Aprendizaje: El equilibrio perfecto entre generar datos útiles y proteger la memoria del sistema.
Requisitos e Instalación
Este proyecto fue desarrollado en Python 3.
Dependencias
La interfaz gráfica usa tkinter (incluido en Python), pero las gráficas de rendimiento requieren matplotlib.
Clonar el repositorio:
git clone https://github.com/TU_USUARIO/decode-ways-genetic.git
cd decode-ways-genetic
Instalar librerías:
pip install matplotlib

Ejecución
Cada módulo funciona como un programa independiente. Ejecuta el que quieras probar desde la terminal:
# Para probar Fuerza Bruta
python fuerza_bruta_dashboard.py

# Para probar Divide y Vencerás
python divide_y_venceras_dashboard.py

# Para probar el Algoritmo Voraz
python algoritmo_voraz.py

# Para probar Programación Dinámica (Recomendado)
python dinamica_generativa.py

Algoritmo	Velocidad	Memoria	Caso de Uso
Fuerza Bruta	 Lento	Alta	Secuencias cortas y demostración de complejidad.
Divide y Vencerás	 Medio	 Media	Paralelización teórica y recursividad.
Voraz (Greedy)	Instantáneo	 Mínima	Encontrar 1 solución óptima en Big Data.
Dinámica	Rápido	Controlada	Generación masiva segura con límites (Topes).

Autores
Oliver Geovanni Melendrez Gutierrez
Ricardo Gregorio Rivera
•	Materia: Análisis de Algoritmos
•	Universidad: Universidad de Guadalajara


