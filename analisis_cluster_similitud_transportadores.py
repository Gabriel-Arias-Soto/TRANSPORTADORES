#!/usr/bin/env python3

"""
Análisis de clústers por similitud de secuencia para transportadores
Genera mapas de calor y dendrogramas basados en distancias de alineamiento por familia.

Mejoras:
- Limpieza de IDs largos (deja solo el primer campo antes de '|').
- Títulos informativos y correctamente ubicados.
- Ajuste estético en ejes y leyenda.

Requisitos:
- Biopython
- pandas
- seaborn
- matplotlib
- tqdm
"""

import os
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")  # Evita error Qt en WSL
import matplotlib.pyplot as plt

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator

# -----------------------------
# Configuración de rutas
# -----------------------------
BASE_DIR = Path("/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas")
ALIGNMENTS_DIR = BASE_DIR / "analisis_filogenia" / "alineamientos"
SALIDA_DIR = BASE_DIR / "analisis_filogenia" / "clustering"

FAMILIAS = ["HMA", "PIP", "TIP", "NIP", "NRAMP", "PCS", "MT-L"]

SALIDA_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Función para procesar una familia
# -----------------------------
def procesar_familia(familia):
    aln_file = ALIGNMENTS_DIR / f"{familia}_alineado.faa"
    
    if not aln_file.exists():
        print(f"⚠️  No se encontró el alineamiento para {familia}")
        return

    print(f"🔬 Procesando {familia}...")

    alignment = AlignIO.read(aln_file, "fasta")

    # Limpieza opcional de caracteres no estándar
    valid_aminoacids = set("ACDEFGHIKLMNPQRSTVWY-")
    for record in alignment:
        cleaned_seq = ''.join(aa if aa in valid_aminoacids else '-' for aa in str(record.seq))
        record.seq = Seq(cleaned_seq)

    # IDs simplificados
    names = [record.id.split("|")[0] for record in alignment]

    # Calcular matriz de distancias
    calculator = DistanceCalculator("blosum62")
    matrix = pd.DataFrame(index=names, columns=names, dtype=float)

    print(f"🔄 Calculando distancias para {familia}...")

    for i, record1 in enumerate(tqdm(alignment, desc="Comparando")):
        for j, record2 in enumerate(alignment):
            if pd.isna(matrix.iat[i, j]):
                dist = calculator._pairwise(record1, record2)
                matrix.iat[i, j] = dist
                matrix.iat[j, i] = dist

    # Visualizar clustermap
    g = sns.clustermap(matrix, cmap="viridis", figsize=(14, 12))
    
    # Limpiar ejes
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=6)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)

    # Leyenda (colorbar) ajustada
    cbar = g.ax_heatmap.collections[0].colorbar
    cbar.ax.set_position([0.97, 0.3, 0.02, 0.4])

    # Título centrado arriba
    total_seq = len(names)
    g.fig.suptitle(f"Mapa de calor de similitud - Familia {familia}\nTotal de secuencias: {total_seq}", fontsize=14, y=1.05)

    output_path = SALIDA_DIR / f"{familia}_cluster_heatmap.pdf"
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close()

    print(f"✅ Clustermap guardado en: {output_path}")

# -----------------------------
# Loop sobre familias
# -----------------------------
for familia in FAMILIAS:
    try:
        procesar_familia(familia)
    except Exception as e:
        print(f"⚠️  Error en {familia}: {e}")

print("\n🎉 Análisis de clúster finalizado.")
