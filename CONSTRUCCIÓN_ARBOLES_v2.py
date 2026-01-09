#!/usr/bin/env python3
"""
Combina árboles filogenéticos generados por IQ-TREE y genera metadatos asociados
Salida lista para visualización o análisis en R
"""

import os
import re
from pathlib import Path
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
import pandas as pd

# ---------------- CONFIGURACIÓN ----------------
BASE_DIR = Path("/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas")
INPUT_DIR = BASE_DIR / "analisis_filogenia" / "arboles"
OUTPUT_DIR = INPUT_DIR
FAMILIAS = ["HMA", "PIP", "TIP", "NIP", "NRAMP", "PCS", "MT-L"]

# ---------------- COMBINACIÓN DE ÁRBOLES ----------------
arboles = []
metadatos = []

for familia in FAMILIAS:
    tree_path = INPUT_DIR / f"{familia}.treefile"
    if not tree_path.exists():
        print(f"⚠️ No se encontró el árbol para {familia}")
        continue

    tree = Phylo.read(tree_path, "newick")

    for clade in tree.get_terminals():
        original_name = clade.name
        clade.name = f"{familia}-{original_name}"

        especie_match = re.search(r"-(\w{2,3})-", original_name)
        codigo_especie = especie_match.group(1) if especie_match else "NA"

        especie = {
            "ATH": "A. thaliana",
            "EGL": "E. globulus",
            "OSA": "O. sativa",
            "PNI": "P. nigra",
            "SPU": "S. purpurea"
        }.get(codigo_especie, "Desconocida")

        metadatos.append({
            "new_id": clade.name,
            "familia": familia,
            "codigo_especie": codigo_especie,
            "especie": especie
        })

    arboles.append(tree)

if not arboles:
    raise ValueError("❌ No se encontraron árboles para combinar.")

# ---------------- CONSTRUCCIÓN DE ÁRBOL COMBINADO ----------------
raiz = Clade(branch_length=0.1, name="ROOT")
raiz.clades = [t.root for t in arboles]
arbol_combinado = Tree(root=raiz)

# ---------------- GUARDADO DE RESULTADOS ----------------
# Guardar árbol combinado
tree_out = OUTPUT_DIR / "arbol_consolidado.nwk"
Phylo.write(arbol_combinado, tree_out, "newick")
print(f"✅ Árbol combinado guardado en: {tree_out}")

# Guardar metadatos
df_meta = pd.DataFrame(metadatos)
meta_out = OUTPUT_DIR / "arbol_consolidado_metadatos.tsv"
df_meta.to_csv(meta_out, sep="\t", index=False)
print(f"✅ Metadatos guardados en: {meta_out}")
