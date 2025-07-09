#!/usr/bin/env python3

"""
Ejecuta hmmscan sobre los FASTA combinados por familia de transportadores.
Genera un .tblout por cada familia.

Requisitos:
- pfam-A.hmm previamente construido
- hmmscan en el PATH
"""

import subprocess
from pathlib import Path

# =======================
# Configuración de rutas
# =======================
BASE_DIR = Path("/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia")
FASTAS_COMBINADOS = BASE_DIR / "fastas_combinados"
PFAM_DB = Path("/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/ELIPS-MICROBIOTA/pfam/Pfam-A.hmm")
SALIDA_DIR = BASE_DIR / "hmmscan_resultados"
SALIDA_DIR.mkdir(exist_ok=True)

FAMILIAS = ["TIP", "NIP", "NRAMP", "PCS", "MT-L"] #retiré HMA y PIP porque ya estaban procesados

# =======================
# Loop por familia
# =======================
for familia in FAMILIAS:
    fasta_file = FASTAS_COMBINADOS / f"{familia}_combinado.faa"
    
    if not fasta_file.exists():
        print(f"⚠️  No se encontró {fasta_file.name}, saltando.")
        continue

    salida_tblout = SALIDA_DIR / f"{familia}_pfam.tblout"

    cmd = [
        "hmmscan",
        "--tblout", str(salida_tblout),
        str(PFAM_DB),
        str(fasta_file)
    ]

    print(f"🔎 Ejecutando hmmscan para {familia}...")
    subprocess.run(cmd, check=True)
    print(f"✅ Resultados guardados en: {salida_tblout}")

print("\n🎉 hmmscan finalizado para todas las familias.")
