#!/usr/bin/env python3
"""
SCRIPT BLAST CON ETIQUETAS DE ESPECIE INTEGRADAS
- Ejecuta BLASTp contra múltiples bases de datos
- Añade automáticamente el código de especie al ID de secuencia
- Genera archivos TSV listos para análisis filogenético
"""

import subprocess
import sys
from pathlib import Path
import time
import os

def run_command(cmd):
    """Ejecuta un comando mostrando tiempo y detectando errores."""
    start_time = time.time()
    try:
        result = subprocess.run(
            cmd, 
            shell=True, 
            check=True, 
            capture_output=True, 
            text=True,
            executable='/bin/bash'
        )
        elapsed = time.time() - start_time
        print(f"✅ [{(elapsed):.1f}s] {cmd[:80]}...")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error en: {cmd}\n{e.stderr}")
        return False

def verify_db(db_path):
    """Verifica que la base de datos BLAST esté completa."""
    required = ['.phr', '.pin', '.psq']
    return all(Path(f"{db_path}{ext}").exists() for ext in required)

# ===== CONFIGURACIÓN =====
BASE_DIR = Path("/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas")
REF_DIR = BASE_DIR / "referencias"
SPECIES = [
    "arabidopsis_thaliana", "eucalyptus_globulus", "oryza_sativa",
    "populus_nigra", "salix_purpurea"
]
SPECIES_CODES = {
    "arabidopsis_thaliana": "ATH",
    "eucalyptus_globulus": "EGL",
    "oryza_sativa": "OSA",
    "populus_nigra": "PNI",
    "salix_purpurea": "SPU"
}
NUM_THREADS = max(6, os.cpu_count() - 2)  # Aprovecha los núcleos disponibles

# ===== PASO 1: CREAR BASES DE DATOS =====
print("\n" + "="*40 + "\nCREANDO BASES DE DATOS\n" + "="*40)
for ref_file in REF_DIR.glob("referencia_*.faa"):
    db_name = ref_file.stem
    db_path = REF_DIR / db_name
    
    if verify_db(db_path):
        print(f"⚡ DB {db_name} ya existe")
        continue
        
    cmd = f"makeblastdb -in {ref_file} -dbtype prot -out {db_path} -parse_seqids"
    if not run_command(cmd) or not verify_db(db_path):
        sys.exit(f"Error crítico con {db_name}")

# ===== PASO 2: EJECUTAR BLASTP CON ETIQUETAS =====
print("\n" + "="*40 + "\nEJECUTANDO BLASTP CON ETIQUETAS\n" + "="*40)
for species in SPECIES:
    species_dir = BASE_DIR / species
    protein_file = species_dir / "protein.faa"
    
    if not protein_file.exists():
        print(f"⚠️ {species}: Archivo no encontrado")
        continue
        
    blast_dir = species_dir / "blast_results"
    blast_dir.mkdir(exist_ok=True)
    
    species_code = SPECIES_CODES[species]
    print(f"\n🌱 Procesando {species} ({species_code}) ({NUM_THREADS} hilos)...")
    
    for ref_db in REF_DIR.glob("referencia_*.phr"):
        ref_name = ref_db.stem.split('.')[0]
        output_file = blast_dir / f"blast_{ref_name}.tsv"
        temp_file = blast_dir / f"temp_{ref_name}.tsv"
        
        if output_file.exists():
            print(f"⏩ {ref_name} ya procesado")
            continue
            
        # Comando BLAST con formato extendido SIN headers
        cmd = (
            f"blastp -query {protein_file} -db {REF_DIR/ref_name} "
            f"-out {temp_file} -outfmt '6 qseqid sseqid pident length "
            f"mismatch gapopen qstart qend sstart send evalue bitscore "
            f"qlen slen qcovs scovs' -evalue 1e-10 -qcov_hsp_perc 50 "
            f"-num_threads {NUM_THREADS}"
        )
        
        if run_command(cmd):
            # Procesar archivo temporal para añadir etiquetas
            with open(temp_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line in infile:
                    parts = line.strip().split('\t')
                    if len(parts) >= 12:
                        # Añadir código de especie al qseqid
                        parts[0] = f"{parts[0]}_{species_code}"
                        # Escribir línea modificada
                        outfile.write('\t'.join(parts[:12]) + '\n')
            
            # Eliminar archivo temporal
            temp_file.unlink()
            
            # Contar hits y tamaño del archivo
            with open(output_file, 'r') as f:
                line_count = sum(1 for _ in f)
            size_mb = output_file.stat().st_size / (1024 * 1024)
            
            print(f"📊 {ref_name}: {size_mb:.1f} MB | {line_count} hits | Etiquetas: {species_code}")

# ===== RESUMEN =====
print("\n" + "="*40 + "\nRESUMEN FINAL\n" + "="*40)
subprocess.run(
    f"find {BASE_DIR} -name 'blast_*.tsv' -size +0k -exec ls -lh {{}} +", 
    shell=True
)
print("\n✅ ¡Proceso completado! Archivos BLAST listos para filogenia:")
print(f"- IDs de secuencia con formato: ID_ESPECIE (ej: NP_001321606.1_ATH)")
print(f"- Formato TSV estándar (12 columnas) compatible con herramientas filogenéticas")