# === SINCRONIZAR SINTEGIA AUTOMATIZADA ENTRE ESPECIES ===
import os
import subprocess
from pathlib import Path

# === CONFIGURACIÓN ===
carpeta_base = Path("/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas")
carpeta_salida = carpeta_base / "analisis_sintenia"
carpeta_salida.mkdir(exist_ok=True)

# === ESPECIES DISPONIBLES Y ARCHIVOS RELACIONADOS ===
especies = {
    "salix_purpurea": "Spurpurea_519_v5.1.gene.gff3",
    "oryza_sativa": "genomic.gff",
    "populus_nigra": "genomic.gff",
    "eucalyptus_globulus": "genomic.gff",
    "arabidopsis_thaliana": "GCF_000001735.4/genomic.gff"  # usar GCF por ser referencia
}

# === FUNCIONES ===
def convertir_a_bed(nombre_especie, ruta_gff):
    """Convierte un archivo GFF a formato BED compatible con MCScan (Python)."""
    archivo_bed = carpeta_salida / f"{nombre_especie}.bed"
    comando = f"python -m jcvi.formats.gff bed --type=gene {ruta_gff} -o {archivo_bed}"
    subprocess.run(comando, shell=True, check=True)
    return archivo_bed

def ejecutar_blastp(nombre_especie):
    """Ejecuta BLASTP entre proteínas de la misma especie."""
    fasta = carpeta_base / nombre_especie / "protein.faa"
    db = carpeta_salida / f"{nombre_especie}.blastdb"
    salida_blast = carpeta_salida / f"{nombre_especie}.blast"
    
    subprocess.run(f"makeblastdb -in {fasta} -dbtype prot -out {db}", shell=True, check=True)
    subprocess.run(
        f"blastp -query {fasta} -db {db} -evalue 1e-5 -outfmt 6 -num_threads 4 -out {salida_blast}",
        shell=True, check=True
    )
    return salida_blast

def correr_mcscan(nombre_especie):
    """Ejecuta jcvi.compara.synteny para duplicaciones segmentales."""
    archivo_bed = carpeta_salida / f"{nombre_especie}.bed"
    archivo_blast = carpeta_salida / f"{nombre_especie}.blast"
    
    comando = f"python -m jcvi.compara.synteny {archivo_bed} {archivo_bed} --blastfile={archivo_blast} -o {nombre_especie}.anchors"
    subprocess.run(comando, shell=True, check=True)

    comando_blocks = f"python -m jcvi.compara.synteny screen --minspan=30 --simple {nombre_especie}.anchors"
    subprocess.run(comando, shell=True, check=True)

    # Graficar dotplot
    dotplot = carpeta_salida / f"{nombre_especie}_dotplot.pdf"
    comando_plot = f"python -m jcvi.graphics.dotplot {nombre_especie}.anchors.filtered --title={nombre_especie} -o {dotplot}"
    subprocess.run(comando_plot, shell=True, check=True)

# === LOOP PRINCIPAL ===
for especie, archivo_gff in especies.items():
    print(f"\n🧬 Procesando especie: {especie}")
    ruta_gff = carpeta_base / especie / archivo_gff
    
    try:
        convertir_a_bed(especie, ruta_gff)
        ejecutar_blastp(especie)
        correr_mcscan(especie)
        print(f"✅ Finalizado para {especie}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error procesando {especie}: {e}")

print("\n🎉 Análisis de duplicaciones segmentales completado")
