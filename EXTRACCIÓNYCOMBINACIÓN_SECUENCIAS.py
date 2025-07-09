#!/usr/bin/env python3
"""
SCRIPT FINAL OPTIMIZADO PARA EXTRACCIÓN Y COMBINACIÓN DE SECUENCIAS
Versión que no depende de pandas para guardar correspondencias
"""

from Bio import SeqIO
from pathlib import Path
import sys

# ================= CONFIGURACIÓN =================
BASE_DIR = Path("/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas")
SPECIES_DIRS = [
    "arabidopsis_thaliana",
    "eucalyptus_globulus",
    "oryza_sativa",
    "populus_nigra",
    "salix_purpurea"
]
FAMILIAS = ["HMA", "PIP", "TIP", "NIP", "NRAMP", "PCS", "MT-L"]
RESULTS_DIR = BASE_DIR / "analisis_filogenia"

# Mapeo de códigos de especie
SPECIES_CODES = {
    "arabidopsis_thaliana": "ATH",
    "eucalyptus_globulus": "EGL",
    "oryza_sativa": "OSA",
    "populus_nigra": "PNI",
    "salix_purpurea": "SPU"
}

# Parámetros de filtrado
EVALUE_CUTOFF = 1e-10
COVERAGE_CUTOFF = 50  # %

# ================= FUNCIONES OPTIMIZADAS =================
def crear_directorios():
    """Crea los directorios necesarios para los resultados"""
    (RESULTS_DIR / "fastas_combinados").mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "correspondencias").mkdir(exist_ok=True)

def leer_blast(archivo):
    """Lee archivo BLAST de manera robusta"""
    datos = []
    with open(archivo, 'r') as f:
        for linea in f:
            if linea.startswith('#'):
                continue
            partes = linea.strip().split('\t')
            if len(partes) >= 12:
                datos.append(partes[:12])
    return datos

def guardar_correspondencias(corr_output, correspondencias):
    """Guarda las correspondencias sin usar pandas"""
    try:
        with open(corr_output, 'w') as f:
            # Encabezados
            f.write("new_id\toriginal_id\tspecies\tfamily\tevalue\tcoverage\tpident\n")
            # Datos
            for corr in correspondencias:
                linea = (
                    f"{corr['new_id']}\t{corr['original_id']}\t{corr['species']}\t"
                    f"{corr['family']}\t{corr['evalue']}\t{corr['coverage']}\t"
                    f"{corr['pident']}\n"
                )
                f.write(linea)
        return True
    except Exception as e:
        print(f"  ❌ Error crítico guardando {corr_output}: {str(e)}")
        return False

def procesar_familia(familia):
    """Procesa una familia completa de transportadores"""
    print(f"\n{'='*60}\nPROCESANDO FAMILIA: {familia}\n{'='*60}")
    
    # Paso 1: Combinar resultados BLAST
    todos_datos = []
    for especie in SPECIES_DIRS:
        blast_file = BASE_DIR / especie / "blast_results" / f"blast_referencia_{familia}.tsv"
        
        if not blast_file.exists():
            print(f"  ⚠️ Archivo no encontrado: {blast_file}")
            continue
            
        try:
            datos = leer_blast(blast_file)
            if not datos:
                print(f"  ⚠️ Archivo vacío: {blast_file.name}")
                continue
                
            for hit in datos:
                try:
                    qseqid, sseqid = hit[0], hit[1]
                    original_id, especie_codigo = qseqid.rsplit('_', 1)
                    
                    # Conversión y cálculos
                    pident = float(hit[2])
                    length = int(hit[3])
                    evalue = float(hit[10])
                    qstart, qend = int(hit[6]), int(hit[7])
                    coverage = (length / (qend - qstart + 1)) * 100
                    
                    # Filtrado
                    if evalue < EVALUE_CUTOFF and coverage > COVERAGE_CUTOFF:
                        todos_datos.append({
                            'original_id': original_id,
                            'especie': especie_codigo,
                            'sseqid': sseqid,
                            'pident': pident,
                            'evalue': evalue,
                            'coverage': coverage
                        })
                except (ValueError, ZeroDivisionError, IndexError) as e:
                    continue
            
            print(f"  ✅ {especie}: {len(datos)} hits leídos")
            
        except Exception as e:
            print(f"  ❌ Error procesando {blast_file.name}: {str(e)}")
    
    if not todos_datos:
        print("  ⚠️ No hay datos BLAST válidos para esta familia")
        return False
    
    # Paso 2: Extraer secuencias
    secuencias = []
    correspondencias = []
    
    for especie in SPECIES_DIRS:
        especie_codigo = SPECIES_CODES[especie]
        protein_file = BASE_DIR / especie / "protein.faa"
        
        if not protein_file.exists():
            print(f"  ⚠️ Archivo FASTA no encontrado: {protein_file}")
            continue
            
        # Obtener IDs para esta especie
        ids_especie = {d['original_id'] for d in todos_datos if d['especie'] == especie_codigo}
        
        if not ids_especie:
            continue
            
        # Extraer secuencias
        try:
            with open(protein_file) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id in ids_especie:
                        # Crear nuevo ID
                        new_id = f"{familia}-{especie_codigo}-{record.id}"
                        record.id = new_id
                        record.description = f"family={familia}|species={especie_codigo}"
                        secuencias.append(record)
                        
                        # Encontrar el mejor hit para esta secuencia
                        best_hit = None
                        for hit in [d for d in todos_datos if d['original_id'] == record.id and d['especie'] == especie_codigo]:
                            if best_hit is None or hit['evalue'] < best_hit['evalue']:
                                best_hit = hit
                        
                        if best_hit:
                            correspondencias.append({
                                'new_id': new_id,
                                'original_id': record.id,
                                'species': especie_codigo,
                                'family': familia,
                                'evalue': best_hit['evalue'],
                                'coverage': best_hit['coverage'],
                                'pident': best_hit['pident']
                            })
            
            print(f"  ➡️ {especie_codigo}: {len([s for s in secuencias if s.id.split('-')[1] == especie_codigo])} secuencias extraídas")
            
        except Exception as e:
            print(f"  ❌ Error extrayendo secuencias de {protein_file.name}: {str(e)}")
    
    # Guardar resultados
    if secuencias:
        # Guardar FASTA combinado
        fasta_output = RESULTS_DIR / "fastas_combinados" / f"{familia}_combinado.faa"
        with open(fasta_output, "w") as handle:
            SeqIO.write(secuencias, handle, "fasta")
        
        # Guardar correspondencias
        corr_output = RESULTS_DIR / "correspondencias" / f"{familia}_correspondencia.tsv"
        if guardar_correspondencias(corr_output, correspondencias):
            print(f"\n📊 RESULTADOS PARA {familia}:")
            print(f"  - Secuencias: {len(secuencias)}")
            print(f"  - Archivo FASTA: {fasta_output}")
            print(f"  - Correspondencias: {corr_output}")
            return True
    
    print("  ⚠️ No se generaron resultados válidos")
    return False

# ================= EJECUCIÓN PRINCIPAL =================
def main():
    print("\n" + "="*60)
    print(" EXTRACCIÓN Y COMBINACIÓN DE SECUENCIAS - VERSIÓN OPTIMIZADA ")
    print("="*60)
    
    crear_directorios()
    
    # Procesar cada familia
    exitosos = 0
    for familia in FAMILIAS:
        if procesar_familia(familia):
            exitosos += 1
    
    # Resumen final
    print("\n" + "="*60)
    print(" RESUMEN FINAL ")
    print("="*60)
    print(f"Familias procesadas: {len(FAMILIAS)}")
    print(f"Procesamientos exitosos: {exitosos}")
    print(f"Procesamientos fallidos: {len(FAMILIAS) - exitosos}")
    
    if exitosos > 0:
        print("\n✅ Directorios de resultados:")
        print(f"  - FASTA combinados: {RESULTS_DIR / 'fastas_combinados'}")
        print(f"  - Tablas correspondencia: {RESULTS_DIR / 'correspondencias'}")
    else:
        print("\n❌ No se generaron resultados válidos")

if __name__ == "__main__":
    main()