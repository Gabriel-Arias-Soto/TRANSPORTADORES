#!/usr/bin/env python3
"""
SCRIPT PARA ALINEAMIENTO MÚLTIPLE CON MAFFT (VERSIÓN PARALELIZADA)
Realiza alineamientos paralelos de todas las familias de transportadores
"""

import subprocess
from pathlib import Path
import sys
from concurrent.futures import ThreadPoolExecutor
import os

# ================= CONFIGURACIÓN =================
BASE_DIR = Path("/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas")
INPUT_DIR = BASE_DIR / "analisis_filogenia" / "fastas_combinados"
OUTPUT_DIR = BASE_DIR / "analisis_filogenia" / "alineamientos"
FAMILIAS = ["HMA", "PIP", "TIP", "NIP", "NRAMP", "PCS", "MT-L"]

# Configuración MAFFT optimizada
MAFFT_CMD = "mafft --auto --thread 1 {input} > {output}"  # Cada alineamiento usa 1 thread
THREADS = max(1, os.cpu_count() - 2)  # Usa todos los cores menos 2

# ================= FUNCIONES OPTIMIZADAS =================
def crear_directorios():
    """Crea los directorios necesarios"""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def ejecutar_mafft(familia):
    """Ejecuta MAFFT para una familia específica"""
    input_file = INPUT_DIR / f"{familia}_combinado.faa"
    output_file = OUTPUT_DIR / f"{familia}_alineado.faa"
    
    if not input_file.exists():
        return (familia, False, f"⚠️ Archivo de entrada no encontrado: {input_file.name}")
    
    cmd = MAFFT_CMD.format(input=input_file, output=output_file)
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            check=True,
            capture_output=True,
            text=True,
            executable='/bin/bash'
        )
        
        # Verificar salida
        if output_file.exists() and output_file.stat().st_size > 0:
            return (familia, True, f"✅ {familia}: Alineamiento completado ({output_file.stat().st_size/1024:.1f} KB)")
        else:
            return (familia, False, f"❌ {familia}: Archivo de salida vacío")
            
    except subprocess.CalledProcessError as e:
        error_msg = f"❌ {familia}: Error en MAFFT (code {e.returncode})"
        if "memory" in e.stderr.lower():
            error_msg += " - Insuficiente memoria"
        return (familia, False, error_msg)
    except Exception as e:
        return (familia, False, f"❌ {familia}: Error inesperado - {str(e)}")

# ================= EJECUCIÓN PARALELA =================
def main():
    print("\n" + "="*60)
    print(f" PIPELINE DE ALINEAMIENTO PARALELO (Usando {THREADS} hilos) ")
    print("="*60)
    
    crear_directorios()
    
    # Verificar archivos de entrada primero
    familias_a_procesar = []
    for familia in FAMILIAS:
        input_file = INPUT_DIR / f"{familia}_combinado.faa"
        if input_file.exists():
            familias_a_procesar.append(familia)
        else:
            print(f"⚠️ Saltando {familia}: archivo de entrada no encontrado")
    
    if not familias_a_procesar:
        print("\n❌ No hay archivos de entrada válidos para procesar")
        sys.exit(1)
    
    # Ejecución paralela
    print(f"\n🔍 Iniciando alineamiento paralelo para {len(familias_a_procesar)} familias...")
    exitosos = 0
    
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        resultados = list(executor.map(ejecutar_mafft, familias_a_procesar))
    
    # Mostrar resultados
    print("\n" + "="*40)
    print(" RESULTADOS POR FAMILIA ")
    print("="*40)
    for familia, status, mensaje in resultados:
        print(mensaje)
        if status:
            exitosos += 1
    
    # Resumen final
    print("\n" + "="*60)
    print(" RESUMEN FINAL ")
    print("="*60)
    print(f"Familias procesadas: {len(familias_a_procesar)}")
    print(f"Alineamientos exitosos: {exitosos}")
    print(f"Alineamientos fallidos: {len(familias_a_procesar) - exitosos}")
    
    if exitosos > 0:
        print(f"\n✅ Alineamientos guardados en: {OUTPUT_DIR}")
        print("\n🛠️ Pasos siguientes:")
        print("  1. Verificar alineamientos:")
        print(f"     less {OUTPUT_DIR}/HMA_alineado.faa")
        print("  2. Construir árboles con IQ-TREE:")
        print(f"     iqtree -s {OUTPUT_DIR}/HMA_alineado.faa -m MFP -bb 1000 -nt AUTO")
    else:
        print("\n❌ No se generaron alineamientos válidos")

if __name__ == "__main__":
    main()