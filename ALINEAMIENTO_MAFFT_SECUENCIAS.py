#!/usr/bin/env python3
"""
SCRIPT PARA ALINEAMIENTO MÚLTIPLE CON MAFFT
Realiza alineamientos de todas las familias de transportadores
"""

import subprocess
from pathlib import Path
import sys

# ================= CONFIGURACIÓN =================
BASE_DIR = Path("/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas")
INPUT_DIR = BASE_DIR / "analisis_filogenia" / "fastas_combinados"
OUTPUT_DIR = BASE_DIR / "analisis_filogenia" / "alineamientos"
FAMILIAS = ["HMA", "PIP", "TIP", "NIP", "NRAMP", "PCS", "MT-L"]

# Configuración MAFFT
MAFFT_CMD = "mafft --auto {input} > {output}"
THREADS = 10  # Ajusta según tu CPU

# ================= FUNCIONES =================
def crear_directorios():
    """Crea los directorios necesarios"""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def ejecutar_mafft(input_file, output_file):
    """Ejecuta MAFFT y maneja errores"""
    cmd = MAFFT_CMD.format(input=input_file, output=output_file)
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            check=True,
            executable='/bin/bash'
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error al alinear {input_file.name}: {str(e)}")
        return False

def procesar_familia(familia):
    """Procesa el alineamiento para una familia"""
    input_file = INPUT_DIR / f"{familia}_combinado.faa"
    output_file = OUTPUT_DIR / f"{familia}_alineado.faa"
    
    if not input_file.exists():
        print(f"⚠️ Archivo de entrada no encontrado: {input_file}")
        return False
    
    print(f"\n🔍 Alineando familia {familia}...")
    print(f"   Input: {input_file}")
    print(f"   Output: {output_file}")
    
    if ejecutar_mafft(input_file, output_file):
        # Verificar que el archivo de salida se creó correctamente
        if output_file.exists() and output_file.stat().st_size > 0:
            print(f"✅ Alineamiento completado: {output_file.name}")
            return True
        else:
            print(f"❌ El archivo de salida no se generó correctamente")
            return False
    return False

# ================= EJECUCIÓN PRINCIPAL =================
def main():
    print("\n" + "="*60)
    print(" PIPELINE DE ALINEAMIENTO CON MAFFT ")
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
    print(f"Alineamientos exitosos: {exitosos}")
    print(f"Alineamientos fallidos: {len(FAMILIAS) - exitosos}")
    
    if exitosos > 0:
        print(f"\n✅ Alineamientos guardados en: {OUTPUT_DIR}")
        print("\n🛠️ Pasos siguientes:")
        print("  1. Verificar alineamientos:")
        print("     less HMA_alineado.faa")
        print("  2. Construir árboles con IQ-TREE:")
        print("     iqtree -s HMA_alineado.faa -m MFP -bb 1000 -nt AUTO")
    else:
        print("\n❌ No se generaron alineamientos válidos")

if __name__ == "__main__":
    main()