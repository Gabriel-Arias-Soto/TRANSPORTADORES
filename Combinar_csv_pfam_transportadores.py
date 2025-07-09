import pandas as pd
import glob
import os

# Carpeta con los .csv generados
carpeta = "/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/hmmscan_resultados"

# Buscar todos los csv que terminan en _pfam_completo.csv
patron_busqueda = os.path.join(carpeta, "*_pfam_completo.csv")
archivos_csv = glob.glob(patron_busqueda)

print(f"Archivos encontrados:\n{archivos_csv}\n")

df_list = []

for archivo in archivos_csv:
    df = pd.read_csv(archivo)
    
    # Extraer familia desde el nombre del archivo (antes del primer guion bajo)
    nombre_archivo = os.path.basename(archivo)
    familia = nombre_archivo.split("_")[0]
    
    df["Familia"] = familia
    df_list.append(df)

# Verifica si hay archivos para combinar
if df_list:
    df_final = pd.concat(df_list, ignore_index=True)
    salida = os.path.join(carpeta, "pfam_transportadores.csv")
    df_final.to_csv(salida, index=False)
    print(f"\n✅ Archivo combinado guardado como:\n{salida}")
else:
    print("⚠️ No se encontraron archivos para combinar. Revisa la ruta o el patrón de búsqueda.")
