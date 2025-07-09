import pandas as pd
import os

# Ruta donde están los .tblout
ruta_tblouts = r"/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/hmmscan_resultados"

# Columnas esperadas según salida hmmscan --tblout
columnas = [
    "target_name", "accession", "query_name", "accession2",
    "E-value_full", "score_full", "bias_full",
    "E-value_best", "score_best", "bias_best",
    "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc",
    "description"
]

# Procesar todos los .tblout en la carpeta
for archivo in os.listdir(ruta_tblouts):
    if archivo.endswith(".tblout"):
        ruta_archivo = os.path.join(ruta_tblouts, archivo)
        ruta_csv = ruta_archivo.replace(".tblout", "_completo.csv")
        
        datos = []
        with open(ruta_archivo) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split(maxsplit=18)
                if len(parts) == 19:
                    datos.append(parts)

        if datos:
            df = pd.DataFrame(datos, columns=columnas)
            df.to_csv(ruta_csv, index=False)
            print(f"✅ Procesado: {archivo} --> {os.path.basename(ruta_csv)}")
        else:
            print(f"⚠️ Archivo vacío o sin datos válidos: {archivo}")

print("\n✅ Todos los .tblout han sido procesados.")
