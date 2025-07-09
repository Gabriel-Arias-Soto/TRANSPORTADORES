import pandas as pd
import os
import re
from pathlib import Path

# ==================== CONFIGURACIÓN ====================
ruta_csv = Path("/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/hmmscan_resultados/pfam_transportadores.csv")
ruta_salida = ruta_csv.parent

# ==================== LECTURA DE DATOS ====================
df = pd.read_csv(ruta_csv)

# Extraer familia y especie desde query_name
df["Familia"] = df["query_name"].str.extract(r"^([^-]+)")
df["Especie"] = df["query_name"].str.extract(r"-(.*?)-")

familias = df["Familia"].unique()

print(f"Procesando metadatos para {len(familias)} familias: {', '.join(familias)}")

# ==================== LOOP POR FAMILIA ====================
for familia in familias:
    print(f"\n➡ Procesando familia: {familia}")

    df_fam = df[df["Familia"] == familia].copy()

    # Obtener todos los dominios Pfam encontrados en la familia
    dominios_unicos = df_fam["target_name"].unique()
    print(f"  Dominios identificados: {len(dominios_unicos)}")

    # Pivot para obtener matriz binaria completa (0/1)
    df_fam["presente"] = 1
    matriz_binaria = df_fam.pivot_table(
        index="query_name",
        columns="target_name",
        values="presente",
        fill_value=0
    ).reset_index()

    # Agregar columna num_dominios y especie
    matriz_binaria["num_dominios"] = matriz_binaria.drop(columns=["query_name"]).sum(axis=1)
    matriz_binaria["Especie"] = matriz_binaria["query_name"].str.extract(r"-(.*?)-")

    # Guardar resultados
    carpeta_salida = ruta_salida / f"graficos_{familia}"
    carpeta_salida.mkdir(exist_ok=True)

    matriz_binaria.to_csv(carpeta_salida / "matriz_binaria.csv", index=False)
    matriz_binaria[["query_name", "num_dominios", "Especie"]].to_csv(
        carpeta_salida / f"metadatos_dominios_{familia}.tsv", sep="\t", index=False
    )

    print(f"  ✅ Guardado: matriz_binaria.csv y metadatos_dominios_{familia}.tsv")

print("\n✨ Generación de metadatos completa.")
