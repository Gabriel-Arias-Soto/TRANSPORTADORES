import os
from Bio import Phylo
import pandas as pd
import re

# ---------------- CONFIGURACIÓN ----------------
ruta_arboles = "/mnt/c/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/arboles"
ruta_salida = ruta_arboles

familias = ["HMA", "PIP", "TIP", "NIP", "NRAMP", "PCS", "MT-L"]

# ---------------- COMBINAR ÁRBOLES ----------------
arboles = []
for familia in familias:
    patron = re.compile(f"{familia.lower()}\\.(treefile|nwk)$", re.IGNORECASE)
    archivos = [f for f in os.listdir(ruta_arboles) if patron.search(f)]
    
    if not archivos:
        print(f"⚠ No se encontró árbol para {familia}")
        continue
    
    archivo = archivos[0]
    path = os.path.join(ruta_arboles, archivo)
    
    tree = Phylo.read(path, "newick")
    
    # Asegurar que todos los terminales tengan un nombre
    for clade in tree.get_terminals():
        if clade.name is None:
            raise ValueError(f"Clado sin nombre en árbol {familia}")
    
    arboles.append(tree)

# Validación básica
if not arboles:
    raise ValueError("No se encontraron árboles para combinar")

# Unir los árboles
from Bio.Phylo.BaseTree import Clade, Tree

def unir_arboles(lista_arboles):
    """Une múltiples árboles como subclados de una raíz común."""
    raiz = Clade(branch_length=0.1, name="ROOT")
    raiz.clades = [t.root for t in lista_arboles]
    return Tree(root=raiz)

arbol_combinado = unir_arboles(arboles)

# Guardar árbol combinado
output_tree = os.path.join(ruta_salida, "arbol_combinado.nwk")
Phylo.write(arbol_combinado, output_tree, "newick")
print(f"✅ Árbol combinado guardado en: {output_tree}")

# ---------------- EXTRAER METADATOS ----------------
datos = []
for clade in arbol_combinado.get_terminals():
    if clade.name is None:
        print("⚠ Terminal sin nombre, se omite")
        continue

    new_id = clade.name
    familia = new_id.split("-")[0] if "-" in new_id else "Desconocida"
    
    especie_cod = re.search(r"-(\w{2,3})-", new_id)
    codigo_especie = especie_cod.group(1) if especie_cod else "NA"
    
    especie = {
        "ATH": "A. thaliana",
        "EGL": "E. globulus",
        "OSA": "O. sativa",
        "PNI": "P. nigra",
        "SPU": "S. purpurea"
    }.get(codigo_especie, "Desconocida")
    
    datos.append({
        "new_id": new_id,
        "familia": familia,
        "codigo_especie": codigo_especie,
        "especie": especie
    })

if not datos:
    raise ValueError("No se pudieron extraer metadatos, revisar los nombres de los terminales")

df_meta = pd.DataFrame(datos)
output_meta = os.path.join(ruta_salida, "arbol_combinado_metadatos.tsv")
df_meta.to_csv(output_meta, sep="\t", index=False)
print(f"✅ Metadatos guardados en: {output_meta}")
