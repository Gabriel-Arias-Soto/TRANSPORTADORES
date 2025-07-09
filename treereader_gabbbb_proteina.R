# ================= LIMPIEZA Y PAQUETES =================
rm(list = ls())
cat("\014")

library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# ================= CONFIGURACIÓN DE RUTAS =================
ruta_base <- "C:/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/arboles"
archivo_arbol <- file.path(ruta_base, "arbol_consolidado.nwk")
archivo_metadatos <- file.path(ruta_base, "arbol_consolidado_metadatos.tsv")
ruta_salida <- file.path(ruta_base, "arbol_consolidado_visual_V02.png")

# ================= CARGA DE DATOS =================
arbol <- read.tree(archivo_arbol)
metadatos <- read_tsv(archivo_metadatos, show_col_types = FALSE)

# ================= CONFIGURACIÓN VISUAL =================
colores_familias <- c(
  "HMA" = "#1F77B4", "PIP" = "#FF7F0E", "TIP" = "#2CA02C",
  "NIP" = "#D62728", "NRAMP" = "#9467BD", "PCS" = "#8C564B",
  "MT-L" = "#E377C2", "Desconocida" = "grey70"
)

formas_especies <- c(
  "A. thaliana" = 15, "E. globulus" = 16, "O. sativa" = 17,
  "P. nigra" = 18, "S. purpurea" = 8, "Desconocida" = 4
)

# ================= VISUALIZACIÓN =================
p <- ggtree(arbol, layout = "fan", open.angle = 5, size = 0.05) %<+% metadatos +
  geom_tippoint(
    aes(color = familia, shape = especie),
    size = 0.3, alpha = 0.4  # Puntos muy pequeños y semitransparentes
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    shape = guide_legend(override.aes = list(size = 3))
  ) +
  scale_color_manual(values = colores_familias, na.value = "grey70") +
  scale_shape_manual(values = formas_especies, na.value = 4) +
  labs(
    title = "Árbol Filogenético Consolidado de Transportadores",
    subtitle = paste("Total de secuencias:", length(arbol$tip.label)),
    color = "Familia",
    shape = "Especie"
  ) +
  theme_tree2() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 25, face = "bold", hjust= 0.5),
    plot.subtitle = element_text(hjust=0.4)
  )

# ================= GUARDAR Y MOSTRAR =================
ggsave(
  filename = ruta_salida,
  plot = p,
  width = 20, 
  height = 16, 
  dpi = 1200,
  bg = "white"
)

cat("\n✅ Visualización guardada en:", ruta_salida, "\n")
print(p)

