# ================= LIMPIEZA Y PAQUETES =================
rm(list = ls())
cat("\014")

library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# ================= CONFIGURACIÓN DE FAMILIA =================
familia <- "MT"
nombre_familia_completo <- "MT-L"  # si quieres que aparezca con el sufijo

# ================= CONFIGURACIÓN DE RUTAS =================
ruta_base <- "C:/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/arboles_por_familia"
archivo_arbol <- file.path(ruta_base, paste0("arbol_", familia, ".nwk"))
archivo_metadatos <- file.path(ruta_base, paste0("metadatos_dominios_", familia, ".tsv"))
ruta_salida <- file.path(ruta_base, paste0("arbol_", familia, "_visual.png"))

# ================= CARGA DE DATOS =================
arbol <- read.tree(archivo_arbol)
metadatos <- read_tsv(archivo_metadatos, show_col_types = FALSE)

# ================= CONFIGURACIÓN VISUAL =================
color_familia <- "#E377C2"
formas_especies <- c(
  "A. thaliana" = 15, "E. globulus" = 16, "O. sativa" = 17,
  "P. nigra" = 18, "S. purpurea" = 8, "Desconocida" = 4
)

# ================= VISUALIZACIÓN =================
p <- ggtree(arbol, layout = "fan", open.angle = 5, size = 0.05) %<+% metadatos +
  geom_tippoint(
    aes(color = familia, shape = especie),
    size = 0.6, alpha = 0.5
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    shape = guide_legend(override.aes = list(size = 3))
  ) +
  scale_color_manual(values = setNames(color_familia, nombre_familia_completo), na.value = "grey70") +
  scale_shape_manual(values = formas_especies, na.value = 4) +
  labs(
    title = paste("Árbol Filogenético - Familia", nombre_familia_completo),
    subtitle = paste("Total de secuencias:", length(arbol$tip.label)),
    color = "Familia",
    shape = "Especie"
  ) +
  theme_tree2() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.4)
  )

# ================= GUARDAR Y MOSTRAR =================
ggsave(
  filename = ruta_salida,
  plot = p,
  width = 16,
  height = 14,
  dpi = 1000,
  bg = "white"
)

cat("\n✅ Árbol de", familia, "guardado en:", ruta_salida, "\n")
print(p)
