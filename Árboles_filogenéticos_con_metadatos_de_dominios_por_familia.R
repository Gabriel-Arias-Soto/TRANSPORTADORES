# ================= LIMPIEZA Y PAQUETES =================
rm(list = ls())
cat("\014")

library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(RColorBrewer)
library(ape)

# ================= CONFIGURACIÓN DE RUTAS =================
ruta_arboles <- "C:/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/arboles"
ruta_metadatos <- "C:/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/hmmscan_resultados"

# Lista de familias a procesar
familias <- c("HMA", "PIP", "NIP", "TIP", "PCS", "NRAMP", "MT-L")

# ================= LOOP POR FAMILIA =================
for (familia in familias) {
  
  cat("\n🌲 Procesando familia:", familia, "\n")
  
  # ---------------- Carga de árbol ----------------
  archivo_arbol <- file.path(ruta_arboles, paste0(familia, ".treefile"))
  if (!file.exists(archivo_arbol)) {
    cat("⚠️ Árbol no encontrado para", familia, "\n")
    next
  }
  arbol <- read.tree(archivo_arbol)
  
  # ---------------- Carga de metadatos ----------------
  archivo_matriz <- file.path(ruta_metadatos, paste0("graficos_", familia, "/matriz_binaria.csv"))
  if (!file.exists(archivo_matriz)) {
    cat("⚠️ Matriz binaria no encontrada para", familia, "\n")
    next
  }
  
  mat <- read_csv(archivo_matriz, show_col_types = FALSE)
  
  if (!"query_name" %in% colnames(mat)) {
    cat("⚠️ No se encontró columna 'query_name' en matriz para", familia, "\n")
    next
  }
  
  metadatos <- mat %>%
    mutate(
      num_dominios = rowSums(select(., -query_name)),
      Especie = str_extract(query_name, "(?<=-)[^-]+(?=-)")
    ) %>%
    select(query_name, num_dominios, Especie)
  
  # ---------------- Paleta de colores y formas ----------------
  max_dom <- max(metadatos$num_dominios, na.rm = TRUE)
  paleta <- colorRampPalette(brewer.pal(min(9, max(3, max_dom)), "YlOrRd"))
  colores_dominios <- paleta(max_dom)
  
  especies <- unique(metadatos$Especie)
  formas <- seq(21, 21 + length(especies) - 1)
  names(formas) <- especies
  
  # ---------------- Visualización ----------------
  p <- ggtree(arbol) %<+% metadatos +
    geom_tippoint(
      aes(color = as.factor(num_dominios), shape = Especie),
      size = 3
    ) +
    scale_color_manual(values = colores_dominios, name = "N° Dominios") +
    scale_shape_manual(values = formas, name = "Especie") +
    labs(
      title = paste("Árbol Filogenético -", familia),
      subtitle = paste("Total de secuencias:", length(arbol$tip.label)),
      color = "N° de Dominios",
      shape = "Especie"
    ) +
    theme_tree2() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # ---------------- Guardado de resultados ----------------
  carpeta_salida <- file.path(ruta_metadatos, paste0("graficos_", familia))
  
  archivo_pdf <- file.path(carpeta_salida, paste0("arbol_dominios_", familia, ".pdf"))
  archivo_png <- file.path(carpeta_salida, paste0("arbol_dominios_", familia, ".png"))
  archivo_metadatos <- file.path(carpeta_salida, paste0("metadatos_dominios_", familia, ".tsv"))
  
  write_tsv(metadatos, archivo_metadatos)
  ggsave(archivo_pdf, p, width = 10, height = 8)
  ggsave(archivo_png, p, width = 10, height = 8, dpi = 300)
  
  cat("✅ Guardados: PDF, PNG y metadatos para", familia, "\n")
}

cat("\n🎯 Procesamiento finalizado para", length(familias), " familias.\n")
