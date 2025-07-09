# ================= LIMPIEZA Y PAQUETES =================
rm(list = ls())
cat("\014")

library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# ================= CONFIGURACIÓN =================
familias <- c("HMA", "MT-L", "NIP", "NRAMP", "PCS", "PIP", "TIP")
ruta_arboles <- "C:/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/arboles"
ruta_resultados <- "C:/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/hmmscan_resultados"

colores_dominios <- scales::viridis_pal(option = "D", direction = -1)(100)
formas_especies <- c(
  "A. thaliana" = 15, "E. globulus" = 16, "O. sativa" = 17,
  "P. nigra" = 18, "S. purpurea" = 8, "Desconocida" = 4
)

# ================= BUCLE POR FAMILIA =================
for (familia in familias) {
  cat("\n🌿 Procesando familia:", familia, "\n")
  
  archivo_arbol <- file.path(ruta_arboles, paste0(familia, ".treefile"))
  ruta_graficos <- file.path(ruta_resultados, paste0("graficos_", familia))
  archivo_metadatos <- file.path(ruta_graficos, paste0("metadatos_dominios_", familia, ".tsv"))
  
  if (!file.exists(archivo_arbol)) {
    warning("⚠️ Árbol no encontrado para ", familia)
    next
  }
  
  if (!file.exists(archivo_metadatos)) {
    warning("⚠️ Metadatos no encontrados para ", familia)
    next
  }
  
  # Leer árbol
  arbol <- read.tree(archivo_arbol)
  
  # Leer metadatos
  metadatos <- read_tsv(archivo_metadatos, show_col_types = FALSE)
  
  # Relacionar directamente con el nombre completo del tip
  tips_df <- tibble(
    tip_label = arbol$tip.label,
    tip_id = arbol$tip.label
  )
  
  # Alinear con metadatos
  metadatos <- metadatos %>% 
    mutate(
      tip_id = query_name,
      especie = case_when(
        str_detect(query_name, "\\bATH\\b") ~ "A. thaliana",
        str_detect(query_name, "\\bEGL\\b") ~ "E. globulus",
        str_detect(query_name, "\\bOSA\\b") ~ "O. sativa",
        str_detect(query_name, "\\bPNI\\b") ~ "P. nigra",
        str_detect(query_name, "\\bSPU\\b") ~ "S. purpurea",
        TRUE ~ "Desconocida"
      )
    )
  
  datos_plot <- tips_df %>% 
    left_join(metadatos, by = "tip_id") %>% 
    mutate(num_dominios = as.numeric(num_dominios))
  
  if (all(is.na(datos_plot$num_dominios))) {
    warning("⚠️ No se pudo mapear dominios para ", familia)
    next
  }
  
  arbol$tip.label <- datos_plot$tip_label
  
  # Visualización
  p <- ggtree(arbol, layout = "rectangular", size = 0.3) %<+% datos_plot +
    geom_tippoint(
      aes(color = num_dominios, shape = especie),
      size = 2, alpha = 0.7
    ) +
    scale_color_gradientn(colors = colores_dominios, na.value = "grey80") +
    scale_shape_manual(values = formas_especies, na.value = 4) +
    theme_tree2() +
    labs(
      title = paste("Árbol Filogenético -", familia),
      subtitle = paste("Secuencias:", length(arbol$tip.label)),
      color = "# Dominios",
      shape = "Especie"
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # Guardar
  ggsave(filename = file.path(ruta_graficos, paste0("arbol_", familia, ".png")), plot = p, width = 12, height = 10, dpi = 1600)
  ggsave(filename = file.path(ruta_graficos, paste0("arbol_", familia, ".pdf")), plot = p, width = 12, height = 10, dpi = 1600)
  
  cat("✅ Guardado para", familia, "\n")
}
