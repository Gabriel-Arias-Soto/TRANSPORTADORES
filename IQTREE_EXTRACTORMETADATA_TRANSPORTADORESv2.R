rm(list = ls()) #%
cat("\014") #%

# Instalar paquetes si no están instalados
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio"))
install.packages(c("ape", "phytools", "ggplot2", "dplyr", "stringr", "tidyr", "purrr"))

# Cargar librerías
library(ape)
library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

# Definir directorios (WSL)
tree_dir <- "C:/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/arboles/"
output_dir <- "C:/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/metadatos_arboles/"

# Crear directorio de salida si no existe
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Definir colores y formas
colores_familias <- c(
  "HMA" = "#1F77B4", "PIP" = "#FF7F0E", "TIP" = "#2CA02C",
  "NIP" = "#D62728", "NRAMP" = "#9467BD", "PCS" = "#8C564B",
  "MT-L" = "#E377C2", "Desconocida" = "grey70"
)

aformas_especies <- c(
  "A. thaliana" = 15, "E. globulus" = 16, "O. sativa" = 17,
  "P. nigra" = 18, "S. purpurea" = 8, "Desconocida" = 4
)

# Función mejorada para procesar nombres de hojas
parse_leaf_labels <- function(labels) {
  data.frame(
    nodo = labels,
    stringsAsFactors = FALSE
  ) %>%
    separate(nodo, into = c("familia", "especie_code", "id"), 
             sep = "-", remove = FALSE, extra = "merge") %>%
    mutate(
      familia = str_extract(familia, "^[A-Za-z]+"), # Extraer solo letras
      especie = case_when(
        str_detect(especie_code, "ATH") ~ "A. thaliana",
        str_detect(especie_code, "EGL") ~ "E. globulus",
        str_detect(especie_code, "OSA") ~ "O. sativa",
        str_detect(especie_code, "PNI") ~ "P. nigra",
        str_detect(especie_code, "SPU") ~ "S. purpurea",
        TRUE ~ "Desconocida"
      )
    ) %>%
    left_join(
      data.frame(
        familia = names(colores_familias),
        color = colores_familias,
        stringsAsFactors = FALSE
      ),
      by = "familia"
    ) %>%
    left_join(
      data.frame(
        especie = names(formas_especies),
        shape = formas_especies,
        stringsAsFactors = FALSE
      ),
      by = "especie"
    ) %>%
    mutate(
      familia = ifelse(familia %in% names(colores_familias), familia, "Desconocida"),
      color = ifelse(is.na(color), colores_familias["Desconocida"], color),
      shape = ifelse(is.na(shape), formas_especies["Desconocida"], shape)
    )
}

# Función para procesar y visualizar cada árbol
process_tree <- function(tree_file) {
  # Extraer nombre de la familia del nombre del archivo
  familia_nombre <- tools::file_path_sans_ext(basename(tree_file))
  
  # Leer árbol (manejando ambos formatos)
  tree <- tryCatch(
    {
      read.iqtree(tree_file) # Para archivos IQ-TREE con soporte
    },
    error = function(e) {
      read.tree(tree_file) # Para árboles básicos
    }
  )
  
  # Crear metadatos
  if (class(tree) == "phylo") {
    metadatos <- parse_leaf_labels(tree$tip.label)
    support <- NULL
  } else {
    metadatos <- parse_leaf_labels(tree@phylo$tip.label)
    support <- as.numeric(tree@node.label)
    support[support < 50] <- NA # Filtrar soportes bajos
  }
  
  # Guardar metadatos
  write.csv(metadatos, 
            file = paste0(output_dir, "metadatos_", familia_nombre, ".csv"), 
            row.names = FALSE)
  
  # Crear visualización circular
  p <- ggtree(tree, layout = "circular") %<+% metadatos +
    geom_tippoint(
      aes(shape = especie, color = familia), 
      size = 3,
      show.legend = TRUE
    ) +
    geom_tiplab(
      aes(color = familia), 
      size = 2.5, 
      offset = 0.01, 
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = colores_familias, 
      name = "Familia Proteica",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = formas_especies, 
      name = "Especie",
      drop = FALSE
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    ggtitle(paste("Árbol Filogenético", familia_nombre)) +
    guides(
      color = guide_legend(override.aes = list(size = 5)),
      shape = guide_legend(override.aes = list(size = 5))
    )
  
  # Añadir soporte bootstrap si existe
  if (!is.null(support) && !all(is.na(support))) {
    p <- p + geom_nodelab(
      aes(label = ifelse(is.na(support), "", paste0(support, "%"))),
      size = 2.5,
      color = "black",
      hjust = 1.3,
      vjust = -0.5
    )
  }
  
  # Guardar gráfico
  ggsave(
    paste0(output_dir, "arbol_circular_", familia_nombre, ".png"),
    plot = p,
    width = 12,
    height = 12,
    dpi = 300
  )
  
  # Guardar gráfico en PDF (mejor calidad para publicaciones)
  ggsave(
    paste0(output_dir, "arbol_circular_", familia_nombre, ".pdf"),
    plot = p,
    width = 12,
    height = 12,
    device = cairo_pdf
  )
  
  return(p)
}

# Procesar todos los árboles en el directorio
arboles <- list.files(tree_dir, pattern = "\\.contree$|\\.treefile$", full.names = TRUE)

# Aplicar la función a todos los árboles (puede tardar)
walk(arboles, safely(process_tree))

# Mensaje de finalización
message("\n¡Procesamiento completado!\n",
        "Se procesaron ", length(arboles), " árboles.\n",
        "Resultados guardados en: ", output_dir)

