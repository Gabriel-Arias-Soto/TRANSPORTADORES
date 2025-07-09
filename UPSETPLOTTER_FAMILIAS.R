rm(list = ls())
cat("\014")

# ----------------------------
# Cargar librerías
# ----------------------------
required_packages <- c("readr", "dplyr", "tidyr", "UpSetR", "RColorBrewer", "grid", "stringr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, library, character.only = TRUE)

# ----------------------------
# Leer datos
# ----------------------------
file_path <- "C:/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/hmmscan_resultados/pfam_transportadores.csv"
pfam_data <- read_csv(file_path, show_col_types = FALSE) %>%
  mutate(
    Familia = sub("-.*", "", query_name),
    Especie = sub(".*-(.*?)-.*", "\\1", query_name)
  )

# ----------------------------
# Parámetros configurables
# ----------------------------
umbral_score <- 80              
min_proteinas_por_combinacion <- 3 
max_combinaciones <- 10         
n_dominios_top <- 7              # Cantidad de dominios más soportados a mostrar en gráfico

# ----------------------------
# Función para preparar objetos y gráfico de una familia
# ----------------------------
preparar_objetos_upset <- function(familia) {
  pfam_sub <- pfam_data %>% filter(Familia == familia)
  
  if(nrow(pfam_sub) == 0) {
    message("⚠️ No hay datos para la familia: ", familia)
    return(NULL)
  }
  
  # Determinar los N dominios mejor soportados (por score promedio)
  dominios_top <- pfam_sub %>%
    group_by(target_name) %>%
    summarise(score_medio = mean(score_full, na.rm = TRUE)) %>%
    arrange(desc(score_medio)) %>%
    slice_head(n = n_dominios_top) %>%
    pull(target_name)
  
  resultados <- tryCatch({
    pfam_sub %>%
      filter(target_name %in% dominios_top) %>%
      group_by(query_name) %>%
      filter(mean(score_full) >= umbral_score) %>%
      ungroup() %>%
      distinct(query_name, target_name) %>%
      mutate(presente = 1) %>%
      pivot_wider(names_from = target_name, values_from = presente, values_fill = 0)
  }, error = function(e) {
    message("Error creando matriz binaria: ", e$message)
    return(NULL)
  })
  
  if(is.null(resultados) || nrow(resultados) == 0) return(NULL)
  
  mat <- as.data.frame(resultados)
  rownames(mat) <- mat$query_name
  mat$query_name <- NULL
  
  combinaciones <- mat %>%
    mutate(
      combi = apply(., 1, function(x) paste(sort(names(.)[x == 1]), collapse = " + ")),
      n_dominios = rowSums(.)
    ) %>%
    count(combi, n_dominios, sort = TRUE) %>%
    filter(n >= min_proteinas_por_combinacion) %>%
    head(max_combinaciones)
  
  mat_filtrado <- mat %>%
    mutate(combi = apply(., 1, function(x) paste(sort(names(.)[x == 1]), collapse = " + "))) %>%
    filter(combi %in% combinaciones$combi) %>%
    select(-combi)
  
  dominios <- colnames(mat_filtrado)
  
  colores <- tryCatch({
    paleta <- colorRampPalette(brewer.pal(min(12, length(dominios)), "Set3"))
    cols <- paleta(length(dominios))
    names(cols) <- dominios
    cols
  }, error = function(e) {
    message("Error generando colores: ", e$message)
    rep("#4B0082", length(dominios))
  })
  
  # Agregar explícitamente query_name como columna
  mat_filtrado$query_name <- rownames(mat_filtrado)
  
  return(list(
    matriz_filtrada = mat_filtrado,
    combinaciones = combinaciones,
    colores = colores,
    dominios = dominios,
    total_proteinas = length(unique(pfam_sub$query_name)),
    proteinas_filtradas = nrow(mat_filtrado)
  ))
}

# ----------------------------
# Loop por familia con guardado
# ----------------------------
familias <- unique(pfam_data$Familia)
output_dir_base <- dirname(file_path)

cat("\n🔧 Iniciando procesamiento de familias...\n")

for (fam in familias) {
  cat(paste0("🛠️ Procesando familia: ", fam, "...\n"))
  
  obj <- preparar_objetos_upset(fam)
  
  if (is.null(obj)) {
    cat("⚠️  Sin resultados para ", fam, "\n")
    next
  }
  
  fam_dir <- file.path(output_dir_base, paste0("graficos_", fam))
  if(!dir.exists(fam_dir)) dir.create(fam_dir)
  
  # Guardar archivos
  write_csv(obj$matriz_filtrada, file.path(fam_dir, "matriz_binaria.csv"))
  write_csv(obj$combinaciones, file.path(fam_dir, "combinaciones_frecuentes.csv"))
  save(obj, file = file.path(fam_dir, paste0("objeto_", fam, ".RData")))
  
  # Generar y guardar gráfico en PDF
  tryCatch({
    pdf(file.path(fam_dir, paste0("upset_frecuentes_", fam, ".pdf")), width = 12, height = 7)
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 8), "inches"))))
    
    # Título
    pushViewport(viewport(layout.pos.row = 1))
    grid.text(
      paste0("Top ", n_dominios_top, " dominios en ", fam), 
      x = 0.5, y = 0.7,
      gp = gpar(fontsize = 16, fontface = "bold", col = "#333333")
    )
    grid.text(
      paste0("Total proteínas: ", obj$total_proteinas, 
             " | Proteínas mostradas: ", obj$proteinas_filtradas,
             " | Dominios visualizados: ", length(obj$dominios)),
      x = 0.5, y = 0.3,
      gp = gpar(fontsize = 10, col = "#555555")
    )
    popViewport()
    
    # Gráfico UpSet
    pushViewport(viewport(layout.pos.row = 2))
    print(
      upset(
        obj$matriz_filtrada %>% select(-query_name),
        sets = obj$dominios,
        sets.bar.color = obj$colores,
        matrix.color = "gray40",
        main.bar.color = "#4B0082",
        order.by = "freq",
        text.scale = c(1.3, 1.2, 1, 1, 1.3, 1),
        mb.ratio = c(0.7, 0.3),
        nintersects = nrow(obj$combinaciones),
        point.size = 3.5,
        line.size = 1.2
      )
    )
    popViewport()
    
    dev.off()
    
    cat("✅ Guardados: matriz, combinaciones, objeto y gráfico PDF para ", fam, "\n\n")
    
  }, error = function(e) {
    message("Error generando gráfico para ", fam, ": ", e$message)
    cat("❌ Error generando gráfico\n")
  })
}

cat("🎯 Procesamiento finalizado para ", length(familias), " familias.\n")
