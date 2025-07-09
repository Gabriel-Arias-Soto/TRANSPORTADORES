rm(list = ls())
cat("\014")

# ----------------------------
# Cargar librerías
# ----------------------------
required_packages <- c("readr", "dplyr", "ggplot2", "tidyr", "forcats", "UpSetR", "RColorBrewer")
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
# Detectar familias
# ----------------------------
familias <- unique(pfam_data$Familia)
output_dir_base <- dirname(file_path)

# ----------------------------
# Loop por familia
# ----------------------------
for (fam in familias) {
  
  cat(paste0("\n🔧 Procesando familia: ", fam, "...\n"))
  
  pfam_sub <- pfam_data %>% filter(Familia == fam)
  fam_dir <- file.path(output_dir_base, paste0("graficos_", fam))
  if (!dir.exists(fam_dir)) dir.create(fam_dir)
  
  # 1. Histograma de dominios por proteína
  dominios_por_proteina <- pfam_sub %>%
    group_by(query_name) %>%
    summarise(n_dominios = n_distinct(target_name), .groups = "drop") %>%
    arrange(desc(n_dominios))
  
  g1 <- ggplot(dominios_por_proteina, aes(x = n_dominios)) +
    geom_histogram(binwidth = 1, fill = "#0072B2", alpha = 0.8, color = "black") +
    geom_vline(aes(xintercept = mean(n_dominios)), color = "orange", linetype = "dashed", size = 1) +
    labs(
      title = paste("Distribución de dominios Pfam en proteínas de la familia", fam),
      subtitle = "Cada barra representa la cantidad de proteínas con un número específico de dominios únicos",
      x = "N° de dominios únicos por proteína",
      y = "Frecuencia"
    ) +
    theme_bw()
  
  ggsave(file.path(fam_dir, "histograma_dominios.png"), g1, width = 8, height = 6, bg = "white")
  
  # 2. Top dominios
  top_dominios <- pfam_sub %>%
    count(target_name, sort = TRUE, name = "n_proteinas") %>%
    mutate(porcentaje = n_proteinas / nrow(dominios_por_proteina) * 100)
  
  g2 <- ggplot(top_dominios %>% head(15), aes(x = fct_reorder(target_name, n_proteinas), y = n_proteinas, fill = porcentaje)) +
    geom_col(color = "black") +
    coord_flip() +
    scale_fill_gradient(low = "skyblue", high = "darkblue") +
    labs(
      title = paste("Top 15 dominios Pfam en la familia", fam),
      subtitle = "Número de proteínas que presentan cada dominio",
      x = "Dominio Pfam",
      y = "N° de proteínas",
      fill = "Porcentaje (%)"
    ) +
    theme_bw()
  
  ggsave(file.path(fam_dir, "top_dominios.png"), g2, width = 8, height = 6, bg = "white")
  
  # 3. Distribución de scores HMM
  top_10_dominios <- top_dominios$target_name[1:10]
  
  g3 <- pfam_sub %>%
    filter(target_name %in% top_10_dominios) %>%
    ggplot(aes(x = fct_reorder(target_name, score_full, .fun = median), y = score_full)) +
    geom_boxplot(aes(fill = target_name), alpha = 0.7, show.legend = FALSE, color = "black") +
    stat_summary(fun = median, geom = "point", color = "red", size = 2) +
    labs(
      title = paste("Distribución de scores HMM en la familia", fam),
      subtitle = "Top 10 dominios Pfam con mayor presencia",
      x = "Dominio Pfam",
      y = "Score HMM (full)"
    ) +
    coord_flip() +
    theme_bw()
  
  ggsave(file.path(fam_dir, "scores_hmm.png"), g3, width = 8, height = 6, bg = "white")
  
  # 4. UpSetPlot con filtros configurables
  # ----------------------------
  # FILTROS DISPONIBLES - ELIJE SOLO UNO (DESCOMENTA EL DESEADO)
  # ----------------------------
  
  ## 1. Top 10 proteínas por score promedio (riguroso, enfoque en alta confianza)
  proteinas_filtradas <- pfam_sub %>%
    group_by(query_name) %>%
    summarise(mean_score = mean(score_full), .groups = "drop") %>%
    arrange(desc(mean_score)) %>%
    slice_head(n = 10) %>%
    pull(query_name)
  
  ## 2. Alternativa: Proteínas con ≥3 dominios únicos (estructuralmente complejas)
  # proteinas_filtradas <- dominios_por_proteina %>%
  #   filter(n_dominios >= 3) %>%
  #   pull(query_name)
  
  ## 3. Alternativa: Proteínas con score promedio ≥ 200
  # proteinas_filtradas <- pfam_sub %>%
  #   group_by(query_name) %>%
  #   summarise(mean_score = mean(score_full), .groups = "drop") %>%
  #   filter(mean_score >= 200) %>%
  #   pull(query_name)
  
  if (length(proteinas_filtradas) > 1) {
    
    cat("🔬 Proteínas seleccionadas para UpSetPlot:\n")
    print(proteinas_filtradas)
    
    binaria <- pfam_sub %>%
      filter(query_name %in% proteinas_filtradas) %>%
      distinct(query_name, target_name) %>%
      mutate(presente = 1) %>%
      pivot_wider(names_from = target_name, values_from = presente, values_fill = 0)
    
    if (nrow(binaria) > 1 && ncol(binaria) > 2) {
      mat <- as.data.frame(binaria)
      rownames(mat) <- mat$query_name
      mat$query_name <- NULL
      
      dominios <- colnames(mat)
      n_dominios <- length(dominios)
      colores <- colorRampPalette(brewer.pal(12, "Set3"))(n_dominios)
      names(colores) <- dominios
      
      pdf(file.path(fam_dir, "upsetplot_filtrado.pdf"), width = 14, height = 8)
      par(bg = "white")
      upset(mat,
            sets = dominios,
            sets.bar.color = colores,
            matrix.color = "gray20",
            main.bar.color = "#4B0082",
            order.by = "freq",
            text.scale = 1.4,
            keep.order = TRUE)
      dev.off()
      
      if (file.exists(file.path(fam_dir, "upsetplot_filtrado.pdf"))) {
        cat("📊 UpSetPlot generado exitosamente para familia ", fam, "\n")
      } else {
        cat("⚠️ Hubo un problema generando el UpSetPlot para familia ", fam, "\n")
      }
      
    } else {
      cat(paste0("⚠️ No hay suficientes datos para UpSetPlot en familia ", fam, "\n"))
    }
    
  } else {
    cat(paste0("⚠️ No hay suficientes proteínas para UpSetPlot en familia ", fam, "\n"))
  }
  
  cat(paste0("✅ Análisis completado para familia: ", fam, "\n"))
}

cat("\n🎯 Todos los análisis por familia han sido generados.\n")
