rm(list = ls())
cat("\014")

# ----------------------------
# Cargar librerías
# ----------------------------
required_packages <- c("readr", "dplyr", "ggplot2", "tidyr", "forcats", "RColorBrewer")
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
      title = paste("Distribución de dominios Pfam en proteínas - Familia", fam),
      subtitle = paste0("Total de proteínas: ", nrow(dominios_por_proteina)),
      x = "N° de dominios únicos por proteína",
      y = "Frecuencia"
    ) +
    theme_bw()
  
  ggsave(file.path(fam_dir, "histograma_dominios.png"), g1, width = 8, height = 6, bg = "white")
  
  # 2. Top dominios
  top_n <- min(15, nrow(pfam_sub %>% count(target_name)))
  
  top_dominios <- pfam_sub %>%
    count(target_name, sort = TRUE, name = "n_proteinas") %>%
    mutate(porcentaje = n_proteinas / nrow(dominios_por_proteina) * 100)
  
  g2 <- ggplot(top_dominios %>% head(top_n), aes(x = fct_reorder(target_name, n_proteinas), y = n_proteinas, fill = porcentaje)) +
    geom_col(color = "black") +
    coord_flip() +
    scale_fill_gradient(low = "skyblue", high = "darkblue") +
    labs(
      title = paste0("Top ", top_n, " dominios Pfam detectados - Familia ", fam),
      subtitle = "Cantidad de proteínas que presentan cada dominio (escala proporcional al porcentaje de aparición)",
      x = "Dominio Pfam",
      y = "N° de proteínas",
      fill = "Porcentaje (%)"
    ) +
    theme_bw()+
    theme(plot.subtitle = element_text(size = 8.4))
  
  ggsave(file.path(fam_dir, "top_dominios.png"), g2, width = 8, height = 6, bg = "white")
  
  # 3. Distribución de scores HMM
  top_n_scores <- min(10, nrow(top_dominios))
  
  top_10_dominios <- top_dominios$target_name[1:top_n_scores]
  
  g3 <- pfam_sub %>%
    filter(target_name %in% top_10_dominios) %>%
    ggplot(aes(x = fct_reorder(target_name, score_full, .fun = median), y = score_full)) +
    geom_boxplot(aes(fill = target_name), alpha = 0.7, show.legend = FALSE, color = "black") +
    stat_summary(fun = median, geom = "point", color = "red", size = 2) +
    labs(
      title = paste0("Distribución de scores HMM - Top ", top_n_scores, " dominios en familia ", fam),
      subtitle = "Se muestra la variabilidad de los scores HMM para los dominios más frecuentes",
      x = "Dominio Pfam",
      y = "Score HMM (full)"
    ) +
    coord_flip() +
    theme_bw()
  
  ggsave(file.path(fam_dir, "scores_hmm.png"), g3, width = 8, height = 6, bg = "white")
  
  cat(paste0("✅ Análisis exploratorio completado para familia: ", fam, "\n"))
}

cat("\n🎯 Todos los análisis exploratorios por familia han sido generados.\n")

