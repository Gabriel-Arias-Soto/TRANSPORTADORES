rm(list = ls())
cat("\014")

# -----------------------
# Cargar librerías
# -----------------------
required_packages <- c("readr", "dplyr", "ggplot2", "tidyr", "forcats", "UpSetR", "RColorBrewer")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, library, character.only = TRUE)

# -----------------------
# Leer y preparar datos
# -----------------------
file_path <- "C:/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/hmmscan_resultados/pfam_transportadores.csv"

pfam_data <- read_csv(file_path, show_col_types = FALSE) %>%
  mutate(
    Especie = sub("_.*", "", query_name),
    Familia = sub(".*_(.*)_.*", "\\1", query_name)
  )

# -----------------------
# 1. Número de dominios únicos por proteína
# -----------------------
dominios_por_proteina <- pfam_data %>%
  group_by(query_name) %>%
  summarise(n_dominios = n_distinct(target_name), .groups = "drop") %>%
  arrange(desc(n_dominios))

ggplot(dominios_por_proteina, aes(x = n_dominios)) +
  geom_histogram(binwidth = 1, fill = "#0072B2", alpha = 0.8) +
  geom_vline(aes(xintercept = mean(n_dominios)), color = "orange", linetype = "dashed") +
  labs(title = "Distribución de dominios Pfam por proteína", x = "N° de dominios únicos", y = "Frecuencia") +
  theme_minimal()

# -----------------------
# 2. Conteo de dominios por especie
# -----------------------
conteo_dominios <- pfam_data %>%
  group_by(Especie, target_name) %>%
  summarise(n_veces = n(), .groups = "drop") %>%
  arrange(Especie, desc(n_veces))

print(conteo_dominios)

# -----------------------
# 3. Top dominios globales
# -----------------------
top_dominios <- pfam_data %>%
  count(target_name, sort = TRUE, name = "n_proteinas") %>%
  mutate(porcentaje = n_proteinas / nrow(dominios_por_proteina) * 100)

ggplot(top_dominios %>% head(15), aes(x = fct_reorder(target_name, n_proteinas), y = n_proteinas, fill = porcentaje)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "skyblue", high = "darkblue") +
  labs(title = "Top 15 dominios más frecuentes", x = "Dominio", y = "N° de proteínas", fill = "%") +
  theme_minimal()

# -----------------------
# 4. Distribución de scores HMM por dominio (global)
# -----------------------
top_10_dominios <- top_dominios$target_name[1:10]

pfam_data %>%
  filter(target_name %in% top_10_dominios) %>%
  ggplot(aes(x = fct_reorder(target_name, score_full, .fun = median), y = score_full)) +
  geom_boxplot(aes(fill = target_name), alpha = 0.7, show.legend = FALSE) +
  stat_summary(fun = median, geom = "point", color = "red", size = 2) +
  labs(
    title = "Distribución de scores HMM por dominio Pfam",
    subtitle = "Top 10 dominios más frecuentes",
    x = "Dominio Pfam",
    y = "Score HMM (full)"
  ) +
  coord_flip() +
  theme_minimal()

# ----------------------------
# 5. UpSetPlot filtrado y coloreado
# ----------------------------

# DOMINIOS DE INTERÉS (ajustable según análisis previo)
dominios_top5 <- top_dominios$target_name[1:5]

# ----------------------------
# FILTROS DISPONIBLES - ELIGE SOLO UNO
# ----------------------------

## 1. Proteínas que contienen al menos un dominio del Top 5 (moderadamente amplio)
#proteinas_filtradas <- pfam_data %>%
#  filter(target_name %in% dominios_top5) %>%
#  distinct(query_name) %>%
#  pull(query_name)

## 2. Proteínas con ≥3 dominios únicos (más restringido, proteínas complejas)
# proteinas_filtradas <- dominios_por_proteina %>%
# filter(n_dominios >= 3) %>%
# pull(query_name)

# 3. Proteínas cuyo score HMM promedio supera 100 (ejemplo de filtro por score de confianza)
 proteinas_filtradas <- pfam_data %>%
   group_by(query_name) %>%
   summarise(mean_score = mean(score_full), .groups = "drop") %>%
   filter(mean_score >= 200) %>%
   pull(query_name)

## 4. Solo proteínas de una especie específica, por ejemplo Arabidopsis (ajustar especie)
# proteinas_filtradas <- pfam_data %>%
#   filter(Especie == "HMA-ATH-NP") %>%
#   distinct(query_name) %>%
#   pull(query_name)

# ----------------------------
# CREACIÓN MATRIZ Y GRAFICADO
# ----------------------------

binaria <- pfam_data %>%
  filter(query_name %in% proteinas_filtradas) %>%
  distinct(query_name, target_name) %>%
  mutate(presente = 1) %>%
  pivot_wider(names_from = target_name, values_from = presente, values_fill = 0)

mat <- as.data.frame(binaria)
rownames(mat) <- mat$query_name
mat$query_name <- NULL

# Paleta de colores proporcional al número de dominios
dominios <- colnames(mat)
n_dominios <- length(dominios)
colores <- colorRampPalette(brewer.pal(12, "Set3"))(n_dominios)
names(colores) <- dominios

# Guardar gráfico
pdf("C:/Users/gabri/Documents/SCRIPTS_ESTUDIO_PYTHON/TRANSPORTADORES/data_proteomas/analisis_filogenia/hmmscan_resultados/upset_colores_filtrado.pdf", width = 14, height = 8)
upset(mat,
      sets = dominios,
      sets.bar.color = colores,
      matrix.color = "gray20",
      main.bar.color = "#4B0082",
      order.by = "freq",
      text.scale = 1.4,
      keep.order = TRUE)
dev.off()

cat("\n✅ UpSetPlot filtrado completado.\n")

