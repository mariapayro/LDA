# Limpiar el entorno de trabajo
rm(list=ls())

# Cargar paquetes necesarios
library(MASS)  
# Para LDA
library(ggplot2)  # Para visualizar
library(readr)  # Para leer archivos CSV y TSV
library(dplyr)  # Para manipulación de datos
library(plotly)
library(htmlwidgets)

# --- Lectura de datos ---
snps <- read_delim("snps_babA_730.tsv", delim = "\t")
metadatos <- read_csv("metadatos.csv")

# --- Ordenar datos por ID (si es necesario) ---
metadatos <- metadatos %>% arrange(ID)

# --- Extraer variables predictoras (SNPs) y variable de clasificación ---

# script sin eliminar primera columna

# Conservar toda la tabla, incluyendo la primera columna
X <- snps  # No eliminamos la primera columna,  antes era snps[, -1]
# Conserva la primera columna como identificador en los resultados
rownames(X) <- snps$POS 
# Ahora sí se elimina primera columna
X <- snps[, -1]  # Excluir la primera columna si es identificador
# volver a correr esto? solo asi funcionó
rownames(X) <- snps$POS 

# Convertir la variable de clasificación a factor
y <- as.factor(metadatos$Phenotype)  

# Transponer los datos, ahora X incluye todas las columnas
X2 <- t(X)  # Excluimos solo la columna de identificadores para la transposición



# --- Realizar LDA ---
lda_model <- lda(X2, grouping = y)

# --- Transformar los datos usando LDA ---
lda_values <- predict(lda_model)$x

lda_df <- as.data.frame(lda_values)
lda_df$Phenotype <- y  # Agregar la variable de clasificación

# --- Guardar resultados en un archivo CSV ---
write.csv(lda_df, "LDA_result_snps_babA.csv", row.names = FALSE)


# colores anteriores #8BE2E5 NAG, #f49c1d IM, #AC2007 GC; '#AC2007', '#f49c1d', '#8BE2E5'

# Definir los colores por fenotipos
colors <- c('#c9270a', '#f59407', '#44bec2')

# Crear el gráfico LDA con ggplot2
ldaplot <- ggplot(lda_df, aes(x = LD1, y = LD2, color = Phenotype)) +
  geom_point(size = 3) +
  labs(title = "LDA Plot babA by SNPs",
       x = "LD1",
       y = "LD2") +
  scale_color_manual(values = colors) +  # Asignar los colores manualmente
  theme_minimal()

# Convertir el gráfico ggplot a plotly
gg_plotly <- ggplotly(ldaplot)
print(gg_plotly)

# Guardar el gráfico ggplot como imagen
ggsave(filename = 'LDA_snps_babA.png', plot = ldaplot, width = 8, height = 6, units = 'in', bg = 'white')




# Script para extraer los SNPs relevantes que permiten la separación de los grupos en el LDA
  
# Extraer los coeficientes del LDA ---
coeficientes_lda <- as.data.frame(lda_model$scaling)

# Calcular la magnitud de cada SNP ---
coeficientes_lda$Importancia <- apply(coeficientes_lda, 1, function(x) sqrt(sum(x^2)))

# Agregar nombres de los SNPs ---
coeficientes_lda$SNP <- rownames(coeficientes_lda)

# Ordenar de mayor a menor importancia ---
coeficientes_lda <- coeficientes_lda %>%
  arrange(desc(Importancia))

head(coeficientes_lda)

# Guardar en un archivo CSV todos los LDAs, solo usar si no se usa el umbral (abajo)
# si no, dejar comentado para que umbral funcione
#write_csv(coeficientes_lda, "snps_importancia_lda.csv")



# Gráficos para visualizar la importancia (extras)

# Gráfico de dispersión de la importancia de los snps
ggplot(coeficientes_lda, aes(x = Importancia, y = SNP)) +
  geom_point(color = "red", size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Importancia de los SNPs en LDA",
    x = "Importancia",
    y = "SNP"
  )


# Curva de distribución de importancia por densidad 
# Eje Y = probabilidad relativa, Suaviza la distribución de los valores sin depender de bins (intervalos).
ggplot(coeficientes_lda, aes(x = Importancia)) +
  geom_density(fill = "blue", alpha = 0.5) +  # Curva de densidad
  theme_minimal() +
  labs(
    title = "Distribución de Importancia de SNPs",
    x = "Importancia",
    y = "Densidad"
  )

# Histograma por conteo, Muestra la frecuencia de valores en intervalos específicos.
ggplot(coeficientes_lda, aes(x = Importancia)) +
  geom_histogram(binwidth = 1, fill = "purple", color = "black", alpha = 0.7) + 
  theme_minimal() +
  labs(
    title = "Distribución de Importancia de SNPs",
    x = "Importancia",
    y = "Frecuencia"
  )

# Histograma + Curva de densidad
ggplot(coeficientes_lda, aes(x = Importancia)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, fill = "purple", color = "black", alpha = 0.5) +
  geom_density(color = "blue", linewidth = 1) +  
  theme_minimal() +
  labs(
    title = "Histograma + Densidad de Importancia",
    x = "Importancia",
    y = "Densidad"
  )


# Gráfico relevante

# Calcular el umbral del 5% superior (percentil 95)
umbral_5 <- quantile(coeficientes_lda$Importancia, 0.95)

# Graficar la densidad con una línea de corte en el percentil 95
distribucion_snps_plot <- ggplot(coeficientes_lda, aes(x = Importancia)) +
  geom_density(fill = "blue", alpha = 0.5, color = "black") +  
  geom_vline(xintercept = umbral_5, color = "red", linetype = "dashed", linewidth = 1) +  
  theme_minimal() +
  labs(
    title = "Distribución de Importancia de SNPs",
    x = "Importancia",
    y = "Densidad"
  ) +
  annotate("text", x = umbral_5, y = 0.275, 
           label = paste("Valor de corte (5%) =", round(umbral_5, 2)), 
           color = "red", hjust = -0.1)
print(distribucion_snps_plot)

# Guardar el gráfico como imagen
ggsave(filename = 'distribucion_importancia_snps_corte.png', plot = distribucion_snps_plot, width = 9, height = 6, units = 'in', bg = 'white')




# Filtrar los SNPs más importantes 5% (0.05) ya calculados en 'umbral_5'
snps_top5 <- coeficientes_lda %>%
  filter(Importancia >= umbral_5)

# Guardar en un archivo CSV
write.csv(snps_top5, "snps_mas_importantes_lda.csv", row.names = FALSE)

# Ver los primeros registros
head(snps_top5)

