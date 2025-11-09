# LDA

## Script de R para hacer un **Análisis de Discriminación lineal (LDA)** a partir de variantes y metadatos del gen *babA* de *Helicobacter pylori*.

El script lleva a cabo el LDA, el cual es un análisis supervisado de reducción de dimensión, a partir de Polimorfismos de un Solo Nucleótido (SNPs) de genes *babA*, donde las etiquetas utilizadas son los fenotipos de enfermedad asociados a las cepas de *Helicobacter pylori* que pertenecen los genes *babA*. 

Los input son los siguientes:
- `snps_babA_730.csv` : SNPs de los genes.
- `metadatos.csv` : archivo donde vienen los fenotipos de enfermedad para etiquetar a los SNPs.


Los resultados del LDA se guardan en el archivo `LDA_result_snps_babA.csv` y se construye la gráfica, la cual es la siguiente:

<img width="2400" height="1800" alt="LDA_snps_babA" src="https://github.com/user-attachments/assets/9808850d-552e-4254-86ff-f362086d52a5" />


Se procede a extraer los coeficientes del LDA con la finalidad de extraer aquellos SNPs relevantes en la separación de los grupos en el LDA. Se grafica la distribución de los coeficientes de importancia de los SNPs con un umbral de corte de los 5% más importantes:

<img width="2700" height="1800" alt="distribucion_importancia_snps_corte" src="https://github.com/user-attachments/assets/a56da6d2-fb21-4fe9-8c34-3a35b55f94bf" />

Por último, se guardan en el archivo `snps_mas_importantes_lda.csv` el 5% de los SNPs más importantes.

