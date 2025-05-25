## Redes de co-ocurrencia con y sin materia obscura, tomando en cuenta todos los
## datos de SprockettTHData, los cuales, viven como objeto phyloseq en physeq_especie

## Procedimiento: 

# CARGA DE LIBRERÍAS NECESARIAS
library(igraph)
library(igraphdata)
library(networkD3)
library(RCy3)
library(tidyverse)  
library(microbiomeDataSets)
library(mia)
library(phyloseq)
library(SpiecEasi)

# CARGA DE DATOS
datos_sprockett <- SprockettTHData()

ncol(datos_sprockett)

# MATRICES
matriz_otu <- assays(datos_sprockett)$counts %>% as.matrix()
matriz_tax <- rowData(datos_sprockett) %>% as.data.frame() %>% as.matrix()
datos_muestras <- colData(datos_sprockett) %>% as.data.frame()

# CONSTRUCCIÓN DE OBJETO PHYLOSEQ
OTU <- otu_table(matriz_otu, taxa_are_rows = TRUE)
TAX <- tax_table(matriz_tax)
MUESTRAS <- sample_data(datos_muestras)
physeq <- phyloseq(OTU, TAX, MUESTRAS)

# AGRUPACIÓN A NIVEL DE ESPECIE (INCLUYENDO NA)
physeq_especie <- tax_glom(physeq, taxrank = "Species", NArm = FALSE)

# DIVISIÓN SEGÚN PRESENCIA DE MATERIA OSCURA
physeq_conocido <- subset_taxa(physeq_especie, !is.na(Species))  # sin materia obscura
physeq_todos <- physeq_especie  # con materia obscura

# FUNCIÓN DE FILTRADO POR PREVALENCIA
filtrar_prevalencia <- function(physeq_objeto, umbral = 0.1) {
  prevalencia <- apply(otu_table(physeq_objeto), 1, function(x) mean(x > 0))
  conservar <- names(prevalencia[prevalencia >= umbral])
  prune_taxa(conservar, physeq_objeto)
}

physeq_conocido_filtrado <- filtrar_prevalencia(physeq_conocido, 0.1)
physeq_todos_filtrado    <- filtrar_prevalencia(physeq_todos, 0.1)

# Renombrar los taxa en physeq_especie (aplicable a physeq_todos y physeq_conocido)

renombrar_especies <- function(physeq_objeto) {
  tax <- tax_table(physeq_objeto)
  especies <- as.character(tax[, "Species"])
  
  # Reemplazar NA con "Desconocido"
  especies_limpias <- ifelse(is.na(especies), "Desconocido", especies)
  
  # Asegurar unicidad en los nombres
  nombres_unicos <- make.unique(especies_limpias)
  
  # Asignar los nuevos nombres como nombres de los taxones
  taxa_names(physeq_objeto) <- nombres_unicos
  return(physeq_objeto)
}

# Aplicar renombramiento
physeq_conocido_filtrado <- renombrar_especies(physeq_conocido_filtrado)
physeq_todos_filtrado    <- renombrar_especies(physeq_todos_filtrado)

# CONSTRUCCIÓN DE REDES CON SPIEC-EASI
red_conocido <- spiec.easi(
  physeq_conocido_filtrado,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

red_todos <- spiec.easi(
  physeq_todos_filtrado,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

# CONVERTIR A OBJETOS igraph
g_conocido <- adj2igraph(getRefit(red_conocido), vertex.attr = list(name = taxa_names(physeq_conocido_filtrado)))
g_todos    <- adj2igraph(getRefit(red_todos),    vertex.attr = list(name = taxa_names(physeq_todos_filtrado)))

# FUNCIÓN PARA CALCULAR MÉTRICAS DE LA RED
calcular_metricas <- function(grafo) {
  list(
    nodos = vcount(grafo),
    aristas = ecount(grafo),
    grado_medio = mean(degree(grafo)),
    densidad = edge_density(grafo),
    clustering = transitivity(grafo, type = "global"),
    modularidad = modularity(cluster_louvain(grafo)),
    closeness_media = mean(closeness(grafo)),
    betweenness_media = mean(betweenness(grafo))
  )
}

metricas_conocido <- calcular_metricas(g_conocido)
metricas_todos    <- calcular_metricas(g_todos)

# COMPARACIÓN DE MÉTRICAS
tabla_comparativa <- tibble(
  Métrica              = names(metricas_todos),
  Con_materia_obscura  = unlist(metricas_todos),
  Sin_materia_obscura  = unlist(metricas_conocido)
)
print(tabla_comparativa)

# VISUALIZACIÓN DE LAS REDES
plot(g_todos,
     vertex.size = degree(g_todos)*2,
     vertex.color = cluster_louvain(g_todos)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red con Materia Obscura")

plot(g_conocido,
     vertex.size = degree(g_conocido)*2,
     vertex.color = cluster_louvain(g_conocido)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red sin Materia Obscura")

# EXPORTAR A CYTOSCAPE
createNetworkFromIgraph(g_todos,    title = "Red con materia obscura")
createNetworkFromIgraph(g_conocido, title = "Red sin materia obscura")

# GRAFICAR COMPARACIÓN DE MÉTRICAS CON ggplot2
library(ggplot2)
library(tidyr)

df_metricas <- tabla_comparativa %>%
  pivot_longer(cols = -Métrica, names_to = "Red", values_to = "Valor")

ggplot(df_metricas, aes(x = Métrica, y = Valor, fill = Red)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Comparación de métricas de las redes",
       x = "Métrica",
       y = "Valor",
       fill = "Tipo de red") +
  coord_flip()
