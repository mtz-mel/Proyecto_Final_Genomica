## Redes de co-ocurrencia con y sin materia obscura, tomando en cuenta todos los
## datos de SprockettTHData, los cuales, viven como objeto phyloseq en physeq_spp

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
dataspo <- SprockettTHData()

# MATRICES
otu <- assays(dataspo)$counts %>% as.matrix()
tax <- rowData(dataspo) %>% as.data.frame() %>% as.matrix()
sample_data_df <- colData(dataspo) %>% as.data.frame()

# CONSTRUCCIÓN DE PHYLOSEQ
OTU     <- otu_table(otu, taxa_are_rows = TRUE)
TAX     <- tax_table(tax)
SAMPLES <- sample_data(sample_data_df)
physeq <- phyloseq(OTU, TAX, SAMPLES)

# AGRUPACIÓN A NIVEL DE ESPECIE (INCLUYENDO NA)
physeq_spp <- tax_glom(physeq, taxrank = "Species", NArm = FALSE)

# DIVIDIR physeq_spp

physeq_known     <- subset_taxa(physeq_spp, !is.na(Species))  # sin tomar en cuenta la materia oscura
physeq_all       <- physeq_spp  # tomando en cuenta la materia oscura


# FILTRADO DE PREVALENCIA

prevalence_filter <- function(phy, threshold = 0.2) {
  prev <- apply(otu_table(phy), 1, function(x) mean(x > 0))
  keep <- names(prev[prev >= threshold])
  prune_taxa(keep, phy)
}

physeq_known_filt <- prevalence_filter(physeq_known, 0.2)
physeq_all_filt   <- prevalence_filter(physeq_all, 0.2)

# REDES SPIEC-EASI

se_known <- spiec.easi(
  physeq_known_filt,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

se_all <- spiec.easi(
  physeq_all_filt,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

g_known <- adj2igraph(getRefit(se_known), vertex.attr = list(name = taxa_names(physeq_known_filt)))
g_all   <- adj2igraph(getRefit(se_all),   vertex.attr = list(name = taxa_names(physeq_all_filt)))


# MÉTRICAS COMPARATIVAS DE LA RED

metricas <- function(g) {
  list(
    nodos = vcount(g),  # cuantos nodos
    aristas = ecount(g),  # Número de aristas
    grado_medio = mean(degree(g)),  # degree
    densidad = edge_density(g),  # Densidad de la red
    clustering = transitivity(g, type = "global"),  # coeficiente de clustering global
    modularidad = modularity(cluster_louvain(g)),  # Modularidad
    closeness_media = mean(closeness(g)),  # Promedio de closeness
    betweenness_media = mean(betweenness(g))  # Promedio de betweenness
  )
}


m_known <- metricas(g_known)
m_known

m_all   <- metricas(g_all)
m_all

# Mostrar comparación
tibble(
  Métrica        = names(m_all),
  Con_materia_obscura = unlist(m_all),
  Sin_materia_obscura = unlist(m_known)
)



# VISUALIZACIÓN DE g_all (con materia obscura)
plot(g_all,
     vertex.size = degree(g_all)*2,
     vertex.color = cluster_louvain(g_all)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red con Materia Obscura")

# VISUALIZACIÓN DE g_known (sin materia obscura)
plot(g_known,
     vertex.size = degree(g_known)*2,
     vertex.color = cluster_louvain(g_known)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red sin Materia Obscura")



createNetworkFromIgraph(g_all,   title = "Red con materia obscura")
createNetworkFromIgraph(g_known, title = "Red sin materia obscura")



####################
#comparacion

library(ggplot2)
library(tidyr)

# Crear tabla de métricas en formato largo
df_metricas <- tibble(
  Métrica        = names(m_all),
  Con_materia_obscura = unlist(m_all),
  Sin_materia_obscura = unlist(m_known)
) %>% pivot_longer(cols = -Métrica, names_to = "Red", values_to = "Valor")

# Ver estructura de los datos
print(df_metricas)

#hacer la grafica con ggplot
ggplot(df_metricas, aes(x = Métrica, y = Valor, fill = Red)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Comparación de métricas de las redes",
       x = "Métrica",
       y = "Valor",
       fill = "Tipo de red") +
  coord_flip() 