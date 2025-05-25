###Filtrado por edades:
## CARGA DE LIBRERÍAS NECESARIAS
library(igraph)
library(igraphdata)
library(networkD3)
library(RCy3)
library(tidyverse)  
library(microbiomeDataSets)
library(mia)
library(phyloseq)
library(SpiecEasi)

#-------------------------------------------------------------------------------
#Funciones y datos:
load("DATOS/physeq_adultos.RData")
load("DATOS/physeq_infantes.RData")
load("DATOS/physeq_niños.RData")

source("FUNCIONES.R")
#-------------------------------------------------------------------------------
#NIÑOS
physeq_spp_niños <- tax_glom(physeq_niños, taxrank = "Species", NArm = FALSE)

physeq_niños.sin     <- subset_taxa(physeq_spp_niños, !is.na(Species))  # sin materia oscura
physeq_niños.com       <- physeq_spp_niños  # con materia oscura


physeq_conocido_filtrado_niños <- filtrar_prevalencia(physeq_niños.sin, 0.1)
physeq_todos_filtrado_niños    <- filtrar_prevalencia(physeq_niños.com, 0.1)

# Renombrar los taxa en physeq_especie (aplicable a physeq_todos y physeq_conocido)
physeq_conocido_filtrado_niños <- renombrar_especies(physeq_conocido_filtrado_niños)
physeq_todos_filtrado_niños    <- renombrar_especies(physeq_todos_filtrado_niños)

# CONSTRUCCIÓN DE REDES CON SPIEC-EASI
red_conocido_niños <- spiec.easi(
  physeq_conocido_filtrado_niños,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

red_todos_niños <- spiec.easi(
  physeq_todos_filtrado_niños,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

g_conocido.niños <- adj2igraph(getRefit(red_conocido_niños), vertex.attr = list(name = taxa_names(physeq_conocido_filtrado_niños)))
g_todos.niños    <- adj2igraph(getRefit(red_todos_niños),    vertex.attr = list(name = taxa_names(physeq_todos_filtrado_niños)))

metricas_conocido_niños <- calcular_metricas(g_conocido.niños)
metricas_todos_niños    <- calcular_metricas(g_todos.niños)

tabla_comparativa.niños <- tibble(
  Métrica              = names(metricas_todos_niños),
  Con_materia_obscura  = unlist(metricas_todos_niños),
  Sin_materia_obscura  = unlist(metricas_conocido_niños)
)
print(tabla_comparativa.niños)

plot(g_todos.niños,
     vertex.size = degree(g_todos.niños)*2,
     vertex.color = cluster_louvain(g_todos.niños)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red con Materia Obscura")

plot(g_conocido.niños,
     vertex.size = degree(g_conocido.niños)*2,
     vertex.color = cluster_louvain(g_conocido.niños)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red sin Materia Obscura")


# GRAFICAR COMPARACIÓN DE MÉTRICAS CON ggplot2
library(ggplot2)
library(tidyr)

df_metricas.niños <- tabla_comparativa.niños %>%
  pivot_longer(cols = -Métrica, names_to = "Red", values_to = "Valor")

ggplot(df_metricas.niños, aes(x = Métrica, y = Valor, fill = Red)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Comparación de métricas de las redes",
       x = "Métrica",
       y = "Valor",
       fill = "Tipo de red") +
  coord_flip()

#--------------------------------------------------------------------------------------------------
#INFANTES:
