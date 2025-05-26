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

source("FUNCIONES/FUNCIONES.R")
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
     main = "Red con Materia Obscura en niños")

plot(g_conocido.niños,
     vertex.size = degree(g_conocido.niños)*2,
     vertex.color = cluster_louvain(g_conocido.niños)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red sin Materia Obscura en niños")

library(RCy3)
cytoscapePing()

createNetworkFromIgraph(g_todos.niños,   title = "Red con materia oscura en niños")
createNetworkFromIgraph(g_conocido.niños, title = "Red sin materia oscura en niños")

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


# REDES CON BOOSTRAP

# función para eliminar un porcentaje aleatorio de taxones
eliminar_taxones_azar <- function(physeq, porcentaje = 0.3) {
  # extraemos los nombres de taxones
  taxones <- taxa_names(physeq)
  
  # Seleccionar 30% de taxones al azar
  eliminar <- sample(taxones, size = round(length(taxones) * porcentaje))
  
  # filtrar el phyloseq sin esos taxones
  prune_taxa(setdiff(taxones, eliminar), physeq)
}



#una vez Hecha la funcion la aplicamos para eliminar 30% de taxones

physeq_reducido.niños <- eliminar_taxones_azar(physeq_conocido_filtrado_niños)

# filtrado de prevalencia en la nueva red
physeq_reducido_filt.niños <- prevalence_filter(physeq_reducido.niños, 0.1)

# hacemos la red con el mismo SPIEC-EASI pero con el phyloseq reducido
se_reducido.niños <- spiec.easi(
  physeq_reducido_filt.niños,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

# convertimos a objeto igraph
g_reducido.niños <- adj2igraph(getRefit(se_reducido.niños), vertex.attr = list(name = taxa_names(physeq_reducido_filt.niños)))

# se alcular métricas de la nueva red
m_reducido.niños <- metricas(g_reducido.niños)

# Visualización de la red reducida
plot(g_reducido.niños,
     vertex.size = degree(g_reducido.niños)*2,
     vertex.color = cluster_louvain(g_reducido.niños)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red con 30% de taxones eliminados")

# Crear red en Cytoscape pa verla más bonita 
  #createNetworkFromIgraph(g_reducido.niños, title = "Red con 30% de taxones eliminados")

#prueba t de Student 

str(metricas_conocido_niños)
str(metricas_todos_niños)


t_test_resultado.niños <- t.test(unlist(metricas_conocido_niños), unlist(metricas_todos_niños), paired = TRUE)

# Mostrar resultados
tibble(
  Métrica = names(metricas_todos_niños),
  P_value_t_Test = t_test_resultado.niños$p.value
) %>% arrange(P_value_t_Test)

#NO es significativo
unlist(metricas_conocido_niños$grado_medio)
unlist(metricas_todos_niños$grado_medio)

t.test(2.35483871, 6.555556)#, paired = TRUE)

#--------------------------------------------------------------------------------------------------
#INFANTES:
physeq_spp_infantes <- tax_glom(physeq_infantes, taxrank = "Species", NArm = FALSE)

physeq_infantes.sin     <- subset_taxa(physeq_spp_infantes, !is.na(Species))  # sin materia oscura
physeq_infante.com       <- physeq_spp_infantes  # con materia oscura


physeq_conocido_filtrado.infantes <- filtrar_prevalencia(physeq_infantes.sin, 0.1)
physeq_todos_filtrado.infantes   <- filtrar_prevalencia(physeq_infante.com, 0.1)

# Renombrar los taxa en physeq_especie (aplicable a physeq_todos y physeq_conocido)
physeq_conocido_filtrado.infantes <- renombrar_especies(physeq_conocido_filtrado.infantes)
physeq_todos_filtrado.infantes    <- renombrar_especies(physeq_todos_filtrado.infantes)

# CONSTRUCCIÓN DE REDES CON SPIEC-EASI
red_conocido_infantes <- spiec.easi(
  physeq_conocido_filtrado.infantes,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

red_todos_infantes <- spiec.easi(
  physeq_todos_filtrado.infantes,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

g_conocido.infantes <- adj2igraph(getRefit(red_conocido_infantes), vertex.attr = list(name = taxa_names(physeq_conocido_filtrado.infantes)))
g_todos.infantes  <- adj2igraph(getRefit(red_todos_infantes),    vertex.attr = list(name = taxa_names(physeq_todos_filtrado.infantes)))

metricas_conocido_infantes <- calcular_metricas(g_conocido.infantes)
metricas_todos_infantes    <- calcular_metricas(g_todos.infantes)

tabla_comparativa.infantes <- tibble(
  Métrica              = names(metricas_todos_infantes),
  Con_materia_obscura  = unlist(metricas_todos_infantes),
  Sin_materia_obscura  = unlist(metricas_conocido_infantes)
)
print(tabla_comparativa.infantes)

plot(g_todos.infantes,
     vertex.size = degree(g_todos.infantes)*2,
     vertex.color = cluster_louvain(g_todos.infantes)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red con Materia Obscura en infantes")

plot(g_conocido.infantes,
     vertex.size = degree(g_conocido.infantes)*2,
     vertex.color = cluster_louvain(g_conocido.infantes)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red sin Materia Obscura en infantes")

library(RCy3)
cytoscapePing()

createNetworkFromIgraph(g_todos.infantes,   title = "Red con materia oscura en infantes")
createNetworkFromIgraph(g_conocido.infantes, title = "Red sin materia oscura en infantes")

# GRAFICAR COMPARACIÓN DE MÉTRICAS CON ggplot2
library(ggplot2)
library(tidyr)

df_metricas.infantes <- tabla_comparativa.infantes %>%
  pivot_longer(cols = -Métrica, names_to = "Red", values_to = "Valor")

ggplot(df_metricas.infantes, aes(x = Métrica, y = Valor, fill = Red)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Comparación de métricas de las redes",
       x = "Métrica",
       y = "Valor",
       fill = "Tipo de red") +
  coord_flip()

#-----------------------------------------------------------------------------------------
#ADULTOS
physeq_spp_adultos <- tax_glom(physeq_adultos, taxrank = "Species", NArm = FALSE)

physeq_adultos.sin     <- subset_taxa(physeq_spp_adultos, !is.na(Species))  # sin materia oscura
physeq_adultos.com       <- physeq_spp_adultos  # con materia oscura

physeq_conocido_filtrado.adultos <- filtrar_prevalencia(physeq_adultos.sin, 0.1)
physeq_todos_filtrado.adultos   <- filtrar_prevalencia(physeq_adultos.com, 0.1)

# Renombrar los taxa en physeq_especie (aplicable a physeq_todos y physeq_conocido)
physeq_conocido_filtrado.adultos <- renombrar_especies(physeq_conocido_filtrado.adultos)
physeq_todos_filtrado.adultos    <- renombrar_especies(physeq_todos_filtrado.adultos)

# CONSTRUCCIÓN DE REDES CON SPIEC-EASI
red_conocido_adultos <- spiec.easi(
  physeq_conocido_filtrado.adultos,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

red_todos_adultos <- spiec.easi(
  physeq_todos_filtrado.adultos,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)
g_conocido.adultos <- adj2igraph(getRefit(red_conocido_adultos), 
                                  vertex.attr = list(name = taxa_names(physeq_conocido_filtrado.adultos)))

g_todos.adultos  <- adj2igraph(getRefit(red_todos_adultos),    
                               vertex.attr = list(name = taxa_names(physeq_todos_filtrado.adultos)))

metricas_conocido_adultos <- calcular_metricas(g_conocido.adultos)
metricas_todos_adultos    <- calcular_metricas(g_todos.adultos)

tabla_comparativa.adultos <- tibble(
  Métrica              = names(metricas_todos_adultos),
  Con_materia_obscura  = unlist(metricas_todos_adultos),
  Sin_materia_obscura  = unlist(metricas_conocido_adultos)
)
print(tabla_comparativa.adultos)

plot(g_todos.adultos,
     vertex.size = degree(g_todos.adultos)*2,
     vertex.color = cluster_louvain(g_todos.adultos)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red con Materia Obscura en adultos")

plot(g_conocido.adultos,
     vertex.size = degree(g_conocido.adultos)*2,
     vertex.color = cluster_louvain(g_conocido.adultos)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red sin Materia Obscura en infantes")


#PARA CARGARLAR EN CYTOSCAPE
library(RCy3)
cytoscapePing()

createNetworkFromIgraph(g_todos.adultos,   title = "Red con materia obscura adultos")
createNetworkFromIgraph(g_conocido.adultos, title = "Red sin materia obscura adultos")

# GRAFICAR COMPARACIÓN DE MÉTRICAS CON ggplot2
library(ggplot2)
library(tidyr)

df_metricas.adultos <- tabla_comparativa.adultos %>%
  pivot_longer(cols = -Métrica, names_to = "Red", values_to = "Valor")

ggplot(df_metricas.adultos, aes(x = Métrica, y = Valor, fill = Red)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Comparación de métricas de las redes",
       x = "Métrica",
       y = "Valor",
       fill = "Tipo de red") +
  coord_flip()
