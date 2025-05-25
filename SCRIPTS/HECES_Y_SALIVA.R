library(dada2)
library(phyloseq)
library(Biostrings)
  #remotes::install_github("jfq3/RDPutils")
library(RDPutils)
library(igraph)
library(igraphdata)
library(networkD3)
library(RCy3)
library(tidyverse)  # Incluye dplyr
library(microbiomeDataSets)
library(mia)
library(SpiecEasi)

#-------------------------------------------------------------------------------
#Funciones y datos:
load("DATOS/physeq_data.RData")
load("DATOS/physeq_feces.RData")
load("DATOS/physeq_saliva.RData")

source("FUNCIONES/FUNCIONES.R")
#--------------------------------------------------------------------------------
#HECES:
physeq_spp_f <- tax_glom(physeq_feces, taxrank = "Species", NArm = FALSE) #agrupacion a nivel de especie
  #esta incluye las NA

physeq_sindm     <- subset_taxa(physeq_spp_f, !is.na(Species))  # sin materia oscura
physeq_conmd      <- physeq_spp_f  # con materia oscura

#FILTRADO DE PREVALENCIA:
#-------------------------------------------------------------------------------------------------------
prevalence_filter <- function(phy, threshold = 0.1) {
  prev <- apply(otu_table(phy), 1, function(x) mean(x > 0))
  keep <- names(prev[prev >= threshold])
  prune_taxa(keep, phy)
}
#--------------------------------------------------------------------------------------------------------


physeq_known_filt.feces <- prevalence_filter(physeq_sindm, 0.2) #aplciar la funcion y sacar nuevos objetos
physeq_all_filt.feces   <- prevalence_filter(physeq_conmd, 0.2)

#Renomabrar: si tienen NA como desconocido para formar las redes.
#Funcion: renombrar_especies:

# Aplicar renombramiento
physeq_known_filt.feces <- renombrar_especies(physeq_known_filt.feces)
physeq_all_filt.feces    <- renombrar_especies(physeq_all_filt.feces)


# REDES SPIEC-EASI
#analsiis de redes de co-ocurrencia 
feces_known <- spiec.easi( #NO DM
  physeq_known_filt.feces,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

feces_all <- spiec.easi(
  physeq_all_filt.feces,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

feces_known <- adj2igraph(getRefit(feces_known), 
                          vertex.attr = list(name = taxa_names(physeq_known_filt.feces)))
feces_all   <- adj2igraph(getRefit(feces_all),   
                          vertex.attr = list(name = taxa_names(physeq_all_filt.feces)))

#-----------------------------------------------------------------------------------------------------
metricas <- function(g) {
  list(
    nodos = vcount(g),
    aristas = ecount(g),
    grado_medio = mean(degree(g)),
    densidad = edge_density(g),
    clustering = transitivity(g, type = "global"),
    modularidad = modularity(cluster_louvain(g))
  )
}
#----------------------------------------------------------------------------------------------------------

feces.m_known <- metricas(feces_known)
feces.m_all   <- metricas(feces_all)

#COMPARACION:
tibble(
  Métrica        = names(feces.m_all),
  Con_materia_obscura = unlist(feces.m_all),
  Sin_materia_obscura = unlist(feces.m_known)
)


# VISUALIZACIÓN con materia obscura)
plot(feces_all,
     vertex.size = degree(feces_all)*2,
     vertex.color = cluster_louvain(feces_all)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red con Materia Obscura HECES")

# VISUALIZACIÓN DE g_known (sin materia obscura)
plot(feces_known,
     vertex.size = degree(feces_known)*2,
     vertex.color = cluster_louvain(feces_known)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red sin Materia Obscura HECES")


#PARA CARGARLAR EN CYTOSCAPE
library(RCy3)
cytoscapePing()

createNetworkFromIgraph(feces_all,   title = "Red con materia obscura solo heces")
createNetworkFromIgraph(feces_known, title = "Red sin materia obscura solo heces")

#-------------------------------------------------------------------------------
#SALIVA muestras:
physeq_spp_s <- tax_glom(physeq_saliva, taxrank = "Species", NArm = FALSE)

physeq_saliva.sin     <- subset_taxa(physeq_spp_s, !is.na(Species))  # sin materia oscura
physeq_saliva.com       <- physeq_spp_s  # con materia oscura

#FILTRADO DE PREVALENCIA:

physeq_known_filt.saliva <- filtrar_prevalencia(physeq_saliva.sin, 0.2) #aplciar la funcion y sacar nuevos objetos
physeq_all_filt.saliva   <- filtrar_prevalencia(physeq_saliva.com, 0.2)

# Aplicar renombramiento
physeq_known_filt.saliva <- renombrar_especies(physeq_known_filt.saliva)
physeq_all_filt.saliva    <- renombrar_especies(physeq_all_filt.saliva)

# REDES SPIEC-EASI
#analsiis de redes de co-ocurrencia 
saliva_known <- spiec.easi( #NO DM
  physeq_known_filt.saliva,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

saliva_all <- spiec.easi(
  physeq_all_filt.saliva,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

saliva_known <- adj2igraph(getRefit(saliva_known), 
                          vertex.attr = list(name = taxa_names(physeq_known_filt.saliva)))
saliva_all   <- adj2igraph(getRefit(saliva_all),   
                          vertex.attr = list(name = taxa_names(physeq_all_filt.saliva)))


saliva.m_known <- metricas(saliva_known)
saliva.m_all   <- metricas(saliva_all)

#COMPARACION:
tibble(
  Métrica        = names(saliva.m_all),
  Con_materia_obscura = unlist(saliva.m_all),
  Sin_materia_obscura = unlist(saliva.m_known)
)

# VISUALIZACIÓN con materia obscura)
plot(saliva_all,
     vertex.size = degree(saliva_all)*2,
     vertex.color = cluster_louvain(saliva_all)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red con Materia Obscura saliva")

# VISUALIZACIÓN DE g_known (sin materia obscura)
plot(saliva_known,
     vertex.size = degree(saliva_known)*2,
     vertex.color = cluster_louvain(saliva_known)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red sin Materia Obscura saliva")


createNetworkFromIgraph(saliva_all,   title = "Red con materia obscura solo saliva")
createNetworkFromIgraph(saliva_known, title = "Red sin materia obscura solo saliva")