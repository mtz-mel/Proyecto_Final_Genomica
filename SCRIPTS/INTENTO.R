library(dada2)
library(phyloseq)
library(Biostrings)

remotes::install_github("jfq3/RDPutils")
library(RDPutils)

#Con el objeto physeq
load("DATOS/physeq_data.RData")

library(igraph)
library(igraphdata)
library(networkD3)
library(RCy3)
library(tidyverse)  # Incluye dplyr
library(microbiomeDataSets)
library(mia)
library(phyloseq)
library(SpiecEasi)
#--------------------------------------------------------------------------------
#Encontrar la prevalencia de cada taxa por cada muestra.
#Une al informacion taxonomicay al informacion de prevalencia de cada taxa en un data frame.
prevdf <- apply(X = otu_table(physeq),MARGIN = ifelse(taxa_are_rows(physeq), yes=1, no=2), 
               FUN=function(x){sum(x>0)}) 

prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(physeq), tax_table(physeq))
prevdf <- prevdf[prevdf$Prevalence > 0,] #remove taxa not present in samples

#Ahora, utilizando el marco de datos prevdf resultante como entrada, 
#se puede crear un gráfico de líneas de la prevalencia de taxones conocidos vs desconocidos a
#nivel de género usando la función get_back_counts_for_line_plotf.
  #prev_lineplot <- get_back_counts_for_line_plot(prevdf, "AMBIENTE")

phylo_for_net <- get_back_res_meeting_min_occ(physeq, filter_val_percent=0.4)
#resulting data will be list of 2- phyloseq, graph
#phyloseq with only taxa meeting sample prevalence kept
#graph created using SpiecEasi neighborhood algorithm on phyloseq

################################################################################
physeq_spp_f <- tax_glom(physeq_feces, taxrank = "Species", NArm = FALSE)

physeq_sindm     <- subset_taxa(physeq_spp_f, !is.na(Species))  # sin materia oscura
physeq_sindm       <- physeq_spp_f  # con materia oscura

#FILTRADO DE PREVALENCIA:

prevalence.filter.feces <- function(phy_1, threshold = 0.2) { 
  prev <- apply(otu_table(phy_1), 1, function(x) mean(x > 0)) #calcular prevalencia de cda taxon
  keep <- names(prev[prev >= threshold]) #seleccioanr los que son mayores
  prune_taxa(keep, phy_1) #filtrado
}

physeq_known_filt.feces <- prevalence_filter(physeq_sindm, 0.2) #aplciar la funcion y sacar nuevos objetos
physeq_all_filt.feces   <- prevalence_filter(physeq_sindm, 0.2)

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
createNetworkFromIgraph(feces_all,   title = "Red con materia obscura solo heces")
createNetworkFromIgraph(feces_known, title = "Red sin materia obscura solo heces")