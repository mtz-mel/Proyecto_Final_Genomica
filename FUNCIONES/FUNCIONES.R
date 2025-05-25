##Funciones###

## FUNCIÓN DE FILTRADO POR PREVALENCIA
prevalence_filter <- function(physeq_objeto, umbral = 0.1) {
  prevalencia <- apply(otu_table(physeq_objeto), 1, function(x) mean(x > 0))
  conservar <- names(prevalencia[prevalencia >= umbral])
  prune_taxa(conservar, physeq_objeto)
}

### Renombrar los taxa en physeq_especie (aplicable a physeq_todos y physeq_conocido)
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

### FUNCIÓN PARA CALCULAR MÉTRICAS DE LA RED
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

#
filtrar_prevalencia <- function(physeq_objeto, umbral = 0.1) {
  prevalencia <- apply(otu_table(physeq_objeto), 1, function(x) mean(x > 0))
  conservar <- names(prevalencia[prevalencia >= umbral])
  prune_taxa(conservar, physeq_objeto)
}