
# Juntarw todas las OTUs al nivel de Family, preservando NAs
physeq_fam <- tax_glom(physeq, taxrank = "Family", NArm = FALSE)

# Función de filtrado por prevalencia
prevalence_filter <- function(phy, threshold = 0.1) {
  prev <- apply(otu_table(phy), 1, function(x) mean(x > 0))
  keep <- names(prev[prev >= threshold])
  prune_taxa(keep, phy)
}

# Aplica, por ejemplo, filtrado: familias presentes en ≥ 20% de las muestras
physeq_fam_filt <- prevalence_filter(physeq_fam, threshold = 0.2)

library(SpiecEasi)

se_fam <- spiec.easi(
  physeq_fam_filt,
  method           = "mb",            # neighborhood selection (Meinshausen–Bühlmann)
  lambda.min.ratio = 1e-2,            # rango de penalización
  nlambda          = 20,              # pasos de lambda
  sel.criterion    = "bstars",        # criterio de selección de modelo
  pulsar.params    = list(thresh = 0.05)  # opcional, ajuste de estabilidad
)

# Extraer la red como objeto igraph
g_fam_se <- adj2igraph(
  getRefit(se_fam),
  vertex.attr = list(name = taxa_names(physeq_fam_filt))
)

# phyloseq listo para SPIEC-EASI
physeq_fam_filt


library(igraph)

# Número de nodos y aristas
nodos  <- vcount(g_fam_se)
aristas <- ecount(g_fam_se)

# Grado medio
grado_medio <- mean(degree(g_fam_se))

# Densidad
densidad <- graph.density(g_fam_se)

# Coeficiente de clustering global
clustering_global <- transitivity(g_fam_se, type = "global")

# Distribución de grado (puedes plotear un histograma)
grado   <- degree(g_fam_se)
hist(grado, main="Distribución de grado", xlab="Grado", ylab="Frecuencia")

# Centralidades (ejemplo: intermediación “betweenness”)
betweenness_vals <- betweenness(g_fam_se)
top10_betw       <- sort(betweenness_vals, decreasing=TRUE)[1:10]

# Detectar comunidades con el método de Louvain
coms <- cluster_louvain(g_fam_se)

# Añadir al grafo como atributo
V(g_fam_se)$community <- membership(coms)

# Calcular modularidad
modularity_score <- modularity(coms)

# Filtrar aristas débiles (por ejemplo, peso < 0.2)
E(g_fam_se)$weight <- E(g_fam_se)$weight  # ya viene con pesos de SpiecEasi

g_fam_strong <- delete_edges(g_fam_se, E(g_fam_se)[weight < 0.2])


# Color por comunidad y tamaño de nodo según grado
windows()
plot(g_fam_strong,
     vertex.size = degree(g_fam_strong)*2,
     vertex.color = V(g_fam_strong)$community,
     vertex.label.cex = 0.7,
     edge.width = E(g_fam_strong)$weight * 2,
     main = "Red de Familias (SPIEC-EASI)")


