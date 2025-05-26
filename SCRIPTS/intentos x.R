######
#meee   #INTENTOS ISA

################

# Agrupar OTUs a nivel de género
physeq_genus <- tax_glom(physeq, taxrank = "Genus", NArm = FALSE)

# Agrupar OTUs a nivel de especie
physeq_species <- tax_glom(physeq, taxrank = "Species", NArm = FALSE)



# Función de filtrado por prevalencia
prevalence_filter <- function(phy, threshold = 0.1) {
  prev <- apply(otu_table(phy), 1, function(x) mean(x > 0))
  keep <- names(prev[prev >= threshold])
  prune_taxa(keep, phy)
}

# Filtrar géneros y especies presentes en ≥ 20% de las muestras
physeq_genus_filt <- prevalence_filter(physeq_genus, threshold = 0.2)
physeq_species_filt <- prevalence_filter(physeq_species, threshold = 0.2)


library(SpiecEasi)

# Inferir red de géneros
se_genus <- spiec.easi(
  physeq_genus_filt,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

# Inferir red de especies
se_species <- spiec.easi(
  physeq_species_filt,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)


library(igraph)

# Red de géneros
g_genus_se <- adj2igraph(getRefit(se_genus), vertex.attr = list(name = taxa_names(physeq_genus_filt)))

# Red de especies
g_species_se <- adj2igraph(getRefit(se_species), vertex.attr = list(name = taxa_names(physeq_species_filt)))


# Número de nodos y conexiones
nodos_genus  <- vcount(g_genus_se)
aristas_genus <- ecount(g_genus_se)

nodos_species  <- vcount(g_species_se)
aristas_species <- ecount(g_species_se)

# Grado medio
grado_medio_genus <- mean(degree(g_genus_se))
grado_medio_species <- mean(degree(g_species_se))

# Densidad de la red
densidad_genus <- graph.density(g_genus_se)
densidad_species <- graph.density(g_species_se)

# Coeficiente de clustering global
clustering_genus <- transitivity(g_genus_se, type = "global")
clustering_species <- transitivity(g_species_se, type = "global")

# Distribución de grado
hist(degree(g_genus_se), main="Distribución de grado (Géneros)", xlab="Grado", ylab="Frecuencia")
hist(degree(g_species_se), main="Distribución de grado (Especies)", xlab="Grado", ylab="Frecuencia")

# Centralidades (ejemplo: intermediación “betweenness”)
betweenness_genus <- betweenness(g_genus_se)
betweenness_species <- betweenness(g_species_se)

# Detectar comunidades con el método de Louvain
coms_genus <- cluster_louvain(g_genus_se)
coms_species <- cluster_louvain(g_species_se)

# Añadir comunidades al grafo
V(g_genus_se)$community <- membership(coms_genus)
V(g_species_se)$community <- membership(coms_species)

# Calcular modularidad
modularity_genus <- modularity(coms_genus)
modularity_species <- modularity(coms_species)



E(g_genus_se)$weight <- E(g_genus_se)$weight
E(g_species_se)$weight <- E(g_species_se)$weight

g_genus_strong <- delete_edges(g_genus_se, E(g_genus_se)[weight < 0.2])
g_species_strong <- delete_edges(g_species_se, E(g_species_se)[weight < 0.2])



windows()
plot(g_genus_strong,
     vertex.size = degree(g_genus_strong)*2,
     vertex.color = V(g_genus_strong)$community,
     vertex.label.cex = 0.7,
     edge.width = E(g_genus_strong)$weight * 2,
     main = "Red de Géneros (SPIEC-EASI)")

plot(g_species_strong,
     vertex.size = degree(g_species_strong)*2,
     vertex.color = V(g_species_strong)$community,
     vertex.label.cex = 0.7,
     edge.width = E(g_species_strong)$weight * 2,
     main = "Red de Especies (SPIEC-EASI)")



#####################################

# Convertir tax_table a data frame
tax <- as.data.frame(tax_table(physeq))

# Identificar taxones no determinados en género y especie
tax$Identificado_Genero <- ifelse(is.na(tax$Genus) | tax$Genus == "", "No Identificado", "Identificado")
tax$Identificado_Especie <- ifelse(is.na(tax$Species) | tax$Species == "", "No Identificado", "Identificado")

otu <- as.data.frame(otu_table(physeq))
otu$Identificado_Genero <- tax$Identificado_Genero[rownames(otu)]
otu$Identificado_Especie <- tax$Identificado_Especie[rownames(otu)]


# Agrupar OTUs a nivel de género
physeq_genus <- tax_glom(physeq, taxrank = "Genus", NArm = FALSE)

# Agrupar OTUs a nivel de especie
physeq_species <- tax_glom(physeq, taxrank = "Species", NArm = FALSE)


# Función de filtrado por prevalencia
prevalence_filter <- function(phy, threshold = 0.1) {
  prev <- apply(otu_table(phy), 1, function(x) mean(x > 0))
  keep <- names(prev[prev >= threshold])
  prune_taxa(keep, phy)
}

# Filtrar géneros y especies presentes en ≥ 20% de las muestras
physeq_genus_filt <- prevalence_filter(physeq_genus, threshold = 0.2)
physeq_species_filt <- prevalence_filter(physeq_species, threshold = 0.2)



# Inferir red de géneros
se_genus <- spiec.easi(
  physeq_genus_filt,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

# Inferir red de especies
se_species <- spiec.easi(
  physeq_species_filt,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)



# Red de géneros
g_genus_se <- adj2igraph(getRefit(se_genus), vertex.attr = list(name = taxa_names(physeq_genus_filt)))

# Red de especies
g_species_se <- adj2igraph(getRefit(se_species), vertex.attr = list(name = taxa_names(physeq_species_filt)))



# Número de nodos y conexiones
nodos_genus  <- vcount(g_genus_se)
aristas_genus <- ecount(g_genus_se)

nodos_species  <- vcount(g_species_se)
aristas_species <- ecount(g_species_se)

# Grado medio
grado_medio_genus <- mean(degree(g_genus_se))
grado_medio_species <- mean(degree(g_species_se))

# Densidad de la red
densidad_genus <- graph.density(g_genus_se)
densidad_species <- graph.density(g_species_se)

# Coeficiente de clustering global
clustering_genus <- transitivity(g_genus_se, type = "global")
clustering_species <- transitivity(g_species_se, type = "global")

# Distribución de grado
hist(degree(g_genus_se), main="Distribución de grado (Géneros)", xlab="Grado", ylab="Frecuencia")
hist(degree(g_species_se), main="Distribución de grado (Especies)", xlab="Grado", ylab="Frecuencia")

# Centralidades (ejemplo: intermediación “betweenness”)
betweenness_genus <- betweenness(g_genus_se)
betweenness_species <- betweenness(g_species_se)

# Detectar comunidades con el método de Louvain
coms_genus <- cluster_louvain(g_genus_se)
coms_species <- cluster_louvain(g_species_se)

# Añadir comunidades al grafo
V(g_genus_se)$community <- membership(coms_genus)
V(g_species_se)$community <- membership(coms_species)

# Calcular modularidad
modularity_genus <- modularity(coms_genus)
modularity_species <- modularity(coms_species)


E(g_genus_se)$weight <- E(g_genus_se)$weight
E(g_species_se)$weight <- E(g_species_se)$weight

g_genus_strong <- delete_edges(g_genus_se, E(g_genus_se)[weight < 0.2])
g_species_strong <- delete_edges(g_species_se, E(g_species_se)[weight < 0.2])


windows()
plot(g_genus_strong,
     vertex.size = degree(g_genus_strong)*2,
     vertex.color = V(g_genus_strong)$community,
     vertex.label.cex = 0.7,
     edge.width = E(g_genus_strong)$weight * 2,
     main = "Red de Géneros (SPIEC-EASI)")

plot(g_species_strong,
     vertex.size = degree(g_species_strong)*2,
     vertex.color = V(g_species_strong)$community,
     vertex.label.cex = 0.7,
     edge.width = E(g_species_strong)$weight * 2,
     main = "Red de Especies (SPIEC-EASI)")




# Ver nombres únicos en la clasificación de género
unique(tax_table(physeq)[, "Genus"])

# Ver nombres únicos en la clasificación de especie
unique(tax_table(physeq)[, "Species"])


table(tax$Identificado_Genero)  # Debe mostrar cuántos son identificados vs no identificados
table(tax$Identificado_Especie)


all(rownames(otu) %in% rownames(tax))  # Debe devolver TRUE si están alineados correctamente



V(g_species_se)$name
V(g_genus_se)$name


# Filtrar organismos no identificados en la red de género
no_identificados_genus <- V(g_genus_se)$name[grepl("Unclassified|Unknown|NA", V(g_genus_se)$name)]

# Filtrar organismos no identificados en la red de especie
no_identificados_species <- V(g_species_se)$name[grepl("Unclassified|Unknown|NA", V(g_species_se)$name)]

# Mostrar resultados
no_identificados_genus
no_identificados_species



unique(tax$Genus)
unique(tax$Species)


rownames(otu)[tax$Identificado_Genero == "No Identificado"]
rownames(otu)[tax$Identificado_Especie == "No Identificado"]


# Incluir taxones no identificados antes del filtro de prevalencia
physeq_genus_filt <- prune_taxa(rownames(tax)[tax$Identificado_Genero == "Identificado" | tax$Identificado_Genero == "No Identificado"], physeq_genus)
physeq_species_filt <- prune_taxa(rownames(tax)[tax$Identificado_Especie == "Identificado" | tax$Identificado_Especie == "No Identificado"], physeq_species)


####no sirvio a ver otra vez
##################


# Convertir tax_table a data frame
tax <- as.data.frame(tax_table(physeq))

# Asignar "Unclassified_Species" a especies sin identificar
tax$Species <- ifelse(is.na(tax$Species) | tax$Species == "", "Unclassified_Species", tax$Species)

# Actualizar tax_table en physeq
tax_table(physeq) <- as.matrix(tax)


rownames(otu_table(physeq))[tax$Species == "Unclassified_Species"]

all(rownames(otu_table(physeq_fam_filt)) %in% V(g_fam_se)$name)
V(g_fam_se)$name[grepl("Unclassified_Species", V(g_fam_se)$name)]

tax_no_identificados <- rownames(tax)[tax$Species == "Unclassified_Species"]
physeq_fam_filt <- prune_taxa(c(tax_no_identificados, taxa_names(physeq_fam_filt)), physeq_fam)


# Función de filtrado por prevalencia que **no excluye** especies no identificadas
prevalence_filter <- function(phy, threshold = 0.1) {
  prev <- apply(otu_table(phy), 1, function(x) mean(x > 0))
  keep <- rownames(tax)[prev >= threshold | tax$Species == "Unclassified_Species"]  # Incluye especies sin clasificar
  prune_taxa(keep, phy)
}

# Aplicar filtrado sin excluir especies no clasificadas
physeq_fam_filt <- prevalence_filter(physeq_fam, threshold = 0.2)




se_fam <- spiec.easi(
  physeq_fam_filt,
  method           = "mb",
  lambda.min.ratio = 1e-2,
  nlambda          = 20,
  sel.criterion    = "bstars",
  pulsar.params    = list(thresh = 0.05)
)





g_fam_se <- adj2igraph(
  getRefit(se_fam),
  vertex.attr = list(name = taxa_names(physeq_fam_filt))
)


V(g_fam_se)$name[grepl("Unclassified_Species", V(g_fam_se)$name)]















##################


any(tax$Species == "Unclassified_Species")  # Debería ser TRUE

sum(taxa_names(physeq_fam_filt) %in% tax_no_identificados)  # ¿Cuántos no clasificados quedaron tras el filtrado?


# Nombres de todos los OTUs tras SpiecEasi
all_taxa <- taxa_names(physeq_fam_filt)

# OTUs que sí están en la red
present_in_graph <- V(g_fam_se)$name

# OTUs que faltan (nodos sin conexiones)
missing_nodes <- setdiff(all_taxa, present_in_graph)

# Añadir nodos aislados
for (node in missing_nodes) {
  g_fam_se <- igraph::add_vertices(g_fam_se, 1, name = node)
}

setdiff(tax_no_identificados, V(g_fam_se)$name)

V(g_fam_se)$name[grepl("Unclassified_Species", V(g_fam_se)$name)]

