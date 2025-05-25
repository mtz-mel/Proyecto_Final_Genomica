
# Juntarw todas las OTUs al nivel de Family, preservando NAs
physeq_fam <- tax_glom(physeq, taxrank = "Species", NArm = FALSE)

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
  lambda.min.ratio = 1e-1,            # rango de penalización
  nlambda          = 20,              # pasos de lambda
  sel.criterion    = "bstars",        # criterio de selección de modelo
  pulsar.params    = list(thresh = 0.1)  # opcional, ajuste de estabilidad
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

plot(g_fam_strong,
     vertex.size = degree(g_fam_strong)*2,
     vertex.color = V(g_fam_strong)$community,
     vertex.label.cex = 0.7,
     edge.width = E(g_fam_strong)$weight * 2,
     main = "Red de Especies (SPIEC-EASI)")

library(RCy3)
createNetworkFromIgraph(g_fam_strong, title = "Especies SPIEC-EASI")



###############################################################################
### Creo que asi deberia de quedar 

### Esto es tomando en cuenta lo que ya habia hech Alo, pero cambiando cosas para 
## poder hacer las dos redes a la vez y lo de SPIEC-EASI

# CARGA DE LIBRERÍAS NECESARIAS
library(igraph)
library(igraphdata)
library(networkD3)
library(RCy3)
library(tidyverse)  # Incluye dplyr
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

physeq_known     <- subset_taxa(physeq_spp, !is.na(Species))  # sin materia oscura
physeq_all       <- physeq_spp  # con materia oscura


# BLOQUE NUEVO: FILTRADO DE PREVALENCIA

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
    betweenness_media = mean(betweenness(g))  # Promedio de betweenness
  )
}


m_known <- metricas(g_known)
m_all   <- metricas(g_all)

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





#################

#analisis
#creamos listas para las metricas de cada red 

metricas_known <- unlist (m_known)

metricas_all <- unlist (m_all)

#hacemos la prueba de wilcoxon para ver si es significativo el cambio 

wilcoxon_com<- mapply (function (x,y) wilcox.test (x, y, paired =TRUE ), metricas_known, metricas_all, SIMPLIFY = FALSE)

#vemos los valores de P

P_values<- sapply( wilcoxon_com, function (x) x$p.value)

#ahora mostramos los resultados nen forma de tabla 

tibble( metrica= names(m_all), P_value = P_values) %>% arrange(P_value)




# 

hist(metricas_known, main = "Métricas sin materia oscura")
hist(metricas_all, main = "Métricas con materia oscura")

################

# analisis tstudent
library(stats)

# Clistas con las metricas
metricas_known <- unlist(m_known)
metricas_all   <- unlist(m_all)

length(metricas_known)
length(metricas_all)


var(metricas_known)
var(metricas_all)

# Aplicar prueba t de Student considerando varianzas diferentes
t_test_resultado <- t.test(metricas_known, metricas_all, paired = TRUE, var.equal = FALSE)

# Mostrar resultados
tibble(
  Métrica = names(m_all),
  P_value = t_test_resultado$p.value
) %>% arrange(P_value)
###########



###################

#redes con boostrap

# creamos una función para eliminar un porcentaje aleatorio de taxones
eliminar_taxones_azar <- function(physeq, porcentaje = 0.3) {
  # extraemos los nombres de taxones
  taxones <- taxa_names(physeq)
  
  # Seleccionar 30% de taxones al azar #es un valor arbitrario no se cual ponerle sjsj
  eliminar <- sample(taxones, size = round(length(taxones) * porcentaje))
  
  # filtrar el phyloseq sin esos taxones
  prune_taxa(setdiff(taxones, eliminar), physeq)
}


#una vez echa la funcion la aplicamos para eliminar 30% de taxones

physeq_reducido <- eliminar_taxones_azar(physeq_all)

# filtrado de prevalencia en la nueva red
physeq_reducido_filt <- prevalence_filter(physeq_reducido, 0.1)

# hacemos la red con el mismo SPIEC-EASI pero con el phyloseq reducido
se_reducido <- spiec.easi(
  physeq_reducido_filt,
  method = "mb",
  lambda.min.ratio = 1e-1,
  nlambda = 20,
  sel.criterion = "bstars",
  pulsar.params = list(thresh = 0.1)
)

# convertimos a objeto igraph
g_reducido <- adj2igraph(getRefit(se_reducido), vertex.attr = list(name = taxa_names(physeq_reducido_filt)))

# se alcular métricas de la nueva red
m_reducido <- metricas(g_reducido)

# Visualización de la red reducida
plot(g_reducido,
     vertex.size = degree(g_reducido)*2,
     vertex.color = cluster_louvain(g_reducido)$membership,
     vertex.label.cex = 0.7,
     edge.width = 1,
     main = "Red con 30% de taxones eliminados")

# Crear red en Cytoscape pa verla más bonita 
createNetworkFromIgraph(g_reducido, title = "Red con 30% de taxones eliminados")


############################


#boxplot de las redes 

library(ggplot2)
library(tidyr)

# primera creaamos  una tabla comparativa con los datos de las tres redes
tabla_comparativa_especie <- tibble(
  Métrica = names(metricas_todos),
  Con_materia_obscura = unlist(metricas_todos),
  Sin_materia_obscura = unlist(metricas_conocido),
  Bootstrap = unlist(m_reducido)  
)

# convertir en  formato largo para ggplot
df_metricas_especie <- tabla_comparativa_especie %>%
  pivot_longer(cols = -Métrica, names_to = "Red", values_to = "Valor")

# hacemos los boxplot
ggplot(df_metricas_especie, aes(x = Red, y = Valor, fill = Red)) +
  geom_boxplot() +
  facet_wrap(~ Métrica, scales = "free") +  # para separar cada metrica en su propio box
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(face = "bold", size = 16),
    legend.position = "top",
     ) +
  labs(title = "Comparación de métricas entre redes a nivel de especie y bootstrap",
       x = "Tipo de Red",
       y = "Valor de Métrica") 


#################

#mejor grafica de barras 


library(ggplot2)
library(tidyr)



# Generar la gráfica de barras
ggplot(df_metricas_especie, aes(x = Métrica, y = Valor, fill = Red)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +  #
  scale_fill_manual(values = c("Con_materia_obscura" = "lightblue", 
                               "Sin_materia_obscura" = "pink", 
                               "Bootstrap" = "purple")) +  # colores para cada red
  theme_minimal() +
  coord_flip()+
  labs(title = "Comparación de métricas entre redes",
       x = "Métrica",
       y = "Valor de métrica") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
