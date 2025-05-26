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
OTU <- otu_table(matriz_otu, taxa_are_rows = TRUE) #se crean tablas apartir de las matrices de los datos
TAX <- tax_table(matriz_tax)
MUESTRAS <- sample_data(datos_muestras)
physeq <- phyloseq(OTU, TAX, MUESTRAS) #se crea un objeto phyloseq que une los OTUS, TAXAS y las muestras

# AGRUPACIÓN A NIVEL DE ESPECIE (INCLUYENDO NA)
physeq_especie <- tax_glom(physeq, taxrank = "Species", NArm = FALSE) # tax_glom agruoa las OTUS a nivel de especie, taxrank especifica el nivel taxonomico que queremos agurpar,
                                                                      #NArm=FALSE indica que queresmos quedarnos con los datos aunque no tengan asignacion a nivel de especie

# DIVISIÓN SEGÚN PRESENCIA DE MATERIA OSCURA
physeq_conocido <- subset_taxa(physeq_especie, !is.na(Species))  # sin materia obscura
physeq_todos <- physeq_especie  # con materia obscura

# FUNCIÓN DE FILTRADO POR PREVALENCIA
#se realizó una funcion donde para cada otu calcular las muestras donde haya una abundancia mayor a 0 y se conservan unicamente los nombres de los OTUS que tengan una prevalencia mayor a 10%
#y se crea un objeto phyloseq que contiene los otus antes seleccionados

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

#se usa el metodo meinshausen-Buhlmann neighborhood selection
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
    betweenness_media = mean(betweenness(grafo))
  )
}

metricas_conocido <- calcular_metricas(g_conocido)
metricas_todos    <- calcular_metricas(g_todos)

# COMPARACIÓN DE MÉTRICAS
#se hace una tabla que contenga las metricas tanto de la red con y sin meteria obscura
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

#------------------------------------------------------------------------------#

# REDES CON BOOSTRAP

# se crea una función para eliminar un porcentaje aleatorio de taxones
# extraemos los nombres de taxones
# seleccionar 30% de taxones al azar
# filtrar el phyloseq sin esos taxones
eliminar_taxones_azar <- function(physeq, porcentaje = 0.3) {
  taxones <- taxa_names(physeq)
  eliminar <- sample(taxones, size = round(length(taxones) * porcentaje))
  prune_taxa(setdiff(taxones, eliminar), physeq)
}


#una vez Hecha la funcion la aplicamos para eliminar 30% de taxones 

physeq_reducido <- eliminar_taxones_azar(physeq_conocido_filtrado)


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

# se alcular métricas de la nueva red con la funcion antes creada
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

#------------------------------------------------------------------------------#

# ANÁLISIS ESTADÍSTICOS PARA DETERMINAR SI LOS CAMBIOS SON SIGNIFICATIVOS 
#se realiza la prueba de xilcoxon, se utiliza mapply para aplicarla a cada una de las metricas 

wilcoxon_com <- mapply(function(x, y) wilcox.test(x, y, paired = TRUE), 
                       metricas_conocido, metricas_todos, SIMPLIFY = FALSE)

# Extraer valores p
#se utiliza sapply para que lo de en un vector numerico
P_values_wilcoxon <- sapply(wilcoxon_com, function(x) x$p.value)

# Mostrar resultados en una tabla con los nombres de las metricas y los valores de p
tibble(
  Métrica = names(metricas_todos),
  P_value_Wilcoxon = P_values_wilcoxon
) %>% arrange(P_value_Wilcoxon)



#------------------------------------------------------------------------------#

# BOXPLOT DE LAS REDES

library(ggplot2)
library(tidyr)


# primera creaamos  una tabla comparativa con los datos de las tres redes (metricas)
tabla_comparativa_especie <- tibble(
  Métrica = names(metricas_all),
  Con_materia_obscura = unlist(m_all),
  Sin_materia_obscura = unlist(m_known),
  Bootstrap = unlist(m_reducido)  
)

tabla_comparativa_especie

#boxplot de las redes 

# convertir en  formato largo para poder trabajar con ggplot 
#esto se hace con la funcion pivot_longer, cols hace que la columna metrica sea el identificador
#names_to combina los nombres de columnas en una variable categorica
#values_to combina todo en una sola columna
df_metricas_especie <- tabla_comparativa_especie %>%
  pivot_longer(cols = -Métrica, names_to = "Red", values_to = "Valor")

# hacemos los boxplot utilizando ggplot, el color de relleno es deacuerdo al tipo de red.  

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


#------------------------------------------------------------------------------#

# MEJORAR G´RAFICA DE BARRAS


library(ggplot2)
library(tidyr)


# se crea unerar una grafica de barras con ggplot 2 que muuestra las metricas de las 3 redes
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



