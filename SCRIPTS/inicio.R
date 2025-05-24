###BASES DE DATOS###
#Librerias----------------------------------------------------------------------
#install.packages("igraph")
library(igraph)
library(igraphdata)
library(networkD3)
library(RCy3)

#install.packages("remotes")
#remotes::install_github("twbattaglia/MicrobeDS")
#BiocManager::install("microbiomeDataSets")
library(microbiomeDataSets)

#BiocManager::install("mia")
library(mia)
library(phyloseq)

availableDataSets()

SprockettTHData() #Los datos que vamos a utilizar.

SprockettTHData()->dataspo
class(dataspo)
#TreSummarizedExperiment.

#Hacer el objeto phyloseq: 
otu <- assays(dataspo)$counts #extraer los datos donde esta todo y el counts para los conteos.
otu <- as.matrix(otu) #matriz 
tax <- rowData(dataspo) #lso datos taoxnomicos.
tax <- as.matrix(as.data.frame(tax))  
sample_data <- colData(dataspo) #para sacar los metadatos
sample_data <- as.data.frame(sample_data) #para poder unirlos
#CREAR LOS OBJETOS OTU:
OTU <- otu_table(otu, taxa_are_rows = TRUE) #crear el objeto tipo OTU table
TAX <- tax_table(tax) #
SAMPLES <- sample_data(sample_data)  #OBJETO PHYLOSEQ.

class(SAMPLES)


physeq <- phyloseq(OTU, TAX, SAMPLES)  #objeto phyloseq experiment-level object
#combinacion de todos.
View(sample_data(physeq))

#del combinado:
otu <- as.data.frame(otu_table(physeq)) #extraer la tabla de OTU
tax <- as.data.frame(tax_table(physeq)) #extraer la tabla TAX 




################################################################################
# Agregar la información taxonómica a la tabla de abundancias
otu$tax <- tax$Genus
library(tidyverse)
# Sumar las abundancias por género
otu_summarized <- otu %>%
  group_by(tax) %>%
  summarise(across(where(is.numeric), sum))
View(otu_summarized) #nivel de genero

#a nivel de clase si hay numero mas
################################################################################