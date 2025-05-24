###BASES DE DATOS###
#Librerias----------------------------------------------------------------------
#install.packages("igraph")
library(igraph)
library(igraphdata)
library(networkD3)
library(RCy3)
library(tidyverse)
library(dplyr)

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
View(OTU)
TAX <- tax_table(tax) #
SAMPLES <- sample_data(sample_data)  #OBJETO PHYLOSEQ.

class(SAMPLES)


physeq <- phyloseq(OTU, TAX, SAMPLES)  #objeto phyloseq experiment-level object
#combinacion de todos.
View(sample_data(physeq))
View(tax_table(physeq))
View(otu_table(physeq))

#del combinado:
otu <- as.data.frame(otu_table(physeq)) #extraer la tabla de OTU
str(otu)

tax <- as.data.frame(tax_table(physeq)) #extraer la tabla TAX 
names(tax)
View(tax)

################################################################################
# Agregar la información taxonómica a la tabla de abundancias
otu$tax <- tax$Family
str(otu)
View(otu)

  ################################################################################
# Sumar las abundancias por familia

otu_sumas <- otu %>%
  dplyr::group_by(tax) %>%
  dplyr::summarise(across(where(is.numeric), sum), .groups = "drop")

View(otu_sumas) #nivel de Familia


################################################################################

#Aestos despues los podemos separar por mas, pero para no mezclar el tipo de muestra primero :)
#Base de datos con las muestras de heces:
metadatos <- sample_data(physeq) %>% as ("data.frame") %>%  rownames_to_column("sample_id")
# Sacar los id y separar:
heces_muestra <- metadatos %>% 
  filter(Sample_Type == "Feces") %>% 
  pull(sample_id) #extraer los id
#Abundancias de cada uno:
heces_abunancia<-otu_table(physeq)[,heces_muestra]

#Base de datos con las muestras de saliva:
# Sacar los id y separar:
saliva_muestras <- metadatos %>% 
  filter(Sample_Type == "Saliva") %>% 
  pull(sample_id) #extraer los id
#Abundancias de cada uno:
saliva_abunancia<-otu_table(physeq)[,saliva_muestras]
##########################


#holi alo

#puedes ver mis cambios ???
