library(dada2)
library(phyloseq)
library(Biostrings)

remotes::install_github("jfq3/RDPutils")
library(RDPutils)

#Con el objeto physeq
load("DATOS/physeq_data.RData")

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