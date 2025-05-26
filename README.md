# **Proyecto Final Genómica Funcional: Materia oscura bacteriana y su impacto en la arquitectura de redes ecológicas.**   

**Integrantes:**    
Melissa Martínez.  
Isabel Herrera.  
Alondra Dominguez.  

**Descripción del proyecto:**      
La **materia oscura microbiana** consiste en aquellos taxones que no han podido ser identificados, lo cual dificulta la caracterización del microbioma en diversos ambientes, especialmente en aquellos que se encuentran bajo condiciones extremas. En este proyecto simulamos cómo la eliminación de esa materia oscura afecta la arquitectura de la red ecológica, comparando sus métricas de centralidad. Para ello trabajamos con la base de datos 
**SprockettTHData()** de la librería **microbiomeDataSets** de Bioconductor, considerando como “materia oscura” los taxones identificados únicamente hasta género, y como “materia conocida” los que cuentan con asignación de especie. Esta distinción nos permitió generar redes con una densidad que nos permitiera visualizar con claridad los cambios estructurales.  

**Objetivos:**  

- Simular cómo la eliminación de materia oscura afecta la estructura de redes ecológicas  
- Comparar métricas de centralidad entre redes con y sin materia oscura  

**Librerías utilizadas:**   
- library(igraph)    
- library(igraphdata)    
- library(networkD3)  
- library(RCy3)  
- library(tidyverse)   
- library(microbiomeDataSets)  
- library(mia)  
- library(phyloseq)  
- library(SpiecEasi)    

**Programas utilizados:**  
- R 4.4.1    
- RStudio    
- Cytoscape 3.10.3    

**Referencias:**   
- Sprockett, D. D., Martin, M., Costello, E. K., Burns, A. R., Holmes, S. P., Gurven, M. D., & Relman, D. A. (2020). Microbiota assembly, structure, and dynamics among Tsimane horticulturalists of the Bolivian Amazon. Nature Communications, 11(1). https://doi.org/10.1038/s41467-020-17541-6  
- Zamkovaya, T., Foster, J. S., De Crécy-Lagard, V., & Conesa, A. (2020). A network approach to elucidate and prioritize microbial dark matter in microbial communities. The ISME Journal, 15(1), 228–244. https://doi.org/10.1038/s41396-020-00777-x
- Del Campo-Moreno, R., Alarcón-Cavero, T., D’Auria, G., Delgado-Palacio, S., & Ferrer-Martínez, M. (2018). Microbiota and Human Health: Characterization techniques and transference. Enfermedades Infecciosas y Microbiologia Clinica (English Ed ), 36(4), 241-245. https://doi.org/10.1016/j.eimce.2018.02.016
