# **Proyecto Final Genómica Funcional: Materia oscura bacteriana y su impacto en la arquitectura de redes ecológicas.**   

## **Integrantes:**    
Melissa Martínez, Isabel Herrera y Alondra Dominguez.  

## **Descripción del proyecto:**      
La **materia oscura microbiana** consiste en aquellos taxones que no han podido ser identificados, lo cual dificulta la caracterización del microbioma en diversos ambientes, especialmente en aquellos que se encuentran bajo condiciones extremas. En este proyecto simulamos cómo la eliminación de esa materia oscura afecta la arquitectura de la red ecológica, comparando sus métricas de centralidad. Para ello trabajamos con la base de datos 
**SprockettTHData()** de la librería **microbiomeDataSets** de Bioconductor, considerando como “materia oscura” los taxones identificados únicamente hasta género, y como “materia conocida” los que cuentan con asignación de especie. Esta distinción nos permitió generar redes con una densidad que nos permitiera visualizar con claridad los cambios estructurales. Para la construcción de redes se utilizó el paquete SpiecEasi, el cual emplea el método SPIEC-EASI, el cual fue diseñado para inferir redes ecológicas microbianas a partir de datos de secuenciación. Fue propuesto por Kurtz et al. (2015) para superar limitaciones de métodos anteriores que aplicaban correlaciones simples.  

## **Objetivos:**  

- Simular cómo la eliminación de materia oscura afecta la estructura de redes ecológicas  
- Comparar métricas de centralidad entre redes con y sin materia oscura  

## **Librerías utilizadas:**   
- library(igraph)    
- library(igraphdata)    
- library(networkD3)  
- library(RCy3)  
- library(tidyverse)   
- library(microbiomeDataSets)  
- library(mia)  
- library(phyloseq)  
- library(SpiecEasi)    

## **Programas utilizados:**  
- R 4.4.1    
- RStudio    
- Cytoscape 3.10.3

## **Flujo de trabajo:**  
1.- Carpeta "Datos": Contiene los objetos en formato RData de el filtrado que se le hizo al objeto physeq que contiene todos los metadatos y muestras.  
2.- Carpeta "Funciones": Contiene las funciones utilizadas durante el trabajo, primero se tiene que cargar este archivo para poder correr los scripts.  
3.- Carpeta "Scripts":  
-El primer script es el de PAQUETES, contiene todos los paquetes utilizados en los analisis.
-El segundo script es el manejo de los datos, donde se crea el objeto phyloseq y se filtra en base a los criterios de edad y tipo de muestra.  
-Los script de POR EDADES, REDES_TODO y HECES_Y_SALIVA contienen los analisis realizados a cada objeto phyloseq que se creo.  

## **Interpretación de resultados:**  
El estudio de las comunidades microbianas desconocidas es relevante en el contexto de la biología. A pesar de los múltiples estudios y la aplicación de tecnologías de secuenciación, solo se conoce alrededor del 2% de las especies microbianas presentes en ambientes con condiciones extremas (Zha et al., 2022). El estudio de estas comunidades microbianas, enfocado en el bacterioma, es una herramienta que puede servir para entender la importancia de estos microorganismos en una comunidad y sus funciones ecológicas (Zha et al., 2022).   

Para el estudio del bacterioma con un enfoque en taxones desconocidos, recreamos la metodología usada en Zamkovaya et al., 2020. En este estudio y en otros estudios, el enfoque se centra en ambientes extremos donde se desconoce la mayoría de los taxones. En el caso de la microbiota humana, el enfoque se centra más en estudios médicos y se conoce qué taxones están asociados a condiciones de enfermedad y cuáles a un estado de salud (Del Campo-Moreno et al., 2018).   

En nuestro proyecto observamos que las redes cambian su estructura si se toman en cuenta los taxones no identificados y si no se toman en cuenta. Algunos de estos taxones no identificados son "hubs" y están conectados a taxones como *Holdemanella biformis, Oxalobacter formigenes o Ruminococcus bromii*, las cuales son especies con funciones en el mantenimiento de la homeostasis en el intestino. Al quitar los taxones desconocidos, observamos que la red se fragmenta.  

La red creada con materia oscura tiene presencia de múltiples taxones caracterizados de color morado, los cuales representan nodos centrales con una importancia ecológica significativa, ya que en comparación con la red creada sin la materia oscura, esta última se observa con una densidad menor, al igual que una disminución en el número de nodos con un alto grado de conexiones, disminuyendo así la interconectividad de la comunidad. Esto se ve reflejado en la disminución de las conexiones entre taxones, así como del betweenness y el diámetro, lo cual se ve reflejado en la estructura de la red. Esto señala que la materia oscura microbiana es parte de la estabilidad de la comunidad. Sin embargo, a pesar de que estas diferencias son observables, no se encontraron diferencias significativas entre el degree, diámetro, betweenness y clustering de las redes cuando se realizaron las pruebas estadísticas. Esto sugiere que, debido a que la microbiota es clave para el desarrollo del sistema inmunitario tanto para que un individuo conserve la homeostasis (Álvarez et al., 2021), los taxones importantes ya han sido clasificados, por lo que aunque la materia oscura forma parte importante de la comunidad, no son indispensables para la red, ya que podría tener mecanismos de resiliencia que hagan que la red pueda mantener su estabilidad aun cuando la estructura de la red cambia.   

## **Licencia:**

Este proyecto está bajo la Licencia MIT [LICENSE](LICENSE). Como parte del curso de Genómica Funcional con fines educativos. Copyright (c) 2025 mtz-mel

## **Referencias:**     
- Álvarez, J., Fernández Real, J. M., Guarner, F., Gueimonde, M., Rodríguez, J. M., Saenz de Pipaon, M., & Sanz, Y. (2021). Microbiota intestinal y salud. Gastroenterología y Hepatología, 44(7), 519-535. https://doi.org/10.1016/j.gastrohep.2021.01.009
- Del Campo-Moreno, R., Alarcón-Cavero, T., D’Auria, G., Delgado-Palacio, S., & Ferrer-Martínez, M. (2018). Microbiota and Human Health: Characterization techniques and transference. Enfermedades Infecciosas y Microbiologia Clinica (English Ed ), 36(4), 241-245. https://doi.org/10.1016/j.eimce.2018.02.016  
- Sprockett, D. D., Martin, M., Costello, E. K., Burns, A. R., Holmes, S. P., Gurven, M. D., & Relman, D. A. (2020). Microbiota assembly, structure, and dynamics among Tsimane horticulturalists of the Bolivian Amazon. Nature Communications, 11(1). https://doi.org/10.1038/s41467-020-17541-6
- Kurtz, Z. D., Müller, C. L., Miraldi, E. R., Littman, D. R., Blaser, M. J., & Bonneau, R. A. (2015). Sparse and compositionally robust inference of microbial ecological networks. PLoS Computational Biology, 11(5), e1004226. https://doi.org/10.1371/journal.pcbi.1004226
- Zamkovaya, T., Foster, J. S., De Crécy-Lagard, V., & Conesa, A. (2020). A network approach to elucidate and prioritize microbial dark matter in microbial communities. The ISME Journal, 15(1), 228–244. https://doi.org/10.1038/s41396-020-00777-x
  

