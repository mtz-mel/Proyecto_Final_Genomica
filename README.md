**Proyecto final Genómica Funcional: Materia Obscura bacteriana y su impacto en la arquitectura de redes ecológicas**   

Integrantes:  
Melissa Martínez, Isabel Herrera y Alondra Dominguez.  

Descripción del proyecto:
La **materia oscura microbiana** consiste en aquellos taxones que no han podido ser identificados, lo cual dificulta la caracterización del microbioma en diversos ambientes, especialmente en aquellos que se encuentran bajo condiciones extremas. En este proyecto simulamos cómo la eliminación de esa materia oscura afecta la arquitectura de la red ecológica, comparando sus métricas de centralidad. Para ello trabajamos con la base de datos 
**SprockettTHData()** de la librería **microbiomeDataSets**, considerando como “materia oscura” los taxones identificados únicamente hasta género, y como “materia conocida” los que cuentan con asignación de especie. Esta distinción nos permitió generar redes con una densidad que nos permitiera visualizar con claridad los cambios estructurales.  

Librerías utilizadas: 
library(igraph)  
library(igraphdata)  
library(networkD3)  
library(RCy3)  
library(tidyverse)   
library(microbiomeDataSets)  
library(mia)  
library(phyloseq)  
library(SpiecEasi)  
