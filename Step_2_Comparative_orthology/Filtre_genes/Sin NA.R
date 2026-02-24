######################################################
## Tema: Filtrar filas con NA para identificar      ##
## las filas coincidentes                           ##
## Autor: Olga Andrea Hernandez Miranda, Miranda H  ##
## Fecha: 09/13/2024                                ##
## Nota: Filtrar ortologos coincidentes             ##
######################################################

directorio <- "C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutorial4/Prueba/OrthoGrupos_version2"
setwd(directorio)


# Cargar el archivo CSV
datos <- read.csv("ParaFiltrarNA.csv", header = TRUE)

# Eliminar filas con NA
datos_sin_na <- na.omit(datos)

# Escribir el nuevo archivo CSV
write.csv(datos_sin_na, file = "archivo_sin_na.csv", row.names = FALSE)

