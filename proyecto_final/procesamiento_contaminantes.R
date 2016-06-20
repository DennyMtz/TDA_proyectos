#setwd("/Users/fernanda/Downloads/")
setwd("/home/denny/itam/topologia/proyecto_final")
install.packages('lubridate')
install.packages('dummies')
install.packages('reshape')
library(lubridate)
library(dummies)
library(reshape)

#leemos los datos
############################## ESTA PARTE SE CORRE MANUAL TRES VECES ##############################
#para cada archivo le cambio manualmente 2014 por 2015 y así
datos_raw <- read.csv("contaminantes_2009.csv",skip=10)
#no sirven si vienen sin valor
datos_sin_na <- datos_raw[is.na(datos_raw$value)==FALSE,]

#queremos crear una columna con el numero de semana
datos_sin_na$fecha_real <- strptime(substr(as.character(datos_sin_na$date),1,10),format="%d/%m/%Y")
datos_sin_na$semana <- week(datos_sin_na$fecha_real)
#para cada archivo le cambio manualmente 2014 por 2015 y así
datos_sin_na$id <- paste0(datos_sin_na$id_station,datos_sin_na$semana,"2009")

#eliminamos columnas que no se usaron
datos_sin_na$date <- NULL; datos_sin_na$unit <- NULL; datos_sin_na$id_station <- NULL
datos_sin_na$fecha_real <- NULL;datos_sin_na$semana <- NULL
remove(datos_raw)

#para cada archivo le cambio manualmente 2014 por 2015 y así
data_2009 <- cast(datos_sin_na, id~id_parameter, mean)

remove(datos_sin_na)
###################################################################################################

#ya que estén todos
todos_años <- rbind(data_2014,data_2015,data_2016)
write.csv(todos_años,'contaminantes.csv',row.names=FALSE)
