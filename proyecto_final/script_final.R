library(fpc)
library(dbscan)
library(kernlab)
library(jsonlite)
library(igraph)
library(RColorBrewer)
library(qrage)
library(caret)

#----------------------------- LOAD DATA -----------------------------
rm(list = ls())
rutawork = ('/home/denny/itam/topologia/proyecto_final/')
datos <- read.csv(paste(rutawork,'contaminantes.csv',sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)
datos_origin<-datos
str(datos)
datos<-datos_origin[c("CO","NO","NO2","NOX","O3","PM10","PM2.5","SO2")]
fechas <- datos_origin["id"]

#----------------------------- PCA ----------------------------------
#Esta parte nos ayuda a seleccionar la variable mas significativa que explica la dinamica de los datos

#PCA
#install.packages('FactoMineR')
library('FactoMineR')
res.pca <- PCA(datos, scale.unit=TRUE, ncp=1, graph=T)
res.pca$eig[3] #cumulativa percentaje of variance
datos_prueba <- as.data.frame(c(res.pca$ind[1],res.pca$ind[2]))
names(datos_prueba)[1] <- "V1"
names(datos_prueba)[2] <- "V2"

#grafica:
res.pca = PCA(datos[,1:8], scale.unit=TRUE, ncp=5, graph=T)



# #KERNEL PCA
# variables <- length(datos)
# sigma <- 1
# kres <- kpca(~., data=datos,features=variables,kernel="rbfdot",kpar = list(sigma = sigma))
# data_kpca <- as.data.frame(kres@rotated)
# print("WEEEEEEEEEY  YA QUEDÓ EL PINCHE PCA")
# #ordena de mayor a menor (por la que tiene mayor varianza)
# datos_prueba <- data_kpca[c("V1","V2")]

#----------------------------- GENERATE INTERVALS----------------------------------

df <- datos_prueba      #choose a dataset

#----------------------------- NECESSARY PARAMETERS -----------------------------
#var_o <- data$x1    #variable we will use to make the overlapping subsets
var_o <- datos_prueba$V1   #if we want to use kernel pca variable to cut
n_int <- 6  #number of intervals we want
p <- 0.2  #proportion of each interval that should overlap with the next

#----------------------------- CREATING THE INTERVALS -----------------------------

#this section will create a data frame in which we will construct overlapping intervals
intervals_centers <- seq(min(var_o),max(var_o),length=n_int)  #basic partition = centers
interval_length <- intervals_centers[2]-intervals_centers[1]  #to create the overlaps of p% of this length
intervals <- data.frame(centers=intervals_centers)            #create a data frame
#create the overlapping intervals  
intervals$min <- intervals_centers - (0.5+p)*interval_length                     
intervals$max <- intervals_centers + (0.5+p)*interval_length
intervals$interval <- seq(1,n_int)
intervals$name <- with(intervals, sprintf("[%.2f;%.2f)",min,max))

#function that will split the variable according to the invervals
res <- lapply(split(intervals,intervals$interval), function(x){   
  return(df[var_o> x$min & var_o <= x$max,])     #res will be a list with each element res[i]
})                                               #being the points on the i'th subset

df$Clusters<-numeric(nrow(df))
for(i in 1:length(res)){
  # i<-1
  interval_temp<-res[[i]]
  distmat <- as.matrix(dist(interval_temp))
  neigh <- 3
  clusters_i<-clustITAM(distmat, neigh)
  interval_temp$Cluster<-numeric(nrow(interval_temp))
  
  for(j in 1:length(clusters_i)){
    #j<-1
    
    for(z in 1:length(clusters_i[[j]])){
    #  z<-1
      interval_temp$Cluster[as.numeric(clusters_i[[j]][z])]<-j
      #df$Clusters[as.numeric(rownames(interval_temp)[nrow(interval_temp)])]<-j
    }
  }
  res[[i]]<-interval_temp
}

for(i in 1:length(res)){
  
  df[[ncol(df)+1]]<-numeric(nrow(df))
  interval_temp<-res[[i]]
  interval_temp$filas<-rownames(interval_temp)
  for(j in 1:nrow(df)){
    
    if(rownames(df)[j] %in% rownames(interval_temp)){
      df[[ncol(df)]][j]<-interval_temp$Cluster[which(interval_temp$filas==rownames(df)[j])]
    }#else{
#       df[[ncol(df)+1]][j]<-0
#     }
    
  }
  
}

####---------Empieza reetiquetamiento de Clusters----------####
df_<-df
# df_<-df_salvado

for(i in 1:length(res)){
 
  ncol<-ncol(df_)
  df_[[ncol+1]]<-numeric(nrow(df_))
  ncol<-ncol(df_)
  
  if(i==1){
    
    contador_clusters<-0
    df_[[ncol]]<-df_[[3+i]]
    
  }else{
    
    for(j in 1:nrow(df_)){
      
      if(df_[[3+i]][j]!=0){
        
        df_[[ncol]][j]<-df_[[3+i]][j]+contador_clusters
        
        }
      }
    }
  contador_clusters<-contador_clusters+max(df_[[3+i]])
  
}

df_salvado<-df
df<-df_[,-c(3:9)]

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

#---------------------DISTANCE BETWEEN CLUSTERS AND ADJ_MATRIX -----------------------------

#Bajo el supuesto que leemos una base que tiene columnas para cada intervalo
#Leemos la base de datos con los intervalos
base <- df
#Creamos una columna para los clusters
base$clusters<-0


#Columna en donde se ubica el ultimo intervalo:
int_fin <- ncol(base)-1
#Columna en donde empieza el intervalo 1:
int_ini <- int_fin-n_int+1
#Columna donde se creo la columna de "clusters":
col_cluster <- ncol(base)


for(i in seq(nrow(base[,int_ini:int_fin]))){
  temp<-c()
  for(m in seq(int_ini,int_fin)){
    if (base[i,m] > 0){
      
      temp<-c(paste0(temp,base[i,m],sep = ","))
    }
  }
  if(length((temp))>0){
    aux<-unlist(strsplit(temp,","))
    aux2<-unique(aux)
    aux3<-paste(aux2,collapse=",")
    base[i,col_cluster]<-aux3
  }
}

base <- data.frame(base$clusters)
names(base) <- c("clusters")

#Creamos una variable para enumerar la observacion
base$obs <- paste("obs",seq(1,length(base$clusters)))

#Detectamos los clusters
num_clusters <- sort(as.numeric(unique(strsplit(paste0(base$clusters, collapse=","),",")[[1]])))
a<-strsplit(paste0(base$clusters, collapse=","),",")
clusters <- length((num_clusters))

#Creamos una columna con ceros para cada cluster
for(x in num_clusters){
  base[[paste("c",x,sep="_")]] <- rep(0,nrow(base))
}

#Para cada columna que creamos agregamos un 1 según el cluster al que pertenece la obs
base$clusters<- as.character(base$clusters)


#Ojo es x+3 (int_ini), porque en la columna 3 en adelante es donde va vaciar los "1" de cada cluster
for(i in seq(nrow(base))){
  vector <- strsplit(base$clusters[i], ",")[[1]]
  vector <- sort(as.numeric(vector))
  for(x in vector){
    base[i,(x+int_ini)] <- 1
  }
}

if(sum(base[[3]])==0){
  dummy_mat<-base[,5:ncol(base)-1]
}else{dummy_mat<-base[,4:ncol(base)-1]}

n_clusters<-ncol(dummy_mat)
cluster_list<-lapply(1:n_clusters,function (col) {which(dummy_mat[,col]==1)})

adj_matrix<-matrix(0,nrow=n_clusters,ncol=n_clusters)

for(i in 1:(n_clusters-1)){
  for(j in (i+1):n_clusters){
    distancia<-setdiff(cluster_list[[i]], cluster_list[[j]])
    cercania<-length(cluster_list[[i]])-length(distancia)
    adj_matrix[i,j]<-round(cercania/min(length(cluster_list[[i]]),length(cluster_list[[j]])),2)
    adj_matrix[j,i]<-adj_matrix[i,j]
  }
}

summary_cluster<-matrix(0,nrow=1,ncol=n_clusters)
for(i in 1:n_clusters){
  summary_cluster[1,i]<-length(cluster_list[[i]])
}



#######################
library(dplyr)
master <- cbind(datos,dummy_mat)
ww <- paste0("c_",seq(2,clusters,1))
master <- cbind(master,fechas)

bb <- data.frame()
for(x in 1:length(ww)){
  myCols <- (c("CO","NO","NO2", "NOX", "O3", "PM10", "PM2.5", "SO2", ww[x]))
  colNums <- match(myCols,names(master))
  a  <- (master  %>% select(colNums) %>% filter(master[8+x]=='1') %>% 
           summarise(NOX_prom=mean(NOX),PM10_prom=mean(PM10), PM2.5_prom=mean(PM2.5)))
  bb <- rbind.data.frame(bb,a)
}
#print(bb)



master$anio <- (as.numeric(substr(master$id,7,10)))
for(i in 1:length(master$anio)){
  if((master$anio[i]==7) | (master$anio[i]==8) | (master$anio[i]==9)){
    a <- paste0("200",master$anio[i])
    master$anio[i] <- a
  }
  else{
    a <- paste0("20",master$anio[i])
    master$anio[i] <- a
  }
}

#install.packages('modeest')
library(modeest)

master$anio <- as.integer(master$anio)
aa <- data.frame()
for(x in 1:length(ww)){
  myCols <- (c(ww[x], "anio"))
  colNums <- match(myCols,names(master))
  a  <- (master  %>% select(colNums) %>% filter(master[x+1]=='1') %>% 
           summarise(mode_anio= mlv(anio, method='mfv')[['M']]))
  aa <- rbind.data.frame(aa,a)
}
print(aa)


ss <- data.frame()
for(i in 1:length(ww)){
  s <- subset(master, master[i+8] == 1)
  mode <- mlv(s$anio)
  ss <- rbind.data.frame(ss,mode$M)
}
names(ss) <- "anio"

master$sem <- floor(as.numeric(substring(master$id,4))/10000)
# for(i in 1:length(master$sem)){
#   if(nchar(master$sem[i])==6){
#     a <- substring(master$sem[i],1,2)
#     master$sem[i] <- as.numeric(a)
#   } else{
#     a <- substring(master$sem[i],1,1)
#     master$sem[i] <- as.numeric(a)
#   }
# }

yy <- data.frame()
for(i in 1:length(ww)){
  y <- subset(master, master[i+8] == 1)
  mode <- mlv(y$sem)
  yy <- rbind.data.frame(yy,as.integer(mode$M))
}
names(yy) <- "sem"



###CREAR TOOLTIPLS


bb$tooltip <- paste0("NOX_prom: ",round(bb$NOX_prom,2), " ppb, PM10_prom: ",round(bb$PM10_prom,2), " mm, PM2.5_prom: ",round(bb$PM2.5_prom,2), " mm", " moda año:",ss$anio, ", moda semana:", yy$sem) 
tooltip <- as.array(bb$tooltip)

#install.packages('gplots')
library('gplots')
#vamos a ordenar los colores segun la variable "var_color"
var_color <- bb$PM2.5_prom
a <- col2hex(colorRampPalette(c("red", "blue"))(clusters))
a <- a[order(var_color)]
c <- paste(shQuote(a,type="cmd"), collapse=", ")
b <- paste(shQuote(order(var_color),type="cmd"), collapse=", ")
##############################



#KEPLER
nodes.n <- n_clusters
nodes.size<- as.numeric(summary_cluster)
#nodes.tooltips <- tooltip
nodes.tooltips <- paste("Grupo:", 1:nodes.n, tooltip)
nodes.names <- 1:nodes.n
nodes.color <- as.character(1:nodes.n)

# ------- AHORA TENEMOS QUE CREAR UN JSON DE ESO -----------------------------
adj.matrix <- adj_matrix
aux_mat <- data.frame()
for(i in 1:nodes.n) {
  
  for(j in 1:nodes.n){
    
  if(adj.matrix[i, j]!=0) {

        aux_mat <- rbind(aux_mat, data.frame(source=i-1, target=j-1, value=adj.matrix[i, j]))
    }
  }
}

linksJSON <- toJSON(aux_mat)
nodesJSON <- toJSON(data.frame(color=nodes.color, group=nodes.size, name=nodes.names, tooltip=nodes.tooltips))
graphJSON <- sprintf("{\"nodes\": %s, \"links\": %s}", nodesJSON, linksJSON)
#head(graphJSON)

# ------------  CREAMOS EL HTML ----------------------------------------------------------
htmlFile <- readLines('/home/denny/itam/topologia/proyecto_final/ManifoldLearning/www/index.html')
#htmlFile <- readLines("www/index.html")
graph_def_line <- which(grepl("graph =", htmlFile))
#htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
#writeLines(htmlFile, "www/index.html")
writeLines(htmlFile, '/home/denny/itam/topologia/proyecto_final/ManifoldLearning/www/index.html')

######colores FALTA CAMBIAR COLOR POR VARIABLE
#modificar el index
#setwd('..')
setwd("/home/denny/itam/topologia/proyecto_final/ManifoldLearning/www/")
htmlFile <- readLines('index.html')
numbers_line <- which(grepl("domain", htmlFile))
colors_line <- which(grepl("range", htmlFile))
htmlFile[numbers_line] <- paste0('.domain([',b,'])')
htmlFile[numbers_line+1] <- paste0('.range([',c,']);')
writeLines(htmlFile,'index.html')
###
browseURL("file:////home/denny/itam/topologia/proyecto_final/ManifoldLearning/www/index.html")

