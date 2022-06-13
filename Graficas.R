############# Preparacion ##############
# Librerias
library(RColorBrewer) #brewer.pal
library(ggplot2)
library(gridExtra)
library(knitr)
library(dplyr)

# Datos
# cargar datos y llamarlos Data
# Data <- load()
# Cargar cartografia y llamarla Map
# Map <- load()

# Definicion de numero total de areas y numero total de años
# Area = nº de municipios, barrios,... que se esten considerando como unidad de vecindad repetidas el numero de años que se tengan disponibles. (eg. Si hay 5 municipios y 2 años: 
# Area = (1,2,3,4,5,1,2,3,4,5))

# Anyos = nº de años que se tengan disponibles cada uno de ellos repetidos por el numero total de municipios (barrios,..etc que se esten considerando como unidad de vecindad) (eg. Si hay 2 años y 5 municipios: 
#Anyos = (1,1,1,1,1,2,2,2,2,2)) 

# Siguiendo los ejemplos anteriores:
TotalMunicipios <- 5
TotalAnyos <- 2
Area <- rep(1:TotalMunicipios, TotalAnyos) 
Anyos <- rep(1:TotalAnyos, each = TotalMunicipios)
Data <- data.frame(Data, Area, Anyos)
S <- length(unique(Data$Area))
T <- length(unique(Data$Anyos))


############# Representacion variable respuesta ##############
# Graficos de la suma de los casos observados de cada año (columnas del vector suma del 1 al T, suma total de casos observados de todos los años columna T+1 del vector suma)
# R = variable respuesta (definida dentro de Data)
R <- matrix(Data$R, nrow=S, ncol=T)
suma <- matrix(NA, nrow=1, ncol=T+1)
for (i in 1:T){
  suma[i] <- sum(ordenes[,i])
  for (i in (T+1)){
    suma[9] <- sum(ordenes[,1:8])
  }
}

RTodas <- data.frame(Años = 1:T, 
                           R = suma[1:8])

ggplot(RTodas, aes(x=Años, y=R)) + geom_point(colour = "black", size = 4) + labs(y="R\n", x="\nAños") + theme(axis.title=element_text(size=25), axis.text=element_text(size=20)) +
  scale_x_continuous(breaks= seq(1, T, 1)) + scale_y_continuous(limits = c(600, 850), breaks = seq(600, 850, by = 50))






############# Modelo espacial (2). Distribuciones a posteriori efectos espaciales ##############
# cargar resultados del modelo (y llamarlo ResMod2):
# load(file.path())

# BYM CONJUNTO:
MargSD <- inla.tmarginal(function(x) 1/sqrt(x),
                         ResMod2$marginals.hyperpar$`Precision for ID (iid component)`)

MargSD2 <- inla.tmarginal(function(x) 1/sqrt(x),
                          ResMod2$marginals.hyperpar$`Precision for ID (spatial component)`)


grid.arrange(
  
  ggplot(data.frame(inla.smarginal(MargSD2)), aes(x, y)) +
    geom_line(size = 1.2) + ggtitle("INLA: Desviación típica del efecto \nespacial estructurado") +
    theme_gray() + theme(plot.title = element_text(size=15)) + theme(axis.title=element_text(size=15), axis.text=element_text(size=15)) +
    scale_x_continuous(breaks= seq(0, 0.7, 0.1)) + scale_y_continuous(breaks = seq(0, 7, by = 1)) ,
  
  ggplot(data.frame(inla.smarginal(MargSD)), aes(x, y)) +
    geom_line(size = 1.2) + ggtitle("INLA: Desviación típica del efecto \nespacial no estructurado") +
    theme_gray() + theme(plot.title = element_text(size=15)) + theme(axis.title=element_text(size=15), axis.text=element_text(size=15))    
  ,ncol=2, nrow=1)







############# Modelo espacial con heterogeneidad temporal (3). Distribuciones efectos espaciales a posteriori ##############
# cargar resultados del modelo (y llamarlo ResMod3):
# load(file.path())

MargSD <- inla.tmarginal(function(x) 1/sqrt(x),
                         ResMod3$marginals.hyperpar$`Precision for ID (iid component)`)
MargSD2 <- inla.tmarginal(function(x) 1/sqrt(x),
                          ResMod3$marginals.hyperpar$`Precision for ID (spatial component)`)


grid.arrange(
  
  ggplot(data.frame(inla.smarginal(MargSD2)), aes(x, y)) +
    geom_line(size = 1.2) + ggtitle("INLA: Desviación típica del efecto \nespacial estructurado") +
    theme_grey() + theme(plot.title = element_text(size=9.5)) + theme(plot.title = element_text(size=15)) + theme(axis.title=element_text(size=15), axis.text=element_text(size=15)),

  
  ggplot(data.frame(inla.smarginal(MargSD)), aes(x, y)) +
    geom_line(size = 1.2) + ggtitle("INLA: Desviación típica del efecto \nespacial no estructurado") +
    theme_grey() + theme(plot.title = element_text(size=9.5)) + theme(plot.title = element_text(size=15)) + theme(axis.title=element_text(size=15), axis.text=element_text(size=15))
  
  ,ncol=2, nrow=1)





############# Modelo espacial con heterogeneidad temporal (3). Distribuciones efecto temporal a posteriori ##############
# cargar resultados del modelo (y llamarlo ResMod3):
# load(file.path())

MargSD3 <- inla.tmarginal(function(x) 1/sqrt(x),
                          ResMod3$marginals.hyperpar$`Precision for Anyos`)

ggplot(data.frame(inla.smarginal(MargSD3)), aes(x, y)) +
    geom_line(size = 1.2) + ggtitle("INLA: Desviación típica del efecto \ntemporal no estructurado") +
    theme_grey() + theme(plot.title = element_text(size=9)) + theme(plot.title = element_text(size=15)) + theme(axis.title=element_text(size=15), axis.text=element_text(size=15)) + scale_y_continuous(breaks = seq(0, 50, by = 10))





############# Modelo espacio-temporal (4). Riesgos relativos (Mapa) ##############
# cargar resultados del modelo (y llamarlo ResMod4):
# load(file.path())
# T y S definidas al principio
# Map = cartografia (SpatialPolygonsDataFrame)

par(oma = c(0,0,0,0), mfrow = c(2, 4), mar = c(0, 0, 0, 0))
for (i in 1:T){
  
  RME <- matrix(ResMod4$summary.fitted.values$mean, nrow = S, ncol=T)
  Grupos <- c(0,0.7,0.9,1.1,1.5,2)
  GruposRBYM <- findInterval(RME[,i], Grupos)
  Paleta2 <- brewer.pal(9, "Blues")[c(1,3,4,5,8,9)]
  plot(Map, col=Paleta2[GruposRBYM])
  tit <- c("Año 1", "Año 2", "Año 3", "Año 4", "Año 5", "Año 6", "Año 7", "Año 8")
  title(paste(tit[i]), cex.main=2, line = -3.25)
  
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
Paleta2 <- brewer.pal(9, "Blues")[c(1,3,4,5,8,9)]
legend(0.830654,1.1031937, legend = c("<0.7", "0.7-<0.9", "0.9-<1.1", "1.1-<1.5", "1.5-<2", ">=2"), fill=Paleta2, lwd = 1, xpd = F, horiz = FALSE, cex = 1.4, seg.len=0, x.intersp = 0.3, y.intersp = 0.5, text.width=0.1, bty = "n",)





############# Modelo espacio-temporal (4). Riesgos relativos a través de los anyos ##############
# Grafico de los riesgos relativos que se han mantenido altos o bajos de forma constante con el paso de los años.
# cargar resultados del modelo (y llamarlo ResMod4):
# load(file.path())
# T y S definidas al principio
# Map = cartografia (SpatialPolygonsDataFrame)


RME <- matrix(ResMod4$summary.fitted.values$mean, nrow = S, ncol=T)

library(dplyr)
a <- as.data.frame(RME) %>% filter_all(all_vars(.> 1.69))
a <- as.matrix(a) # riesgos altos constantes
b <- as.data.frame(RME) %>% filter_all(all_vars(.< 0.5))
b <- as.matrix(b) # riesgos bajos constantes


library(ggplot2)
xy <- data.frame(on_x = 1:T, on_y = a[1,])
p <- ggplot(xy, aes(x = on_x, y = on_y)) + 
  scale_x_continuous(name="Años",
                     limits=c(1,T),
                     breaks=seq(1,T,1)) + 
  scale_y_continuous(name="Riesgo relativo",
                     limits=c(0,4.2),
                     breaks=seq(0,4.2,1)) +
  theme_minimal() + 
  theme(axis.title=element_text(size=29),
        axis.text=element_text(size=22)) + theme(axis.title.x = element_text(vjust = -0.1)) + theme(axis.title.y = element_text(vjust = 1.1))

for (i in 1:dim(a)[[1]]) {
  
  xy <- data.frame(on_x = 1:T, on_y = a[i,])  
  p <- p + geom_line(data = xy, aes(x = on_x, y = on_y), color="firebrick4")
  
}

for (i in 1:dim(b)[[1]]) {
  
  xy <- data.frame(on_x = 1:T, on_y = b[i,])  
  p <- p + geom_line(data = xy, aes(x = on_x, y = on_y), color="dodgerblue4")
  
}
p + geom_hline(yintercept=1, linetype="dashed")





############# Modelo 4. Riesgos relativos a través de los anyos MAPA ##############
# Grafico de los riesgos relativos que se han mantenido altos o bajos de forma constante con el paso de los años. VERSION MAPA
# cargar resultados del modelo (y llamarlo ResMod4):
# load(file.path())
# T y S definidas al principio
# Map = cartografia (SpatialPolygonsDataFrame)

library(dplyr)
RME <- matrix(ResMod4$summary.fitted.values$mean, nrow = S, ncol=T)


a <- as.data.frame(RME) %>% filter_all(all_vars(.> 1.69))
a <- as.matrix(a) # riesgos altos constantes
b <- as.data.frame(RME) %>% filter_all(all_vars(.< 0.5))
b <- as.matrix(b) # riesgos bajos constantes

PosAlto <- matrix(NA, nrow=1, ncol=dim(a)[[1]])
for (i in 1:dim(a)[[1]]) {
  PosAlto[,i] <- which(RME[]==a[i,1])
}
PosBajo <- matrix(NA, nrow=1, ncol=dim(b)[[1]])
for (i in 1:dim(b)[[1]]) {
  PosBajo[,i] <- which(RME[]==b[i,1])
}

par(mar=c(0.5,0.5,0.5,0.5))
plot(Map)
for (i in PosAlto) {
  polygon(Map@polygons[[i]]@Polygons[[1]]@coords[,1], 
          Map@polygons[[i]]@Polygons[[1]]@coords[,2], col="firebrick4")
}
for (i in PosBajo) {
  polygon(Map@polygons[[i]]@Polygons[[1]]@coords[,1], 
          Map@polygons[[i]]@Polygons[[1]]@coords[,2], col="dodgerblue4")
}
legend(727374.3,4376627, legend = c("Riesgo bajo", "Riesgo alto"), fill=c("dodgerblue4", "firebrick4"), lwd = 1, xpd = F,  cex = 1.7, seg.len=0, x.intersp = 0.5, y.intersp = 0.7, text.width=1, bty = "n",)







############# Modelo espacio-temporal (4). Riesgos relativos que mas han variado a traves de los anyos ##############
# Riesgos que mas han cambiado con el paso de los años.
# cargar resultados del modelo (y llamarlo ResMod4):
# load(file.path())
# T y S definidas al principio

RME <- matrix(ResMod4$summary.fitted.values$mean, nrow = S, ncol=T)
RME_1 <- RME

# 1.)
MenosAMas <- RME_1
## Elegir las secciones en las que el ultimo año tenga un riesgo mas alto que el primer año
MenosAMasDF <- data.frame("Año1" = MenosAMas[,1],
                          "AñoFin" = MenosAMas[,T],
                          Posicion = 1:S)
posMenosAMas <- NULL
for (i in 1:dim(RME_1)[[1]]){
  if (MenosAMasDF[i,1] < MenosAMasDF[i,2])
    posMenosAMas <- rbind(posMenosAMas,
                          data.frame(Pos = MenosAMasDF$Posicion[i],
                                     Dif = (MenosAMasDF[i,2]-MenosAMasDF[i,1])))
}
View(posMenosAMas)

posMenosAMas[tail(order(posMenosAMas$Dif), 9),] # las secciones con las diferencias mas grandes entre el ultimo año y el primero
posMenosAMas$Pos[tail(order(posMenosAMas$Dif), 9)] # el numero de las secciones con las mayores diferencias
a <- MenosAMas[posMenosAMas$Pos[tail(order(posMenosAMas$Dif), 9)], ]



# 2.)
MasAMenos <- RME_1
## Elegir las secciones en las que el ultimo año tenga un riesgo mas bajo que el primer año
MasAMenosDF <- data.frame("Año1" = MasAMenos[,1],
                          "AñoFin" = MasAMenos[,T],
                          Posicion = 1:S)
posMasAMenos <- NULL
for (i in 1:dim(RME_1)[[1]]){
  if (MasAMenosDF[i,1] > MasAMenosDF[i,2])
    posMasAMenos <- rbind(posMasAMenos,
                          data.frame(Pos = MasAMenosDF$Posicion[i],
                                     Dif = (MasAMenosDF[i,2]-MasAMenos[i,1])))
}
View(posMasAMenos)

posMasAMenos[head(order(posMasAMenos$Dif), 9), ]
posMasAMenos$Pos[head(order(posMasAMenos$Dif), 9)]
b <- MasAMenos[posMasAMenos$Pos[head(order(posMasAMenos$Dif), 9)], ]



# Grafica
xy <- data.frame(on_x = 1:T, on_y = a[1,])
p <- ggplot(xy, aes(x = on_x, y = on_y)) + 
  scale_x_continuous(name="Años",
                     limits=c(1,T),
                     breaks=seq(1,T,1)) + 
  scale_y_continuous(name="Riesgo relativo",
                     limits=c(0,4.2),
                     breaks=seq(0,4.2,1)) +
  theme_minimal() + 
  theme(axis.title=element_text(size=29),
        axis.text=element_text(size=22)) +
  theme(axis.title.x = element_text(vjust = -0.1)) + 
  theme(axis.title.y = element_text(vjust = 1.1)) 


for (i in 1:dim(a)[[1]]) {
  
  xy <- data.frame(on_x = 1:T, on_y = a[i,])  
  p <- p + geom_line(data = xy, aes(x = on_x, y = on_y), color="firebrick4")
  
}

for (i in 1:dim(b)[[1]]) {
  
  xy <- data.frame(on_x = 1:T, on_y = b[i,])  
  p <- p + geom_line(data = xy, aes(x = on_x, y = on_y), color="dodgerblue4")
  
}
p + geom_hline(yintercept=1, linetype="dashed")





############# Modelo espacio-temporal (4). Riesgos relativos que mas han variado a traves de los anyos MAPA ##############
# Riesgos que mas han cambiado con el paso de los años. VERSION MAPA
# cargar resultados del modelo (y llamarlo ResMod4):
# load(file.path())
# T y S definidas al principio
# Map = cartografia (SpatialPolygonsDataFrame)

# posMenosAMas y posMasAMenos calculadas en la grafica anterior.

PosMenosAMas <- posMenosAMas$Pos[tail(order(posMenosAMas$Dif), 9)]
PosMasAMenos <- posMasAMenos$Pos[head(order(posMasAMenos$Dif), 9)]


par(mar=c(0.5,0.5,0.5,0.5))
plot(Map)
for (i in PosMenosAMas) {
  polygon(Map@polygons[[i]]@Polygons[[1]]@coords[,1], 
          Map@polygons[[i]]@Polygons[[1]]@coords[,2], col="firebrick4")
}
for (i in PosMasAMenos) {
  polygon(Map@polygons[[i]]@Polygons[[1]]@coords[,1], 
          Map@polygons[[i]]@Polygons[[1]]@coords[,2], col="dodgerblue4")
}
legend(726724.3,4376627, legend = c("Riesgo creciente", "Riesgo decreciente"), fill=c("firebrick", "dodgerblue4"), lwd = 1, xpd = F,  cex = 1.7, seg.len=0, x.intersp = 0.5, y.intersp = 0.7, text.width=1, bty = "n",)

