# Modelo espacial + heterogeneidad temporal
############# Preparacion ##############
# Librerias
library(rgdal) #readOGR
library(spdep) #poly2nb
library(INLA)
library(RColorBrewer) #brewer.pal
library(ggplot2)
library(gridExtra)
library(knitr)



# Datos
# cargar datos y llamarlos Data
# Data <- load()



# Estructura de vecindad
# Map es la cartografia, de tipo SpatialToPolygonsDataFrama
# una vez cargada la cartografia se construye la estructura de vecindad. Para ello se utilizan las funciones poly2nb y nb2INLA, que determinan como vecinos aquellos poligonos que son contiguos y transforma el sistema de vecindad en listas de vecinos, respectivamente.
nb <- poly2nb(Map, snap = 0.000001)
head(nb)
nb2INLA("MapGraph", nb)
g <- inla.read.graph("MapGraph")
plot(g)
image(inla.graph2matrix(g),xlab="",ylab="")





############# Formula y Modelo ##############
# Efectos aleatorios espaciales
# U = sin estructura espacial
# S = con estructura espacial
# H = efecto temporal (Anyos)
# length(nb) = numero de municipios, barrios,... que se esten usando como unidad de vecindad
U <- 1:length(length(nb)) 
S <- 1:length(length(nb))
Anyos <- rep(1:years, each = length(nb)) # years = nº de años (eg: si hay informacion desde 2011 hasta 2015, years = 5)
Data2 <- Data
Data2 <- data.frame(Data2, U, S, Anyos)


# ESPECIFICANCO LAS DOS COMPONENTES BYM DE FORMA SEPARADA
# V1,..., Vn son las n covariables (definidas dentro de Data)
# R es la variable respuesta (definida dentro de Data) (casos observados de la enfermedad)
# Rellenar el campo de formula con la variable respuesta y las variables que se quieran estudiar
# Formula
Formula <- R ~ V1 + V2 + V3 + V4 + V5 + V6 + 
  f(U, 
    model = "iid") + 
  f(S, 
    model = "besag", 
    graph = g, 
    scale.model = TRUE) + 
  f(Anyos, 
    model = "iid")

# Esperados = Casos esperados (variable definida dentro de Data)
# Modelo
ResMod3 <- inla(Formula, 
                data=Data2,
                E=Esperados,
                family="poisson",
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                control.predictor = list(compute=TRUE, cdf=c(log(1))),
                control.fixed = list(mean.intercept = 0, prec.intercept = 0, mean = 0, prec = 0))
# verbose = TRUE



# ESPECIFICANDO LAS COMPONENTES DEL BYM A LA VEZ (MEJOR, EL QUE ESTÁ GUARDADO)
ID <- 1:length(length(nb)) 
Anyos <- rep(1:years, each = length(nb)) # years = nº de años (eg: si hay informacion desde 2011 hasta 2015, years = 5)
Data2 <- Data
Data2 <- data.frame(Data2, ID, Anyos)

# Formula
Formula <- R ~ V1 + V2 + V3 + V4 + V5 + V6 + 
  f(ID, 
    model="bym",
    graph=g, 
    hyper  =
      list(prec.unstruct  =  list(prior="pc.prec",param=c(1,0.01)),
           prec.spatial = list(prior="pc.prec",param=c(1,0.01)))) + 
  f(Anyos, 
    model = "iid",
    hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01))))


# Modelo
ResMod3 <- inla(Formula, 
                data=Data2,
                E=Esperados,
                family="poisson",
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                control.predictor = list(compute=TRUE, cdf=c(log(1))),
                control.fixed = list(mean.intercept = 0, prec.intercept = 0, mean = 0, prec = 0),
                verbose = TRUE)


# Modelo BYM2
# Formula
ID <- 1:length(length(nb))
Anyos <- rep(1:years, each = length(nb)) # years = nº de años (eg: si hay informacion desde 2011 hasta 2015, years = 5)
Data2 <- Data
Data2 <- data.frame(Data2, ID, Anyos)

Formula <- R ~ V1 + V2 + V3 + V4 + V5 + V6 +
  f(ID, 
    model="bym2",
    graph=g, 
    hyper  =
      list(prec  =  list(prior="pc.prec", param=c(1,0.01)),
           phi = list(prior="pc.prec", param=c(1,0.01)))) + 
  f(Anyos, 
    model = "iid",
    hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01))))


# Modelo
ResMod3 <- inla(Formula, 
                data=Data2,
                E=Esperados,
                family="poisson",
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                control.predictor = list(compute=TRUE, cdf=c(log(1))),
                control.fixed = list(mean.intercept = 0, prec.intercept = 0, mean = 0, prec = 0),
                verbose = TRUE)

# Convergencia del modelo con verbose = TRUE


ResMod3$cpu.used

ResMod3$dic$dic


############# Resultados ##############
# Resultados
round(ResMod3$summary.fixed, 5)

round(ResMod3$summary.hyperpar, 5)

# Proporción de varianza explicada (Blangiardo p.12).
# Cogemos el efecto estructurado y calculamos la varianza empírica. Creamos matriz con filas=secciones y columnas=1000.
# Para cada seccion extraemos 1000 valores de la distribucion marginal correspondiente y calculamos la varianza empirica. 
# Extraemos el valor esperado de la varianza del efecto no estructurado y construimos  la varianza espacial.

Nsecciones <- length(length(nb))

M <- ResMod3$marginals.random$ID
MatMarg <- matrix(NA, nrow=Nsecciones, ncol=1000)
for (i in 1:Nsecciones){
  u <- M[[Nsecciones+i]]
  s <- inla.rmarginal(1000, u)
  MatMarg[i,] <- s
}

var.RRspatial <- mean(apply(MatMarg, 2, sd))^2
var.RRhet <- inla.emarginal(function(x) 1/x,
                            ResMod3$marginals.hyper$"Precision for ID (iid component)")
var.RRspatial/(var.RRspatial + var.RRhet)




############# Graficas de las distribuciones a posteriori del modelo ##############

# intercepto
ggplot(data.frame(data.frame(inla.smarginal(ResMod3$marginals.fixed[[1]]))), aes(x,y)) +
    geom_line() + ggtitle("INLA: Intercept") +
    theme_bw() + theme(plot.title = element_text(size=9.5))

kable(data.frame(INLAalpha = inla.emarginal(function(x) x, ResMod3$marginals.fixed$`(Intercept)`)))


# betas
# entre corchetes se indica el numero de cada covariable
# 1 = intercept, por eso no se pone
BetasINLA <- c(ResMod3$marginals.fixed[2], ResMod3$marginals.fixed[3], ResMod3$marginals.fixed[4], ResMod3$marginals.fixed[5], ResMod3$marginals.fixed[6], ResMod3$marginals.fixed[7])



# n es el numero de covariables
n <- 6
par(mfrow=c(2,4))
par(mar=c(3,3,3,3))
for (i in 1:n){
  plot(inla.smarginal(BetasINLA[[i]]), type="l", xlab="", ylab="")
  tit <- c("INLA: V1", "INLA: V2", "INLA: V3", "INLA: V4", "INLA V5", "INLA V6")
  title(tit[i], cex.main=1)
}  





# Transformar las precisiones en desviaciones tipicas
# COMPONENTES BYM DE FORMA SEPARADA:
MargSD <- inla.tmarginal(function(x) 1/sqrt(x),
                         ResMod3$marginals.hyperpar$`Precision for U`)

MargSD2 <- inla.tmarginal(function(x) 1/sqrt(x),
                          ResMod3$marginals.hyperpar$`Precision for S`)

MargSD3 <- inla.tmarginal(function(x) 1/sqrt(x),
                          ResMod3$marginals.hyperpar$`Precision for Anyos`)


# BYM CONJUNTO:
MargSD <- inla.tmarginal(function(x) 1/sqrt(x),
                         ResMod3$marginals.hyperpar$`Precision for ID (iid component)`)

MargSD2 <- inla.tmarginal(function(x) 1/sqrt(x),
                          ResMod3$marginals.hyperpar$`Precision for ID (spatial component)`)

MargSD3 <- inla.tmarginal(function(x) 1/sqrt(x),
                          ResMod3$marginals.hyperpar$`Precision for Anyos`)


# Obtener la media a posteriori de la desviacion tipica de los efectos espaciales de INLA
MeanMargSD <- inla.emarginal(function(x) x, MargSD)
MeanMargSD2 <- inla.emarginal(function(x) x, MargSD2)
MeanMargSD3 <- inla.emarginal(function(x) x, MargSD3)

# Estadisticos basicos de las desviaciones tipicas
SumMargSD <- inla.zmarginal(MargSD)
SumMargSD2 <- inla.zmarginal(MargSD2)
SumMargSD3 <- inla.zmarginal(MargSD3)

kable(data.frame(rbind(sigma.u=c(SumMargSD$mean, SumMargSD$sd, SumMargSD$quant0.025, SumMargSD$quant0.975), sigma.s=c(SumMargSD2$mean, SumMargSD2$sd, SumMargSD2$quant0.025, SumMargSD2$quant0.975), sigma.h=c(SumMargSD3$mean, SumMargSD3$sd, SumMargSD3$quant0.025, SumMargSD3$quant0.975))),col.names=c("Mean", "SD", "2.5%", "97.5%"))




# Summary
round(ResMod3$summary.fixed, 5) 




# Representacion de las distribuciones a posteriori de las desviaciones tipicas de los efectos aleatorios espaciales del modelo de INLA
grid.arrange(
  
  ggplot(data.frame(inla.smarginal(MargSD)), aes(x, y)) +
    geom_line() + ggtitle("INLA: Desviación de U") +
    theme_bw() + theme(plot.title = element_text(size=10.5)),
  
  ggplot(data.frame(inla.smarginal(MargSD2)), aes(x, y)) +
    geom_line() + ggtitle("INLA: Desviación de S") +
    theme_bw() + theme(plot.title = element_text(size=10.5))
  
  ,ncol=2)


# Representacion de la distribucion a posteriori de la desviacion tipica del efecto temporal del modelo de INLA
ggplot(data.frame(inla.smarginal(MargSD3)), aes(x, y)) +
  geom_line() + ggtitle("INLA: Desviación de H") +
  theme_bw() + theme(plot.title = element_text(size=10.5))








############# Graficas RME ##############
# Dibujo RME del primer año
# years = numero de años (eg: si hay informacion desde 2011 hasta 2015, years = 5)
# length(nb) = numero de municipios, barrios,... que se esten usando como unidad de vecindad
# Map = cartografia (SpatialPolygonsDataFrame)
RME <- matrix(ResMod3$summary.fitted.values$mean, nrow = length(nb), ncol=years)
par(mar=c(0,0,2,0))
par(mfrow=c(1,1))
quantile(RME[,1], probs=seq(0,1,0.1))*100
Grupos <- c(0,35,70,105,140,175,210,245,280)
GruposRBYM <- findInterval(RME[,1]*100, Grupos)
table(GruposRBYM)
Paleta2 <- brewer.pal(9, "BrBG")[9:1]
plot(Map, col=Paleta2[GruposRBYM], main = "RMEs INLA Año 1")
legend(727759.3,4376852, c("<35","35-<75","75-<105", "105-<140","140-<175","175-<210", "210-<245", "245-280",">280"), 
       fill=Paleta2, border=NULL, bty="n", cex = 0.6)



# Mostrar regiones con RME > 280
par(mar=c(0,0,2,0))
par(mfrow=c(2,4))
for (i in 1:years){
  
  RME <- matrix(ResMod3$summary.fitted.values$mean, nrow = length(nb), ncol=years)
  pos <- which((RME*100)[,i] > 280)
  
  Grupos <- c(0,35,70,105,140,175,210,245,280)
  GruposRBYM <- findInterval(RME[,i]*100, Grupos)
  Paleta2 <- brewer.pal(9, "BrBG")[9:1]
  plot(Map, col=Paleta2[GruposRBYM])
  tit <- c("RME Año 1", "RME Año 2", "RME Año 3", "RME Año 4", "RME Año 5", "RME Año 6", "RME Año 7", "RME Año 8")
  title(paste(tit[i]), cex.main=2)
  
  
  for (j in 1:length(pos)) {
    lines(Map@polygons[[pos[j]]]@Polygons[[1]]@coords, col="red")
  }
  
}



# p(RME) > 1 (exceso de riesgo)
# Dibujo p(RME) > 1 del primer año
# years = numero de años (eg: si hay informacion desde 2011 hasta 2015, years = 5)
# length(nb) = numero de municipios, barrios,... que se esten usando como unidad de vecindad
# Map = cartografia (SpatialPolygonsDataFrame)
ProbRME <- matrix(1 - ResMod3$summary.fitted.values$`1 cdf`, nrow = length(nb), ncol=years)

GruposProb <- c(0,0.50,0.75,0.90,0.95)
GruposProbBYM <- findInterval(ProbRME[,1], GruposProb)
Paleta2 <- brewer.pal(5, "BrBG")[5:1]
par(mar=c(0,0,2,0))
par(mfrow=c(1,1))
plot(Map, col=Paleta2[GruposProbBYM], main = "p(RME > 1) INLA 2011")
legend(727878,4376546,c("<0.50", "0.51-0.75","0.76-0.90", "0.91-0.95",">0.95"), 
       fill=Paleta2, border=NULL, bty="n", cex=0.8)




# Mostrar regiones con p(RME) > 0.95
par(mar=c(0,0,2,0))
par(mfrow=c(2,4))
for (i in 1:years){
  
  ProbRME <- matrix(1 - ResMod3$summary.fitted.values$`1 cdf`, nrow = length(nb), ncol=years)
  pos <- which((ProbRME)[,i] > 0.95)
  
  GruposProb <- c(0,0.50,0.75,0.90,0.95)
  GruposRBYM <- findInterval(ProbRME[,i], GruposProb)
  Paleta2 <- brewer.pal(5, "BrBG")[5:1]
  plot(Map, col=Paleta2[GruposRBYM])
  tit <- c("p(RME) Año 1", "p(RME) Año 2", "p(RME) Año 3", "p(RME) Año 4", "p(RME) Año 5", "p(RME) Año 6", "p(RME) Año 7", "p(RME) Año 8")
  title(paste(tit[i]), cex.main=2)
  
  
  for (j in 1:length(pos)) {
    lines(Map@polygons[[pos[j]]]@Polygons[[1]]@coords, col="red")
  }
}



