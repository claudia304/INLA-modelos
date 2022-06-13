# Modelo 4 (Artículo 1)
# Interaccion Tipo IV y RW1 prior para el tiempo



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

# Hay que definir las variables Area y Anyos dentro de Data

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
t.from <- min(Data$Anyos)
t.to <- max(Data$Anyos)



# Ordenar los datos para hacer el producto de kronecker de la interaccion espacio-temporal
V <- Data
V <- V[order(V[,"Anyos"],V[,"Area"]),]
head(V,30)

# N = numero de variables (casos observados (variable respuesta), casos esperados y las n covariables)
Data.INLA <- data.frame(V[1:N], 
                        ID.area = rep(1:S, T),
                        ID.year = rep(1:T, each=S),
                        ID.area.year = seq(1, T*S))
head(Data.INLA)


# Definir la matriz de estructura espacial de un LCAR
# Map es la cartografia, de tipo SpatialToPolygonsDataFrama
# una vez cargada la cartografia se construye la estructura de vecindad. Para ello se utilizan las funciones poly2nb y nb2INLA, que determinan como vecinos aquellos poligonos que son contiguos y transforma el sistema de vecindad en listas de vecinos, respectivamente.
nb <- poly2nb(Map, snap = 0.000001)
head(nb)
nb2INLA("MapGraph", nb)
g <- inla.read.graph("MapGraph")

Q.xi <- matrix(0, g$n, g$n)
for (i in 1:g$n) {
  Q.xi[i, i] = g$nnbs[[i]]
  Q.xi[i, g$nbs[[i]]] = -1
}

Q.Leroux <- diag(S) - Q.xi


# Definir la estructura temporal de la matriz de un RW1
D1 <- diff(diag(T), differences=1)
Q.gammaRW1 <- t(D1)%*%D1


# Definir hiperpriors apropiadas
# Unif(0, Inf) para las desviaciones estándar
# Unid(0,1) para el parámetro de suavizado espacial

sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

lunif = "expression:
    a = 1;
    b = 1;
    beta = exp(theta)/(1+exp(theta));
    logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
    log_jacobian = log(beta*(1-beta));
    return(logdens+log_jacobian)"

compute.patterns <- TRUE


strategy <- "gaussian"
# strategy <- "simplified.laplace"


# Restricciones necesarias
R <- kronecker(Q.gammaRW1, Q.xi)
r.def <- S+T-1
A1 <- kronecker(matrix(1,1,T), diag(S))
A2 <- kronecker(diag(T), matrix(1,1,S))
A.constr <- rbind(A1[-1,], A2[-1,])




############# Formula y Modelo ##############
# V1,..., Vn son las N covariables (definidas dentro de Data)
# R es la variable respuesta (definida dentro de Data) (casos observados de la enfermedad)
# Rellenar el campo de formula con la variable respuesta y las variables que se quieran estudiar

# Formula
Formula <- R ~ V1 + V2 + V3 + V4 + V5 + V6 + 
  
  f(ID.area,
    model="generic1",
    Cmatrix=Q.Leroux,
    constr = TRUE,
    hyper = list(prec=list(prior=sdunif), beta=list(prior=lunif))) +
  
  f(ID.year,
    model="rw1",
    constr = TRUE,
    hyper = list(prec=list(prior=sdunif))) +
  
  f(ID.area.year,
    model="generic0",
    Cmatrix = R,
    rankdef = r.def,
    constr = TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr = list(A=A.constr, e=rep(0, S+T-2)))


# Esperados = Casos esperados (variable definida dentro de Data)
# Modelo
ResMod4 <- inla(Formula,
                family="poisson",
                data=Data.INLA,
                E=Esperados,
                control.predictor = list(compute=TRUE, cdf=c(log(1))),
                control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE),
                control.inla = list(strategy=strategy))



ResMod4$cpu.used

ResMod4$dic$dic






############# Resultados ##############
# Resultados
round(ResMod4$summary.fixed, 4)

round(ResMod4$summary.hyperpar, 5)







############# Graficas de las distribuciones a posteriori del modelo ##############
# intercepto
ggplot(data.frame(data.frame(inla.smarginal(ResMod4$marginals.fixed[[1]]))), aes(x,y)) +
    geom_line() + ggtitle("INLA: Intercept") +
    theme_bw() + theme(plot.title = element_text(size=9.5))

kable(data.frame(INLAalpha = inla.emarginal(function(x) x, ResMod4$marginals.fixed$`(Intercept)`)))


# betas
# entre corchetes se indica el numero de cada covariable
# 1 = intercept, por eso no se pone
BetasINLA <- c(ResMod4$marginals.fixed[2], ResMod4$marginals.fixed[3], ResMod4$marginals.fixed[4], ResMod4$marginals.fixed[5], ResMod4$marginals.fixed[6], ResMod4$marginals.fixed[7])


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
MargSD <- inla.tmarginal(function(x) 1/sqrt(x),
                         ResMod4$marginals.hyperpar$`Precision for ID.area`)

MargSD2 <- inla.tmarginal(function(x) 1/sqrt(x),
                          ResMod4$marginals.hyperpar$`Precision for ID.year`)

MargSD3 <- inla.tmarginal(function(x) 1/sqrt(x),
                          ResMod4$marginals.hyperpar$`Precision for ID.area.year`)




# Obtener la media a posteriori de la desviacion tipica de los efectos espaciales de INLA
MeanMargSD <- inla.emarginal(function(x) x, MargSD)
MeanMargSD2 <- inla.emarginal(function(x) x, MargSD2)
MeanMargSD3 <- inla.emarginal(function(x) x, MargSD3)

# Comparar estadisticos basicos de las desviaciones tipicas de ambos modelos
SumMargSD <- inla.zmarginal(MargSD)
SumMargSD2 <- inla.zmarginal(MargSD2)
SumMargSD3 <- inla.zmarginal(MargSD3)

kable(data.frame(rbind(sigma.u=c(SumMargSD$mean, SumMargSD$sd, SumMargSD$quant0.025, SumMargSD$quant0.975), sigma.s=c(SumMargSD2$mean, SumMargSD2$sd, SumMargSD2$quant0.025, SumMargSD2$quant0.975), sigma.h=c(SumMargSD3$mean, SumMargSD3$sd, SumMargSD3$quant0.025, SumMargSD3$quant0.975))),col.names=c("Mean", "SD", "2.5%", "97.5%"))




# Summary
round(ResMod4$summary.fixed, 5)





# Representacion de las distribuciones a posteriori de las desviaciones tipicas de los efectos aleatorios espaciales del modelo de INLA
ggplot(data.frame(inla.smarginal(MargSD)), aes(x, y)) +
  geom_line() + ggtitle("INLA: Desviación de U") +
  theme_bw() + theme(plot.title = element_text(size=10.5))








############# Graficas RME INLA ##############
# Dibujo RME de un año en particular (i)
# S = numero de años (eg: si hay informacion desde 2011 hasta 2015, years = 5) (definido al principio)
# S = numero de municipios, barrios,... que se esten usando como unidad de vecindad (definido al principio)
# Map = cartografia (SpatialPolygonsDataFrame)
i <- 1 # (el numero del año que queramos representar)
RME <- matrix(ResMod4$summary.fitted.values$mean, nrow = S, ncol=T)
par(mar=c(0,0,2,0))
par(mfrow=c(1,1))
quantile(RME[,1], probs=seq(0,1,0.1))
Grupos <- c(0,0.7,0.9,1.1,1.5,2)
GruposRBYM <- findInterval(RME[,i], Grupos)
table(GruposRBYM)
Paleta2 <- brewer.pal(9, "BrBG")[9:1]
plot(Map, col=Paleta2[GruposRBYM], main = "RMEs INLA Año 1")
legend(727759.3,4376852, c("<0.7", "0.7-<0.9", "0.9-<1.1", "1.1-<1.5", "1.5-<2", ">=2"), 
       fill=Paleta2, border=NULL, bty="n", cex = 0.6)



# Mostrar regiones con RME >= 2
par(mar=c(0,0,2,0))
par(mfrow=c(2,4))
for (i in 1:T){
  
  RME <- matrix(ResMod4$summary.fitted.values$mean, nrow = S, ncol=T)
  pos <- which((RME)[,i] > 2)
  
  Grupos <- c(0,0.7,0.9,1.1,1.5,2)
  GruposRBYM <- findInterval(RME[,i], Grupos)
  Paleta2 <- brewer.pal(9, "Blues")[c(1,3,4,5,8,9)]
  plot(Map, col=Paleta2[GruposRBYM])
  legend(726659.3,4377922, c("<0.7", "0.7-<0.9", "0.9-<1.1", "1.1-<1.5", "1.5-<2", ">=2"), 
         fill=Paleta2, border=NULL, bty="n", cex = 0.9, y.intersp = 0.7)
  tit <- c("RME Año 1", "RME Año 2", "RME Año 3", "RME Año 4", "RME Año 5", "RME Año 6", "RME Año 7", "RME Año 8")
  title(paste(tit[i]), cex.main=2)
  
  for (j in 1:length(pos)) {
    lines(Map@polygons[[pos[j]]]@Polygons[[1]]@coords, col="red")
  }
  
}




############# Graficas p(RME) > 1 INLA ##############
# Dibujo p(RME) > 1 del primer año
# T = numero de años (eg: si hay informacion desde 2011 hasta 2015, years = 5) (definido al principio)
# S = numero de municipios, barrios,... que se esten usando como unidad de vecindad (definido al principio)
# Map = cartografia (SpatialPolygonsDataFrame)
ProbRME <- matrix(1 - ResMod4$summary.fitted.values$`1 cdf`, nrow = S, ncol=T)

GruposProb <- c(0,0.50,0.75,0.90,0.95)
GruposProbBYM <- findInterval(ProbRME[,1], GruposProb)
Paleta2 <- brewer.pal(5, "BrBG")[5:1]
par(mar=c(0,0,2,0))
par(mfrow=c(1,1))
plot(Map, col=Paleta2[GruposProbBYM], main = "p(RME > 1) INLA Año 1")
legend(727878,4376546,c("<0.50", "0.51-0.75","0.76-0.90", "0.91-0.95",">0.95"), 
       fill=Paleta2, border=NULL, bty="n", cex=0.8)




# Mostrar regiones con p(RME) > 0.95
par(mar=c(0,0,2,0))
par(mfrow=c(2,4))
for (i in 1:T){
  
  ProbRME <- matrix(1 - ResMod4$summary.fitted.values$`1 cdf`, nrow = S, ncol=T)
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


