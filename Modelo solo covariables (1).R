# Modelo solo con covariables
############# Preparacion ##############
# Librerias
library(spdep) #poly2nb
library(INLA)
library(RColorBrewer) #brewer.pal
library(ggplot2)
library(gridExtra)
library(knitr)

# Datos
# cargar datos y llamarlos Data
# Data <- load()



############# Formula y Modelo ##############
# V1,..., Vn son las N covariables (definidas dentro de Data)
# R es la variable respuesta (definida dentro de Data) (casos observados de la enfermedad)
# Rellenar el campo de formula con la variable respuesta y las variables que se quieran estudiar

# Formula
Formula <- R ~ V1 + V2 + V3 + V4 + V5 + V6 

# Esperados = Casos esperados (variable definida dentro de Data)
# Modelo
ResMod1 <- inla(Formula, 
                data=Data,
                E=Esperados, 
                family="poisson",
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                control.predictor = list(compute=TRUE, cdf=c(log(1))),
                control.fixed = list(mean.intercept = 0, prec.intercept = 0, mean = 0, prec = 0), 
                verbose = TRUE)




############# Resultados ##############

# Resultados
round(ResMod1$summary.fixed, 5)

ResMod1$dic$dic

ResMod1$cpu.used




############# Graficas de las distribuciones a posteriori del modelo ##############

# intercepto
ggplot(data.frame(data.frame(inla.smarginal(ResMod1$marginals.fixed[[1]]))), aes(x,y)) +
    geom_line() + ggtitle("INLA: Intercept") +
    theme_bw() + theme(plot.title = element_text(size=9.5))

kable(data.frame(INLAalpha = inla.emarginal(function(x) x, ResMod1$marginals.fixed$`(Intercept)`)))



# betas 
# entre corchetes se indica el numero de cada covariable
# 1 = intercept, por eso no se pone
BetasINLA <- c(ResMod1$marginals.fixed[2], ResMod1$marginals.fixed[3], ResMod1$marginals.fixed[4], ResMod1$marginals.fixed[5], ResMod1$marginals.fixed[6], ResMod1$marginals.fixed[7])


# n es el numero de covariables
n <- 6
par(mfrow=c(2,4))
par(mar=c(3,3,3,3))
for (i in 1:n){
  plot(inla.smarginal(BetasINLA[[i]]), type="l", xlab="", ylab="")
  tit <- c("INLA: V1", "INLA: V2", "INLA: V3", "INLA: V4", "INLA V5", "INLA V6")
  title(tit[i], cex.main=1)
}  



