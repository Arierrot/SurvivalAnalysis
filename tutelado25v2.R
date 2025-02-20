#################################
###### PAQUETES NECESARIOS ######
#################################

library(survival)  # Para realizar Análisis de Supervivencia
library(survminer) # Para gráficos avanzados de Análisis de Supervivencia
library(corrplot) # Para matriz de correlaciones
library(car) # Para el cálculo de VIF
library(RcmdrMisc) # Para función numSummary

#################################
######### CARGAR DATOS ##########
#################################

setwd("D:/MUBICS/Estadistica/Tutelado")

data <- read.csv("D:/MUBICS/Estadistica/Tutelado/data/BC_cardiotox_clinical_variables.csv", header = TRUE, sep = ";", dec = ",")
#head(data)

##### VARIABLES DE ESTUDIO ######
#'* CTRCD: 1 = cardiotoxicidad, 0 = no (censura) *
#'* time: tiempo en días desde el inicio del tratamiento hasta el evento o fin del seguimiento *
#'* LVEF: fracción de eyección ventricular izquierda (%) *
#'* heart_rate: frecuencia cardíaca (lpm) *
#'* PWT: grosor de la pared posterior del ventrículo izquierdo (cm) *
#'* LAd: diámetro de la aurícula izquierda (cm) *
#'* LVDd: diámetro diastólico del ventrículo izquierdo (cm) *
#'* LVSd: diámetro sistólico del ventrículo izquierdo (cm) *
#'* AC: tratamiento con antraciclinas (1 = sí, 0 = no) *


##########################
##### PREPARAR DATOS #####
##########################

# nos quedamos con las columnas heart_rate, CTRCD, time, LVEF, PWT, LAd, LVDd, LVSd y AC
data = data[,c(1, 5:6, 8:13)]

# transformamos el tipo de dato las variables cuantitativas a numéricas
for (i in (1:length(data))){
  data[,i] <- as.numeric(data[,i])
}

# ordenamos el dataset en función del tiempo
data <- data[order(data[, 3]),]
#data

# summary de nuestros datos cuantitativos
summary(data) #'*Observamos que la columna 'AC' presenta un buen número de valores faltantes (NA)*

porcentaje_NA <- sum(is.na(data)) / nrow(data) * 100
print(porcentaje_NA) #'* 10.73% de valores faltantes -> inputarlos o eliminarlos?*

#Decido eliminar los NAs, ya que no suponen un gran porcentaje del dataset
data <- na.omit(data)

summary(data) #'* ya no hay valores faltantes*


#################################
##### ANÁLISIS EXPLORATORIO #####
#################################

head(data)
str(data)

table(data$CTRCD) # tabla de frecuencias para CTRCD
#'* TENEMOS MUY POCOS CASOS DE CTRCD=1 (47) EN COMPARACIÓN CON CTRCD=0 (430) *

# Medidas descriptivas para variables cuantitativas:
var_num = data[,c(1, 4:8)]
var_num
numSummary(var_num,
           statistics = c("mean", "sd", "IQR", "quantiles", "cv", "skewness", "kurtosis"),
           quantiles = c(0, 0.25, 0.5, 0.75, 1), type = "2")

# Histogramas para observar distribución de variables cuantitativas
windows()
par(mfrow=c(2, 3))
hist(data$heart_rate, main = "Distribución de Heart Rate", xlab = "Heart Rate")
#'* heart_rate: Observamos cierta asimetría hacias la derecha*
hist(data$LVEF, main = "Distribución de LVEF", xlab = "LVEF")
#'* LVEF: Observamos cierta asimetría hacias la izquierda*
hist(data$PWT, main = "Distribución de PWT", xlab = "PWT")
#'* PWT: Observamos algún valor extremo*
hist(data$LAd, main = "Distribución de LAd", xlab = "LAd")
hist(data$LVDd, main = "Distribución de LVDd", xlab = "LVDd")
hist(data$LVSd, main = "Distribución de LVSd", xlab = "LVSd")
#'* El resto parecen presentar una distribución normal*


# Posibles valores atípicos?
boxplot(data$heart_rate, main = "Boxplot de Heart Rate")
boxplot(data$LVEF, main = "Boxplot de LVEF")
boxplot(data$PWT, main = "Boxplot de PWT")
boxplot(data$LAd, main = "Boxplot de LAd")
boxplot(data$LVDd, main = "Boxplot de LVDd")
boxplot(data$LVSd, main = "Boxplot de LVSd")
#'* Parecen presentar valores atípicos *


# Comprobamos que no haya problemas de multicolinealidad entre las variables independientes
# Matriz de correlación
windows()
corr_matrix <- cor(var_num[, sapply(var_num, is.numeric)], use = "complete.obs")
corrplot.mixed(corr_matrix)
#'* Observamos cierta correlación negativa entre LVSd y LVEF*
#'* Cierta correlación positiva entre LVSd y LVDd*

# Como se comporta cada variable independiente frente a la respuesta?
windows()
par(mfrow=c(2, 3))
boxplot(heart_rate ~ CTRCD, data = data, main = "heart_rate vs CTRCD", xlab = "CTRCD", ylab = "heart_rate")
boxplot(LVEF ~ CTRCD, data = data, main = "LVEF vs CTRCD", xlab = "CTRCD", ylab = "LVEF")
boxplot(PWT ~ CTRCD, data = data, main = "PWT vs CTRCD", xlab = "CTRCD", ylab = "PWT")
boxplot(LAd ~ CTRCD, data = data, main = "LAd vs CTRCD", xlab = "CTRCD", ylab = "LAd")
boxplot(LVDd ~ CTRCD, data = data, main = "LVDd vs CTRCD", xlab = "CTRCD", ylab = "LVDd")
boxplot(LVSd ~ CTRCD, data = data, main = "LVSd vs CTRCD", xlab = "CTRCD", ylab = "LVSd")
#'* Las variables más influyentes sobre CTRCD parecen ser heart_rate y LVEF*


# Pie-chart para 'AC'
# Crear una tabla de frecuencias cruzadas entre AC y CTRCD
tabla_AC_CTRCD <- table(data$AC, data$CTRCD)
tabla_AC_CTRCD <- frec(tabla_AC_CTRCD)
windows()
par(mfrow = c(1,2))
for (i in 1:length(unique(data$AC))) {
  categoria <- unique(data$AC)[i]
  counts <- tabla_AC_CTRCD[i, ]
  pie(counts, main = paste("AC =", categoria), labels = c("CTRCD = 0", "CTRCD = 1"), col = c("lightgreen", "pink"))
}
#'* Parece que existe una ligera influencia de AC*

table(data$CTRCD)
#'* Observamos que algunas pacientes no han sufrido el evento de interés (CTRCD) en el momento *
#'* de finalización del estudio (Por lo que no sabemos si en algún momento lo sufrirán), esto es *
#'* lo que se conoce como una censura por la derecha *


######################################################
###### AJUSTAMOS UN MODELO DE REGRESIÓN DE COX #######
######       ( ANÁLISIS DE SUPERVICENCIA )     #######             
######################################################

#'* Debido a que existe una censura, la mejor aproximación es ajustar un modelo de regresión de Cox*

# GRÁFICO DE FUNCIÓN DE SUPERVIVENCIA
plot(survfit(Surv(data$time, data$CTRCD) ~ 1), xlab = "Tiempo de Seguimiento", ylab = "Probabilidad de Supervivencia", main = "Curva de supervivencia", conf.int = 0.95)
#'* Linea continua: probabilidad de supervivencia acumulativa en función del tiempo *
#'* Lineas punteadas: IC a nivel de confianza = 0.95 *
# otra opción
windows()
ggsurvplot(survfit(Surv(data$time, data$CTRCD) ~ 1), data = data, conf.int = TRUE,
           xlab = "Tiempo de Seguimiento",
           ylab = "Probabilidad de Supervivencia",
)

###############################
##### MODELOS UNIVARIADOS #####
###############################

##### VARIABLE AC #####
#'* VARIABLE CATEGORICA, el resto son cuantitativas *

# Función de supervivencia (Kaplan-Meier) para el tiempo hasta la aparicion de cardiotoxicidad
AC_survfit <- survfit(Surv(time, CTRCD) ~ AC, data = data, conf.int = 0.95)
AC_survfit

plot(AC_survfit, col = c("red", "blue"), lty = 1, xlab = "Time", ylab ="Fracción de pacientes sin cardiotoxicidad")
legend("bottomleft", legend = c("Ausencia tratamiento AC","Presencia Tratamiento AC"), lty = 1, col = c("red", "blue"))
#'* Las pacientes tratadas con antraciclina parecen presentar menor supervivencia a corto plazo *
#'* luego parece que se estabiliza * 

# Función de riesgos
ggsurvplot(
  AC_survfit,
  fun = "cumhaz",              # Función de riesgo acumulado
  conf.int = TRUE,             # Intervalo de confianza
  legend.title = "Grupo AC",   # Título de la leyenda
  legend.labs = c("AC = 0", "AC = 1"), # Etiquetas de los grupos
  xlab = "Tiempo", 
  ylab = "Riesgo acumulado (CTRCD)", 
  ggtheme = theme_minimal(),
  data = data
)

# Modelo univariado para AC
modelo_AC <- coxph(Surv(time, CTRCD) ~ AC, data = data)
summary(modelo_AC, conf.int=0.95) #'* p-value = 0.055 *
plot(survfit(modelo_AC), xlab = "Tiempo de Seguimiento", ylab = "Probabilidad de Supervivencia", main = "Curvas de Supervivencia según AC", conf.int = 0.95)

survdiff(Surv(time, CTRCD) ~ AC, data = data) # chi-cuadrado
#'* p-value = 0.05 -> existen diferencias significativas entre los dos niveles del factor *


##### RESTO DE VARIABLES #####
#'* Ajustamos curvas de Kaplan-Meier (exploratorio) y modelos de Cox para cada variable *
# Curvas de Kaplan-Meier para var.cuantitativas: dividiendo en alto y bajo basado en la mediana

### HEART_RATE
# Funcion de supervivencia para heart_rate
data$heart_rate_group <- ifelse(data$heart_rate > median(data$heart_rate, na.rm = TRUE), "Alto", "Bajo")
heart_rate_survfit <- survfit(Surv(time, CTRCD) ~ heart_rate_group, data = data)
plot(heart_rate_survfit, col = c("red", "blue"), lty = 1, xlab = "Time", ylab ="Fracción de pacientes sin cardiotoxicidad")
legend("bottomleft", legend = c("Alto heart_rate","Bajo heart_rate"), lty = 1, col = c("red", "blue"))
# Modelo univariado para heart_rate
modelo_heart_rate <- coxph(Surv(time, CTRCD) ~ heart_rate, data = data)
summary(modelo_heart_rate, conf.int=0.95) #'* p-value = 0.000589 *
# Función de riesgos
ggsurvplot(
  heart_rate_survfit,
  fun = "cumhaz",              # Función de riesgo acumulado
  conf.int = TRUE,             # Intervalo de confianza
  legend.title = "heart_rate estratificado",   # Título de la leyenda
  legend.labs = c("alto", "bajo"), # Etiquetas de los grupos
  xlab = "Tiempo", 
  ylab = "Riesgo acumulado (CTRCD)", 
  ggtheme = theme_minimal(),
  data = data
)


### LVEF
# Funcion de supervivencia para LVEF
data$LVEF_group <- ifelse(data$LVEF > median(data$LVEF, na.rm = TRUE), "Alto", "Bajo")
LVEF_survfit <- survfit(Surv(time, CTRCD) ~ LVEF_group, data = data)
plot(LVEF_survfit, col = c("red", "blue"), lty = 1, xlab = "Time", ylab ="Fracción de pacientes sin cardiotoxicidad")
legend("bottomleft", legend = c("Alto LVEF","Bajo LVEF"), lty = 1, col = c("red", "blue"))
# Modelo univariado para LVEF
modelo_LVEF <- coxph(Surv(time, CTRCD) ~ LVEF, data = data)
summary(modelo_LVEF, conf.int=0.95) #'* p-value = 0.0192 *
# Función de riesgos
ggsurvplot(
  LVEF_survfit,
  fun = "cumhaz",              # Función de riesgo acumulado
  conf.int = TRUE,             # Intervalo de confianza
  legend.title = "LVEF estratificado",   # Título de la leyenda
  legend.labs = c("alto", "bajo"), # Etiquetas de los grupos
  xlab = "Tiempo", 
  ylab = "Riesgo acumulado (CTRCD)", 
  ggtheme = theme_minimal(),
  data = data
)


### PWT
# Funcion de supervivencia para PWT
data$PWT_group <- ifelse(data$PWT > median(data$PWT, na.rm = TRUE), "Alto", "Bajo")
PWT_survfit <- survfit(Surv(time, CTRCD) ~ PWT_group, data = data)
plot(PWT_survfit, col = c("red", "blue"), lty = 1, xlab = "Time", ylab ="Fracción de pacientes sin cardiotoxicidad")
legend("bottomleft", legend = c("Alto PWT","Bajo PWT"), lty = 1, col = c("red", "blue"))
# Modelo univariado para PWT
modelo_PWT <- coxph(Surv(time, CTRCD) ~ PWT, data = data)
summary(modelo_PWT, conf.int=0.95) #'* p-value = 0.889 *
# Función de riesgos
ggsurvplot(
  PWT_survfit,
  fun = "cumhaz",              # Función de riesgo acumulado
  conf.int = TRUE,             # Intervalo de confianza
  legend.title = "PWT estratificado",   # Título de la leyenda
  legend.labs = c("alto", "bajo"), # Etiquetas de los grupos
  xlab = "Tiempo", 
  ylab = "Riesgo acumulado (CTRCD)", 
  ggtheme = theme_minimal(),
  data = data
)


### LAd
# Funcion de supervivencia para LAd
data$LAd_group <- ifelse(data$LAd > median(data$LAd, na.rm = TRUE), "Alto", "Bajo")
LAd_survfit <- survfit(Surv(time, CTRCD) ~ LAd_group, data = data)
plot(LAd_survfit, col = c("red", "blue") , lty = 1, xlab = "Time", ylab ="Fracción de pacientes sin cardiotoxicidad")
legend("bottomleft", legend = c("Alto LAd","Bajo LAd"), lty = 1, col = c("red", "blue"))
# Modelo univariado para LAd
modelo_LAd <- coxph(Surv(time, CTRCD) ~ LAd, data = data)
summary(modelo_LAd, conf.int=0.95) #'* p-value = 0.855 *
# Función de riesgos
ggsurvplot(
  LAd_survfit,
  fun = "cumhaz",              # Función de riesgo acumulado
  conf.int = TRUE,             # Intervalo de confianza
  legend.title = "LAd estratificado",   # Título de la leyenda
  legend.labs = c("alto", "bajo"), # Etiquetas de los grupos
  xlab = "Tiempo", 
  ylab = "Riesgo acumulado (CTRCD)", 
  ggtheme = theme_minimal(),
  data = data
)


### LVDd
# Funcion de supervivencia para LVDd
data$LVDd_group <- ifelse(data$LVDd > median(data$LVDd, na.rm = TRUE), "Alto", "Bajo")
LVDd_survfit <- survfit(Surv(time, CTRCD) ~ LVDd_group, data = data)
plot(LVDd_survfit, col = c("red", "blue"), lty = 1, xlab = "Time", ylab ="Fracción de pacientes sin cardiotoxicidad")
legend("bottomleft", legend = c("Alto LVDd","Bajo LVDd"), lty = 1, col = c("red", "blue"))
# Modelo univariado para LVDd
modelo_LVDd <- coxph(Surv(time, CTRCD) ~ LVDd, data = data)
summary(modelo_LVDd, conf.int=0.95) #'* p-value = 0.866 *
# Función de riesgos
ggsurvplot(
  LVDd_survfit,
  fun = "cumhaz",              # Función de riesgo acumulado
  conf.int = TRUE,             # Intervalo de confianza
  legend.title = "LVDd estratificado",   # Título de la leyenda
  legend.labs = c("alto", "bajo"), # Etiquetas de los grupos
  xlab = "Tiempo", 
  ylab = "Riesgo acumulado (CTRCD)", 
  ggtheme = theme_minimal(),
  data = data
)


### LVSd
# Funcion de supervivencia para LVSd
data$LVSd_group <- ifelse(data$LVSd > median(data$LVSd, na.rm = TRUE), "Alto", "Bajo")
LVSd_survfit <- survfit(Surv(time, CTRCD) ~ LVSd_group, data = data)
plot(LVSd_survfit, col = c("red", "blue"), lty = 1, xlab = "Time", ylab ="Fracción de pacientes sin cardiotoxicidad")
legend("bottomleft", legend = c("Alto LVSd","Bajo LVSd"), lty = 1, col = c("red", "blue"))
# Modelo univariado para LVSd
modelo_LVSd <- coxph(Surv(time, CTRCD) ~ LVSd, data = data)
summary(modelo_LVSd, conf.int=0.95) #'* p-value = 0.208 *
# Función de riesgos
ggsurvplot(
  LVSd_survfit,
  fun = "cumhaz",              # Función de riesgo acumulado
  conf.int = TRUE,             # Intervalo de confianza
  legend.title = "LVSd estratificado",   # Título de la leyenda
  legend.labs = c("alto", "bajo"), # Etiquetas de los grupos
  xlab = "Tiempo", 
  ylab = "Riesgo acumulado (CTRCD)", 
  ggtheme = theme_minimal(),
  data = data
)


###############################
#### MODELOS MULTIVARIADOS ####
###############################

#'* MODELO CON TODAS LAS VARIABLES *
#'* decido ajustar un modelo con todas las variables y voy excluyendo aquellas no significativas con mayor p-valor *

Cox1 <- coxph(Surv(time,CTRCD) ~ AC + heart_rate + LVEF + PWT + LAd + LVDd + LVSd, method = "efron", data = data)
summary(Cox1, conf.int=0.95)
#                 coef exp(coef)  se(coef)      z Pr(>|z|)    
# AC          0.532281  1.702813  0.363408  1.465 0.143005    
# heart_rate  0.033023  1.033574  0.009691  3.408 0.000655 ***
# LVEF       -0.028754  0.971655  0.036764 -0.782 0.434140    
# PWT         0.025202  1.025522  1.098419  0.023 0.981695    
# LAd         0.107475  1.113463  0.319500  0.336 0.736581    
# LVDd       -0.328986  0.719653  0.666962 -0.493 0.621828    
# LVSd        0.525686  1.691620  0.954193  0.551 0.581687    
Cox2 <- coxph(Surv(time,CTRCD) ~ AC + heart_rate + LVEF + LAd + LVDd + LVSd, method = "efron", data = data)
summary(Cox2, conf.int=0.95)
#                 coef exp(coef)  se(coef)      z Pr(>|z|)    
# AC          0.532307  1.702856  0.363417  1.465 0.142995    
# heart_rate  0.033011  1.033562  0.009677  3.411 0.000647 ***
# LVEF       -0.028725  0.971684  0.036755 -0.782 0.434491    
# LAd         0.109833  1.116092  0.302470  0.363 0.716515    
# LVDd       -0.329916  0.718984  0.665980 -0.495 0.620328    
# LVSd        0.525098  1.690625  0.954350  0.550 0.582172 
Cox3 <- coxph(Surv(time,CTRCD) ~ AC + heart_rate + LVEF + LVDd + LVSd, method = "efron", data = data)
summary(Cox3, conf.int=0.95)
#                 coef exp(coef)  se(coef)      z Pr(>|z|)    
# AC          0.524127  1.688984  0.362367  1.446 0.148065    
# heart_rate  0.033171  1.033727  0.009681  3.426 0.000612 ***
# LVEF       -0.027830  0.972554  0.036999 -0.752 0.451942    
# LVDd       -0.325469  0.722189  0.675753 -0.482 0.630063    
# LVSd        0.556184  1.744004  0.965742  0.576 0.564674    
Cox4 <- coxph(Surv(time,CTRCD) ~ AC + heart_rate + LVEF + LVSd, method = "efron", data = data)
summary(Cox4, conf.int=0.95)
#                 coef exp(coef)  se(coef)      z Pr(>|z|)    
# AC          0.526713  1.693357  0.362409  1.453 0.146122    
# heart_rate  0.033288  1.033848  0.009646  3.451 0.000558 ***
# LVEF       -0.040733  0.960086  0.024980 -1.631 0.102969    
# LVSd        0.134139  1.143552  0.405904  0.330 0.741044    
Cox5 <- coxph(Surv(time,CTRCD) ~ AC + heart_rate + LVEF, method = "efron", data = data)
summary(Cox5, conf.int=0.95)
#                 coef exp(coef) se(coef)      z Pr(>|z|)    
# AC          0.52355   1.68802  0.36218  1.446 0.148297    
# heart_rate  0.03275   1.03329  0.00951  3.444 0.000573 ***
# LVEF       -0.04587   0.95516  0.01946 -2.358 0.018382 *  
Cox6 <- coxph(Surv(time,CTRCD) ~ heart_rate + LVEF, method = "efron", data = data)
summary(Cox6, conf.int=0.95)
#                 coef exp(coef)  se(coef)      z Pr(>|z|)    
# heart_rate  0.034789  1.035401  0.009546  3.644 0.000268 ***
# LVEF       -0.049595  0.951615  0.019360 -2.562 0.010417 *  
  
#'* INTERVALOS DE CONFIANZA *

#               exp(coef) exp(-coef) lower .95 upper .95
# heart_rate    1.0354     0.9658    1.0162    1.0550
# LVEF          0.9516     1.0508    0.9162    0.9884



## MODELO CON INTERACCION
Cox7 <- coxph(Surv(time,CTRCD) ~ heart_rate * LVEF, method = "efron", data = data)
summary(Cox7, conf.int=0.95)
#'* interaccion no significativa *


######################################################
################# BONDAD DE AJUSTE ###################
######################################################

# CONCORDANCIA -> mide qué tan bien el modelo predice el orden de los eventos
# C = 0.5  significa nula capacidad predictiva (equivalente al azar)
# C = 1    significa que el modelo predice completamente el orden de los eventos
summary(Cox6)$concordance
#'* C = 0.6937 -> concordancia aceptable, el modelo tiene una precisión del 69.37% *
summary(Cox5)$concordance
#'* C = 0.7166 -> mejora la concordancia al añadir AC, PERO ESTA NO ES SIGNIFICATIVA *


######################################################
######## COMPROBACION HIPOTESIS ESTRUCTURALES ########
######################################################

#'* SON LOS RIESGOS PROPORCIONALES? *

# Grafico de la estimación de los coeficientes beta(t)
windows()
par(mfrow = c(1, 2)) # Dividir el área gráfica en dos columnas
# Gráfico para heart_rate
plot(cox.zph(Cox6)[1], main = "Beta(t) for heart_rate")
abline(h = Cox6$coef["heart_rate"], col = "blue", lwd = 2)
# Gráfico para LVEF
plot(cox.zph(Cox6)[2], main = "Beta(t) for LVEF")
abline(h = Cox6$coef["LVEF"], col = "blue", lwd = 2)
par(mfrow = c(1, 1)) # Restaurar configuración gráfica por defecto

#'* Ambos gráficos muestran que los coeficientes beta(t) para heart_rate y LVEF son aproximadamente  *
#'* constantes en el tiempo, cumpliendo el supuesto de riesgos proporcionales.                       *
#'* Esto indica que el modelo de Cox ajustado es adecuado para estas variables sin necesidad de      *
#'* ajustes adicionales.                                                                             *

cox.zph(Cox6)
#         chisq df    p
#heart_rate 0.114  1 0.74
#LVEF       1.875  1 0.17
#GLOBAL     1.987  2 0.37
#'* pvalor=0.37 Se puede asumir riesgos proporcionales *


#'*##### MULTICOLINEALIDAD #####*

# Calcular el VIF basado en un modelo de regresión lineal
lm6 <- lm(CTRCD ~ heart_rate + LVEF, data = data)
vif_values <- vif(lm6)
print(vif_values)
# heart_rate    LVEF 
# 1.001606      1.001606 

#'* valores muy bajos, no hay multicolinealidad entre las variables *


#'*##### detectar VALORES ATIPICOS #####*

# residuos de Martingale -> útiles para detectar valores atípicos en el tiempo de supervivencia
#'* Observaciones con residuos muy cercanos a -1 pueden indicar casos atípicos*
martingale_res <- residuals(Cox6, type = "martingale")
plot(martingale_res, main = "Residuos de Martingale", ylab = "Residuos de Martingale", xlab = "Observaciones")
abline(h = c(-3, 3), col = "red", lty = 2)

# Residuos deviance (transformados a partir de los residuos de Martingale) -> mas robusto para detectar outliers
deviance_res <- residuals(Cox6, type = "deviance")
plot(deviance_res, main = "Residuos de Deviance", ylab = "Residuos de Deviance", xlab = "Observaciones")
abline(h = c(-3, 3), col = "red", lty = 2)

outliers <- which(abs(deviance_res) > 3)
print(outliers)
# named integer(0)
#'* ningún residuo clasificado como outlier "residuo deviance > 3"

#'* no se observan valores problematicos *
