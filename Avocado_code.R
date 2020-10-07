
# ========================================================================================================== #
# Title: Soil moisture and vegetation productivity across the green gold beltin central Mexico (2001-2018);
# Producers: Mario Guevara & Sebastian Gutierrez;
# Last uptade: October 2020;
# ========================================================================================================== #

rm(list=ls())


# =================== #
# Required libraries
# =================== #
pckg <- c('raster',         #Geographic data analysis and modelling
          'sp',             #Classes and methods for spatial data 
          'rgdal',          #Bindings for the 'Geospatial' data
          'soilassessment', #Assessment moel for agriculture soil conditions and crop suitability
          'magrittr',       #A forward-Pipe operator for R
          'rasterVis',      #Visualizaton methods for raster data
          'dichromat',      #Color schemes for dichromats
          'RColorBrewer')   #ColorBrewer palettes

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

lapply(pckg,usePackage)
getwd()
# ================== #
# Working directory
# ================== #
setwd("E:/IGAC2020/AVOC_MG/Avocado_Mexico_2020")

# ================== #
# Functions
# ================== #

# Function to normalize data
normalize <- function(x) {
  return ((x - min(x,na.rm=T)) / (max(x,na.rm=T) - min(x,na.rm=T)))
}


# ================== #
# Load files
# ================== #

## DEM derivatives
demder <- stack("E:/IGAC2020/AVOC_MG/terrain_parameters.tif")
names(demder) <- c("ASPECT","CAREA","CHNL_BASE",
                   "CHNK_DIST","CONVERGENCE","HCURV",
                   "ELEVATION","LSFACTOR","RSP","SHADE",
                   "SINKS","SLOPE","VALL_DEPTH","VCURV","TWI")
x11()
myTheme <- rasterTheme(region = brewer.pal(9,"Reds"))
levelplot(demder[[7]],par.settings = myTheme,main="Altitude")

## Study area boundary and municipalities
lim <- readOGR("muni_2018gw.shp")
plot(lim)

## Study area boundary
lim1 <- readOGR(".\\para_Sebastian\\muni_2018gw_diss.shp")
plot(lim1)

## Study area - avocado orchards
lim2 <- readOGR("Huertas_clipped.shp")
plot(lim2)

## Study area - water wells- polygon centroids
lim3 <- readOGR("points_wells.shp")
plot(lim3)

## Temperature
tmp <- raster("TEMPERATURE.tif")
myTheme <- rasterTheme(region = brewer.pal(9,"Reds"))
x11()
levelplot(tmp,par.settings = myTheme,main="Annual Mean Temperature")

## Precipitation
ppt <- raster("RAINFALL.tif")
myTheme <- rasterTheme(region = brewer.pal(9,"Reds"))
x11()
levelplot(ppt,par.settings = myTheme,main="Annual rainfall")

## Ecoregions level 3
ecor <- readOGR("ECOREGIONS_LEVELIII.shp")
sp::spplot(ecor,zcol="NA_L2CODE",names.attr="NA_L2CODE")

## Landsat 8 image
LS8 <- stack("E:\\IGAC2020\\AVOC_MG\\LS8.tif")
names(LS8)<- c("B1", "BLUE","GREEN","RED","NIR","SWIR1",
               "SWIR2","B8","B9","B10","B11","B12")

## NDVI
NDVI_clip <-raster("NDVI_2018_31072020.tif")

## BSI
BSI2_clip <- raster("BSI_2018_31072020.tif")

# ================== #
# Data processing
# ================== #

## Computing NDVI
NDVI <- (LS8$NIR - LS8$RED)/(LS8$NIR + LS8$RED) %>% 
  as("SpatialGridDataFrame") %>% 
  as("data.frame")
{
  NDVI$layer <- ifelse(NDVI$layer>1,NA,NDVI$layer)
  NDVI$layer <- ifelse(NDVI$layer< -1,NA,NDVI$layer)
  NDVI$layer <- ifelse(NDVI$layer< 0,NA,NDVI$layer)
  }

NDVI <- data.frame(x = NDVI$s1, y = NDVI$s2, 
                   NDVI = NDVI$layer)
gridded(NDVI) <- ~x+y
proj4string(NDVI) <- CRS("+init=epsg:4326")
NDVI <- raster(NDVI) %>% mask(lim1)
plot(NDVI)


## Computing bare soil index (BSI)
BSI2 <- (LS8$RED+LS8$BLUE)-(LS8$GREEN)/(LS8$RED+LS8$BLUE)+(LS8$GREEN)%>% 
  as("SpatialGridDataFrame") %>% 
  as("data.frame")
BSI2$layer <- ifelse(BSI2$layer=="-Inf",NA,BSI2$layer)
BSI2$layer <- ifelse(BSI2$layer<0,NA,BSI2$layer)
BSI2$layer <- ifelse(BSI2$layer>5000,NA,BSI2$layer)
# BSI2$layer<-normalize(BSI2$layer)
summary(BSI2$layer)
head(BSI2)
BSI2 <- data.frame(x = BSI2$s1, y = BSI2$s2, 
                   BSI2 = BSI2$layer)
gridded(BSI2) <- ~x+y
proj4string(BSI2) <- CRS("+init=epsg:4326")
BSI2 <- raster(BSI2)
plot(BSI2)

## Clipping indices
NDVI_clip <- mask(NDVI,lim1)
x11()
myTheme <- rasterTheme(region = dichromat(brewer.pal(9,"Blues")))
levelplot(NDVI_clip,par.settings = myTheme,main="NDVI")
writeRaster(NDVI_clip,"NDVI_2018_31072020.tif")
NDVI_clip <-raster("NDVI_2018_31072020.tif")

BSI2_clip <- mask(BSI2,lim1)
x11()
myTheme1 <- rasterTheme(region = dichromat(brewer.pal(9,"Blues")))
levelplot(BSI2_clip,par.settings = myTheme1,main="BSI")
writeRaster(BSI2_clip,"BSI_2018_31072020.tif")
BSI2_clip <- raster("BSI_2018_31072020.tif")


## Stacking BSI2_clip
INDC <- stack(NDVI_clip,BSI2)

## sOIL MOISTURE TIME SERIES
sm <- stack("SM_1km_kknn_annual_1979_2018.tif")
names(sm)
sm <- sm[[23:40]]
names(sm) <- paste0("SM_",2001:2018)
sm
sm <- raster::mask(sm,lim1)
sm_av <- calc(sm, mean)
sm_sd <- calc(sm, sd)
sm_desc <- stack(sm_av,sm_sd)
names(sm_desc) <- c("Mean_SM_time_series","Sd_SM_time_series")
x11()
levelplot(sm_desc)
myTheme <- rasterTheme(region = brewer.pal(9,"Blues"))
x11()
levelplot(sm_av,par.settings = myTheme,main="Mean SM time series")
x11()
levelplot(sm_sd,par.settings = myTheme,main="Sd SM time series")
writeRaster(sm_av,"Mean_SM_time_series.tif")
writeRaster(sm_sd,"Sd_SM_time_series.tif")

spplot(gpp)
an.sm <- c()
for (i in 1:dim(sm)[3]) {
  temp <- cellStats(sm[[i]],mean)
  an.sm <- c(an.sm,temp)
}
an.sm


## GPP TIME SERIES
gpp <- stack("GPP_TIMESERIES.tif")
names(gpp)
names(gpp) <- paste0("SM_",2001:2018)
gpp <- resample(gpp,sm,method="bilinear")
gpp <- raster::mask(gpp,lim1)
levelplot(gpp,layout=c(3,6))
gpp_av <- calc(gpp, mean)
gpp_sd <- calc(gpp, sd)
gpp_desc <- stack(gpp_av,gpp_sd)
names(gpp_desc) <- c("Mean_GPP_time_series","Sd_GPP_time_series")
myTheme <- rasterTheme(region = dichromat(brewer.pal(9,"Blues")))
x11()
levelplot(gpp_av,par.settings = myTheme,main="Mean GPP time series")
x11()
levelplot(gpp_sd,par.settings = myTheme,main="Sd GPP time series")
writeRaster(gpp_av,"Mean_GPP_time_series.tif")
writeRaster(gpp_sd,"Sd_GPP_time_series.tif")

spplot(gpp)
an.gpp <- c()
for (i in 1:dim(gpp)[3]) {
  temp <- cellStats(gpp[[i]],mean)
  an.gpp <- c(an.gpp,temp)
}
an.gpp

alma.data <- data.frame(Anio=2001:2018,SM.mean=an.sm,GP.mean=an.gpp)
View(alma.data)
write.csv(alma.data,"Datos_Dra_Alma_huertas.csv",row.names=F)

library(spatialEco)
cols <- colorRampPalette(brewer.pal(9,"RdGy"))
library(raster)

### Tendencia de productividad primaria bruta
t_gpp <- stack("stack_prob_gpp.tif")[[2]]
meanTendencias <- raster::cellStats(t_gpp, mean)
t_gpp_c <- (t_gpp - meanTendencias)
t_gpp_c <- resample(t_gpp_c,sm,method="bilinear")
values(t_gpp_c) <- ifelse(values(t_gpp_c)>220,220,values(t_gpp_c))
x11()
levelplot(t_gpp_c, par.settings = RdBuTheme,main="GPP trend")
writeRaster(t_gpp_c,"t_gpp_c.tif",overwrite=T)


## Probabilidad GPP
p_gpp <- stack("stack_prob_gpp.tif")[[3]]
p_gpp <- resample(p_gpp,sm,method="bilinear")
meanTendencias <- raster::cellStats(p_gpp, mean)
p_gpp_c <- p_gpp - meanTendencias
x11()
levelplot(p_gpp_c, par.settings = RdBuTheme)


## Probabilidad SM
p_sm <- stack("stack_SM.tif")[[3]]
p_sm <- mask(p_sm,lim1)
meanTendencias <- raster::cellStats(p_sm, mean)
p_sm_c <- p_sm - meanTendencias
x11()
levelplot(p_sm_c, par.settings = RdBuTheme)


## Tendencias SM
t_sm <- stack("stack_SM.tif")[[2]]
t_sm <- mask(t_sm,lim1)
meanTendencias <- raster::cellStats(t_sm, mean)
t_sm <- t_sm - meanTendencias
t_sm[t_sm>0.001] <- 0.001
t_sm[t_sm< -0.001] <- -0.001

x11()
levelplot(t_sm , par.settings = RdBuTheme,main="SM trend")
writeRaster(t_sm,"t_sm.tif")
t_sm
r_NDVI_clip <- resample(NDVI_clip,t_sm,method="bilinear")
r_BSI2_clip <- resample(BSI2_clip,t_sm,method="bilinear")
x11()
levelplot(p_gpp_c, par.settings = RdBuTheme)

#### Clasificacion de capas raster
library(raster)

### Tendencias de humedad t_sm
t_sm 
reclass_tsm <- c(-0.003020063, 0, 1,
                 0, 0.0025, 0)
reclass_tsm
reclass_tsm <- matrix(reclass_tsm,
                      ncol = 3,
                      byrow = TRUE)
reclass_tsm
tsm_classified <- reclassify(t_sm,
                             reclass_tsm)
spplot(tsm_classified)
tsm_classified<-as.factor(tsm_classified)
levels(tsm_classified)[[1]]




### Tendencias de proudctividad primaria t_gpp
(t_gpp <- stack("stack_prob_gpp.tif")[[2]])
reclass_tgpp <- c(-135.0083, 0, 1,
                  0, 681, 0)
reclass_tgpp
reclass_tgpp <- matrix(reclass_tgpp,
                       ncol = 3,
                       byrow = TRUE)
reclass_tgpp
tgpp_classified <- reclassify(t_gpp,
                              reclass_tgpp)
tgpp_classified <- resample(tgpp_classified,t_sm,method="ngb")
spplot(tgpp_classified)
tgpp_classified<-as.factor(tgpp_classified)
levels(tgpp_classified)[[1]]


### ïndice de vegetacion NDVI
NDVI_clip
cellStats(NDVI_clip,mean)
reclass_ndvi <- c(0,cellStats(NDVI_clip,mean), 1,
                  cellStats(NDVI_clip,mean), 1, 0)
reclass_ndvi
reclass_ndvi <- matrix(reclass_ndvi,
                       ncol = 3,
                       byrow = TRUE)
reclass_ndvi
ndvi_classified <- reclassify(NDVI_clip,
                              reclass_ndvi)
ndvi_classified <- resample(ndvi_classified,t_sm,method="ngb")
spplot(ndvi_classified)
ndvi_classified<-as.factor(ndvi_classified)
levels(ndvi_classified)[[1]]


### Indice de suelo desnudo BSI
BSI2_clip
cellStats(BSI2_clip,mean)
reclass_bsi <- c(0,cellStats(BSI2_clip,mean), 0,
                 cellStats(BSI2_clip,mean), 5000, 1)
reclass_bsi
reclass_bsi <- matrix(reclass_bsi,
                      ncol = 3,
                      byrow = TRUE)
reclass_bsi
bsi_classified <- reclassify(BSI2_clip,
                             reclass_bsi)
bsi_classified <- resample(bsi_classified,t_sm,method="ngb")
spplot(bsi_classified)
bsi_classified<-as.factor(bsi_classified)
levels(bsi_classified)[[1]]


## Final risk index computation
final <- bsi_classified + ndvi_classified + tgpp_classified + tsm_classified
final<-ratify(final)
(catg <- levels(final)[[1]])
catg$index <- c('Very Low','Low', 'Medium', 'High','Very High')
(levels(final) <- catg)
final
levelplot(final)

writeRaster(final,"CLASS_FINAL.tif",overwrite=T)

myTheme <- rasterTheme(region = dichromat(brewer.pal(5,"Oranges")))
x11()
levelplot(final,par.settings = myTheme,main="Mean GPP time series")
dev.off()




## Masking out urban zones and water bodies and zonal stats by municipality

## final --> raster
## lim --> municipalities polygons
## mask.ca --> water bodies mask 
## mask.zu --> urban zones mask
final <- raster("CLASS_FINAL.tif")
mask.ca <- readOGR("E:\\IGAC2020\\AVOC_MG\\mapa_base_cuerpos_agua_areas_urbanas\\C_A_Etapa_I_a_VI_2016actual.shp")
mask.zu <- readOGR("E:\\IGAC2020\\AVOC_MG\\mapa_base_cuerpos_agua_areas_urbanas\\agebur15gw.shp")
final1 <- final
final1 <- mask(final,mask.ca,inverse=T)
final1 <- mask(final1, mask.zu,inverse=T)
plot(final1)
myTheme <- rasterTheme(region = rev(brewer.pal(5,"RdYlGn")))
x11()
levelplot(final1,par.settings =myTheme,main="Degradation index")

#FIGURE INDEX + ORCHARDS (POLYGONS) + WELLS (POINTS)
wells <- readOGR("Ollas_centroides.shp")
myTheme <- rasterTheme(region = rev(brewer.pal(5,"RdYlGn")))
x11()
p1
p1 <- levelplot(final1,par.settings =myTheme,main="Indice de degradación")+
  layer(sp.polygons(lim2))+
  layer(sp.points(wells))
#end

final <- final1

(e <- raster::extract(final, lim))

class(e)
( class.counts <- lapply(e, table) ) 
( class.prop <- lapply(e, FUN = function(x) { prop.table(table(x)) }) ) 


class(class.prop)


rbind.fill <- function(x) {
  nam <- sapply(x, names)
  unam <- unique(unlist(nam))
  len <- sapply(x, length)
  out <- vector("list", length(len))
  for (i in seq_along(len)) {
    out[[i]] <- unname(x[[i]])[match(unam, nam[[i]])]
  }
  setNames(as.data.frame(do.call(rbind, out), stringsAsFactors=FALSE), unam)
}

( p.prop <- rbind.fill(class.prop) )
lim@data <- cbind(lim@data, p.prop) #add to polygons
lim@data #display data.frame in polygon object

names(lim)
lim@data$class_max <- colnames(lim@data)[10:14][apply(lim@data[,10:14],1,which.max)]
writeOGR(lim,layer="Indice_final", dsn="E:/IGAC2020/AVOC_MG",driver="ESRI Shapefile",overwrite_layer=T)
dim(lim@data)
write.csv(lim@data,"E:/IGAC2020/AVOC_MG/zonales_final.csv",row.names = F)


## Stacked barplot percentage
library(reshape)
library(readxl)
datos <- read_excel("E:\\IGAC2020\\AVOC_MG\\ZONAL_IND_FINAL.xlsx",sheet="DATA") %>% data.frame
head(datos)
datos <- datos[order(datos$High),]
datos1 <- as.matrix(datos[1:34,2:6])
colnames(datos1) <- c("Very Low","Low","Medium","High","Very High")
rownames(datos1) <- datos$NOM_MUN  
datos1 <- t(datos1)

x11()
par(oma = c(5,0,4,0))
p <- barplot(datos1, ylim=c(0,1),col = c("green","yellow","orange","red","darkred"),
             legend.text = rownames(datos1),las=2,pch=1,
             args.legend = list(x = "top", bty="n", inset=c(0,-0.17), xpd = T,ncol=5),
             names.arg = datos$NOM_MUN)
dev.off()




## Indicators validation boxplots
NDVI_clip
BSI2_clip
ind <- stack(NDVI_clip,
             BSI2_clip)
names(ind) <- c("NDVI","BSI")

final <- raster("CLASS_FINAL.tif")
final <- rasterToPoints(final) 
final <- data.frame(x=final[,1],y=final[,2],risk=final[,3])
final1 <- final
coordinates(final1) <- ~ x + y
str(final1)

final1
proy <- CRS("+proj=longlat +datum=WGS84")
final1@proj4string <- proy
final1 <- spTransform (final1, CRS=projection(t_sm))
final <- cbind(final, extract(ind, final1))
summary(final)


final1 <- final
coordinates(final1) <- ~ x + y
str(final1)
final1
proy <- CRS("+proj=longlat +datum=WGS84")
final1@proj4string <- proy
final1 <- spTransform (final1, CRS=projection(t_sm))
final <- cbind(final, extract(t_sm, final1))
summary(final)

final1 <- final
coordinates(final1) <- ~ x + y
str(final1)
final1
proy <- CRS("+proj=longlat +datum=WGS84")
final1@proj4string <- proy
final1 <- spTransform (final1, CRS=projection(t_sm))
final <- cbind(final, extract(t_gpp, final1))
summary(final)
names(final)[c(6,7)]<- c("t_SM","t_GPP")


final
str(final)
final$risk <- as.factor(final$risk)
x11()
par(mfrow=c(2,2))
boxplot(final$NDVI~final$risk,main="NDVI",xlab="Risk category",ylab="NDVI")
boxplot(final$BSI~final$risk,main="BSI",xlab="Risk category",ylab="BSI")
final$t_SM_norm <- normalize(final$t_SM)
boxplot(final$t_SM_norm~final$risk,main="t_SM",xlab="Risk category",ylab="Soil moisture trend")
boxplot(final$t_GPP~final$risk,main="t_GPP",xlab="Risk category",ylab="GPP trend")



### Zonal stats for avocado orchards areas
final <- raster("CLASS_FINAL.tif")
levelplot(final)
orch <- readOGR("Huertas_clipped.shp")
##final class raster only for avocado orchards area
final1 <- crop(final,orch)
final1 <- mask(final,orch)
levelplot(final1)
myTheme <- rasterTheme(region = rev(brewer.pal(5,"RdYlGn")))
x11()
levelplot(final1,par.settings =myTheme,main="Indice de degradación zona aguacatera")
writeRaster(final1,"IndicexHuertas.tif")
## Study area boundary and municipalities
myTheme <- rasterTheme(region = rev(brewer.pal(5,"RdYlGn")))
p1 <- levelplot(final1,par.settings =myTheme)+
  layer(sp.polygons(lim))
x11()
p1+ layer(sp.text(coordinates(lim),txt = lim$NAME_MUN,
                  pos = 1,font=2,cex=0.9,col="black"))


###zonal stats
#final <- raster("CLASS_FINAL.tif")
lim <- readOGR("E:\\IGAC2020\\AVOC_MG\\Huertas_Mpios.shp")
names(lim)
final <- final1

(e <- raster::extract(final, lim))

class(e)
( class.counts <- lapply(e, table) ) 
( class.prop <- lapply(e, FUN = function(x) { prop.table(table(x)) }) ) 


class(class.prop)


rbind.fill <- function(x) {
  nam <- sapply(x, names)
  unam <- unique(unlist(nam))
  len <- sapply(x, length)
  out <- vector("list", length(len))
  for (i in seq_along(len)) {
    out[[i]] <- unname(x[[i]])[match(unam, nam[[i]])]
  }
  setNames(as.data.frame(do.call(rbind, out), stringsAsFactors=FALSE), unam)
}

( p.prop <- rbind.fill(class.prop) )
lim@data <- cbind(lim@data, p.prop) #add to polygons
lim@data #display data.frame in polygon object

names(lim)
lim@data$class_max <- colnames(lim@data)[29:33][apply(lim@data[,29:33],1,which.max)]
writeOGR(lim,layer="Indice_final", dsn="E:/IGAC2020/AVOC_MG",driver="ESRI Shapefile",overwrite_layer=T)
head(lim@data)
write.csv(lim@data,"E:/IGAC2020/AVOC_MG/zonaleshuertas_final.csv",row.names = F)
plot(final)
writeRaster(final,"CLASS_FINAL_HUERTOS.tif",overwrite=T)

## Barplot percentage
library(reshape)
library(readxl)
library(magrittr)
datos <- read_excel("E:\\IGAC2020\\AVOC_MG\\ZONALES_HUERTAS_VF.xlsx",sheet="DATA") %>% data.frame
head(datos)
datos <- datos[order(datos$High),]
dim(datos)
datos1 <- as.matrix(datos[1:34,2:6])*100
colnames(datos1) <- c("Very Low","Low","Medium","High","Very High")
rownames(datos1) <- datos$Name 
datos1 <- t(datos1)

x11()
par(oma = c(5,0,4,0))
p <- barplot(datos1, ylim=c(0,100),col = c("green","yellow","orange","red","darkred"),
             legend.text = rownames(datos1),las=2,pch=1,ylab="Percentage",
             args.legend = list(x = "top", bty="n", inset=c(0,-0.17), xpd = T,ncol=5),
             names.arg = datos$Name)
dev.off()


install.packages("highcharter")
library(highcharter)
head(datos)
x11()
highchart() %>% 
  hc_yAxis_multiples(
    list(lineWidth = 3, lineColor='blue', title=list(text="Total area (ha)")),
    list(lineWidth = 3, lineColor="green", title=list(text="Avocado area (ha)"))
  ) %>% 
  hc_add_series(data = datos$Tot_Area, color='blue', type = "column") %>% 
  hc_add_series(data = datos$Avoc_area, color='green', type = "column", yAxis = 1) %>%
  hc_xAxis(categories = datos$Name, title = list(text = "Municipalities"))



### Zonal stats for water wells "Ollas de agua"
library(raster)
library(sp)
library(spatialEco)
library(rgdal)
final <- raster("CLASS_FINAL.tif")
values(final) <- round(values(final),0)

wells <- readOGR("Ollas_centroides.shp")
wells@data$point <- paste0("P",1:dim(wells@data)[1])
data <- data.frame(point=wells@data$point,
                   Longitude=wells@data$x_cent,
                   Latitude=wells@data$y_cent)
dim(data)
final <- projectRaster(final, crs= projection(wells))
plot(wells,add=T)
plot(final)

dat_subset_sp <- data
coordinates(dat_subset_sp) <- ~ Longitude + Latitude
str(dat_subset_sp)
dat_subset_sp
proy <- CRS("+proj=longlat +datum=WGS84")
dat_subset_sp@proj4string <- proy
dat_subset_sp <- spTransform (dat_subset_sp, CRS=projection(final))

library(spatialEco)
dat_subset_sp = point.in.poly(dat_subset_sp, lim)

start <- Sys.time()
#dat <- extract(COV84, dat_subset_sp, sp = TRUE)
data <- cbind(dat_subset_sp@data, class=extract(final, dat_subset_sp))
print(Sys.time() - start)
names(data)
head(data)
data$class <- factor(data$class)
levels(data$class) <- c("Very Low","Low","Medium","High","Very High") 
plot(data$class,ylab="Number of Wells",xlab="Category",ylim=c(0,5000),main="Wells per degradation level")
write.csv(data,"indx_wells.csv",row.names=F)

### Stacked barplot percentage showing municipalities degradation
library(reshape)
library(readxl)
datos <- read_excel("E:\\IGAC2020\\AVOC_MG\\indx_wells.xlsx",sheet="DATA") %>% data.frame
head(datos)
datos <- datos[order(datos$High),]
dim(datos)
datos1 <- as.matrix(datos[1:34,2:6])*100
colnames(datos1) <- c("Very Low","Low","Medium","High","Very High")
rownames(datos1) <- datos$Name  
datos1 <- t(datos1)

x11()
par(oma = c(5,0,4,0))
p <- barplot(datos1, ylim=c(0,100),col = c("green","yellow","orange","red","darkred"),
             legend.text = rownames(datos1),las=2,pch=1,ylab="Percentage",
             args.legend = list(x = "top", bty="n", inset=c(0,-0.17), xpd = T,ncol=5),
             names.arg = datos$Name)
dev.off()




### Zonal stats for water wells "Ollas de agua" inside avocado orchards
library(raster)
library(sp)
library(spatialEco)
library(rgdal)
final <- raster("CLASS_FINAL.tif")
values(final) <- round(values(final),0)

wells <- readOGR("Wells_IN_Orchards.shp")
names(wells)
levels(as.factor(wells$NOM_MUN))
library(spatialEco)

dat_subset_sp = point.in.poly(wells, lim2)

start <- Sys.time()
#dat <- extract(COV84, dat_subset_sp, sp = TRUE)
data <- cbind(wells@data, class=extract(final, wells))
print(Sys.time() - start)
names(data)
head(data)
data$class <- factor(data$class)
levels(data$class) <- c("Very Low","Low","Medium","High","Very High") 
plot(data$class,ylab="Number of Wells",xlab="Category",ylim=c(0,5000),main="Wells per degradation level")
write.csv(data,"indx_wells_orchards.csv",row.names=F)

### Stacked barplot percentage showing municipalities degradation
library(reshape)
library(readxl)
datos <- read_excel("E:\\IGAC2020\\AVOC_MG\\ZONALES_WELLSxORCHARDS.xlsx",sheet="DATA") %>% data.frame
head(datos)
datos <- datos[order(datos$High),]
dim(datos)
datos1 <- as.matrix(datos[1:30,2:6])*100
colnames(datos1) <- c("Very Low","Low","Medium","High","Very High")
rownames(datos1) <- datos$Name  
datos1 <- t(datos1)

{
  x11()
  par(oma = c(5,0,4,0))
  p <- barplot(datos1, ylim=c(0,100),col = c("green","yellow","orange","red","darkred"),
               legend.text = rownames(datos1),las=2,pch=1,ylab="Percentage",
               args.legend = list(x = "top", bty="n", inset=c(0,-0.17), xpd = T,ncol=5),
               names.arg = datos$Name)
}
dev.off()
install.packages("highcharter")
library(highcharter)
head(datos)
x11()
highchart() %>% 
  hc_yAxis_multiples(
    list(lineWidth = 3, lineColor='blue', title=list(text="Total area (ha)")),
    list(lineWidth = 3, lineColor="green", title=list(text="Avocado area (ha)"))
  ) %>% 
  hc_add_series(data = datos$Tot_Area, color='blue', type = "column") %>% 
  hc_add_series(data = datos$Avoc_area, color='green', type = "column", yAxis = 1) %>%
  hc_xAxis(categories = datos$Name, title = list(text = "Municipalities"))


