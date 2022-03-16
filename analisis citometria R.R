##
#Analisis de datos de citometria de flujo de celulas S2
#fecha de creacion: 2021-08-18
#autor: Daniel Prieto (dprieto@fcien.edu.uy) y Gonzalez Ana Clara (anaclgonzalez92@gmail.com)
#institucion: Instituto de Investigaciones Biologicas Clemente Estable, Uruguay
##
library(flowCore)
library(flowDensity)
library(flowViz)
library(flowStats)
library(ggcyto)
#Abro archivos
files.test1 = list.files("D:/Maestria neurodesarrollo/fotos y resultados ensayos/archivos citometro nov/2021-11-17 S2 con y sin ind/", all.files = F, full.names = TRUE)#Lee todos los archivos del directorio
fs <- read.flowSet(files.test1[1:8], transformation = F, alter.names =T, emptyValue=FALSE)#[3:25])#Elimino las dos primeras filas que son /. y /..
#Transformo logicle
bx <- logicleTransform()
bxlist <- transformList(c("FSC.A", "SSC.A", "FSC.H", "BL1.A", "VL2.A", "VL1.A", "RL1.A"), bx)
datostrans <- transform(fs, bxlist)
##
#Limpiamos el dataset de debris (fcs/ssc)
clean.dt <- datostrans
for (i in 1:8) { # Loop over the length of the flowSet
  f <- datostrans[[i]]
  # First restrict the FSC - A values :
  fsc.indices <- intersect (which (exprs (f)[, "FSC.A"]<5.5) , which (exprs (f) [, "FSC.A"] > 3.5))
  # Then restrict SSC - A values and intersect with FSC - A restriction above :
  ssc.indices <- intersect (which ( exprs ( f)[, "SSC.A"] > 3) ,
                            which ( exprs (f)[, "SSC.A"] < 5.5) )
  non.debris.indices <- intersect ( fsc.indices , ssc.indices )
  # Now only select the non debris cells and place the cleaned up flowFrame into the flowSet :
  f.clean <- f[non.debris.indices]
  clean.dt [[i]] <- f.clean
}
##
#Muestro los plots
plotDens (clean.dt[[4]], c("FSC.A" ,"SSC.A"))
plotDens (clean.dt[[5]], c ("VL2.A", "SSC.A"))
plotDens (clean.dt[[3]], c ("FSC.A", "FSC.H"))
autoplot(clean.dt[[4]]) #+ labs_cyto("marker")
##
#Matriz de compensacion
#library('flowStats')
#comp.files1 = list.files("~/Documentos/Posdoc_IIBCE/Citometro/", all.files = T, full.names = TRUE)#Lee todos los archivos del directorio
#comp.fs <- read.flowSet(comp.files1[3:5], transformation = F, alter.names = T)#Elimino las dos primeras filas que son /. y /..
##
#cx <- logicleTransform()
#cxlist <- transformList(c("FSC.A", "SSC.A", "FL1.A", "FL2.A"), cx)
#comp.fs <- transform(comp.fs, cxlist)
##
#namepatt <- "FL2.A|FL1.A"
#comp.mat <- spillover(comp.fs, unstained = sampleNames(comp.fs[3]),fsc = "FSC.A",ssc = "SSC.A",patt = namepatt, method = "mean",stain_match = "ordered",useNormFilt = TRUE,prematched = FALSE,exact_match = FALSE)
#write.table(comp.mat, file = "~/Documentos/Posdoc_IIBCE/Citometro/comp_mat_median.csv")#Mando la matriz a un archivo
##
#Compensacion
#comp.mat <- read.csv(file = "~/Documentos/Posdoc_IIBCE/Citometro/comp_mat_median.csv", sep = " ")#uso la matriz guardada para no calcular cada vez
## create a compensation object
#comp <- compensation(comp.mat,compensationId="comp1")
#compensados <- compensate(clean.dt, comp)
##
#Ahora un poco de gating y filtrado sobre todo el flowSet.oO
#Gate elipsoide tipo linfocitos
linfos <- lymphGate(clean.dt, channels=c("FSC.A", "SSC.A"), preselection=NULL,scale=3, plot = T)
#Aplico el gate
linf <- Subset(clean.dt, linfos)
#Gate elipsoide para filtrar singuletes
singfilt1 <- lymphGate(linf, channels = c("FSC.H", "FSC.A"), preselection = NULL, scale=3, bwFac = 2, 
                       filterId = "singGate", evaluate = T, plot =T)
#Aplico el gate
linfossing <- Subset(linf, singfilt1)
#Gate para seleccionar MG negativas (celulas viables)
viablesfilt <- rectangleGate("RL1.A" = c(0, 2), "SSC.A" = c(4, 5))
#Aplico el gate
cellviables <- Subset(linfossing, viablesfilt)
#Cambio nombre de los archivos
pData(cellviables)$name <- c("Ind 2" ,  "Ind 3",  "Ind 4", "Ind 1", "S/Ind 2", "S/Ind 3", "S/Ind 4" , "S/Ind 1")
#Hago el calculo para ver en que punto se separan las poblaciones + de -
test1 <- (cellviables[1:4])
test2 <- (cellviables[1:8])
#hago un marker para cada canal
marcadorFRET <- rangeGate(test1, "VL2.A", plot = T, alpha = "min", filterId = "FRETfilt1")
marcadorCFP <- rangeGate(test1, "VL1.A", plot = T, alpha = "min", filterId = "CFPfilt1", sd=0.5) #agrego sd xq las poblaciones de CFP no se separan
marcadorYFP <- rangeGate(test2, "BL1.A", plot = T, alpha= "min", filterId = "YFPfilt1")
##
#Filtros para conteo de regiones positivas
YFPposfilt <- rectangleGate("BL1.A" = c(2.2, 4.5), "SSC.A" = c(4, 5))
CFPposfilt <- rectangleGate("VL1.A" = c(marcadorCFP@min, 3.5), "SSC.A" = c(4, 5))
FRETfilt <- rectangleGate("VL2.A" = c(marcadorFRET@min, 4), "SSC.A" = c(4, 5))
library(flowStats)
#FRETfilt <- curv2Filter("SSC.A", "VL2.A", bwFac = 1.7)
#FRETfilt.results <- filter(cellviables, FRETfilt)
#xyplot(`SSC.A` ~ `VL2.A`, data=cellviables, filter=FRETfilt.results, smooth=F, ylim=c(4,5)) 
       #names=TRUE, 
       #par.settings=list(gate=list(fill="black", alpha=0.2, 
        #                           col="transparent"),
         #                gate.text=list(col="darkred", alpha=0.7, cex=0.6)))

#FRETwise <- Subset(cellviables, FRETfilt.results)
#pMTviablesfilt <- rectangleGate("FL1.A" = c(3.5, 6), "SSC.A" = c(3.6, 5.4))
#Aplico el gate de Neuroepitelio
CFPpos <- Subset(cellviables, CFPposfilt)
YFPpos <- Subset(cellviables, YFPposfilt)
FretTotales <- Subset(linfossing, FRETfilt)
FretNB <- Subset(YFPpos, FRETfilt)
##
#Defino los nombres
#pData(linfossing)$name <- c("c855a>CUTie2 control", "c855a>CUTie2 +IBMX", "c855a>CUTie2 +IBMX+8BrcGMP 5 min", "c855a>CUTie2 +IMBX+8BrcGMP 10 min", "c855a>CUTie2 +IBMX+8BrcGMP 20 min")
#library(flowStats)`
#Trellis plots para todo el flowSet
xyplot(`SSC.A` ~ `FSC.A`, data = linfossing, smooth=F, xlim=c(3,5.5), ylim=c(3,5.5), stats = T)#, layout = c(3,6))#scatter pĺot
xyplot(`SSC.A` ~ `BL1.A`, data = linfossing, filter=YFPposfilt, smooth=F, xlim=c(0,6), ylim=c(3,6), stats = T)#scatter pĺot
xyplot(`SSC.A` ~ `VL1.A`, data = linfossing, filter=CFPposfilt, smooth=F, xlim=c(0,6), ylim=c(3.5,6), stats = T)#scatter pĺot
xyplot(`SSC.A` ~ `VL2.A`, data = linfossing, filter=FRETfilt, smooth=F, xlim=c(0,6), ylim=c(3.5,6), stats = T, xlab=("FRET-A"), ylab=("SSC-A") )#scatter pĺot
xyplot(`SSC.A` ~ `RL1.A`, data = linfossing, filter=MGviables, smooth=F, xlim=c(0,6), ylim=c(3,6), stats = T)#scatter pĺot
xyplot(`BL1.A` ~ `VL1.A`, data = linfossing, smooth=F, xlim=c(0,6), ylim=c(0,6), stats = T, layout = c(4,3))#scatter pĺot
##
#c2f <- curv2Filter("SSC.A", "BL1.A", bwFac=4)
#c2f.results <- filter(linfossing, c2f)
#xyplot(`SSC.A` ~ `BL1.A`, data = linfossing, filter=c2f.results, smooth=F, xlim=c(0,4), ylim=c(2,6), stats = T)#scatter pĺot
#Los Trellis por separado y enviandolos a un elemento organizable 
#plots para YFP
xy1YFP<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables[1], filter=YFPposfilt, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy2YFP<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables[2], filter=YFPposfilt, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot#xyplot(`SSC.A` ~ `FL1.A`, data=,  xyplot(`FSC.A` ~ `FL1.A`, data = viables, filter=posfilt, smooth=F, xlim=c(0,6), ylim=c(5,6), stats = T)#scatter pĺot
xy3YFP<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables[3], filter=YFPposfilt, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy4YFP<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables[4], filter=YFPposfilt, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy5YFP<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables[5], filter=YFPposfilt, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy6YFP<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables[6], filter=YFPposfilt, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy7YFP<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables[7], filter=YFPposfilt, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy8YFP<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables[8], filter=YFPposfilt, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
library('ggpubr')#cargo biblioteca de funciones para publicar 
plotYFP <- ggarrange (xy4YFP, xy1YFP, xy2YFP, xy3YFP, xy8YFP, xy5YFP, xy6YFP, xy7YFP, ncol = 4, nrow = 2)
plotYFP
#plots para CFP
xy1CFP<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables[1], filter=CFPposfilt, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy2CFP<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables[2], filter=CFPposfilt, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot#xyplot(`SSC.A` ~ `FL1.A`, data=,  xyplot(`FSC.A` ~ `FL1.A`, data = viables, filter=posfilt, smooth=F, xlim=c(0,6), ylim=c(5,6), stats = T)#scatter pĺot
xy3CFP<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables[3], filter=CFPposfilt, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy4CFP<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables[4], filter=CFPposfilt, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy5CFP<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables[5], filter=CFPposfilt, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy6CFP<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables[6], filter=CFPposfilt, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy7CFP<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables[7], filter=CFPposfilt, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy8CFP<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables[8], filter=CFPposfilt, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
plotCFP <- ggarrange (xy4CFP, xy1CFP, xy2CFP, xy3CFP, xy8CFP, xy5CFP, xy6CFP, xy7CFP, ncol = 4, nrow = 2)
plotCFP
#plots para FRET
xy1FRET<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables[1], filter=FRETfilt, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="FRET", ylab="SSC-A")#scatter pĺot
xy2FRET<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables[2], filter=FRETfilt, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="FRET", ylab="SSC-A")#scatter pĺot#xyplot(`SSC.A` ~ `FL1.A`, data=,  xyplot(`FSC.A` ~ `FL1.A`, data = viables, filter=posfilt, smooth=F, xlim=c(0,6), ylim=c(5,6), stats = T)#scatter pĺot
xy3FRET<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables[3], filter=FRETfilt, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="FRET", ylab="SSC-A")#scatter pĺot
xy4FRET<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables[4], filter=FRETfilt, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="FRET", ylab="SSC-A")#scatter pĺot
xy5FRET<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables[5], filter=FRETfilt, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="FRET", ylab="SSC-A")#scatter pĺot
xy6FRET<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables[6], filter=FRETfilt, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="FRET", ylab="SSC-A")#scatter pĺot
xy7FRET<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables[7], filter=FRETfilt, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="FRET", ylab="SSC-A")#scatter pĺot
xy8FRET<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables[8], filter=FRETfilt, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="FRET", ylab="SSC-A")#scatter pĺot
plotFRET <- ggarrange (xy4FRET, xy1FRET, xy2FRET, xy3FRET, xy8FRET, xy5FRET, xy6FRET, xy7FRET, ncol = 4, nrow = 2)
plotFRET
#Multiplot ordenado en formato publicacion
#library('ggpubr')#cargo biblioteca de funciones para publicar 
#Multipanel 1
#ggarrange(xy1, xy2, xy3, xy4, 
 #         labels = c("A", "B", "C", "D"),
  #        ncol = 2, nrow = 2)
#Multipanel 2
#ggarrange(xy5, xy6, xy7, xy8, 
 #         labels = c("E", "F", "G", "H"),
  #        ncol = 2, nrow = 2)
##
#Extraigo valores numericos
#library(flowStats)
#c2f <- curv2Filter("FL2.A", "FL4.A", bwFac = 2.1)
#c2f.results <- filter(compensados, c2f)
#xyplot(`FL2.A` ~ `FL4.A`, data = compensados, filter = c2f.results, ylab="IgM-PE", xlab="S100-APC", smooth=F, xlim=c(0,4), ylim=c(0,4), names=F)
##
#
#parallel(~., clean.dt, filter = c2f.results, alpha = 0.01)
##
#Filtro linfos (en flowSet)
#lymph <- lymphGate(compensados, channels = c("SSC-A", "FSC-A"), preselection = NULL, scale=1.5, bwFac = 2, 
#                  filterId="defaultLymphGate", evaluate = TRUE, plot = T)
##
#Cargo los FMO
#fmo.files = list.files("~/Documentos/citom/emmprin14-6-16/fmo/", all.files = T, full.names = TRUE)#Lee todos los archivos del directorio
#fmo.fs <- read.flowSet(fmo.files[3:5], transformation = F, alter.names = T)#[3:25])#Elimino las dos primeras filas que son /. y /..
#fmo.fs <- fmo.fs[, 7:12]
#Los transformo logicle
#cx <- logicleTransform()
#cxlist <- transformList(c("FL.1.Log","FL.2.Log","FL.8.Log"), cx)
#fmo.trans <- transform(fmo.fs, cxlist)
#fmo.comp <- compensate(fmo.trans, comp.mat)3.437
##
#Filtro linfos en un flowFrame del set
#lymph <- flowDensity(obj = compensados[[2]], channels = c("FSC.A", "SSC.A"), position = c(F, F), debris.gate = c(F,F))#Filtro linfos F/S
#sglt <- flowDensity(obj = lymph, singlet.gate = T)#Filtro singuletes
#IgMpos <- flowDensity(obj = sglt, channels = c(3, 2), position = c(T, NA), upper = c(F, NA))#, use.percentile = T, percentile = c(0.01, NA), use.control = c(T, F), control = c(fmo.comp[[3]], NA))#Filtro las cd19+ dentro de los singuletes
#cd3pos <- flowDensity(obj = sglt, channels = c(5, 2), position = c(T, NA), use.percentile=c(T, F), percentile=c(0.35, NA))#, use.control = c(T, F), control = c(clean.dt[[2]], NA))
#S100pos <- flowDensity(obj = sglt, channels = c(6, 2), position = c(T, NA),use.percentile=c(T, F), percentile=c(0.35, NA))#, use.control = c(T, F), control = c(fmo.comp[[2]], NA))
#IgMemmprinpos <- flowDensity(obj = sglt, channels = c(13, 15), position = c(T, T), upper = c(F, F), use.percentile = c(T, F), percentile = c(0.01, NA), use.control = c(T, T), control = c(fmo.comp[[3]], fmo.comp[[2]]))
#cd3emmprinpos <- flowDensity(obj = sglt, channels = c(11, 15), position = c(T, T), upper = c(T, F), use.percentile = c(F, T), percentile = c(NA, 0.01), use.control = c(T, T), control = c(fmo.comp[[1]], fmo.comp[[2]]))
#IgMnegcd3pos <- flowDensity(obj = sglt, channels = c(11, 13), position = c(T, T), upper = c(T, F),  use.control = c(T, T), control = c(fmo.comp[[1]], fmo.comp[[3]]))
#Ki67pos <- flowDensity(obj = sglt, channels = c(18, 9), position = c(T, NA), upper = c(T, NA), use.control = c(T, F), control = c(fmo.comp[[2]], NA))

#Histogramas comparados
#autoplot(linfossing, "VL2.A") + 
 #xlab("FRET")+
  #ylab("Densidad")#
#comparativo <- densityplot (~ `VL2.A`, cellviables, xlim = c(0,4))#, FRETfilt)#Compara los histogramas de datos transformados en FL1, 3 y4
#comparativo
#Hago un marker con rangegate
#library(flowStats)`
#Hago el calculo para ver en que punto se separan las poblaciones + de -
test <- (cellviables[1:4])
#hago un marker
marcador <- rangeGate(test, "VL1.A", plot = T, alpha = "min", filterId = "FRETfilt1")
#                      axis.title.x=element_text(size=14), axis.title.y=element_text(size=18), axis.text.x = element_text(size=14),
#                      fig.lab = "C�lculo de punto de corte", fig.lab.face = "bold", fig.lab.size=14)
conMarker <- densityplot (~ `VL2.A`,  cellviables, filter=marcador, refline = marcador@min, xlim = c(0,4), xlab=("FRET-A Log"), stack = T, stats=F)
conMarker
#marcador1 <- rangeGate(CFPpos, "VL1.A", plot = T, alpha = "min", filterId = "CFPposfilt") #hago un marker
#conMarker1 <- densityplot (~ `VL1.A`,  CFPpos, layout = c(1,1), filter=CFPposfilt, refline = FRETfilt@min, xlim = c(0,4), xlab=("CFP-A Log"), stack = T, stats=F)
#conMarker1
#calculo de densidades con poblaciones separadas
c1f <- curv1Filter("VL1.A", bwFac=1.5)
c1f.results <- filter(cellviables, c1f)
g <- densityplot (~ `VL1.A`,  cellviables, filter = c1f.results, layout = c(1,1), xlim = c(0,4), xlab=("FRET-A Log"), stack = T)
g
#Pasamos txt info de cada muestra
#pData(linfossing)$name <- c("pCoPuro s/ind", "pCoPuro Cu", "pMT_YC s/ind", "pMT_YC Cu", "pCoPuro s/ind PI", "pCoPuro ind PI", "pMT_YC s/ind PI", "pMT_YC ind PI")
#ggcyto(linfossing, aes(x = `FL1.A`, fill = name)) + geom_density(alpha = 0.2)
#Construyo plots de densidad monoparametricos y los envio a un objeto ordenable
#d1 <- ggplot(linfossing[1:2], aes(x = `FL1.A`, fill = name )) + geom_density(alpha = 0.2) + xlim (2,5) + xlab("YFP") + theme(legend.position = c(0.75, 0.9), legend.title = element_blank(), legend.text = element_text(size=8), legend.key.size = unit(10, "pt"), legend.background = element_blank())
#d2 <- ggplot(linfossing[5:6], aes(x = `FL2.A`, fill = name )) + geom_density(alpha = 0.2) + xlim (2,5) + xlab("PI") + theme(legend.position = c(0.75, 0.9), legend.title = element_blank(), legend.text = element_text(size=8), legend.key.size = unit(10, "pt"), legend.background = element_blank())
#d3 <- ggplot(linfossing[3:4], aes(x = `FL1.A`, fill = name )) + geom_density(alpha = 0.2) + xlim (2,5) + xlab("YFP") + theme(legend.position = c(0.75, 0.9), legend.title = element_blank(), legend.text = element_text(size=8), legend.key.size = unit(10, "pt"), legend.background = element_blank())
#d4 <- ggplot(linfossing[7:8], aes(x = `FL2.A`, fill = name )) + geom_density(alpha = 0.2) + xlim (2,5) + xlab("PI") + theme(legend.position = c(0.75, 0.9), legend.title = element_blank(), legend.text = element_text(size=8), legend.key.size = unit(10, "pt"), legend.background = element_blank())
#Multipanel 3
#ggarrange(d1, d2, d3, d4 + rremove("x.text"), 
 #         labels = c("I", "J", "K", "L"),
  #        ncol = 2, nrow = 2)
#EOF
#Agrego copia del script pero para abrir archivos a dif tiempo para comparar despues
#Abro archivos de 72 hs de inducci�n
files.test2 = list.files("D:/Maestria neurodesarrollo/fotos y resultados ensayos/archivos citometro nov/2021-11-18 S2 con y sin ind 72 hs/", all.files = F, full.names = TRUE)#Lee todos los archivos del directorio
fs2 <- read.flowSet(files.test2[1:8], transformation = F, alter.names =T, emptyValue=FALSE)#[3:25])#Elimino las dos primeras filas que son /. y /..
#Transformo logicle
bx2 <- logicleTransform()
bxlist2 <- transformList(c("FSC.A", "SSC.A", "FSC.H", "BL1.A", "VL2.A", "VL1.A", "RL1.A"), bx2)
datostrans2 <- transform(fs2, bxlist2)
##
#Limpiamos el dataset de debris (fcs/ssc)
clean2.dt <- datostrans2
for (i in 1:8) { # Loop over the length of the flowSet
  f <- datostrans2[[i]]
  # First restrict the FSC - A values :
  fsc.indices2 <- intersect (which (exprs (f)[, "FSC.A"]<5.5) , which (exprs (f) [, "FSC.A"] > 3.5))
  # Then restrict SSC - A values and intersect with FSC - A restriction above :
  ssc.indices2 <- intersect (which ( exprs ( f)[, "SSC.A"] > 3) ,
                            which ( exprs (f)[, "SSC.A"] < 5.5) )
  non.debris.indices2 <- intersect ( fsc.indices2 , ssc.indices2 )
  # Now only select the non debris cells and place the cleaned up flowFrame into the flowSet :
  f2.clean <- f[non.debris.indices2]
  clean2.dt [[i]] <- f2.clean
}
##
#Muestro los plots
plotDens (clean2.dt[[4]], c("FSC.A" ,"SSC.A"))
plotDens (clean2.dt[[5]], c ("VL2.A", "SSC.A"))
plotDens (clean2.dt[[3]], c ("FSC.A", "FSC.H"))
autoplot(clean2.dt[[4]]) #+ labs_cyto("marker")
#Ahora un poco de gating y filtrado sobre todo el flowSet.oO
#Gate elipsoide tipo linfocitos
linfos2 <- lymphGate(clean2.dt, channels=c("FSC.A", "SSC.A"), preselection=NULL,scale=3, plot = T)
#Aplico el gate
linf2 <- Subset(clean2.dt, linfos2)
#Gate elipsoide para filtrar singuletes
singfilt12 <- lymphGate(linf2, channels = c("FSC.H", "FSC.A"), preselection = NULL, scale=3, bwFac = 2, 
                       filterId = "singGate", evaluate = T, plot =T)
#Aplico el gate
linfossing2 <- Subset(linf2, singfilt12)
#Gate para seleccionar MG negativas (celulas viables)
viablesfilt2 <- rectangleGate("RL1.A" = c(0, 2.2), "SSC.A" = c(4, 5))
#Aplico el gate
cellviables2 <- Subset(linfossing2, viablesfilt2)
#Cambio nombre de los archivos
pData(cellviables2)$name <- c("Ind 2" ,  "Ind 3",  "Ind 4", "Ind 1", "S/Ind 2", "S/Ind 3", "S/Ind 4" , "S/Ind 1")
#Hago el calculo para ver en que punto se separan las poblaciones + de -
test12 <- (cellviables2[1:4])
test22 <- (cellviables2[1:8])
#hago un marker para cada canal
marcadorFRET2 <- rangeGate(test12, "VL2.A", plot = T, alpha = "min", filterId = "FRETfilt1")
marcadorCFP2 <- rangeGate(test22, "VL1.A", plot = T, alpha = "min", filterId = "CFPfilt1", sd=0.8) #agrego sd xq las poblaciones de CFP no se separan
marcadorYFP2 <- rangeGate(test22, "BL1.A", plot = T, alpha = "min", filterId = "YFPfilt1")
##
#Filtros para conteo de regiones positivas
YFPposfilt2 <- rectangleGate("BL1.A" = c(marcadorYFP2@min, 4.5), "SSC.A" = c(4, 5))
CFPposfilt2 <- rectangleGate("VL1.A" = c(marcadorCFP2@min, 3.5), "SSC.A" = c(4, 5))
FRETfilt2 <- rectangleGate("VL2.A" = c(marcadorFRET2@min, 4), "SSC.A" = c(4, 5))
library(flowStats)
#FRETfilt <- curv2Filter("SSC.A", "VL2.A", bwFac = 1.7)
#FRETfilt.results <- filter(cellviables, FRETfilt)
#xyplot(`SSC.A` ~ `VL2.A`, data=cellviables, filter=FRETfilt.results, smooth=F, ylim=c(4,5)) 
#names=TRUE, 
#par.settings=list(gate=list(fill="black", alpha=0.2, 
#                           col="transparent"),
#                gate.text=list(col="darkred", alpha=0.7, cex=0.6)))

#FRETwise <- Subset(cellviables, FRETfilt.results)
#pMTviablesfilt <- rectangleGate("FL1.A" = c(3.5, 6), "SSC.A" = c(3.6, 5.4))
#Aplico el gate de Neuroepitelio
#CFPpos <- Subset(cellviables, CFPposfilt)
#YFPpos <- Subset(cellviables, YFPposfilt)
#FretTotales <- Subset(linfossing, FRETfilt)
#FretNB <- Subset(YFPpos, FRETfilt)
##
#Defino los nombres
#pData(linfossing)$name <- c("c855a>CUTie2 control", "c855a>CUTie2 +IBMX", "c855a>CUTie2 +IBMX+8BrcGMP 5 min", "c855a>CUTie2 +IMBX+8BrcGMP 10 min", "c855a>CUTie2 +IBMX+8BrcGMP 20 min")
#library(flowStats)`
#Trellis plots para todo el flowSet
xyplot(`SSC.A` ~ `FSC.A`, data = linfossing2, smooth=F, xlim=c(3,5.5), ylim=c(3,5.5), stats = T)#, layout = c(3,6))#scatter pĺot
xyplot(`SSC.A` ~ `BL1.A`, data = linfossing2, filter=YFPposfilt2, smooth=F, xlim=c(0,6), ylim=c(3,6), stats = T)#scatter pĺot
xyplot(`SSC.A` ~ `VL1.A`, data = linfossing2, filter=CFPposfilt2, smooth=F, xlim=c(0,6), ylim=c(3.5,6), stats = T)#scatter pĺot
xyplot(`SSC.A` ~ `VL2.A`, data = linfossing2, filter=FRETfilt2, smooth=F, xlim=c(0,6), ylim=c(3.5,6), stats = T, xlab=("FRET-A"), ylab=("SSC-A") )#scatter pĺot
xyplot(`SSC.A` ~ `RL1.A`, data = linfossing2, filter=viablesfilt2, smooth=F, xlim=c(0,6), ylim=c(3,6), stats = T)#scatter pĺot
xyplot(`BL1.A` ~ `VL1.A`, data = linfossing2, smooth=F, xlim=c(0,6), ylim=c(0,6), stats = T, layout = c(4,3))#scatter pĺot
##
#c2f <- curv2Filter("SSC.A", "BL1.A", bwFac=4)
#c2f.results <- filter(linfossing, c2f)
#xyplot(`SSC.A` ~ `BL1.A`, data = linfossing, filter=c2f.results, smooth=F, xlim=c(0,4), ylim=c(2,6), stats = T)#scatter pĺot
#Los Trellis por separado y enviandolos a un elemento organizable 
#plots para YFP
xy1YFP2<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables2[1], filter=YFPposfilt2, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy2YFP2<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables2[2], filter=YFPposfilt2, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot#xyplot(`SSC.A` ~ `FL1.A`, data=,  xyplot(`FSC.A` ~ `FL1.A`, data = viables, filter=posfilt, smooth=F, xlim=c(0,6), ylim=c(5,6), stats = T)#scatter pĺot
xy3YFP2<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables2[3], filter=YFPposfilt2, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy4YFP2<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables2[4], filter=YFPposfilt2, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy5YFP2<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables2[5], filter=YFPposfilt2, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy6YFP2<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables2[6], filter=YFPposfilt2, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy7YFP2<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables2[7], filter=YFPposfilt2, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy8YFP2<-xyplot(`SSC.A` ~ `BL1.A`, data = cellviables2[8], filter=YFPposfilt2, smooth=F, xlim=c(0,4.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
#library('ggpubr')#cargo biblioteca de funciones para publicar 
plotYFP2 <- ggarrange (xy4YFP2, xy1YFP2, xy2YFP2, xy3YFP2, xy8YFP2, xy5YFP2, xy6YFP2, xy7YFP2, ncol = 4, nrow = 2)
plotYFP2
#plots para CFP
xy1CFP2<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables2[1], filter=CFPposfilt2, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy2CFP2<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables2[2], filter=CFPposfilt2, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot#xyplot(`SSC.A` ~ `FL1.A`, data=,  xyplot(`FSC.A` ~ `FL1.A`, data = viables, filter=posfilt, smooth=F, xlim=c(0,6), ylim=c(5,6), stats = T)#scatter pĺot
xy3CFP2<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables2[3], filter=CFPposfilt2, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy4CFP2<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables2[4], filter=CFPposfilt2, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy5CFP2<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables2[5], filter=CFPposfilt2, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy6CFP2<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables2[6], filter=CFPposfilt2, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy7CFP2<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables2[7], filter=CFPposfilt2, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy8CFP2<-xyplot(`SSC.A` ~ `VL1.A`, data = cellviables2[8], filter=CFPposfilt2, smooth=F, xlim=c(0,3.5), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
plotCFP2 <- ggarrange (xy4CFP2, xy1CFP2, xy2CFP2, xy3CFP2, xy8CFP2, xy5CFP2, xy6CFP2, xy7CFP2, ncol = 4, nrow = 2)
plotCFP2 
#plots para FRET
xy1FRET2<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables2[1], filter=FRETfilt2, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="", ylab="", )#scatter pĺot
xy2FRET2<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables2[2], filter=FRETfilt2, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot#xyplot(`SSC.A` ~ `FL1.A`, data=,  xyplot(`FSC.A` ~ `FL1.A`, data = viables, filter=posfilt, smooth=F, xlim=c(0,6), ylim=c(5,6), stats = T)#scatter pĺot
xy3FRET2<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables2[3], filter=FRETfilt2, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy4FRET2<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables2[4], filter=FRETfilt2, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy5FRET2<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables2[5], filter=FRETfilt2, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy6FRET2<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables2[6], filter=FRETfilt2, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy7FRET2<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables2[7], filter=FRETfilt2, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
xy8FRET2<-xyplot(`SSC.A` ~ `VL2.A`, data = cellviables2[8], filter=FRETfilt2, smooth=F, xlim=c(0,4), ylim=c(4,5), stats = T, xlab="", ylab="")#scatter pĺot
plotFRET2 <- ggarrange (xy4FRET2, xy1FRET2, xy2FRET2, xy3FRET2, xy8FRET2, xy5FRET2, xy6FRET2, xy7FRET2, ncol = 4, nrow = 2)
plotFRET2
#////////////////////////////////////////////////////////////////////////
#Goznalez Ana Clara (anaclgonzalez92@gmail.com)
#A�exo para hacer plots en conjunto para comparar lso 3 canales en 48 vs 72 hs
#
#
#Extraemos datos para ggplot con ggcyto
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(palmerpenguins)
#tengo que cambiar el nombre de los archivos para que sean todos iguales los de 48 para que me los junte
pData(cellviables)$name <- c("48 hs" ,  "48 hs",  "48 hs", "48 hs", "48 hs", "48 hs", "48 hs" , "48 hs")
pData(cellviables2)$name <- c("72 hs" , "72 hs",  "72 hs", "72 hs", "72 hs", "72 hs", "72 hs" , "72 hs")
#genero los rbind con todos los datos de replicas de 48 por un lado y de 72 por el otro
fs48ind1 <- rbind2(cellviables[1],cellviables[2])
fs48ind2 <- rbind2(cellviables[3],cellviables[4])
fs48indfinal <- rbind2(fs48ind1,fs48ind2)
fs72ind1 <- rbind2(cellviables2[1],cellviables2[2])
fs72ind2 <- rbind2(cellviables2[3],cellviables2[4])
fs72indfinal <- rbind2(fs72ind1,fs72ind2)  
#grafico para comparar con forma de violin
comparados1 <- rbind2(fs48indfinal, fs72indfinal)#genero un nuevo flowSet con los datos que me interesan para cada gr�fico
p1 <- ggcyto(comparados1, aes(y = `BL1.A`))
p1 + geom_violin(aes(x=factor(name), fill=name)) +
  scale_fill_manual(values = c("#FADEE1", "#D9F1F1"), name= "Tiempo")+
  coord_cartesian(ylim =c(0, 5))+#por alg�n endiablado motivo si no usas esto no puedes establecer los l�mites en y, y la barra de estadistica queda afuera
  facet_null()+#esto es para que est�n juntas
  stat_signif(aes(x=factor(name)), comparisons = list(c("48 hs", "72 hs")), test = "wilcox.test", map_signif_level=TRUE, color="black", size = 0.5, textsize = 5, y_position =  4.5)+
  xlab("Tiempo de inducci�n") + ylab("YFP") +
  theme_classic(base_size = 14) +
  labs(title="Comparaci�n de tiempos de inducci�n", size = 14) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "")
                  