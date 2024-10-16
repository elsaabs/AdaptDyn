setwd("/Volumes/ElsaSSD2022/4.LSCE/PhD/Abby_chapter_5/2023_Oct_recovered_from_old_harddrive/nc_files")

#Packages
library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
#Open NetCDF file
#Files r1i1p1
ncsoiltemp1 = nc_open("tsl_Lmon_CCSM4_rcp85_r1i1p1_200601-205002.nc")
ncsoiltemp2 = nc_open("tsl_Lmon_CCSM4_rcp85_r1i1p1_205001-210012.nc")
#Files r2i1p1
#ncsoiltemp1 = nc_open("tsl_Lmon_CCSM4_rcp85_r2i1p1_200601-204912.nc")
#ncsoiltemp2 = nc_open("tsl_Lmon_CCSM4_rcp85_r2i1p1_205001-210012.nc")
#Get the dimension variables
lon = ncvar_get(ncsoiltemp1, "lon")-180
lat = ncvar_get(ncsoiltemp1, "lat")
depth = ncvar_get(ncsoiltemp1, "depth")
#depth2 = ncvar_get(ncsoiltemp1, "depth")
time1 = ncvar_get(ncsoiltemp1, "time")
time2 = ncvar_get(ncsoiltemp2, "time")
#Only keep 1st soil layer (surface) and transform the matrix of "soiltemp" in a 2D matrix (1 location per row, 1 monthly date per column)
T1 = matrix(ncvar_get(ncsoiltemp1, "tsl")[,,1,], nrow = 55296, ncol = 530, byrow = FALSE)
T2 = matrix(ncvar_get(ncsoiltemp2, "tsl")[,,1,], nrow = 55296, ncol = 612, byrow = FALSE)
#Mean annual (from January to Decembre included) surface soil temperature (in °C) from 2010 to 2100 every 10 years
C10 = rowMeans(T1[,49:60]-273, na.rm=F, dims = 1)
#mean(c(T1[32000,169],T1[32000,170],T1[32000,171],T1[32000,172],T1[32000,173],
#       T1[32000,174],T1[32000,175],T1[32000,176],T1[32000,177],T1[32000,178],
#       T1[32000,179],T1[32000,180]))-273
C20 = rowMeans(T1[,169:180]-273, na.rm=F, dims = 1)
C30 = rowMeans(T1[,289:300]-273, na.rm=F, dims = 1)
C40 = rowMeans(T1[,409:420]-273, na.rm=F, dims = 1)
C50 = rowMeans(T2[,1:12]-273, na.rm=F, dims = 1)
C60 = rowMeans(T2[,121:132]-273, na.rm=F, dims = 1)
C70 = rowMeans(T2[,241:252]-273, na.rm=F, dims = 1)
C80 = rowMeans(T2[,361:372]-273, na.rm=F, dims = 1)
C90 = rowMeans(T2[,481:492]-273, na.rm=F, dims = 1)
C100 = rowMeans(T2[,601:612]-273, na.rm=F, dims = 1)
#Transform the lists into matrices wih longitude in rows and latitude in columns
C10a = matrix(C10, nrow = 288, ncol = 192, byrow = FALSE)
C10b1 = C10a[1:145,]
C10b2 = C10a[146:288,]
C10c = rbind(C10b2,C10b1)
C20a = matrix(C20, nrow = 288, ncol = 192, byrow = FALSE)
C20b1 = C20a[1:145,]
C20b2 = C20a[146:288,]
C20c = rbind(C20b2,C20b1)
C30a = matrix(C30, nrow = 288, ncol = 192, byrow = FALSE)
C30b1 = C30a[1:145,]
C30b2 = C30a[146:288,]
C30c = rbind(C30b2,C30b1)
C40a = matrix(C40, nrow = 288, ncol = 192, byrow = FALSE)
C40b1 = C40a[1:145,]
C40b2 = C40a[146:288,]
C40c = rbind(C40b2,C40b1)
C50a = matrix(C50, nrow = 288, ncol = 192, byrow = FALSE)
C50b1 = C50a[1:145,]
C50b2 = C50a[146:288,]
C50c = rbind(C50b2,C50b1)
C60a = matrix(C60, nrow = 288, ncol = 192, byrow = FALSE)
C60b1 = C60a[1:145,]
C60b2 = C60a[146:288,]
C60c = rbind(C60b2,C60b1)
C70a = matrix(C70, nrow = 288, ncol = 192, byrow = FALSE)
C70b1 = C70a[1:145,]
C70b2 = C70a[146:288,]
C70c = rbind(C70b2,C70b1)
C80a = matrix(C80, nrow = 288, ncol = 192, byrow = FALSE)
C80b1 = C80a[1:145,]
C80b2 = C80a[146:288,]
C80c = rbind(C80b2,C80b1)
C90a = matrix(C90, nrow = 288, ncol = 192, byrow = FALSE)
C90b1 = C90a[1:145,]
C90b2 = C90a[146:288,]
C90c = rbind(C90b2,C90b1)
C100a = matrix(C100, nrow = 288, ncol = 192, byrow = FALSE)
C100b1 = C100a[1:145,]
C100b2 = C100a[146:288,]
C100c = rbind(C100b2,C100b1)


#Map Soil Temperature
par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
brks = seq(-60,40, by=0.1)
colfunc = colorRampPalette(c(
  "#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf",
  "#fee090","#fdae61","#f46d43","#d73027","#a50026"))(length(brks)-1)
fields::image.plot(x=lon, y=lat, C100c, ann=F, breaks=brks,
                   horizontal=T, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='Temperature (°C)', 
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("Mean surface soil temperature in 2100", outer=F)

#Map Change in Soil T
par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
brks = seq(-2,15, by=0.1)
colfunc = colorRampPalette(c(
  "#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026"))(length(brks)-1)
fields::image.plot(x=lon, y=lat, C100c-C10c, ann=F, breaks=brks,
                   horizontal=T, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='Temperature (°C)', 
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("Mean surface soil temperature in 2100", outer=F)

max(C100c-C10c, na.rm=T)
points(20,60, pch = 16)
(C100c-C10c)[which(lon==20),which(round(lat)==60)]

#Save the matrices as tables
write.table(C10c, file = "soilT2010", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(C20c, file = "soilT2020", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(C30c, file = "soilT2030", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(C40c, file = "soilT2040", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(C50c, file = "soilT2050", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(C60c, file = "soilT2060", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(C70c, file = "soilT2070", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(C80c, file = "soilT2080", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(C90c, file = "soilT2090", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(C100c, file = "soilT2100", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
#Plot global mean surface soil temperature as a function of years
#Removing soil temperatures below -20°C (mainly in Antarctica)
C10d = C10c
C10d[C10d < -20] = NA
C10e = mean(C10d,na.rm=T)
C20d = C20c
C20d[C20d < -20] = NA
C20e = mean(C20d,na.rm=T)
C30d = C30c
C30d[C30d < -20] = NA
C30e = mean(C30d,na.rm=T)
C40d = C40c
C40d[C40d < -20] = NA
C40e = mean(C40d,na.rm=T)
C50d = C50c
C50d[C50d < -20] = NA
C50e = mean(C50d,na.rm=T)
C60d = C60c
C60d[C60d < -20] = NA
C60e = mean(C60d,na.rm=T)
C70d = C70c
C70d[C70d < -20] = NA
C70e = mean(C70d,na.rm=T)
C80d = C80c
C80d[C80d < -20] = NA
C80e = mean(C80d,na.rm=T)
C90d = C90c
C90d[C90d < -20] = NA
C90e = mean(C90d,na.rm=T)
C100d = C100c
C100d[C100d < -20] = NA
C100e = mean(C100d,na.rm=T)
#Plot
years = c(2010,2020,2030,2040,2050,2060,2070,2080,2090,2100)
meansoilT = c(C10e,C20e,C30e,C40e,C50e,C60e,C70e,C80e,C90e,C100e)
par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(0,0,0,0))
plot(years, meansoilT, type="p", pch=16, col="red", 
     xlab="Years", ylab="Surface soil temperature (°C)")
title("Mean annual surface soil temperature as a function 
      of years with RCP8.5 warming scenario")
