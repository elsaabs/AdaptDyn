setwd("/Volumes/ElsaSSD2022/4.LSCE/PhD/Abby_chapter_5/2023_Oct_recovered_from_old_harddrive/nc_files/")
#Packages
library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
#Open NetCDF file
nclitt = nc_open("fVegLitter_Lmon_CCSM4_rcp85_r1i1p1_200601-210012.nc")
#Get Africa in the middle of the map
lon = ncvar_get(nclitt, "lon")-180
lat = ncvar_get(nclitt, "lat")
time = ncvar_get(nclitt, "time")
#Transform the matrix for litter flux in a 2D matrix (1 location per row, 1 monthly date per column)
D = matrix(ncvar_get(nclitt, "fVegLitter")[,,], nrow = 55296, ncol = 1140, byrow = FALSE)
#Remove absurd value
D[D > 1e+32] = NA
#Change unit to mg/cm3/h
D = D*100*3600
#Average per location of litter flux over 15 years (from January 2006 to Decembre 2020 included) in mg/cm3/h 
D10 = rowMeans(D[,49:60], na.rm=F, dims = 1)
D20 = rowMeans(D[,169:180], na.rm=F, dims = 1)
D30 = rowMeans(D[,289:300], na.rm=F, dims = 1)
D40 = rowMeans(D[,409:420], na.rm=F, dims = 1)
D50 = rowMeans(D[,529:540], na.rm=F, dims = 1)
D60 = rowMeans(D[,649:660], na.rm=F, dims = 1)
D70 = rowMeans(D[,769:780], na.rm=F, dims = 1)
D80 = rowMeans(D[,889:900], na.rm=F, dims = 1)
D90 = rowMeans(D[,1009:1020], na.rm=F, dims = 1)
D100 = rowMeans(D[,1129:1140], na.rm=F, dims = 1)
#Transform the list into a matrix to longitude in rows and latitude in columns
D10a = matrix(D10, nrow = 288, ncol = 192, byrow = FALSE)
D10b1 = D10a[1:145,]
D10b2 = D10a[146:288,]
D10c = rbind(D10b2,D10b1)
D20a = matrix(D20, nrow = 288, ncol = 192, byrow = FALSE)
D20b1 = D20a[1:145,]
D20b2 = D20a[146:288,]
D20c = rbind(D20b2,D20b1)
D30a = matrix(D30, nrow = 288, ncol = 192, byrow = FALSE)
D30b1 = D30a[1:145,]
D30b2 = D30a[146:288,]
D30c = rbind(D30b2,D30b1)
D40a = matrix(D40, nrow = 288, ncol = 192, byrow = FALSE)
D40b1 = D40a[1:145,]
D40b2 = D40a[146:288,]
D40c = rbind(D40b2,D40b1)
D50a = matrix(D50, nrow = 288, ncol = 192, byrow = FALSE)
D50b1 = D50a[1:145,]
D50b2 = D50a[146:288,]
D50c = rbind(D50b2,D50b1)
D60a = matrix(D60, nrow = 288, ncol = 192, byrow = FALSE)
D60b1 = D60a[1:145,]
D60b2 = D60a[146:288,]
D60c = rbind(D60b2,D60b1)
D70a = matrix(D70, nrow = 288, ncol = 192, byrow = FALSE)
D70b1 = D70a[1:145,]
D70b2 = D70a[146:288,]
D70c = rbind(D70b2,D70b1)
D80a = matrix(D80, nrow = 288, ncol = 192, byrow = FALSE)
D80b1 = D80a[1:145,]
D80b2 = D80a[146:288,]
D80c = rbind(D80b2,D80b1)
D90a = matrix(D90, nrow = 288, ncol = 192, byrow = FALSE)
D90b1 = D90a[1:145,]
D90b2 = D90a[146:288,]
D90c = rbind(D90b2,D90b1)
D100a = matrix(D100, nrow = 288, ncol = 192, byrow = FALSE)
D100b1 = D100a[1:145,]
D100b2 = D100a[146:288,]
D100c = rbind(D100b2,D100b1)
#Map litter flux matrix
par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
brks = seq(0,0.025, by=1e-4)
colfunc = colorRampPalette(c(
  "#d9f0d3","#a6dba0","#5aae61","#1b7837","#00441b"))(length(brks)-1)
fields::image.plot(x=lon, y=lat, D100c, ann=F, breaks=brks,
                   horizontal=T, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='Litter flux (mg/cm3/y)', 
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("Mean litter flux in 2100", outer=F)
#Save the matrices as tables
write.table(D10c, file = "litter2010", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(D20c, file = "litter2020", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(D30c, file = "litter2030", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(D40c, file = "litter2040", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(D50c, file = "litter2050", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(D60c, file = "litter2060", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(D70c, file = "litter2070", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(D80c, file = "litter2080", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(D90c, file = "litter2090", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
write.table(D100c, file = "litter2100", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
#Plot global mean litter flux as a function of years
D10d = mean(D10c,na.rm=T)
D20d = mean(D20c,na.rm=T)
D30d = mean(D30c,na.rm=T)
D40d = mean(D40c,na.rm=T)
D50d = mean(D50c,na.rm=T)
D60d = mean(D60c,na.rm=T)
D70d = mean(D70c,na.rm=T)
D80d = mean(D80c,na.rm=T)
D90d = mean(D90c,na.rm=T)
D100d = mean(D100c,na.rm=T)
years = c(2010,2020,2030,2040,2050,2060,2070,2080,2090,2100)
meanlitter = c(D10d,D20d,D30d,D40d,D50d,D60d,D70d,D80d,D90d,D100d)
par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(0,0,0,0))
plot(years, meanlitter, type="p", pch=16, col="chartreuse4", 
     xlab="Years", ylab="Litter flux (mg cm-3)")
title("Mean global annual litter flux as a function 
      of years with RCP8.5 warming scenario")

