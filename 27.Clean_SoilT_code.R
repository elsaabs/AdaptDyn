setwd("/Users/elsa/Desktop/Abby/")
#Packages
library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
#Open NetCDF file
ncsoiltemp = nc_open("tsl_Lmon_CCSM4_rcp85_r1i1p1_200601-205002.nc")
print(ncsoiltemp)
#Get the dimension variables
lon = ncvar_get(ncsoiltemp, "lon")-180
lat = ncvar_get(ncsoiltemp, "lat")
depth = ncvar_get(ncsoiltemp, "depth")
time = ncvar_get(ncsoiltemp, "time")
#Only keep 1st soil layer (surface) and transform the matrix of "soiltemp" in a 2D matrix (1 location per row, 1 monthly date per column)
B = matrix(ncvar_get(ncsoiltemp, "tsl")[,,1,], nrow = 55296, ncol = 530, byrow = FALSE)
#Average per location of litter flux over 15 years (from January 2006 to Decembre 2020 included) in °C
B2 = B[,1:180]
B3 = rowMeans(B2, na.rm=F, dims = 1)
B4 = B3-273
#Transform the matrix to longitude in rows and latitude in columns
B5 = matrix(B4, nrow = 288, ncol = 192, byrow = FALSE)
B51 = B5[1:145,]
B52 = B5[146:288,]
B6 = rbind(B52,B51)
#Save the matrix as a table
write.table(B6, file = "soilT", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)




#Map Soil Temperature
par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
brks = seq(-60,40, by=0.1)
colfunc = colorRampPalette(c(
  "#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf",
  "#fee090","#fdae61","#f46d43","#d73027","#a50026"))(length(brks)-1)
fields::image.plot(x=lon, y=lat, B6, ann=F, breaks=brks,
                   horizontal=T, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='Temperature (°C)', 
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("Surface soil temperature averaged over 15 years 
      (from January 2006 to December 2020 included)", outer=F)

