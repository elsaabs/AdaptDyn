setwd("/Volumes/ElsaSSD2022/4.LSCE/PhD/Abby_chapter_5/2023_Oct_recovered_from_old_harddrive/nc_files/")
#Packages
library(raster)
library(rgdal)
library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
#Create a raster of longitude and latitude in the same dimension as surface soil temperature and litter flux
world <- raster(nrows=192, ncols=288)
s1 = shapefile("./wwf_terr_ecos_oRn/wwf_terr_ecos_oRn.shx")
biomes = rasterize(s1, world, "BIOME", fun='last', background=NA, mask=FALSE, update=FALSE, updateValue='all', filename="", na.rm=TRUE)
biomes[biomes > 15] <- NA
brks = seq(0, 14, by=1)
mycols <- colors()[c(7,8,26,30,31,32,33,37,42,47,51,52,56,62)]
biomesmat <- as.matrix(biomes)
biomesmat2 <- biomesmat[nrow(biomesmat):1, ]
biomesmat3 <- t(biomesmat2)
#Save this 14 biomes matrix
write.table(biomesmat3, file = "14biomes_table", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
