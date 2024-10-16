# Set working directory
setwd("/Volumes/ElsaSSD2022/4.LSCE/PhD/Abby_chapter_5/2023_Oct_recovered_from_old_harddrive/output")

# Load necessary packages
library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
library(rootSolve)
library(deSolve)
library(matrixcalc)

# Function to load matrices
load_matrices <- function(prefix, years, c0_value=NULL, suffix=NULL, special_case=NULL, special_year=NULL) {
  matrices <- list()
  for (year in years) {
    if (!is.null(special_case) && year == special_year) {
      file_path <- special_case
    } else if (!is.null(c0_value) && !is.null(suffix)) {
      file_path <- paste0(prefix, '_c0=', c0_value, '_', suffix, '_', year)
    } else {
      file_path <- paste0(prefix, year)
    }
    df <- read.table(file_path, sep=" ", dec='.', header=F)
    matrices[[as.character(year)]] <- as.matrix(df)
  }
  return(matrices)
}

# Function to multiply each matrix by land area
multiply_by_landarea <- function(matrices, landarea) {
  multiplied_matrices <- lapply(matrices, function(matrix) {
    matrix * landarea
  })
  return(multiplied_matrices)
}

# Function to calculate the sum of each matrix
calculate_sums <- function(matrices) {
  sums <- sapply(matrices, function(matrix) {
    sum(matrix, na.rm=TRUE)
  })
  return(sums)
}

# Function to calculate the difference with the base year (2010)
calculate_differences <- function(matrices, base_year="2010") {
  differences <- lapply(names(matrices), function(year) {
    matrices[[year]] - matrices[[base_year]]
  })
  names(differences) <- names(matrices)
  return(differences)
}

# Function to calculate average global soil temperature
calculate_avg_soilT <- function(soilT_matrices, landarea) {
  avg_soilT <- sapply(soilT_matrices, function(matrix) {
    soilT_multiplied <- matrix * landarea
    valid_landarea_sum <- sum(landarea[!is.na(matrix)], na.rm=TRUE)
    avg_temp <- sum(soilT_multiplied, na.rm=TRUE) / valid_landarea_sum
    return(avg_temp)
  })
  return(avg_soilT)
}

# Calculate average global soil temperature for each decade
special_case <- 'soilT20062020'
years <- seq(2010, 2100, by=10)
soilT_matrices <- load_matrices('soilT', years, special_case=special_case, special_year=2010)
landarea_df <- read.table('landarea', sep= " ", dec='.', header=F)
landarea <- as.matrix(landarea_df)
avg_soilT <- calculate_avg_soilT(soilT_matrices, landarea)

# Ensure avg_soilT is available and correctly computed
years_plot <- as.character(seq(2010, 2100, by=10))
avg_soilT_plot <- avg_soilT[as.character(years_plot)]

# Function to process data for a given c0 value
c0_values <- c(1.1,1.17,1.24)
process_data_for_c0 <- function(c0_value) {
  # Load 'eco' and 'ecoevo' matrices
  soilC_eco_matrices <- load_matrices('soilC', years, c0_value, 'eco')
  soilC_ecoevo_matrices <- load_matrices('soilC', years, c0_value, 'ecoevo')
  
  # Multiply each matrix by landarea
  soilC_eco_multiplied <- multiply_by_landarea(soilC_eco_matrices, landarea)
  soilC_ecoevo_multiplied <- multiply_by_landarea(soilC_ecoevo_matrices, landarea)
  
  # Calculate differences with the base year (2010) for both eco and ecoevo matrices
  soilC_eco_differences <- calculate_differences(soilC_eco_multiplied, base_year="2010")
  soilC_ecoevo_differences <- calculate_differences(soilC_ecoevo_multiplied, base_year="2010")
  
  return(list(eco_differences = soilC_eco_differences, ecoevo_differences = soilC_ecoevo_differences))
}

# Process data for each c0 value
results <- lapply(c0_values, process_data_for_c0)

# Extract the difference matrix for c0 = 1.17 and 'eco' treatment
c0_value <- 1.17
soilC_eco_matrices <- load_matrices('soilC', years, c0_value, 'eco')
soilC_ecoevo_matrices <- load_matrices('soilC', years, c0_value, 'ecoevo')

# Load stable_aat mask
stable_mask <- read.table('stable_c0=1.17', sep=" ", dec='.', header=F)
stable_mask <- as.matrix(stable_mask)

# Directly apply the stable mask to the soilT matrix for 2010 and 2100
masked_soilT_matrix_2010 <- soilT_matrices[['2010']]
masked_soilT_matrix_2010[stable_mask == 0] <- NA
masked_soilT_matrix_2100 <- soilT_matrices[['2100']]
masked_soilT_matrix_2100[stable_mask == 0] <- NA

# Calculate temperature difference for stable locations only
T_change_2010_2100 <- masked_soilT_matrix_2100 - masked_soilT_matrix_2010

# Map change in soil T between 2010 and 2100
old.par <- par(mar = c(0, 0, 0, 0))
par(old.par)
nclitt <- nc_open("fVeglitter_Lmon_CCSM4_rcp85_r1i1p1_200501-210012.nc")
lon <- ncvar_get(nclitt, "lon") - 180
lat <- ncvar_get(nclitt, "lat")
par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
brks = seq(-60,40, by=0.1)
colfunc = colorRampPalette(c(
  "#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf",
  "#fee090","#fdae61","#f46d43","#d73027","#a50026"))(length(brks)-1)
fields::image.plot(x=lon, y=lat, masked_soilT_matrix_2100, ann=F, breaks=brks,
                   horizontal=T, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='Temperature (°C)', 
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("soilT2100", outer=F)

# Map change in soil T between 2010 and 2100
old.par <- par(mar = c(0, 0, 0, 0))
par(old.par)
nclitt <- nc_open("fVeglitter_Lmon_CCSM4_rcp85_r1i1p1_200501-210012.nc")
lon <- ncvar_get(nclitt, "lon") - 180
lat <- ncvar_get(nclitt, "lat")
par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
brks = seq(-8,8, by=0.1)
colfunc = colorRampPalette(c(
  "#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","gray48",
  "#fee090","#fdae61","#f46d43","#d73027","#a50026"))(length(brks)-1)
fields::image.plot(x=lon, y=lat, masked_soilT_matrix_2100-masked_soilT_matrix_2010, ann=F, breaks=brks,
                   horizontal=T, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='Temperature (°C)', 
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("soilT2100 - soilT2010", outer=F)

# Compute the difference matrix between 2100 and 2010 for eco
difference_matrix_eco <- soilC_eco_matrices[['2100']] - soilC_eco_matrices[['2010']]
difference_matrix_ecoevo <- soilC_ecoevo_matrices[['2100']] - soilC_ecoevo_matrices[['2010']]

# Calculate evo effect for 2100
evoeffect_matrix = ((soilC_ecoevo_matrices[['2100']] - soilC_ecoevo_matrices[['2010']]) - (soilC_eco_matrices[['2100']] - soilC_eco_matrices[['2010']])) / ((soilC_eco_matrices[['2100']] - soilC_eco_matrices[['2010']]))

# Calculate evo effect 2100 normalized by temperature change
evoeffect_normalized = evoeffect_matrix / T_change_2010_2100

# Plotting ecoevo differences between 2100 and 2010
old.par <- par(mar = c(0, 0, 0, 0))
par(old.par)
nclitt <- nc_open("fVeglitter_Lmon_CCSM4_rcp85_r1i1p1_200501-210012.nc")
lon <- ncvar_get(nclitt, "lon") - 180
lat <- ncvar_get(nclitt, "lat")
brks_diff <- seq(-15000, 15000, by=1)
colfunc_diff <- colorRampPalette(c(
  "#a50026","#d73027","#f46d43","#fdae61","#fee090","gray48",
  "#e0f3f8","#abd9e9","#74add1","#4575b4","#313695"))(length(brks_diff)-1)
fields::image.plot(x=lon, y=lat, soilC_eco_matrices[['2100']] - soilC_eco_matrices[['2010']], ann=FALSE, breaks=brks_diff,
                   horizontal=TRUE, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='Ecoevo Difference (Pg)',
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc_diff,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("SoilC2100 - SoilC2010, eco, c0=1.17")

# Plotting the evo effect matrix for 2100
old.par <- par(mar = c(0, 0, 0, 0))
par(old.par)
nclitt <- nc_open("fVeglitter_Lmon_CCSM4_rcp85_r1i1p1_200501-210012.nc")
lon <- ncvar_get(nclitt, "lon") - 180
lat <- ncvar_get(nclitt, "lat")
par(xpd=TRUE, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
brks <- seq(-2, 2, by=0.01)
colfunc <- colorRampPalette(c(
  "#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","gray48",
  "#fee090","#fdae61","#f46d43","#d73027","#a50026"))(length(brks)-1)
fields::image.plot(x=lon, y=lat, evoeffect_matrix, ann=FALSE, breaks=brks,
                   horizontal=TRUE, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='Evo Effect',
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("Evo effect, eco, c0=1.17")

# Plotting the normalized evo effect for 2100
old.par <- par(mar = c(0, 0, 0, 0))
par(old.par)
nclitt <- nc_open("fVeglitter_Lmon_CCSM4_rcp85_r1i1p1_200501-210012.nc")
lon <- ncvar_get(nclitt, "lon") - 180
lat <- ncvar_get(nclitt, "lat")
par(xpd=TRUE, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
brks <- seq(-1, 1, by=0.01)
colfunc <- colorRampPalette(c(
  "#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","gray48",
  "#fee090","#fdae61","#f46d43","#d73027","#a50026"))(length(brks)-1)
fields::image.plot(x=lon, y=lat, evoeffect_normalized, ann=FALSE, breaks=brks,
                   horizontal=TRUE, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='Evo Effect',
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("Evo effect normalized, eco, c0=1.17")

max(evoeffect_normalized,na.rm=T)

# # Prepare data for plotting
# years_plot <- as.character(seq(2010, 2100, by=10))
# colors <- c("cornflowerblue", "darkkhaki", "darkgoldenrod1")
# line_types <- c("dashed", "solid")
# 
# # Plot Change in total carbon as a function of years
# # Plot Change in total carbon as a function of years
# old.par = par(mar = c(5, 4, 4, 10), xpd=TRUE)
# par(xpd=TRUE, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1)
# plot(NULL, xlim=c(2010, 2100), ylim=c(-350,10), 
#      xlab="Years", ylab="Change in total soil carbon (Pg)", type="n")
# 
# for (i in 1:length(c0_values)) {
#   lines(years_plot, results[[i]]$eco_differences, type="b", pch=16, col=colors[i], lty=2)  # Dashed lines for 'eco'
#   lines(years_plot, results[[i]]$ecoevo_differences, type="b", pch=16, col=colors[i], lty=1)  # Solid lines for 'ecoevo'
# }
# legend("topright", inset=c(-0.5,0), legend=c("Eco c0=1.1", "Ecoevo c0=1.1", "Eco c0=1.17", "Ecoevo c0=1.17", "Eco c0=1.24", "Ecoevo c0=1.24"), col=rep(colors, each=2), lty=rep(line_types, 3), pch=16)
# 
# # Add a second y-axis for the average global soil temperature
# par(new=TRUE)
# plot(years_plot, avg_soilT_plot, type="b", pch=16, col="black", axes=FALSE, xlab="", ylab="", ylim=range(pretty(avg_soilT_plot)))
# axis(side=4, at=pretty(range(avg_soilT_plot)))
# mtext("Average Global Soil Temperature (°C)", side=4, line=3)

