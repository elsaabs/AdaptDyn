setwd("/Volumes/ElsaSSD2022/4.LSCE/PhD/Abby_chapter_5/2023_Oct_recovered_from_old_harddrive/output")
# Packages
library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
library(rootSolve)
library(deSolve)
library(matrixcalc)

# Load soil temperature data
soilT_df <- read.table('soilT2100', sep=" ", dec='.', header=F)
soilT <- as.matrix(soilT_df)

# Load map of stable points
stable_aat_df <- read.table('stable_c0=1.17', sep=" ", dec='.', header=F)
stable_aat <- as.matrix(stable_aat_df)

# Fixed parameter values
Vsol=0.1
litt = 5e-3*Vsol
V0D = 1.15e+6
EvD = 36.1
K0D = 3.3e+4*Vsol
EkD = 9.7
eC = 1e-6
eD = 1e-2
dM = 2e-4
dZ = 2e-3
V0U = 1e+5
EvU = 38
K0U = 1.6e+3*Vsol
EkU = 21
gM = 0.31
gZ = 0.4
c0 = 1.17

# Initialize matrices
AA = matrix(0L, nrow=288, ncol=192)
VmD = AA
KmD = AA
VmU = AA
KmU = AA
phi = AA
delta = AA
K1 = AA
K2 = AA
K3 = AA
K4 = AA
K5 = AA
K6 = AA
deltaM = AA
deltaZ = AA
C = AA
D = AA
M = AA
Z = AA
param = c(0,0,0,0,0,0,0,0,0,0,0,0)
eig = c(0,0,0,0)
Csol = AA

# Differential equations of the model
eqn <- function (t, state, pars) {
  SOC <- state[1]
  DOC <- state[2]
  MBC <- state[3]
  ENZ <- state[4]
  
  ICp <- pars[1]
  VmDp <- pars[2]
  KmDp <- pars[3]
  eCp <- pars[4]
  VmUp <- pars[5]
  KmUp <- pars[6]
  eDp <- pars[7]
  dZp <- pars[8]
  dMp <- pars[9]
  gMp <- pars[10]
  gZp <- pars[11]
  phip <- pars[12]
  
  dSOC <- ICp - VmDp * SOC * ENZ / (KmDp + SOC) - eCp * SOC
  dDOC <- VmDp * SOC * ENZ / (KmDp + SOC) - VmUp * DOC * MBC / (KmUp + DOC) - eDp * DOC + dZp * ENZ + dMp * MBC
  dMBC <- gMp * (1 - phip) * VmUp * DOC * MBC / (KmUp + DOC) - dMp * MBC
  dENZ <- gZp * phip * VmUp * DOC * MBC / (KmUp + DOC) - dZp * ENZ
  
  list(c(dSOC, dDOC, dMBC, dENZ))
}

# Calculation of SOC stock
for (i in 1:dim(AA)[1]) {
  for (j in 1:dim(AA)[2]) {
    # Calculate only in locations stable across the 20 (eco-evo treatment x year)
    if (stable_aat[i,j] == 0) {
      Csol[i,j] = NA
      phi[i,j] = NA
    }else{
      VmD[i, j] = V0D * exp(-EvD / (8.314e-3 * (soilT[i, j] + 273)))
      KmD[i, j] = K0D * exp(-EkD / (8.314e-3 * (soilT[i, j] + 273)))
      VmU[i, j] = V0U * exp(-EvU / (8.314e-3 * (soilT[i, j] + 273)))
      KmU[i, j] = K0U * exp(-EkU / (8.314e-3 * (soilT[i, j] + 273)))
      phi[i, j] = 1 - dM / (VmU[i, j] * gM) - 1 / c0
      delta[i, j] = 1 - (1 - phi[i, j]) * gM - phi[i, j] * gZ
      K1[i, j] = (1 - phi[i, j]) * gM * VmU[i, j] - dM
      K2[i, j] = eC * KmD[i, j] * K1[i, j] - eD * KmU[i, j] * dM
      K3[i, j] = (phi[i, j] * gZ * VmD[i, j] - dZ * delta[i, j]) * K1[i, j]
      K4[i, j] = K3[i, j] * (litt - eC * KmD[i, j]) + phi[i, j] * gZ * VmD[i, j] * K2[i, j]
      K5[i, j] = K3[i, j] * K1[i, j] * delta[i, j] * dZ * litt * eC * KmD[i, j]
      deltaM[i, j] = 2 * dZ * delta[i, j] * K2[i, j]
      deltaZ[i, j] = -2 * dZ * delta[i, j] * K2[i, j]
  
      if (is.na(K4[i, j]^2 - 4 * K5[i, j])) {
        Csol[i, j] = NA
        phi[i,j] = NA
      } else {
        if (K4[i, j]^2 - 4 * K5[i, j] >= 0) {
          C[i, j] = (K4[i, j] - sqrt(K4[i, j]^2 - 4 * K5[i, j])) / (2 * eC * K3[i, j])
          D[i, j] = dM * KmU[i, j] / K1[i, j]
          M[i, j] = (gM * (1 - phi[i, j])) * (K4[i, j] - deltaM[i, j] + sqrt(K4[i, j]^2 - 4 * K5[i, j])) / (2 * dM * K3[i, j] * delta[i, j])
          Z[i, j] = (gZ * phi[i, j]) * (K4[i, j] + deltaZ[i, j] + sqrt(K4[i, j]^2 - 4 * K5[i, j])) / (2 * dZ * K3[i, j] * delta[i, j])
          
          param <- c(ICp = litt, VmDp = VmD[i, j], KmDp = KmD[i, j], eCp = eC, VmUp = VmU[i, j], KmUp = KmU[i, j], eDp = eD, dZp = dZ, dMp = dM, gMp = gM, gZp = gZ, phip = phi[i, j])
          state <- c(SOC = C[i, j], DOC = D[i, j], MBC = M[i, j], ENZ = Z[i, j])
          
          jacobian_matrix <- jacobian.full(y = state, func = eqn, parms = param)
          eigenvalues <- eigen(jacobian_matrix)$values
          
          if (C[i, j] >= 0 & D[i, j] >= 0 & M[i, j] >= 0 & Z[i, j] >= 0 & all(Re(eigenvalues) <= 0)) {
            Csol[i, j] = C[i, j] + M[i, j]
          } else {
            Csol[i, j] = NA
            phi[i,j] = NA
          }
        } else {
          Csol[i, j] = NA
          phi[i,j] = NA
        }
      }
    }
  }
}

# length(which(!is.na(Csol)))

# Global map of phi
old.par = par(mar = c(0, 0, 0, 0))
par(old.par)
nclitt = nc_open("fVeglitter_Lmon_CCSM4_rcp85_r1i1p1_200501-210012.nc")
lon = ncvar_get(nclitt, "lon")-180
lat = ncvar_get(nclitt, "lat")
par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
brks = seq(0.04, 0.13, by=0.001)
colfunc <- colorRampPalette(c(
  "#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf",
  "#fee090","#fdae61","#f46d43","#d73027","#a50026"))(length(brks)-1)
fields::image.plot(x=lon, y=lat, phi2010, ann=F, breaks=brks,
                   horizontal=T, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='SOCecoevo (g m-2)',
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("phi 2010", outer=F)

# # Save phi
# phi = as.matrix(phi)
# write.table(phi, file = "phiecoevo2100_c0=1.17", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

# Global map of phi2100 - phi2010
phi2010_df <- read.table('phiecoevo2010_c0=1.17',sep= " ",dec='.',header=F)
phi2010 <- as.matrix(phi2010_df)
phi2100_df <- read.table('phiecoevo2100_c0=1.17',sep= " ",dec='.',header=F)
phi2100 <- as.matrix(phi2100_df)
# Plot
old.par = par(mar = c(0, 0, 0, 0))
par(old.par)
nclitt = nc_open("fVeglitter_Lmon_CCSM4_rcp85_r1i1p1_200501-210012.nc")
lon = ncvar_get(nclitt, "lon")-180
lat = ncvar_get(nclitt, "lat")
par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
brks = seq(-0.03, 0.03, by=1e-5)
colfunc <- colorRampPalette(c(
  "#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","gray48",
  "#fee090","#fdae61","#f46d43","#d73027","#a50026"))(length(brks)-1)
fields::image.plot(x=lon, y=lat, phi2100-phi2010, ann=F, breaks=brks,
                   horizontal=T, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='phi',
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("phi2100 - phi2010", outer=F)

max(phi2100-phi2010,na.rm=T)

# # Global map of soil C
# old.par = par(mar = c(0, 0, 0, 0))
# par(old.par)
# nclitt = nc_open("fVeglitter_Lmon_CCSM4_rcp85_r1i1p1_200501-210012.nc")
# lon = ncvar_get(nclitt, "lon")-180
# lat = ncvar_get(nclitt, "lat")
# par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
# brks = seq(0, 18000, by=10)
# colfunc = colorRampPalette(c(
#   "white","blue","green","yellow","orange","red"))(length(brks)-1)
# fields::image.plot(x=lon, y=lat, Csol*1000, ann=F, breaks=brks,
#                    horizontal=T, smallplot=c(0.16, 0.92, 0.01, 0.04),
#                    legend.args=list(text='SOCecoevo (g m-2)',
#                                     side=1, font=1, line=2, cex=1),
#                    col=colfunc,
#                    xlab="Longitude", ylab="Latitude")
# mtext(side=1, text="Longitude", line=2.2, font=1)
# mtext(side=2, text="Latitude", line=2.2, font=1)
# title("ecoevo 2010", outer=F)

# # Save Ctot
# Ctot = Csol*1000
# Ctot = as.matrix(Ctot)
# write.table(Ctot, file = "soilC_c0=1.24_ecoevo_2100", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

# # Map change in soil C without evolution for c0=1.17
# soilCecoevo2010_df <- read.table('soilC_c0=1.17_ecoevo_2010',sep= " ",dec='.',header=F)
# soilCecoevo2010 <- as.matrix(soilCecoevo2010_df)
# soilCecoevo2100_df <- read.table('soilC_c0=1.17_ecoevo_2100',sep= " ",dec='.',header=F)
# soilCecoevo2100 <- as.matrix(soilCecoevo2100_df)
# 
# min(soilCecoevo2100-soilCeco2010,na.rm=T)
# 
# old.par <- par(mar = c(0, 0, 0, 0))
# par(old.par)
# nclitt <- nc_open("fVeglitter_Lmon_CCSM4_rcp85_r1i1p1_200501-210012.nc")
# lon <- ncvar_get(nclitt, "lon") - 180
# lat <- ncvar_get(nclitt, "lat")
# par(xpd=TRUE, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
# brks <- seq(-15000, 15000, by=0.1)
# colfunc <- colorRampPalette(c(
#   "#a50026","#d73027","#f46d43","#fdae61","#fee090","gray48",
#   "#e0f3f8","#abd9e9","#74add1","#4575b4","#313695"))(length(brks)-1)
# fields::image.plot(x=lon, y=lat, soilCecoevo2100-soilCecoevo2010, ann=FALSE, breaks=brks,
#                    horizontal=TRUE, smallplot=c(0.16, 0.92, 0.01, 0.04),
#                    legend.args=list(text='Change in Soil Carbon (Pg)',
#                                     side=1, font=1, line=2, cex=1),
#                    col=colfunc,
#                    xlab="Longitude", ylab="Latitude")
# mtext(side=1, text="Longitude", line=2.2, font=1)
# mtext(side=2, text="Latitude", line=2.2, font=1)

