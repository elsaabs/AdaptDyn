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

#Get phiecoevo2010
phi_df <- read.table('phiecoevo2010_c0=1.24',sep= " ",dec='.',header=F)
phi <- as.matrix(phi_df)

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

# Initialize matrices
AA = matrix(0L, nrow=288, ncol=192)
VmD = AA
KmD = AA
VmU = AA
KmU = AA
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
stable = AA

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
    VmD[i, j] = V0D * exp(-EvD / (8.314e-3 * (soilT[i, j] + 273)))
    KmD[i, j] = K0D * exp(-EkD / (8.314e-3 * (soilT[i, j] + 273)))
    VmU[i, j] = V0U * exp(-EvU / (8.314e-3 * (soilT[i, j] + 273)))
    KmU[i, j] = K0U * exp(-EkU / (8.314e-3 * (soilT[i, j] + 273)))
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
      stable[i, j] = 0
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
          stable[i, j] = 1
        } else {
          Csol[i, j] = NA
          stable[i, j] = 0
        }
      } else {
        Csol[i, j] = NA
        stable[i, j] = 0
      }
    }
  }
}




# # Total soil C map
# old.par = par(mar = c(0, 0, 0, 0))
# par(old.par)
# nclitt = nc_open("fVeglitter_Lmon_CCSM4_rcp85_r1i1p1_200501-210012.nc")
# lon = ncvar_get(nclitt, "lon")-180
# lat = ncvar_get(nclitt, "lat")
# par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
# brks = seq(0, 40000, by=10)
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


# #Save stable
# stable = as.matrix(stable)
# write.table(stable, file = "stable2100_eco_c0=1.24", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

# #Retrieve stable_eco and stable_ecoevo to find points where eco is table but not ecoevo
# stableeco_df <- read.table('stable2010_eco',sep= " ",dec='.',header=F)
# stableeco <- as.matrix(stableeco_df)
# stableecoevo_df <- read.table('stable2010_ecoevo',sep= " ",dec='.',header=F)
# stableecoevo <- as.matrix(stableecoevo_df)
# which(stableecoevo==0 & stableeco==1)
# linear_index=36500
# i_index=(linear_index - 1) %% dim(AA)[1] + 1
# j_index=(linear_index - 1)/dim(AA)[1] + 1

# Create map with points that are always stable (between eco/ecoevo and years)
# Load the 20 stable tables
# Initialize lists to store the matrices
stableeco <- list()
stableecoevo <- list()

# Load stableeco matrices
for (year in seq(2010, 2100, by=10)) {
  filename <- paste0('stable', year, '_eco_c0=1.24')
  stableeco[[as.character(year)]] <- as.matrix(read.table(filename, sep=" ", dec='.', header=F))
}

# Load stableecoevo matrices
for (year in seq(2010, 2100, by=10)) {
  filename <- paste0('stable', year, '_ecoevo_c0=1.24')
  stableecoevo[[as.character(year)]] <- as.matrix(read.table(filename, sep=" ", dec='.', header=F))
}

# Combine all matrices into a single list
stable_matrices <- c(stableeco, stableecoevo)

# Initialize stable_aat matrix with the same dimensions as the input matrices
stable_aat <- matrix(0, nrow=nrow(stableeco[['2010']]), ncol=ncol(stableeco[['2010']]))
# Iterate through each element and check if it is stable across all matrices
for (i in 1:nrow(stable_aat)) {
  for (j in 1:ncol(stable_aat)) {
    stable_aat[i, j] <- all(sapply(stable_matrices, function(mat) mat[i, j] == 1))
  }
}

# Convert logical TRUE/FALSE to 1/0
stable_aat <- as.numeric(stable_aat)
stable_aat <- matrix(stable_aat, nrow=nrow(stableeco[['2010']]), ncol=ncol(stableeco[['2010']]))

# Number of stable points in stable_aat
length(which(stable_aat==1))

# Global map of stable points
old.par = par(mar = c(0, 0, 0, 0))
par(old.par)
nclitt = nc_open("fVeglitter_Lmon_CCSM4_rcp85_r1i1p1_200501-210012.nc")
lon = ncvar_get(nclitt, "lon")-180
lat = ncvar_get(nclitt, "lat")
par(xpd=T, mfrow=c(1,1), lwd=1, font.lab=1, font.axis=1, oma=c(3,0,0,0))
brks = seq(0, 1, by=0.1)
colfunc = colorRampPalette(c(
  "blue","red"))(length(brks)-1)
fields::image.plot(x=lon, y=lat, stable_aat, ann=F, breaks=brks,
                   horizontal=T, smallplot=c(0.16, 0.92, 0.01, 0.04),
                   legend.args=list(text='stability (1: yes, 0:no)',
                                    side=1, font=1, line=2, cex=1),
                   col=colfunc,
                   xlab="Longitude", ylab="Latitude")
mtext(side=1, text="Longitude", line=2.2, font=1)
mtext(side=2, text="Latitude", line=2.2, font=1)
title("c0=1.24")

# Save stable 2010
stable_aat = as.matrix(stable_aat)
write.table(stable_aat, file = "stable_c0=1.24", append = FALSE, quote = TRUE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)

length(which(!is.na(Csol)))







