### LAKE GUS biomarker and PSM plots

setwd('/Users/allisoncluett/Dropbox/LakeGus_Temps/JGR-Biogeosci Submission/Data')

library(dplyr); library(magrittr)
library(corrplot); library(FactoMineR); library(factoextra)
library(ggplot2); library(ggtern); library(ggridges)
library(patchwork); library(cowplot)

GUS_biomarker <- read.csv("GUS_GDGT_Master.csv")
ages <- GUS_biomarker$Age_yr

# -------------- Some GDGT calculations ------------------

brGDGTdf <- GUS_biomarker[,c(21:35)]
brGDGTconc <- ((brGDGTdf/GUS_biomarker$IntStd)* GUS_biomarker$Std_ug)/GUS_biomarker$SedMass_g
brGDGT.Total <- rowSums(brGDGTconc, na.rm = TRUE)
brGDGTfa <- brGDGTdf/rowSums(brGDGTdf)

brGDGT.TPH <- data.frame(tetra = (brGDGTfa$brGDGT.Ia + brGDGTfa$brGDGT.Ib + brGDGTfa$brGDGT.Ic),
                         penta = (brGDGTfa$brGDGT.IIa.5Me + brGDGTfa$brGDGT.IIa.6Me + brGDGTfa$brGDGT.IIb.5Me +
                                    brGDGTfa$brGDGT.IIb.6Me + brGDGTfa$brGDGT.IIc.5Me + brGDGTfa$brGDGT.IIc.6Me),
                         hexa = (brGDGTfa$hexa <- (brGDGTfa$brGDGT.IIIa.5Me + brGDGTfa$brGDGT.IIIa.6Me +
                                                          brGDGTfa$brGDGT.IIIb.5Me + brGDGTfa$brGDGT.IIIb.6Me +
                                                          brGDGTfa$brGDGT.IIIc.5Me + brGDGTfa$brGDGT.IIIc.6Me)))

GDGT0.4 <- GUS_biomarker$GDGT.0/GUS_biomarker$GDGT.4

MBT.5me <- (GUS_biomarker$brGDGT.Ia + GUS_biomarker$brGDGT.Ib + GUS_biomarker$brGDGT.Ic)/
  (GUS_biomarker$brGDGT.Ia + GUS_biomarker$brGDGT.Ib + GUS_biomarker$brGDGT.Ic + GUS_biomarker$brGDGT.IIa.5Me +
     GUS_biomarker$brGDGT.IIb.5Me + GUS_biomarker$brGDGT.IIc.5Me + GUS_biomarker$brGDGT.IIIa.5Me)

GUS_biomarker$UK37 <- (GUS_biomarker$LCA37.2 - GUS_biomarker$LCA37.4)/(GUS_biomarker$LCA37.2 + GUS_biomarker$LCA37.3a + GUS_biomarker$LCA37.3b + GUS_biomarker$LCA37.4)


## -------- FIG 2: brGDGT source discrimination ---------

# Scale brGDGT concentrations such that most abundant compound in each sample = 1
brGDGTscaled <- brGDGTconc

for (i in 1:nrow(brGDGTscaled)) {
  max_val <- max(brGDGTscaled[i,], na.rm = T)
  brGDGTscaled[i,] <- brGDGTconc[i,] / max_val
}

res.pca <- PCA(brGDGTscaled, graph = FALSE)
pca.res <- data.frame(ages, res.pca$ind$coord[,1],res.pca$ind$coord[,2])
eig.val <- get_eigenvalue(res.pca)
var <- get_pca_var(res.pca)
ind <- get_pca_ind(res.pca)

# --- FIG 2A ---

midpt1 = mean(c(max(GDGT0.4), min(GDGT0.4))) # for color scaling

CaldCrenTern <- ggtern(data=brGDGT.TPH,aes(x=tetra,y=penta,z=hexa, color = GDGT0.4 )) + geom_point() +
  scale_colour_gradientn(colours = c("#07beb8", "#f4d06f", "#ee964b"))
CaldCrenTern + theme_minimal()

# --- FIG 2B ---

CaldCrenBiPlot <- fviz_pca_biplot(res.pca, geom.ind = "point", pointsize=3, col.ind = GDGT0.4,  legend.title = "Cald/Cren")+
  scale_color_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1) +
  geom_point(data = data.frame(res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),1]),
             aes(x = res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),1],
                 y = res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),2]), pch = 4, size = 5, col = 'black')

CaldCrenBiPlot + theme_minimal()

# --- FIG 2 C-H ---

set.seed(79)

PC1t <- ggplot() + geom_point(aes(x = res.pca$ind$coord[,1], y = brGDGT.TPH$tetra, col = GDGT0.4)) +
  scale_color_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1) +
  geom_point(aes(x = res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),1], y = brGDGT.TPH$tetra[which(brGDGT.TPH$hexa < 0.3)]), pch = 4, size = 3) +
  theme_linedraw() + theme(legend.position = 'none') + scale_x_continuous(position = 'top') + ylab('% Tetra') + xlab('PC1 Coordinate') + ylim(0.2, 0.35)

PC1p <- ggplot() + geom_point(aes(x = res.pca$ind$coord[,1], y = brGDGT.TPH$penta, col = GDGT0.4)) +
  scale_color_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1) +
  geom_point(aes(x = res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),1], y = brGDGT.TPH$penta[which(brGDGT.TPH$hexa < 0.3)]), pch = 4, size = 3) +
  theme_linedraw() + theme(legend.position = 'none', axis.title.x=element_blank(),
                          axis.text.x=element_blank(), axis.ticks.x=element_blank())  + ylab('% Penta') + ylim(0.39, 0.48)

PC1h <- ggplot() + geom_point(aes(x = res.pca$ind$coord[,1], y = brGDGT.TPH$hexa, col = GDGT0.4)) +
  scale_color_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1) +
  geom_point(aes(x = res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),1], y = brGDGT.TPH$hexa[which(brGDGT.TPH$hexa < 0.3)]), pch = 4, size = 3) +
  theme_linedraw() + theme(legend.position = 'none') + scale_x_continuous(position = 'bottom') + ylab('% Hexa') + xlab('PC1 Coordinate') + ylim(0.22, 0.38)

PC2t <- ggplot() + geom_point(aes(x = res.pca$ind$coord[,2], y = brGDGT.TPH$tetra, col = GDGT0.4)) +
  scale_color_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1) +
  geom_point(aes(x = res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),2], y = brGDGT.TPH$tetra[which(brGDGT.TPH$hexa < 0.3)]), pch = 4, size = 3) +
  theme_linedraw() + theme(legend.position = 'none')  + ylab('% Tetra') + scale_y_continuous(position = 'right', limits = c(0.2, 0.35)) + xlab('PC2 Coordinate') +
  scale_x_continuous(position = 'top')

PC2p <- ggplot() + geom_point(aes(x = res.pca$ind$coord[,2], y = brGDGT.TPH$penta, col = GDGT0.4)) +
  scale_color_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1) +
  geom_point(aes(x = res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),2], y = brGDGT.TPH$penta[which(brGDGT.TPH$hexa < 0.3)]), pch = 4, size = 3) +
  theme_linedraw() + theme(legend.position = 'none') + theme(legend.position = 'none', axis.title.x=element_blank(),
                                                            axis.text.x=element_blank(),
                                                            axis.ticks.x=element_blank(),
                                                            axis.line.x.bottom = element_blank())  + ylab('% Penta') +
  scale_y_continuous(position = 'right', limits = c(0.39, 0.48))

PC2h <- ggplot() + geom_point(aes(x = res.pca$ind$coord[,2], y = brGDGT.TPH$hexa, col = GDGT0.4)) +
  scale_color_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1) +
  geom_point(aes(x = res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),2], y = brGDGT.TPH$hexa[which(brGDGT.TPH$hexa < 0.3)]), pch = 4, size = 3) +
  theme_linedraw() + theme(legend.position = 'none') + ylab('% Hexa') + scale_y_continuous(position = 'right', limits = c(0.22, 0.38)) + xlab('PC2 Coordinate')

patchGrid <- ((PC1t + PC2t)/
  (PC1p + PC2p)/ (PC1h + PC2h))

patchGrid & theme(panel.border = element_rect(size = 2), plot.margin = margin(2,2,2,2, 'pt'),
                  axis.ticks.length = unit(4, 'pt'), axis.ticks = element_line(size = 0.8))

## -------- FIG S5 ---------

# --- FIG S5A ---

AgeTern <- ggtern(data=brGDGT.TPH,aes(x=tetra,y=penta,z=hexa, color =ages ))+ geom_point() + scale_colour_gradientn(colours = c("#07beb8", "#f4d06f", "#ee964b"))
AgeTern + theme_minimal()

# --- FIG S5 B-D ---

hexa <- ggplot() + geom_line(aes(x = ages, y = brGDGT.TPH$hexa)) + geom_point(aes(x = ages,  y = brGDGT.TPH$hexa, fill = GDGT0.4), col = 'black', pch = 21) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1) + scale_x_reverse(limits = c(8500, -100),
                                                                                                               position = 'top') +
  xlab('Age (cal yr BP)') + ylab('Hexa') + geom_hline(yintercept = 0.3, lty = 'dashed') +
  geom_point(data = data.frame(res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),1]),
             aes(x = ages[which(brGDGT.TPH$hexa < 0.3)],
                 y = brGDGT.TPH$hexa[which(brGDGT.TPH$hexa < 0.3)]), pch = 4, size = 3, col = 'black')

pc1 <- ggplot() + geom_line(aes(x = ages,  y = res.pca$ind$coord[,1])) + geom_point(aes(x = ages,  y = res.pca$ind$coord[,1], fill = GDGT0.4), col = 'black', pch = 21) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1) + scale_x_reverse(limits = c(8500, -100)) +
  xlab('Age (cal yr BP)') + ylab('PC1') +  geom_point(data = data.frame(res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),1]),
                                                      aes(x = ages[which(brGDGT.TPH$hexa < 0.3)],
                                                          y = res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),1]), pch = 4, size = 3, col = 'black')

pc2 <- ggplot() + geom_line(aes(x = ages,  y = res.pca$ind$coord[,2])) + geom_point(aes(x = ages,  y = res.pca$ind$coord[,2], fill = GDGT0.4), col = 'black', pch = 21) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1) + scale_x_reverse(limits = c(8500, -100)) +
  xlab('Age (cal yr BP)') + ylab('PC2')+  geom_point(data = data.frame(res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),1]),
                                                     aes(x = ages[which(brGDGT.TPH$hexa < 0.3)],
                                                         y = res.pca$ind$coord[which(brGDGT.TPH$hexa < 0.3),2]), pch = 4, size = 3, col = 'black')

hexa/pc1/pc2 & theme_linedraw() & theme(panel.border = element_rect(size = 2), plot.margin = margin(2,2,2,2, 'pt'),
                                        axis.ticks.length = unit(4, 'pt'), axis.ticks = element_line(size = 0.8),
                                        legend.position = 'NULL')


## -------- FIG 3: Calibrations ---------

source('CalibrationPriors.R')

T2IndexCalib <- function(Temp, slope, int){slope * Temp + int}
Index2TCalib <- function(Index, slope, int){(Index - int)/slope}

GUS_biomarker$MBT.5me <- (GUS_biomarker$brGDGT.Ia + GUS_biomarker$brGDGT.Ib + GUS_biomarker$brGDGT.Ic)/
  (GUS_biomarker$brGDGT.Ia + GUS_biomarker$brGDGT.Ib + GUS_biomarker$brGDGT.Ic + GUS_biomarker$brGDGT.IIa.5Me +
     GUS_biomarker$brGDGT.IIb.5Me + GUS_biomarker$brGDGT.IIc.5Me + GUS_biomarker$brGDGT.IIIa.5Me)

# East African Lakes
GUS_biomarker$Russell2018T <- -1.21 + 32.42 * GUS_biomarker$MBT.5me # RMSE = 2.44

# Lake 578
GUS_biomarker$Zhao2020T <- -1.82 + 56.06 * GUS_biomarker$MBT.5me # RMSE = 0.52

unflagged_subset <- GUS_biomarker[which(brGDGT.TPH$hexa > 0.3),]

GDGTtemp <- ggplot(data = unflagged_subset) + geom_ribbon(aes(x= Age_yr, ymin = Russell2018T - 2.44, ymax = Russell2018T + 2.44), alpha = 0.3, fill = 'gold', col = 'gold', lty = 'dotted') +
  geom_ribbon(aes(x= Age_yr, ymin = Zhao2020T - 0.52, ymax = Zhao2020T + 0.52), alpha = 0.3, fill = 'red', col = 'red', lty = 'dotted') +
  geom_line(aes(x= Age_yr, y = Russell2018T), lwd = 0.8, col = 'gold') + geom_line(aes(x= Age_yr, y = Zhao2020T), lwd = 0.8, col = 'red') +
  xlim(9000, -100)  + xlab('Age (cal yr BP)') + ylab('Calibrated Temp (°C)') + theme_minimal() + theme(text = element_text(size = 14))

MBTscreen <- ggplot(data = unflagged_subset) +
  geom_ribbon(aes(x= Age_yr, ymin = MBT.5me - 0.00322, ymax = MBT.5me + 0.00322), alpha = 0.3, fill = 'black') +
  geom_line(aes(x= Age_yr, y = MBT.5me), lwd = 0.8, col = 'black') + geom_point(aes(x= Age_yr, y = MBT.5me)) +
  xlim(9000, -100) + labs(x="") + guides(x = "none") + ylab('MBT5Me') + theme_minimal() + theme(text = element_text(size = 14)) +
  theme(panel.grid.minor.y = element_blank())

LCAsubset <- GUS_biomarker[which(!is.na(GUS_biomarker$UK37)),]
LCAsubset$minLCApeak <- LCAsubset[,c('LCA37.2', 'LCA37.3a', 'LCA37.3b', 'LCA37.4')] %>%  apply(FUN = min, MARGIN = 1, na.rm = TRUE)

LCAs <- ggplot(data = LCAsubset) +
  geom_ribbon(aes(x= Age_yr, ymin = UK37 - (-0.000118 * minLCApeak + 0.0106), ymax = UK37 + (-0.000118 * minLCApeak + 0.0106)), alpha = 0.3, fill = 'black') +
  geom_line(aes(x= Age_yr, y =UK37), lwd = 0.8, col = 'black') + geom_point(aes(x= Age_yr, y = UK37)) +
  xlim(9000, -100) + labs(x="") + guides(x = "none") + ylab('UK37') + theme_minimal() + theme(text = element_text(size = 14)) +
  theme(panel.grid.minor.y = element_blank())

LCAtemp <- ggplot(data = LCAsubset) +
  geom_ribbon(aes(x= Age_yr, ymin = Index2TCalib(UK37, slope = calibData$Longo16$slopeMean, int= calibData$Longo16$intMean) - calibData$Longo16$RMSEtemp,
                  ymax = Index2TCalib(UK37, slope = calibData$Longo16$slopeMean, int= calibData$Longo16$intMean) + calibData$Longo16$RMSEtemp),
              alpha = 0.3, fill = '#339989', col = '#339989', lty = 'dotted') +
  geom_ribbon(aes(x= Age_yr, ymin = Index2TCalib(UK37, slope = calibData$DAndrea11$slopeMean, int= calibData$DAndrea11$intMean) - calibData$DAndrea11$RMSEtemp,
                  ymax = Index2TCalib(UK37, slope = calibData$DAndrea11$slopeMean, int= calibData$DAndrea11$intMean) + calibData$DAndrea11$RMSEtemp),
              alpha = 0.3, fill = '#0C6291', col = '#0C6291', lty = 'dotted') +
  geom_ribbon(aes(x= Age_yr, ymin = Index2TCalib(UK37, slope = calibData$DAndrea16$slopeMean, int= calibData$DAndrea16$intMean) - calibData$DAndrea16$RMSEtemp,
                  ymax = Index2TCalib(UK37, slope = calibData$DAndrea16$slopeMean, int= calibData$DAndrea16$intMean) + calibData$DAndrea16$RMSEtemp),
              alpha = 0.3, fill = '#28557A', col = '#28557A', lty = 'dotted') +
  geom_line(aes(x= Age_yr, y = Index2TCalib(UK37, slope = calibData$DAndrea16$slopeMean, int= calibData$DAndrea16$intMean)), lwd = 0.8, col = '#28557A') +
  geom_line(aes(x= Age_yr, y = Index2TCalib(UK37, slope = calibData$Longo16$slopeMean, int= calibData$Longo16$intMean)), lwd = 0.8, col = '#339989') +
  geom_line(aes(x= Age_yr, y = Index2TCalib(UK37, slope = calibData$DAndrea11$slopeMean, int= calibData$DAndrea11$intMean)), lwd = 0.8, col = '#0C6291') +
  xlim(9000, -100)  + xlab('Age (cal yr BP)') + ylab('Calibrated Temp (°C)') + theme_minimal() + theme(text = element_text(size = 14)) +
  ylim(2.5, 22.5) +  scale_y_continuous(position = "right")

MBTscreen  + LCAs + GDGTtemp  + ylim(2.5, 22)  + LCAtemp  + ylim(2.5, 22)  +
  plot_layout(heights = c(2,  3)) & theme_linedraw() & theme(panel.border = element_rect(size = 2), plot.margin = margin(2,2,2,2, 'pt'),
                                                             axis.ticks.length = unit(4, 'pt'), axis.ticks = element_line(size = 0.8))


##  ------- FIG 4: Air-water temp sensitivity -----------

# Load model output
load('EnsSummaries-ForcedSurfaceWARM2.Rdata')
load('EnsSummaries-ForcedSurfaceCOLD2.Rdata')
load('EnsSummaries-ForcedSurfaceRunsCOLD2.Rdata')
load('EnsSummaries-ForcedSurfaceRunsWARM2.Rdata')
load('EnsSummaries-ForcedSurfaceMAF3WARM.Rdata')
load('EnsSummaries-ForcedSurfaceMAF3COLD.Rdata')
load('EnsSummaries-ForcedSurfaceJJACOLD2.Rdata')
load('EnsSummaries-ForcedSurfaceJJAWARM2.Rdata')
load('EnsSummaries-Modern.Rdata')

# wrangle the data
modernPt <- colMeans(EnsSummaries) %>% t() %>% as.data.frame()
ForcedRunMAF_JJA <- rbind(ForcedSurfaceEnsMAFWARM[0:20,c('X1_2.5', 'X1_25', 'X1_50', 'X1_75', 'X1_97.5', 'DeltaJJA', 'DeltaAnn')], ForcedSurfaceEnsMAFCOLD[0:20,c('X1_2.5', 'X1_25', 'X1_50', 'X1_75', 'X1_97.5', 'DeltaJJA', 'DeltaAnn')])
ForcedRunMAF_Ann <- rbind(ForcedSurfaceEnsMAFWARM[21:40,], ForcedSurfaceEnsMAFCOLD[21:40,])
ForcedRunJJA_JJA <- rbind(ForcedSurfaceEnsJJAWARM[0:20,c('X1_2.5', 'X1_25', 'X1_50', 'X1_75', 'X1_97.5', 'DeltaJJA', 'DeltaAnn')], ForcedSurfaceEnsJJACOLD[0:20,c('X1_2.5', 'X1_25', 'X1_50', 'X1_75', 'X1_97.5', 'DeltaJJA', 'DeltaAnn')])
ForcedRunJJA_Ann <- rbind(ForcedSurfaceEnsJJAWARM[21:40,c('X1_2.5', 'X1_25', 'X1_50', 'X1_75', 'X1_97.5', 'DeltaJJA', 'DeltaAnn')], ForcedSurfaceEnsJJACOLD[21:40,c('X1_2.5', 'X1_25', 'X1_50', 'X1_75', 'X1_97.5', 'DeltaJJA', 'DeltaAnn')])
ForcedRuns_JJA <- rbind(ForcedSurfaceEnsRunsWARM[0:20,], ForcedSurfaceEnsRunsCOLD[0:20,])
ForcedRuns_Ann <- rbind(ForcedSurfaceEnsRunsWARM[21:40,], ForcedSurfaceEnsRunsCOLD[21:40,])
ForcedRunTs_JJA <- rbind(ForcedSurfaceEnsWARM[0:20,], ForcedSurfaceEnsCOLD[0:20,])
ForcedRunTs_Ann <- rbind(ForcedSurfaceEnsWARM[21:40,], ForcedSurfaceEnsCOLD[21:40,])
ForcedRunTs_JJA$DeltaJJA[0:20] <- ForcedRunTs_JJA$DeltaJJA[0:20] * -1 # this is just to correct a mistake in the ∆T values
ForcedRunTs_Ann$DeltaAnn[0:20] <- ForcedRunTs_Ann$DeltaAnn[0:20] * -1
ForcedRuns_JJA$DeltaJJA[0:20] <- ForcedRuns_JJA$DeltaJJA[0:20] * -1
ForcedRuns_Ann$DeltaAnn[0:20] <- ForcedRuns_Ann$DeltaAnn[0:20] * -1

## ----------- Start plots ---------------

IO2WKtemp <- ggplot() + # temperature in 2 wks following ice out
  geom_ribbon(data = ForcedRunTs_JJA, aes(x= DeltaJJA, ymin = iceout2wkMeanT_2.5, ymax = iceout2wkMeanT_97.5), alpha = 0.3, fill = 'orange') +
  geom_ribbon(data = ForcedRunTs_JJA, aes(x= DeltaJJA, ymin = iceout2wkMeanT_25, ymax = iceout2wkMeanT_75), alpha = 0.3, fill = 'orange') +
  geom_line(data = ForcedRunTs_JJA, aes(x= DeltaJJA, y = iceout2wkMeanT_50), col = 'orange') +
  geom_point(data = ForcedRunTs_JJA, aes(x= DeltaJJA, y = iceout2wkMeanT_50), col = 'orange') +
  geom_ribbon(data = ForcedRunTs_Ann, aes(x= DeltaAnn, ymin = iceout2wkMeanT_2.5, ymax = iceout2wkMeanT_97.5), alpha = 0.3, fill = 'light blue') +
  geom_ribbon(data = ForcedRunTs_Ann, aes(x= DeltaAnn, ymin = iceout2wkMeanT_25, ymax = iceout2wkMeanT_75), alpha = 0.3, fill = 'light blue') +
  geom_line(data = ForcedRunTs_Ann, aes(x= DeltaAnn, y = iceout2wkMeanT_50), col = 'light blue') +
  geom_point(data = ForcedRunTs_Ann, aes(x= DeltaAnn, y = iceout2wkMeanT_50), col = 'light blue') +
  geom_point(aes(x = 0, y = modernPt$iceTs.iceout2wkMeanT_50)) + geom_abline(slope = 1, intercept = modernPt$iceTs.iceout2wkMeanT_50) +
  xlab('Perturbation (°C)') + ylab('2 Wks Following Ice Out (°C)') + ylim(0, 15)

IceFreetemp <- ggplot() + # mean ice free season temperature
  geom_ribbon(data = ForcedRunTs_JJA, aes(x= DeltaJJA, ymin = icefreeMeanT_2.5, ymax = icefreeMeanT_97.5), alpha = 0.3, fill = 'orange') +
  geom_ribbon(data = ForcedRunTs_JJA, aes(x= DeltaJJA, ymin = icefreeMeanT_25, ymax = icefreeMeanT_75), alpha = 0.3, fill = 'orange') +
  geom_line(data = ForcedRunTs_JJA, aes(x= DeltaJJA, y = icefreeMeanT_50), col = 'orange') +
  geom_point(data = ForcedRunTs_JJA, aes(x= DeltaJJA, y = icefreeMeanT_50), col = 'orange') +
  geom_ribbon(data = ForcedRunTs_Ann, aes(x= DeltaAnn, ymin = icefreeMeanT_2.5, ymax = icefreeMeanT_97.5), alpha = 0.3, fill = 'light blue') +
  geom_ribbon(data = ForcedRunTs_Ann, aes(x= DeltaAnn, ymin = icefreeMeanT_25, ymax = icefreeMeanT_75), alpha = 0.3, fill = 'light blue') +
  geom_line(data = ForcedRunTs_Ann, aes(x= DeltaAnn, y = icefreeMeanT_50), col = 'light blue') +
  geom_point(data = ForcedRunTs_Ann, aes(x= DeltaAnn, y = icefreeMeanT_50), col = 'light blue') +
  geom_point(aes(x = 0, y = modernPt$iceTs.icefreeMeanT_50)) + geom_abline(slope = 1, intercept = modernPt$iceTs.icefreeMeanT_50) +
  xlab('Perturbation (°C)') + ylab('Ice Free Season Temp (°C)') + ylim(0, 15)

isothermaltemp <- ggplot() + # temperature during isothermal mixing period
  geom_ribbon(data = ForcedRunTs_JJA, aes(x= DeltaJJA, ymin = isothermalMeanT_2.5, ymax = isothermalMeanT_97.5), alpha = 0.3, fill = 'orange') +
  geom_ribbon(data = ForcedRunTs_JJA, aes(x= DeltaJJA, ymin = isothermalMeanT_25, ymax = isothermalMeanT_75), alpha = 0.3, fill = 'orange') +
  geom_line(data = ForcedRunTs_JJA, aes(x= DeltaJJA, y = isothermalMeanT_50), col = 'orange') +
  geom_point(data = ForcedRunTs_JJA, aes(x= DeltaJJA, y = isothermalMeanT_50), col = 'orange') +
  geom_ribbon(data = ForcedRunTs_Ann, aes(x= DeltaAnn, ymin = isothermalMeanT_2.5, ymax = isothermalMeanT_97.5), alpha = 0.3, fill = 'light blue') +
  geom_ribbon(data = ForcedRunTs_Ann, aes(x= DeltaAnn, ymin = isothermalMeanT_25, ymax = isothermalMeanT_75), alpha = 0.3, fill = 'light blue') +
  geom_line(data = ForcedRunTs_Ann, aes(x= DeltaAnn, y = isothermalMeanT_50), col = 'light blue') +
  geom_point(data = ForcedRunTs_Ann, aes(x= DeltaAnn, y = isothermalMeanT_50), col = 'light blue') +
  geom_point(aes(x = 0, y = modernPt$iceTs.isothermalMeanT_50)) + geom_abline(slope = 1, intercept = modernPt$iceTs.isothermalMeanT_50) +
  xlab('Perturbation (°C)') + ylab('Spring Mixing Temp (°C)') + ylim(0, 15)

MAFtemp <- ggplot() + # average temperature during months above freezing
  geom_ribbon(data = ForcedRunMAF_JJA, aes(x= DeltaJJA, ymin = X1_2.5, ymax = X1_97.5), alpha = 0.3, fill = 'orange') +
  geom_ribbon(data = ForcedRunMAF_JJA, aes(x= DeltaJJA, ymin = X1_25, ymax = X1_75), alpha = 0.3, fill = 'orange') +
  geom_line(data = ForcedRunMAF_JJA, aes(x= DeltaJJA, y = X1_50), col = 'orange') +
  geom_point(data = ForcedRunMAF_JJA, aes(x= DeltaJJA, y = X1_50), col = 'orange') +
  geom_ribbon(data = ForcedRunMAF_Ann, aes(x= DeltaAnn, ymin = X1_2.5, ymax = X1_97.5), alpha = 0.3, fill = 'light blue') +
  geom_ribbon(data = ForcedRunMAF_Ann, aes(x= DeltaAnn, ymin = X1_25, ymax = X1_75), alpha = 0.3, fill = 'light blue') +
  geom_line(data = ForcedRunMAF_Ann, aes(x= DeltaAnn, y = X1_50), col = 'light blue') +
  geom_point(data = ForcedRunMAF_Ann, aes(x= DeltaAnn, y = X1_50), col = 'light blue') +
  geom_point(aes(x = 0, y = modernPt$MAF.X1_50)) + geom_abline(slope = 1, intercept = modernPt$MAF.X1_50) +
  xlab('Perturbation (°C)') + ylab('MAF (°C)') + ylim(0, 15)

JJAtemp <- ggplot() + # average JJA temperatures
  geom_ribbon(data = ForcedRunJJA_JJA, aes(x= DeltaJJA, ymin = X1_2.5, ymax = X1_97.5), alpha = 0.3, fill = 'orange') +
  geom_ribbon(data = ForcedRunJJA_JJA, aes(x= DeltaJJA, ymin = X1_25, ymax = X1_75), alpha = 0.3, fill = 'orange') +
  geom_line(data = ForcedRunJJA_JJA, aes(x= DeltaJJA, y = X1_50), col = 'orange') +
  geom_point(data = ForcedRunJJA_JJA, aes(x= DeltaJJA, y = X1_50), col = 'orange') +
  geom_ribbon(data = ForcedRunJJA_Ann, aes(x= DeltaAnn, ymin = X1_2.5, ymax = X1_97.5), alpha = 0.3, fill = 'light blue') +
  geom_ribbon(data = ForcedRunJJA_Ann, aes(x= DeltaAnn, ymin = X1_25, ymax = X1_75), alpha = 0.3, fill = 'light blue') +
  geom_line(data = ForcedRunJJA_Ann, aes(x= DeltaAnn, y = X1_50), col = 'light blue') +
  geom_point(data = ForcedRunJJA_Ann, aes(x= DeltaAnn, y = X1_50), col = 'light blue') +
  geom_point(aes(x = 0, y = modernPt$JJA.X1_50)) + geom_abline(slope = 1, intercept = modernPt$JJA.X1_50) +
  xlab('Perturbation (°C)') + ylab('JJA (°C)') + ylim(0, 15)

IFD <- ggplot() + # length of ice free season
  geom_ribbon(data = ForcedRuns_JJA, aes(x= DeltaJJA, ymin = icefree_days_2.5, ymax = icefree_days_97.5), alpha = 0.3, fill = 'orange') +
  geom_ribbon(data = ForcedRuns_JJA, aes(x= DeltaJJA, ymin = icefree_days_25, ymax = icefree_days_75), alpha = 0.3, fill = 'orange') +
  geom_line(data = ForcedRuns_JJA, aes(x= DeltaJJA, y = icefree_days_50), col = 'orange') +
  geom_point(data = ForcedRuns_JJA, aes(x= DeltaJJA, y = icefree_days_50), col = 'orange') +
  geom_ribbon(data = ForcedRuns_Ann, aes(x= DeltaAnn, ymin = icefree_days_2.5, ymax = icefree_days_97.5), alpha = 0.3, fill = 'light blue') +
  geom_ribbon(data = ForcedRuns_Ann, aes(x= DeltaAnn, ymin = icefree_days_25, ymax = icefree_days_75), alpha = 0.3, fill = 'light blue') +
  geom_line(data = ForcedRuns_Ann, aes(x= DeltaAnn, y = icefree_days_50), col = 'light blue') +
  geom_point(data = ForcedRuns_Ann, aes(x= DeltaAnn, y = icefree_days_50), col = 'light blue') +
  geom_point(aes(x = 0, y = modernPt$doys.icefree_days_50)) +
   xlab('Perturbation (°C)') + ylab('Ice Free Season Length (days)')

(JJAtemp + MAFtemp + IceFreetemp)/ (IO2WKtemp + isothermaltemp + plot_spacer()) /
  (IFD + plot_spacer() + plot_spacer()) & theme_linedraw()  &
  theme(panel.border = element_rect(size = 1.5), plot.margin = margin(2,2,2,2, 'pt'),
        axis.ticks.length = unit(4, 'pt'), axis.ticks = element_line(size = 0.8))

# --- FIG S10: Summer stratification ---

if_strat <- ggplot() + # Number of stratified days during the ice free season
  geom_ribbon(data = ForcedRuns_JJA, aes(x= DeltaJJA, ymin = icefree_strat_total_2.5, ymax = icefree_strat_total_97.5), alpha = 0.3, fill = 'orange') +
  geom_ribbon(data = ForcedRuns_JJA, aes(x= DeltaJJA, ymin = icefree_strat_total_25, ymax = icefree_strat_total_75), alpha = 0.3, fill = 'orange') +
  geom_line(data = ForcedRuns_JJA, aes(x= DeltaJJA, y = icefree_strat_total_50), col = 'orange') +
  geom_point(data = ForcedRuns_JJA, aes(x= DeltaJJA, y = icefree_strat_total_50), col = 'orange') +
  geom_ribbon(data = ForcedRuns_Ann, aes(x= DeltaAnn, ymin = icefree_strat_total_2.5, ymax = icefree_strat_total_97.5), alpha = 0.3, fill = 'light blue') +
  geom_ribbon(data = ForcedRuns_Ann, aes(x= DeltaAnn, ymin = icefree_strat_total_25, ymax = icefree_strat_total_75), alpha = 0.3, fill = 'light blue') +
  geom_line(data = ForcedRuns_Ann, aes(x= DeltaAnn, y = icefree_strat_total_50), col = 'light blue') +
  geom_point(data = ForcedRuns_Ann, aes(x= DeltaAnn, y = icefree_strat_total_50), col = 'light blue') +
  geom_point(aes(x = 0, y = modernPt$doys.icefree_strat_total_50)) + #geom_abline(slope = 1, intercept = modernPt$iceTs.iceout2wkMeanT_50) +
  theme_cowplot() + xlab('Perturbation (°C)') + ylab('Stratified Ice Free Days')

if_strat + theme_linedraw()  &
  theme(panel.border = element_rect(size = 1.5), plot.margin = margin(2,2,2,2, 'pt'),
        axis.ticks.length = unit(4, 'pt'), axis.ticks = element_line(size = 0.8))

##  ------- FIG 5: Proxy foward modeling -----------

# load brGDGT MBT values
#load("GUS-MBT.Rdata")
load("GUS-MBTs.Rdata")
load("GUS-UK37.Rdata")

# load air temps (PSM input)
load("met_TempSummary.Rdata") # these are prepared and exported in LakeGus_model_simulations_MSversion.R

# load psm ensemble data
load('EnsSummaries-Modern.Rdata')

# load forced psm ensemble data
load('EnsSummaries-ForcedSurface.Rdata')
load('EnsSummaries-ForcedSurfaceJJA.Rdata')
load('EnsSummaries-ForcedSurfaceMAM.Rdata')
load('EnsSummaries-ForcedSurfaceMAF3.Rdata')
load('EnsSummaries-ForcedSurfaceRuns.Rdata')

load('EnsSummaries-ForcedSurfaceCOLD.Rdata')
load('EnsSummaries-ForcedSurfaceJJACOLD.Rdata')
load('EnsSummaries-ForcedSurfaceMAMCOLD.Rdata')
load('EnsSummaries-ForcedSurfaceMAF3COLD.Rdata')
load('EnsSummaries-ForcedSurfaceRunsCOLD.Rdata')

# Functions to forward model MBT and UK37 values from lake temperatures
Russell_T2MBT <- function(Temp){(Temp+1.21)/32.42}
Zhao_T2MBT <- function(Temp){(Temp+1.82)/56.06}

Longo2016_T2UK37 <- function(Temp){0.021*Temp - 0.68}
Longo2018_T2UK37 <- function(Temp){0.029*Temp - 0.49}
Yao2019_T2UK37 <- function(Temp){0.03*Temp - 0.469}
DAndrea2011_T2UK37 <- function(Temp){0.0245*Temp - 0.779} # this is equation from 2011 paper
DAndrea2012_T2UK37 <- function(Temp){0.02811755 *Temp -0.8155494} # this is equation from 2011 paper
DAndrea2016_T2UK37 <- function(Temp){0.0284*Temp - 0.655} # Vikvatnet

# --------------convert all PSM output into MBT and Uk37 space------------------
Ens_RussellMBT <- purrr::map(EnsSummaries, Russell_T2MBT)
Ens_ZhaoMBT <- purrr::map(EnsSummaries, Zhao_T2MBT)
Ens_Longo2016_UK37 <- purrr::map(EnsSummaries, Longo2016_T2UK37)
Ens_DAndrea2011_UK37 <- purrr::map(EnsSummaries, DAndrea2011_T2UK37)
Ens_DAndrea2012_UK37 <- purrr::map(EnsSummaries, DAndrea2012_T2UK37)
Ens_DAndrea2016_UK37 <- purrr::map(EnsSummaries, DAndrea2016_T2UK37)

ForcedEns_Longo2016_UK37 <- purrr::map_df(ForcedSurfaceEns, Longo2016_T2UK37)
ForcedEns_Longo2016_UK37[23:24] <- ForcedSurfaceEns[23:24]
ForcedEns_Longo2016_UK37_JJA <- purrr::map_df(ForcedSurfaceEnsJJA, Longo2016_T2UK37)
ForcedEns_Longo2016_UK37_JJA[23:24] <- ForcedSurfaceEnsJJA[23:24]
ForcedEns_Longo2016_UK37_MAM <- purrr::map_df(ForcedSurfaceEnsMAM, Longo2016_T2UK37)
ForcedEns_Longo2016_UK37_MAM[23:24] <- ForcedSurfaceEnsMAM[23:24]

ForcedEns_Russell_MBT <- purrr::map_df(ForcedSurfaceEns, Russell_T2MBT)
ForcedEns_Russell_MBT[23:24] <- ForcedSurfaceEns[23:24]
ForcedEns_Russell_MBT_JJA <- purrr::map_df(ForcedSurfaceEnsJJA, Russell_T2MBT)
ForcedEns_Russell_MBT_JJA[23:24] <- ForcedSurfaceEnsJJA[23:24]
ForcedEns_Russell_MBT_MAF <- purrr::map_df(ForcedSurfaceEnsMAF, Russell_T2MBT)
ForcedEns_Russell_MBT_MAF[23:24] <- ForcedSurfaceEnsMAF[23:24]

ForcedEnsCOLD_Longo2016_UK37 <- purrr::map_df(ForcedSurfaceEnsCOLD, Longo2016_T2UK37)
ForcedEnsCOLD_Longo2016_UK37[7:8] <- ForcedSurfaceEnsCOLD[7:8]
ForcedEnsCOLD_Longo2016_UK37_JJA <- purrr::map_df(ForcedSurfaceEnsJJACOLD, Longo2016_T2UK37)
ForcedEnsCOLD_Longo2016_UK37_JJA[7:8] <- ForcedSurfaceEnsJJACOLD[7:8]
ForcedEnsCOLD_Longo2016_UK37_MAM <- purrr::map_df(ForcedSurfaceEnsMAMCOLD, Longo2016_T2UK37)
ForcedEnsCOLD_Longo2016_UK37_MAM[7:8] <- ForcedSurfaceEnsMAMCOLD[7:8]

ForcedEnsCOLD_Russell_MBT <- purrr::map_df(ForcedSurfaceEnsCOLD, Russell_T2MBT)
ForcedEnsCOLD_Russell_MBT[7:8] <- ForcedSurfaceEnsCOLD[7:8]
ForcedEnsCOLD_Russell_MBT_JJA <- purrr::map_df(ForcedSurfaceEnsJJACOLD, Russell_T2MBT)
ForcedEnsCOLD_Russell_MBT_JJA[7:8] <- ForcedSurfaceEnsJJACOLD[7:8]
ForcedEnsCOLD_Russell_MBT_MAF <- purrr::map_df(ForcedSurfaceEnsMAFCOLD, Russell_T2MBT)
ForcedEnsCOLD_Russell_MBT_MAF[7:8] <- ForcedSurfaceEnsMAFCOLD[7:8]

# --------------functions to make plots of dists------------------
## Make plots of UK378 and MBT value distributions

#UK37

boxdata_UK37 <- function (df){
  df$iceTs %>%
    dplyr::select(contains("50")& !contains("Delta") & !contains("Min") & !contains("Max")) %>% stack() %>% rbind(
      GUS_UK37s %>% filter(!is.na(UK37ab)) %>%
        dplyr::select(UK37ab) %>% mutate(., ind='Obs') %>%
        set_colnames(c('values', 'ind'))) %>%
    rbind(df$MAM %>% dplyr:: select(dplyr::contains(c("X1")) & dplyr::contains("50")) %>% stack() %>%
            mutate_at(.vars = 'ind', .funs = funs(paste("MAM",.)))) %>%
    rbind(df$JJA %>% dplyr:: select(dplyr::contains(c("X1")) & dplyr::contains("50")) %>% stack() %>%
            mutate_at(.vars = 'ind', .funs = funs(paste("JJA",.))))}

UK37_obsOnly <- function (){
  GUS_UK37s %>% filter(!is.na(UK37ab)) %>%
    dplyr::select(UK37ab) %>% mutate(., ind='Obs') %>%
    set_colnames(c('values', 'ind'))
}

maxmindata_UK37 <- function (df){df$iceTs %>% dplyr::select(c('icefreeMaxT_50', 'icefreeMinT_50'))}

# ==================== Monster UK37 figure =====================================

CombinedLCA <-

  ggplot() +

  # Lines showing full range on ice-free values
  geom_linerange(data = maxmindata_UK37(Ens_Longo2016_UK37), aes(x = icefreeMinT_50, xmin = icefreeMinT_50, xmax = icefreeMaxT_50, y = 'MAM X1_50'),
                 color = '#5995ED', alpha = 1, lwd = 3, position = position_nudge(y = -0.6)) +
  geom_linerange(data = maxmindata_UK37(Ens_DAndrea2011_UK37), aes(x = icefreeMinT_50, xmin = icefreeMinT_50, xmax = icefreeMaxT_50, y = 'MAM X1_50'),
                 color = '#F5C08E', alpha = 1, lwd = 3, position = position_nudge(y = -0.8)) +
  geom_linerange(data = maxmindata_UK37(Ens_DAndrea2016_UK37), aes(x = icefreeMinT_50, xmin = icefreeMinT_50, xmax = icefreeMaxT_50, y = 'MAM X1_50'),
                 color = '#6BAA75', alpha = 1, lwd = 3, position = position_nudge(y = -0.4)) +

  # Core top sample
  geom_vline(xintercept = which(GUS_UK37s$Age_yr == min(GUS_UK37s$Age_yr, na.rm=TRUE)) %>%
               GUS_UK37s$UK37ab[.], lty = 'dashed') +

  # Density ridges for modern simulations
  stat_density_ridges(data = boxdata_UK37(Ens_DAndrea2016_UK37),
                      aes(x = values, y = ind, group = ind),  alpha = 1,
                      fill = '#6BAA75', scale = 0.8, bandwidth = 0.0119, quantile_lines = TRUE,
                      position = position_nudge(y = 0.15)) +
  stat_density_ridges(data = boxdata_UK37(Ens_DAndrea2011_UK37),
                      aes(x = values, y = ind, group = ind), alpha = 1,
                      fill = '#F5C08E', scale = 0.8, bandwidth = 0.0119, quantile_lines = TRUE,
                      position = position_nudge(y = 0.15)) +
  stat_density_ridges(data = boxdata_UK37(Ens_Longo2016_UK37),
                      aes(x = values, y = ind, group = ind), alpha = 1,
                      fill = '#5995ED', scale = 0.8, bandwidth = 0.0119, quantile_lines = TRUE,
                      position = position_nudge(y = 0.15)) +
  stat_density_ridges(data = UK37_obsOnly(),
                      aes(x = values, y = ind, group = ind), bandwidth = 0.0119, alpha = 1,
                      fill = '#313638', scale = 0.025, quantile_lines = TRUE,
                      position = position_nudge(y = 0.15)) +

  geom_point(aes(y = 'Obs', x = UK37_obsOnly() %>% .$values),
             position= position_jitter(width = 0, height = 0.2), #position_nudge(y = 0.4-0.28, x = 0),
             alpha = 1, size = 2.5, fill = 'black', shape = 24) +

  ## Interannual variation jittered points for modern simulations

  geom_point(aes(y = 'isothermalMeanT_50', x = Ens_Longo2016_UK37$iceTs$isothermalMeanT_50),
             position= position_jitter(height = 0.1, width = 0, seed = 116),
             alpha = 1, size = 1.3, color = '#5995ED') +
  geom_errorbar(aes(y = 'isothermalMeanT_50', x = Ens_Longo2016_UK37$iceTs$isothermalMeanT_50,
                    xmin = Ens_Longo2016_UK37$iceTs$isothermalMeanT_2.5,
                    xmax = Ens_Longo2016_UK37$iceTs$isothermalMeanT_97.5),
                position= position_jitter(height = 0.1, width = 0, seed = 116),
                linetype='solid', alpha = 1, size = 0.6,  color = '#5995ED', width = 0) +
  geom_point(aes(y = 'iceout2wkMeanT_50', x = Ens_Longo2016_UK37$iceTs$iceout2wkMeanT_50),
             position= position_jitter(height = 0.1, width = 0, seed = 94),
             alpha = 1, size = 1.3, color = '#5995ED') +
  geom_errorbar(aes(y = 'iceout2wkMeanT_50', x =Ens_Longo2016_UK37$iceTs$iceout2wkMeanT_50,
                    xmin = Ens_Longo2016_UK37$iceTs$iceout2wkMeanT_2.5,
                    xmax = Ens_Longo2016_UK37$iceTs$iceout2wkMeanT_97.5),
                position= position_jitter(height = 0.1, width = 0, seed = 94),
                linetype='solid', alpha = 1, size = 0.6,  color = '#5995ED', width = 0) +
  geom_point(aes(y = 'icefreeMeanT_50', x = Ens_Longo2016_UK37$iceTs$icefreeMeanT_50),
             position= position_jitter(height = 0.1, width = 0, seed = 244),
             alpha = 1, size = 1.3, color = '#5995ED') +
  geom_errorbar(aes(y = 'icefreeMeanT_50', x = Ens_Longo2016_UK37$iceTs$icefreeMeanT_50,
                    xmin = Ens_Longo2016_UK37$iceTs$icefreeMeanT_2.5,
                    xmax = Ens_Longo2016_UK37$iceTs$icefreeMeanT_97.5),
                position= position_jitter(height = 0.1, width = 0, seed = 244),
                linetype='solid', alpha = 1, size = 0.6,  color = '#5995ED', width = 0) +
  geom_point(aes(y = 'JJA X1_50', x = Ens_Longo2016_UK37$JJA$X1_50),
             position= position_jitter(height = 0.1, width = 0, seed = 43),
             alpha = 1, size = 1.3, color = '#5995ED') +
  geom_errorbar(aes(y = 'JJA X1_50', x = Ens_Longo2016_UK37$JJA$X1_50,
                    xmin = Ens_Longo2016_UK37$JJA$X1_2.5,
                    xmax = Ens_Longo2016_UK37$JJA$X1_97.5),
                position= position_jitter(height = 0.1, width = 0, seed = 43),
                linetype='solid', alpha = 1, size = 0.6,  color = '#5995ED', width = 0) +
  geom_point(aes(y = 'MAM X1_50', x = Ens_Longo2016_UK37$MAM$X1_50),
             position= position_jitter(height = 0.1, width = 0, seed = 413),
             alpha = 1, size = 1.3, color = '#5995ED') +
  geom_errorbar(aes(y = 'MAM X1_50', x = Ens_Longo2016_UK37$MAM$X1_50,
                    xmin = Ens_Longo2016_UK37$MAM$X1_2.5,
                    xmax = Ens_Longo2016_UK37$MAM$X1_97.5),
                position= position_jitter(height = 0.1, width = 0, seed = 413),
                linetype='solid', alpha = 1, size = 0.6,  color = '#5995ED', width = 0) +

  ## Triangles showing forced simulation changes to the mean (points and error bars)
  geom_point(data = ForcedEns_Longo2016_UK37[which(ForcedEns_Longo2016_UK37$DeltaAnn == 0 &
                                                     (ForcedEns_Longo2016_UK37$DeltaJJA %in% c(1,2.5,5,10))),],
             aes(y= 'iceout2wkMeanT_50', x = iceout2wkMeanT_mean, size = DeltaJJA),
             fill = '#5995ED', shape = 24, position = position_nudge(x=0, y = 0.4-0.28)) +
  geom_point(data = ForcedEnsCOLD_Longo2016_UK37[which(ForcedEnsCOLD_Longo2016_UK37$DeltaAnn == 0 &
                                                         (ForcedEnsCOLD_Longo2016_UK37$DeltaJJA %in% c(-1,-2.5,-5,-10))),],
             aes(y= 'iceout2wkMeanT_50', x = iceout2wkMeanT_mean, size = abs(DeltaJJA)),
             fill = '#5995ED', shape = 24, position = position_nudge(x=0, y = 0.4-0.28)) +

  geom_point(data = ForcedEns_Longo2016_UK37[which(ForcedEns_Longo2016_UK37$DeltaAnn == 0 &
                                                     (ForcedEns_Longo2016_UK37$DeltaJJA %in% c(1,2.5,5,10))),],
             aes(y= 'isothermalMeanT_50', x = isothermalMeanT_mean,  size = DeltaJJA),
             fill = '#5995ED', shape = 24,  position = position_nudge(x = 0, y = 0.4-0.28))  +
  geom_point(data = ForcedEnsCOLD_Longo2016_UK37[which(ForcedSurfaceEnsCOLD$DeltaAnn == 0 &
                                                         (ForcedSurfaceEnsCOLD$DeltaJJA %in% c(-1,-2.5,-5,-10))),],
             aes(y= 'isothermalMeanT_50', x = isothermalMeanT_mean,  size = abs(DeltaJJA)),
             fill = '#5995ED', shape = 24,  position = position_nudge(x = 0, y = 0.4-0.28))  +

  geom_point(data = ForcedEns_Longo2016_UK37[which(ForcedEns_Longo2016_UK37$DeltaAnn == 0 &
                                                     (ForcedEns_Longo2016_UK37$DeltaJJA %in% c(1,2.5,5,10))),],
             aes(y= 'icefreeMeanT_50', x = icefreeMeanT_mean, size = DeltaJJA),
             fill = '#5995ED', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +
  geom_point(data = ForcedEnsCOLD_Longo2016_UK37[which(ForcedEnsCOLD_Longo2016_UK37$DeltaAnn == 0 &
                                                         (ForcedEnsCOLD_Longo2016_UK37$DeltaJJA %in% c(-1,-2.5,-5,-10))),],
             aes(y= 'icefreeMeanT_50', x = icefreeMeanT_mean, size = abs(DeltaJJA)),
             fill = '#5995ED', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +

  geom_point(data = ForcedEns_Longo2016_UK37_JJA[which(ForcedEns_Longo2016_UK37_JJA$DeltaAnn == 0 &
                                                         (ForcedEns_Longo2016_UK37$DeltaJJA %in% c(1,2.5,5,10))),],
             aes(y= 'JJA X1_50', x = X1_mean, size = DeltaJJA),
             fill = '#5995ED', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +

  geom_point(data = ForcedEnsCOLD_Longo2016_UK37_JJA[which(ForcedEnsCOLD_Longo2016_UK37$DeltaAnn == 0 &
                                                             (ForcedEnsCOLD_Longo2016_UK37$DeltaJJA %in% c(-1,-2.5,-5,-10))),],
             aes(y= 'JJA X1_50', x = X1_mean, size = abs(DeltaJJA)),
             fill = '#5995ED', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +

  geom_point(data = ForcedEns_Longo2016_UK37_MAM[which(ForcedEns_Longo2016_UK37_MAM$DeltaAnn == 0 &
                                                         (ForcedEns_Longo2016_UK37$DeltaJJA %in% c(1,2.5,5,10))),],
             aes(y= 'MAM X1_50', x = X1_mean, size = DeltaJJA),
             fill = '#5995ED', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +

  geom_point(data = ForcedEnsCOLD_Longo2016_UK37_MAM[which(ForcedEnsCOLD_Longo2016_UK37$DeltaAnn == 0 &
                                                             (ForcedEnsCOLD_Longo2016_UK37$DeltaJJA %in% c(-1,-2.5,-5,-10))),],
             aes(y= 'MAM X1_50', x = X1_mean), size = abs(c(-1,-2.5,-5,-10)),
             fill = '#5995ED', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +

  geom_point(data = ForcedEns_Longo2016_UK37[which(ForcedEns_Longo2016_UK37$DeltaAnn %in% c(1,2.5,5,10)),],
             aes(y= 'iceout2wkMeanT_50', x = iceout2wkMeanT_mean, size = DeltaAnn),
             fill = '#5995ED', shape = 2, position = position_nudge(x = 0, y = 0)) +

  geom_point(data = ForcedEnsCOLD_Longo2016_UK37[which(ForcedEnsCOLD_Longo2016_UK37$DeltaAnn %in% c(-1,-2.5,-5,-10)),],
             aes(y= 'iceout2wkMeanT_50', x = iceout2wkMeanT_mean, size = abs(c(-1,-2.5,-5,-10))),
             fill = '#5995ED', shape = 2, position = position_nudge(x = 0, y = 0)) +

  geom_point(data = ForcedEns_Longo2016_UK37[which(ForcedEns_Longo2016_UK37$DeltaAnn %in% c(1,2.5,5,10)),],
             aes(y= 'isothermalMeanT_50', x = isothermalMeanT_mean, size = DeltaAnn),
             fill = '#5995ED', shape = 2,  position = position_nudge(x = 0, y = 0))  +

  geom_point(data = ForcedEnsCOLD_Longo2016_UK37[which(ForcedEnsCOLD_Longo2016_UK37$DeltaAnn %in% c(-1,-2.5,-5,-10)),],
             aes(y= 'isothermalMeanT_50', x = isothermalMeanT_mean, size = abs(DeltaAnn)),
             fill = '#5995ED', shape = 2,  position = position_nudge(x = 0, y = 0))  +

  geom_point(data = ForcedEns_Longo2016_UK37[which(ForcedEns_Longo2016_UK37$DeltaAnn %in% c(1,2.5,5,10)),],
             aes(y= 'icefreeMeanT_50', x = icefreeMeanT_mean, size = DeltaAnn),
             fill = '#5995ED', shape = 2,  position = position_nudge(x = 0, y = 0)) +

  geom_point(data = ForcedEnsCOLD_Longo2016_UK37[which(ForcedEnsCOLD_Longo2016_UK37$DeltaAnn %in% c(-1,-2.5,-5,-10)),],
             aes(y= 'icefreeMeanT_50', x = icefreeMeanT_mean, size = abs(DeltaAnn)),
             fill = '#5995ED', shape = 2,  position = position_nudge(x = 0, y = 0)) +

  geom_point(data = ForcedEns_Longo2016_UK37_JJA[which(ForcedEns_Longo2016_UK37_JJA$DeltaAnn %in% c(1,2.5,5,10)),],
             aes(y= 'JJA X1_50', x = X1_mean, size = DeltaAnn),
             fill = '#5995ED', shape = 2,  position = position_nudge(x = 0, y = 0)) +

  geom_point(data = ForcedEnsCOLD_Longo2016_UK37_JJA[which(ForcedEnsCOLD_Longo2016_UK37$DeltaAnn %in% c(-1,-2.5,-5,-10)),],
             aes(y= 'JJA X1_50', x = X1_mean, size = abs(c(-1,-2.5,-5,-10))),
             fill = '#5995ED', shape = 2,  position = position_nudge(x = 0, y = 0)) +

  geom_point(data = ForcedEns_Longo2016_UK37_MAM[which(ForcedEns_Longo2016_UK37_MAM$DeltaAnn %in% c(1,2.5,5,10)),],
             aes(y= 'MAM X1_50', x = X1_mean, size = DeltaAnn),
             fill = '#5995ED', shape = 2, position = position_nudge(x = 0, y = 0)) +

  # Formatting
  xlab('UK37') + ylab('Period Integrated') + scale_y_discrete(expand=expansion(add=c(1.2)),
                                                              limits=c('MAM X1_50','JJA X1_50', 'icefreeMeanT_50',
                                                                       'isothermalMeanT_50',  'iceout2wkMeanT_50', 'Obs'),
                                                              labels=c('MAM', 'JJA', 'Ice Free\nSeason ', 'Isothermal\nMixing',
                                                                       '2 Weeks\nFollowing\nIce Out', 'Downcore\nObservations')) +
  scale_size_binned(range = c(1.5, 6), breaks = c(0.5, 2, 3, 7, 12, 15))


CombinedLCA + theme_cowplot()#+theme(aspect.ratio = 0.8)

# MBT
boxdata_MBT_X1 <- function (df){
  df$MAF %>% dplyr::select(contains("50")  & contains("X1"))  %>% stack() %>% mutate('ind'= "MAF") %>%
    rbind(df$JJA %>% dplyr:: select(dplyr::contains(c("X1")) & dplyr::contains("50")) %>% stack() %>%
            mutate('ind' = "JJA")) %>%
    rbind(df$iceTs %>% dplyr:: select(icefreeMeanT_50) %>% stack() %>%
            mutate('ind' = "icefree"))}

boxdata_MBT_X6 <- function (df){
  df$MAF %>% dplyr::select(contains("50")  & contains("X6"))  %>% stack() %>% mutate('ind'= "MAF") %>%
    rbind(df$JJA %>% dplyr:: select(dplyr::contains(c("X6")) & dplyr::contains("50")) %>% stack() %>%
            mutate('ind' = "JJA")) %>%
    rbind(df$iceTs %>% dplyr:: select(icefreeMeanT_50) %>% stack() %>%
            mutate('ind' = "icefree"))}

boxdata_MBT_WCT <- function (df){
  df$MAF %>% dplyr::select(contains("50")  & contains("WCT"))  %>% stack() %>% mutate('ind'= "MAF") %>%
    rbind(df$JJA %>% dplyr:: select(dplyr::contains(c("WCT")) & dplyr::contains("50")) %>% stack() %>%
            mutate('ind' = "JJA")) %>%
    rbind(df$iceTs %>% dplyr:: select(icefreeMeanT_50) %>% stack() %>%
            mutate('ind' = "icefree"))}

MBT_obsOnly <- function (){
  GUS_MBTs %>% filter(!is.na(MBT.5me)) %>%
    dplyr::select(MBT.5me) %>% mutate(., ind='Obs') %>%
    set_colnames(c('values', 'ind'))}

MBT_GOODobsOnly <- function (){
  GUS_MBTs %>% filter(!is.na(MBT.5me)) %>%
    dplyr::select(MBT.5me) %>% mutate(., ind='Obs') %>%
    set_colnames(c('values', 'ind'))}

maxmindata_MBT <- function (df){df$iceTs %>% dplyr::select(c('icefreeMaxT_50', 'icefreeMinT_50'))}

GDGT0.4 <- GUS_biomarker$GDGT.0/GUS_biomarker$GDGT.4

# ==================== Monster MBT figure =====================================

CombinedMBT <-

  ggplot() +

  #Lines showing full range on ice-free values
  geom_linerange(data = maxmindata_MBT(Ens_RussellMBT), aes(x = icefreeMinT_50, xmin = icefreeMinT_50, xmax = icefreeMaxT_50, y = 'MAF'),
               color = '#F88A4F', alpha = 1, lwd = 3, position = position_nudge(y = -0.4)) +
  geom_linerange(data = maxmindata_MBT(Ens_ZhaoMBT), aes(x = icefreeMinT_50, xmin = icefreeMinT_50, xmax = icefreeMaxT_50, y = 'MAF'),
                 color = '#1FB8FF', alpha = 1, lwd = 3, position = position_nudge(y = -0.6)) +

  geom_vline(xintercept = GUS_MBTs %>% slice(which.min(Age_yr)) %>% select(MBT.5me) %>% unlist, lty = 'dashed') +

  stat_density_ridges(data = boxdata_MBT_X1(Ens_RussellMBT),
                      aes(x = values, y = ind, group = ind), alpha = 1,
                      fill = '#F88A4F', scale = 0.8,  quantile_lines = TRUE,
                      position = position_nudge(x = 0, y = 0.15)) +
  stat_density_ridges(data = boxdata_MBT_X6(Ens_RussellMBT),
                      aes(x = values, y = ind, group = ind), alpha = 1,
                      fill = NA,  scale = 0.8,  quantile_lines = FALSE,
                      position = position_nudge(x = 0, y = 0.15), linetype = 'dashed') +
  stat_density_ridges(data = boxdata_MBT_WCT(Ens_RussellMBT),
                      aes(x = values, y = ind, group = ind), alpha = 1,
                      fill = NA,  scale = 0.8,  quantile_lines = FALSE,
                      position = position_nudge(x = 0, y = 0.15), linetype = 'dotted') +
  stat_density_ridges(data = boxdata_MBT_X1(Ens_ZhaoMBT),
                      aes(x = values, y = ind, group = ind), alpha = 1,
                      fill = '#1FB8FF', scale = 0.8,  quantile_lines = TRUE,
                      position = position_nudge(x=0, y = 0.15)) +

  stat_density_ridges(aes(x = MBT.5me[which(brGDGT.TPH$hexa/rowSums(brGDGT.TPH) > 0.3)], y = 'Obs'),  alpha = 1,
                      fill = 'dark gray', scale = 0.05, quantile_lines = TRUE,
                      position = position_nudge(x = 0, y = 0.15)) +
  stat_density_ridges(aes(x = MBT.5me, y = 'Obs'),  alpha = 1,
                      fill = 'blue', scale = 0.05, quantile_lines = TRUE,
                      position = position_nudge(x = 0, y = 0.15)) +

# Downcore observations; color scaled to GDGT 0/4 and Xs based on % hexa
geom_point(aes(y = 'Obs', x = MBT.5me , fill = GDGT0.4),
           position= position_jitter(height = 0.2, width = 0, seed = 2022),
           alpha = 1, size = 2.5, shape = 25) +  geom_point(aes(y = 'Obs', x = MBT.5me,
                                                                alpha = as.numeric(ifelse(MBT.5me %in% MBT.5me[which(brGDGT.TPH$hexa/rowSums(brGDGT.TPH) < 0.3)],1,0))), pch = 4, size = 3,
                                                            position= position_jitter(height = 0.2, width = 0, seed = 2022)) +
  scale_colour_gradientn(trans = "log", colors = c('black', '#9E0059'),
                         breaks = c(0,1,2.5,5,10,25,50,100,250)) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1) +


## Interannual variation jittered points for modern simulations
geom_point(aes(y = 'JJA', x = Ens_RussellMBT$JJA$X1_50),
           position= position_jitter(height = 0.1, width = 0, seed = 116),
           alpha = 1, size = 1.3, color = '#F88A4F') +
  geom_errorbar(aes(y = 'JJA', x = Ens_RussellMBT$JJA$X1_50,
                    xmin = Ens_RussellMBT$JJA$X1_2.5,
                    xmax = Ens_RussellMBT$JJA$X1_97.5),
                position= position_jitter(height = 0.1, width = 0, seed = 116),
                linetype='solid', alpha = 1, size = 0.6,  color = '#F88A4F', width = 0) +
  geom_point(aes(y = 'MAF', x = Ens_RussellMBT$MAF$X1_50),
             position= position_jitter(height = 0.1, width = 0, seed = 116),
             alpha = 1, size = 1.3, color = '#F88A4F') +
  geom_errorbar(aes(y = 'MAF', x = Ens_RussellMBT$MAF$X1_50,
                    xmin = Ens_RussellMBT$MAF$X1_2.5,
                    xmax = Ens_RussellMBT$MAF$X1_97.5),
                position= position_jitter(height = 0.1, width = 0, seed = 116),
                linetype='solid', alpha = 1, size = 0.6,  color = '#F88A4F', width = 0) +

  geom_point(aes(y = 'icefree', x = Ens_RussellMBT$iceTs$icefreeMeanT_50),
             position= position_jitter(height = 0.1, width = 0, seed = 244),
             alpha = 1, size = 1.3, color = '#F5C08E') +
  geom_errorbar(aes(y = 'icefree', x = Ens_RussellMBT$iceTs$icefreeMeanT_50,
                    xmin = Ens_RussellMBT$iceTs$icefreeMeanT_2.5,
                    xmax = Ens_RussellMBT$iceTs$icefreeMeanT_97.5),
                position= position_jitter(height = 0.1, width = 0, seed = 244),
                linetype='solid', alpha = 1, size = 0.6,  color = '#F5C08E', width = 0) +

  ## Forced model simulations
  geom_point(data = ForcedEns_Russell_MBT_JJA[which(ForcedEns_Russell_MBT_JJA$DeltaAnn == 0 &
                                                      (ForcedEns_Russell_MBT_JJA$DeltaJJA %in% c(1,2.5,5,10))),],
             aes(y= 'JJA', x = X1_mean, size = DeltaJJA),
             fill = '#F88A4F', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +
  geom_point(data = ForcedEnsCOLD_Russell_MBT_JJA[which(ForcedEns_Russell_MBT_JJA$DeltaAnn == 0 &
                                                          (ForcedEns_Russell_MBT_JJA$DeltaJJA %in% c(1,2.5,5,10))),],
             aes(y= 'JJA', x = X1_mean, size = c(1,2.5,5,10)),
             fill = '#F88A4F', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +
  geom_point(data = ForcedEns_Russell_MBT_JJA[which(ForcedEns_Russell_MBT_JJA$DeltaAnn %in% c(1,2.5,5,10)),],
             aes(y= 'JJA', x = X1_mean, size = DeltaAnn),
             fill = '#F88A4F', shape = 2, position = position_nudge(x = 0, y = 0.4-0.28)) +
  geom_point(data = ForcedEnsCOLD_Russell_MBT_JJA[which(ForcedEns_Russell_MBT_JJA$DeltaAnn %in% c(1,2.5,5,10)),],
             aes(y= 'JJA', x = X1_mean, size = c(1,2.5,5,10)),
             fill = '#F88A4F', shape = 2, position = position_nudge(x = 0, y = 0.4-0.28)) +
  geom_point(data = ForcedEns_Russell_MBT_MAF[which(ForcedEns_Russell_MBT_MAF$DeltaAnn == 0 &
                                                      (ForcedEns_Russell_MBT_MAF$DeltaJJA %in% c(1,2.5,5,10))),],
             aes(y= 'MAF', x = X1_mean, size = DeltaJJA),
             fill = '#F88A4F', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +

  geom_point(data = ForcedEnsCOLD_Russell_MBT_MAF[which(ForcedEns_Russell_MBT_MAF$DeltaAnn == 0 &
                                                          (ForcedEns_Russell_MBT_MAF$DeltaJJA %in% c(1,2.5,5,10))),],
             aes(y= 'MAF', x = X1_50, size = c(1,2.5,5,10)),
             fill = '#F88A4F', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +

geom_point(data = ForcedEns_Russell_MBT_MAF[which(ForcedEns_Russell_MBT_MAF$DeltaAnn %in% c(1,2.5,5,10)),],
             aes(y= 'MAF', x = X1_mean, size = DeltaAnn),
             fill = '#F88A4F', shape = 2, position = position_nudge(x = 0, y = 0.4-0.28)) +
  geom_point(data = ForcedEnsCOLD_Russell_MBT_MAF[which(ForcedEns_Russell_MBT_MAF$DeltaAnn %in% c(1,2.5,5,10)),],
             aes(y= 'MAF', x = X1_50, size = c(1,2.5,5,10)),
             fill = '#F88A4F', shape = 2, position = position_nudge(x = 0, y = 0.4-0.28)) +

  geom_point(data = ForcedEns_Russell_MBT[which(ForcedEns_Russell_MBT$DeltaAnn == 0 &
                                                  (ForcedEns_Russell_MBT$DeltaJJA %in% c(1,2.5,5,10))),],
             aes(y= 'icefree', x = icefreeMeanT_mean, size = DeltaJJA),
             fill = '#F88A4F', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +

  geom_point(data = ForcedEns_Russell_MBT[which(ForcedEns_Russell_MBT$DeltaAnn %in% c(1,2.5,5,10)),],
             aes(y= 'icefree', x = icefreeMeanT_mean, size = DeltaAnn),
             fill = '#F88A4F', shape = 2, position = position_nudge(x = 0, y = 0.4-0.28)) +

  geom_point(data = ForcedEnsCOLD_Russell_MBT[which(ForcedEns_Russell_MBT$DeltaAnn == 0 &
                                                      (ForcedEns_Russell_MBT$DeltaJJA %in% c(1,2.5,5,10))),],
             aes(y= 'icefree', x = icefreeMeanT_mean, size = c(1,2.5,5,10)),
             fill = '#F88A4F', shape = 24, position = position_nudge(x = 0, y = 0.4-0.28)) +

  geom_point(data = ForcedEnsCOLD_Russell_MBT[which(ForcedEns_Russell_MBT$DeltaAnn %in% c(1,2.5,5,10)),],
             aes(y= 'icefree', x = icefreeMeanT_mean, size = c(1,2.5,5,10)),
             fill = '#F88A4F', shape = 2, position = position_nudge(x = 0, y = 0.4-0.28)) +

  xlab('MBT5me') + ylab('Period Integrated') +
  scale_y_discrete(expand=expansion(add=c(1)), limits=c('MAF','JJA', 'icefree', 'Obs'),
                   labels=c('MAF', 'JJA', 'Ice Free Season', 'Downcore\nObservations')) +
  scale_size_binned(range = c(1.5, 6), breaks = c(0.5, 2, 3, 7, 12, 15))

CombinedMBT + theme_cowplot()
