
## Calculate calibration prior distributions

library(distr)
library(dplyr)
library(purrr)

## Make function to extract relevant calibration info

calibVals <- function(calibFile){
  Data <- read.csv(calibFile)
  lmData <- lm(Data[,1] ~ Data[,2])
  intMean <- lmData$coefficients[1]
  slopeMean <- lmData$coefficients[2]
  intSE <- coef(summary(lmData))[1, "Std. Error"]
  slopeSE <- coef(summary(lmData))[2, "Std. Error"]
  lmData2 <- lm(Data[,2] ~ Data[,1])
  RMSEtemp <- sqrt(mean(lmData2$residuals^2))
  res <- list('slopeMean' = slopeMean, 'slopeSE' = slopeSE, 'intMean' = intMean, 'intSE' = intSE, 'RMSEtemp' = RMSEtemp)
  res}

## Read in calibration data

calibFiles = list(Russell18= 'CalibrationData/Russell2018-MAAT-MBT5me.csv', Zhao21='CalibrationData/Zhao2021-IWTjja-MBT5me.csv',
                  Longo16 = 'CalibrationData/Longo2016-insituTemp-UK37.csv', Longo18 = 'CalibrationData/Longo2018-MTSI-UK37.csv',
                  DAndrea11 = 'CalibrationData/DAndrea2011-Braya-UK37.csv', DAndrea16 = 'CalibrationData/DAndrea2016-Vikvat-UK37.csv',
                  DAndrea12 ='CalibrationData/DAndrea2012-Kongress-UK37.csv')

calibData <- calibFiles %>% map(calibVals)

## Define all priors and create r/d sampling functions for each

MBTslopePrior <- distr::UnivarMixingDistribution(Norm(mean=calibData$Russell18$slopeMean, sd=calibData$Russell18$slopeSE),
                                                 Norm(mean=calibData$Zhao21$slopeMean, sd=calibData$Zhao21$slopeSE))

rMBTslopePrior <- r(MBTslopePrior)
dMBTslopePrior <- d(MBTslopePrior)

MBTintPrior <- UnivarMixingDistribution(Norm(mean=calibData$Russell18$intMean, sd=calibData$Russell18$intSE),
                                        Norm(mean=calibData$Zhao21$intMean, sd=calibData$Zhao21$intSE))

rMBTintPrior <- r(MBTintPrior)
dMBTintPrior <- d(MBTintPrior)

UK37slopePrior <- UnivarMixingDistribution(Norm(mean=calibData$Longo16$slopeMean, sd=calibData$Longo16$slopeSE),
                                           Norm(mean=calibData$DAndrea11$slopeMean, sd=calibData$DAndrea11$slopeSE),
                                           Norm(mean=calibData$DAndrea12$slopeMean, sd=calibData$DAndrea12$slopeSE),
                                           Norm(mean=calibData$DAndrea16$slopeMean, sd=calibData$DAndrea16$slopeSE))

rUK37slopePrior <- r(UK37slopePrior)
dUK37slopePrior <- d(UK37slopePrior)

UK37intPrior <- UnivarMixingDistribution(Norm(mean=calibData$Longo16$intMean, sd=calibData$Longo16$intSE),
                                         Norm(mean=calibData$DAndrea11$intMean, sd=calibData$DAndrea11$intSE),
                                         Norm(mean=calibData$DAndrea12$intMean, sd=calibData$DAndrea12$intSE),
                                         Norm(mean=calibData$DAndrea16$intMean, sd=calibData$DAndrea16$intSE))

rUK37intPrior <- r(UK37intPrior)
dUK37intPrior <- d(UK37intPrior)


## ---------------- Plot priors -----------------------------------------------

library(ggplot2)
library(patchwork)

UK37slopeP <- ggplot() + theme_minimal() +
  geom_area(aes(c(0.015, 0.028)), stat = "function", fun = dnorm, na.rm = T, n = 1e4, args = list(mean = calibData$Longo16$slopeMean, sd = calibData$Longo16$slopeSE), fill = "#339989", col = NA, alpha = 0.5,  xlim = c(0, 0.03)) +
  geom_area(aes(c(0, 0.05)), stat = "function", fun = dnorm, na.rm = T, n = 1e4, args = list(mean = calibData$DAndrea11$slopeMean, sd = calibData$DAndrea11$slopeSE), fill = '#0C6291', alpha = 0.5,  col = NA, xlim = c(0, 0.03)) +
  geom_area(aes(c(0, 0.04)), stat = "function", fun = dnorm, na.rm = T, n = 1e4, args = list(mean = calibData$DAndrea12$slopeMean, sd = calibData$DAndrea12$slopeSE), fill = "#ABA9C3", col = NA, alpha = 0.5,  xlim = c(0, 0.04)) +
  geom_area(aes(c(0, 0.05)), stat = "function", fun = dnorm, na.rm = T, n = 1e4, args = list(mean = calibData$DAndrea16$slopeMean, sd = calibData$DAndrea16$slopeSE), fill = '#28557A', alpha = 0.5,  col = NA, xlim = c(0, 0.04)) +
  #geom_density(aes(rUK37slopePrior(1e6)), lwd = 1.2)  +
  xlab('UK37-T slope') + ylab('Probability Density') + xlim(0.01, 0.04) +
  theme(panel.grid.minor.y = element_blank()) + theme(text = element_text(size = 12))

UK37slopeP

UK37intP <- ggplot() + theme_minimal() +
  geom_area(aes(c(-1, -0.5)), stat = "function", fun = dnorm, na.rm = T, n = 1e4, args = list(mean = calibData$Longo16$intMean, sd = calibData$Longo16$intSE), fill = "#339989", col = NA, alpha = 0.5,  xlim = c(-1, -0.5))+
  geom_area(aes(c(-1, -0.5)), stat = "function", fun = dnorm, na.rm = T, n = 1e4, args = list(mean = calibData$DAndrea11$intMean, sd = calibData$DAndrea11$intSE), fill = '#0C6291', alpha = 0.5,  col = NA, xlim = c(-1, -0.5)) +
  geom_area(aes(c(-1, -0.5)), stat = "function", fun = dnorm, na.rm = T, n = 1e4, args = list(mean = calibData$DAndrea12$intMean, sd = calibData$DAndrea12$intSE), fill = "#ABA9C3", col = NA, alpha = 0.5,  xlim = c(-1, -0.5)) +
  geom_area(aes(c(-1, -0.5)), stat = "function", fun = dnorm, na.rm = T, n = 1e4, args = list(mean = calibData$DAndrea16$intMean, sd = calibData$DAndrea16$intSE), fill = '#28557A', alpha = 0.5,  col = NA, xlim = c(-1, -0.5)) +
  #geom_density(aes(rUK37intPrior(1e6)), lwd = 1.2)  +
  xlab('UK37-T Intercept') + ylab('Probability Density') + xlim(-0.95, -0.5) +
  theme(panel.grid.minor.y = element_blank()) + theme(text = element_text(size = 12))

UK37intP

MBTslopeP <- ggplot() + theme_minimal() +
  geom_area(aes(c(0, 0.25)), stat = "function", fun = dnorm, na.rm = T, n = 1e4, args = list(mean = calibData$Russell18$slopeMean, sd = calibData$Russell18$slopeSE), fill = "red", col = NA, alpha = 0.5,  xlim = c(0, 0.25))+
  geom_area(aes(c(0, 0.25)), stat = "function", fun = dnorm, na.rm = T, n = 5e4, args = list(mean = calibData$Zhao21$slopeMean, sd = calibData$Zhao21$slopeSE), fill = 'gold', alpha = 0.5,  col = NA, xlim = c(0, 0.25)) +
  #geom_density(aes(rMBTslopePrior(1e6)), lwd = 1.2)  +
  xlab('MBT5Me-T Slope') + ylab('Probability Density') +  xlim(0, 0.037) +
  theme(panel.grid.minor.y = element_blank()) + theme(text = element_text(size = 12))

MBTslopeP

MBTintP <- ggplot() + theme_minimal() +
  geom_area(aes(c(-0.1, 0.25)), stat = "function", fun = dnorm, na.rm = T, n = 1e4, args = list(mean = calibData$Russell18$intMean, sd = calibData$Russell18$intSE), fill = "red", col = NA, alpha = 0.5,  xlim = c(-0.1, 0.25))+
  geom_area(aes(c(-0.1, 0.25)), stat = "function", fun = dnorm, na.rm = T, n = 1e4, args = list(mean = calibData$Zhao21$intMean, sd = calibData$Zhao21$intSE), fill = 'gold', alpha = 0.5,  col = NA, xlim = c(-0.1, 0.25)) +
  #geom_density(aes(rMBTintPrior(1e6)), lwd = 1.2)  +
  xlab('MBT5Me-T Intercept') + ylab('Probability Density') + xlim(-0.05, 0.2) +
  theme(panel.grid.minor.y = element_blank()) + theme(text = element_text(size = 12))

MBTintP

 (MBTslopeP + MBTintP + UK37slopeP + UK37intP ) + plot_layout(nrow = 1)

UK37seasP <- ggplot() + stat_function(fun = function(x) {dnorm(x, mean = 7, sd = 3)}, lwd = 1.2) +
  xlim(-20,20)   + theme_minimal() + ylab('Probability Density') +
  xlab('UK37 start (days relative to ice out)') +  theme(panel.grid.minor.y = element_blank()) + theme(text = element_text(size = 12))

UK37durP <- ggplot() + stat_function(fun = function(x) {dnorm(x, mean = 10, sd = 3)}, lwd = 1.2) +
  theme_minimal() + ylab('Probability Density') +
  xlab('UK37 duration (days)') + xlim(-10, 30)  +  theme(panel.grid.minor.y = element_blank()) + theme(text = element_text(size = 12))

MBTseasP <- ggplot() + stat_function(fun = function(x) {dnorm(x, mean = 0, sd = 2)}, lwd = 1.2) +
  xlim(-20,20) + theme_minimal() + ylab('Probability Density') +
  xlab('MBT temp threshold (°C)') +  theme(panel.grid.minor.y = element_blank()) + theme(text = element_text(size = 12))

(MBTslopeP + MBTintP + UK37slopeP + UK37intP )  + UK37seasP + UK37durP + MBTseasP + plot_layout(nrow = 2)
