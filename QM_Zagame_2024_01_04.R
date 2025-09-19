#####################################################################
# Project: "R project for Quantitative Methods"                     
# Description:   This script contains the second part of the R project for the QM
#              course, which focuses on the analysis of time series data for 
#             Costa Rica, from 1980 to 2023. 
#                Part a) of the script is devoted
#             to the representation of the selected variables (GDP at constant,
#             prices, Unemployment, and Inflation, measured both with the CPI index
#             and the GDP deflator). 
#               Part b) revolves around the presentation of 
#             correlograms of both the levels, logged, and the growth rates of the indicated
#             variables.
#               In part c) we run an Augmented Dickey-Fuller both on the logged levels and
#             on the first differences of the log to check whether the time series at hand are
#             stationary or not.
#               Part d) has the objective of estimating a Phillips curve in which
#             first differences of logged inflation are explained by their past history and
#             by the past levels of the unemployment rate. We'll also compute the BIC to estimate
#             the number of lags p.
# File_name: QM_Zagame_2024_01_02.R
# Filed by: Davide Zagame
# First filed on: 2024 01 04
#                 Updated: 2025 09 18
#
#####################################################################


# Clear memory
rm(list = ls())

# Download the needed libraries
install.packages("AER")
install.packages("systemfit")
install.packages("moments")
install.packages("stargazer")
install.packages("MASS")
install.packages("lmtest")
install.packages("sandwich")
install.packages("plm")
install.packages("gridExtra")
install.packages("fBasics")
install.packages("ivpack")
install.packages("ggplot2")
install.packages("haven")
install.packages("dplyr")
install.packages("pryr")
install.packages("normtest")
install.packages("het.test")
install.packages("tis")
install.packages("readxl")


# Load the needed libraries
library("AER")
library("readxl")
library("systemfit")
library("moments")
library("stargazer") 
library("MASS")
library("lmtest")
library("sandwich") 
library("plm")
library("gridExtra") 
library("olsrr") 
library("fBasics") 
library("ggplot2")
library("haven") 
library("dplyr")
library("pryr")
library("dynlm")
library("forecast")
library("readxl")
library("scales")
library("quantmod")
library("urca")
library("tis")
library("lubridate")
library("tidyr")
library("readr")
library("purrr")
library("tibble")
library("stringr")
library("forcats")
library("expss")
library("reshape")
library("tidyverse")
library("reshape2")
library("psych")
library("WDI")

#Set work directory
setwd("C:/Users/david/OneDrive/Desktop/Econometrics Project")
#### Preliminary operations

#Load IMF-WEO Database
imfweo <- read_excel("in_data/WEOOct2024all.xlsx", na = "n/a")
# Rename ISO 3 letter code to "ccISO3"
names(imfweo)[names(imfweo) == "ISO"]   <- "ccISO3"
#We also change the name of WEO Subject Code, for handiness
names(imfweo)[names(imfweo) == "WEO Subject Code"]   <- "var" 

# Keep only the variables that are needed: we got rid of everything, except for the country names, and the variables
imfweo <-  select(imfweo, ccISO3, var, as.character(1980:2023))

# Remove thousands separators (commas) and convert all year-columns to numeric
imfweo[as.character(1980:2023)] <- lapply(
  imfweo[as.character(1980:2023)],
  function(col) as.numeric(gsub(",", "", col))
)

# Clean some empty records, which are the last 2 in the dataset: we use the minus symbol and the comma followed by nothing
imfweo <- imfweo[-c(8625),]
imfweo <- imfweo[-c(8626),]

# Reshape using "melt"
imfweom <- melt(imfweo, id=(c("var", "ccISO3")))

# The na.omit R function removes all incomplete cases of the data frame at hand
imfweom <- na.omit(imfweom)

# Rename "year" the time variable
names(imfweom)[names(imfweom) == "variable"] <- "year"

# We want the time variable to be an integer.
imfweom$year <- as.integer(substr(imfweom$year, start = 1, stop = 4))

# Set missing values as "NA", as required by R
imfweom$value[imfweom$value == "n/a"] <- NA
imfweom$value[imfweom$value == ""] <- NA

# Drop observations where "ccISO3" is empty
imfweom <- dplyr::filter(imfweom, ccISO3 != "")

# Reshape again, using "cast"
imfweoc <- cast(imfweom, ccISO3+year~var)

# select subset of dataset with the variables required for our analysis  
imfweoCRI <- dplyr::select(imfweoc, ccISO3, year, NGDP_R, NGDP_D,PCPI, LUR)

# Filter the desired subset of countries. We use ISO3 country codes
#Let's select Costa Rica
imfweoCRI <- dplyr::filter(imfweoCRI, ccISO3 =="CRI")

# Format columns from characters to numeric. Also, scale them correctly.
imfweoCRI$NGDP_R <- as.numeric(as.character(imfweoCRI$NGDP_R))
imfweoCRI$NGDP_R <-  imfweoCRI$NGDP_R / 1000

imfweoCRI$NGDP_D <- as.numeric(as.character(imfweoCRI$NGDP_D))
imfweoCRI$NGDP_D[imfweoCRI$year >= 1982] <-  imfweoCRI$NGDP_D[imfweoCRI$year >= 1982] /1000

imfweoCRI$PCPI <- as.numeric(as.character(imfweoCRI$PCPI))
imfweoCRI$PCPI[imfweoCRI$year >= 1982] <-  imfweoCRI$PCPI[imfweoCRI$year >= 1982] /1000

imfweoCRI$LUR <- as.numeric(as.character(imfweoCRI$LUR))
imfweoCRI$LUR <- imfweoCRI$LUR /1000

# Format the column "year" to be considered as Date, so R can read it as time 
# series
imfweoCRI$year <- lubridate::ymd(imfweoCRI$year, truncated = 2L)

##### a) Represent the data, their growth rates

#a1) NGDP_R (GDP at constant prices or Real GDP) NAME plot NGDP_R_log

# Log GDP
NGDP_R_log <- log(imfweoCRI$NGDP_R)

# Represent GDP log
p1 <- ggplot(imfweoCRI, aes(x = year, y = NGDP_R_log)) + 
  geom_point(alpha = 0.9) +
  geom_line(col = "steelblue") + 
  labs(title = "GDP constant prices (log of NGDP_R) of Costa Rica (1980-2023)",
       y = "Log", x= "Year")

ggsave("out_data/NGDP_R_log.png", plot = p1, width = 8, height = 6, dpi = 300)


#Transform data in time series
NGDP_R_TS <- ts(imfweoCRI$NGDP_R, start = c(1980), end = c(2023),
              frequency = 1)

# Calculate GDP growth rate and plot NAME PLOT NGDP_R_GR

NGDP_R_GR <- growth.rate (NGDP_R_TS, lag = 1, simple = T)

png("out_data/NGDP_R_GR.png", width = 800, height = 600)
plot.ts(NGDP_R_GR, col = "steelblue", 
        main = "Costa Rica GDP Growth Rate - constant LCU (1980-2023)", 
        xlab = "Year", ylab = "% Growth per year", axes = TRUE)
dev.off()


# a2) NGDP_D (GDP deflator, an alternate measure of inflation)

# Log NGDP_D 
NGDP_D_log <- log(imfweoCRI$NGDP_D)

# Represent NGDP_D log name plt NGDP_D_log

p2 <- ggplot(imfweoCRI, aes(x = year, y = NGDP_D_log)) + 
  geom_point(alpha = 0.9) +
  geom_line(col = "darkred") + 
  labs(title = "GDP deflator (log of NGDP_D) of Costa Rica (1980-2023)",
       y = "Log", x= "Year")

ggsave("out_data/NGDP_D_log.png", plot = p2, width = 8, height = 6, dpi = 300)

#Transform data in time series

NGDP_D_TS <- ts(imfweoCRI$NGDP_D, start = c(1980), end = c(2023),
                frequency = 1)

# Computing NGDP_D growth rate and plotting it name NGDP_D_GR
NGDP_D_GR <- growth.rate (NGDP_D_TS, lag = 1, simple = T)


png("out_data/NGDP_D_GR.png", width = 800, height = 600)
plot.ts(NGDP_D_GR, col = "darkred", 
        main = "Costa Rica GDP Deflator Growth Rate (1980-2023)", 
        xlab = "Year", ylab = "% Growth per year", axes = TRUE)
dev.off()


# a3) PCPI (Inflation measured as yearly average CPI)

# Log PCPI
PCPI_log <- log(imfweoCRI$PCPI)

# Represent PCPI log name PCPI_log
p3 <- ggplot(imfweoCRI, aes(x = year, y = PCPI_log)) + 
  geom_point(alpha = 0.9) +
  geom_line(col = "darkgreen") + 
  labs(title = "Consumer Price Index (log of PCPI) of Costa Rica (1980-2023)",
       y = "Log", x= "Year")

ggsave("out_data/PCPI_log.png", plot = p3, width = 8, height = 6, dpi = 300)




#Transform data in time series
PCPI_TS <- ts(imfweoCRI$PCPI, start = c(1980), end = c(2023),
              frequency = 1)

# Compute PCPI growth rate and plot name PCPI_GR
PCPI_GR <- growth.rate (PCPI_TS, lag = 1, simple = T)

png("out_data/PCPI_GR.png", width = 800, height = 600)
plot.ts(PCPI_GR, col = "darkgreen", 
        main = "Costa Rica CPI Growth Rate (1980-2023)", 
        xlab = "Year", ylab = "% Growth per year", axes = TRUE)
dev.off()

# a4) LUR (Unemployment rate)
# Log LUR
LUR_log <- log(imfweoCRI$LUR)

# Represent LUR log name LUR_log
p4 <- ggplot(imfweoCRI, aes(x = year, y = LUR_log)) + 
  geom_point(alpha = 0.9) +
  geom_line(col = "purple") + 
  labs(title = "Unemployment Rate (log of LUR) of Costa Rica (1980-2023)",
       y = "Log", x= "Year")

ggsave("out_data/LUR_log.png", plot = p4, width = 8, height = 6, dpi = 300)



# Transform dataset in time series
LUR_TS <- ts(imfweoCRI$LUR, start = c(1980), end = c(2023),
             frequency = 1)

# Calculate LUR growth rate and plot name LUR_GR
LUR_GR <- growth.rate (LUR_TS, lag = 1, simple = T)

png("out_data/LUR_GR.png", width = 800, height = 600)
plot.ts(LUR_GR, col = "purple", 
        main = "Costa Rica Unemployment Rate Growth (1980-2023)", 
        xlab = "Year", ylab = "% Growth per year", axes = TRUE)
dev.off()



##### b) Present the correlogram of both the (logged) levels, and the growth 
#        rates of the indicated variables. 

#b1) We compute the autocorrelation for the logged levels.

# Plotting the correlogram for logged series
#names: NGDPR_R_corr etc
print(acf(na.omit(NGDP_R_log), plot = T, lag.max = 40, main = "Correlogram of Real GDP"))
print(acf(na.omit(NGDP_D_log), plot = T, lag.max = 40, main = "Correlogram of GDP deflator"))
print(acf(na.omit(PCPI_log), plot = T, lag.max = 40, main = "Correlogram of Inflation rate"))
print(acf(na.omit(LUR_log), plot = T, lag.max = 40, main = "Correlogram of unemployment rate"))

png("out_data/NGDP_R_corr.png", width = 800, height = 600)
acf(NGDP_R_log, main = "Correlogram of Real GDP")
dev.off()

png("out_data/NGDP_D_corr.png", width = 800, height = 600)
acf(NGDP_D_log, main = "Correlogram of GDP Deflator")
dev.off()

png("out_data/PCPI_corr.png", width = 800, height = 600)
acf(PCPI_log, main = "Correlogram of inflation rate")
dev.off()

png("out_data/LUR_corr.png", width = 800, height = 600)
acf(LUR_log, main = "Correlogram of unemployment rate")
dev.off()

# b2) We compute the autocorrelation for the growth rates.

png("out_data/NGDP_R_GR_corr.png", width = 800, height = 600)
acf(NGDP_R_GR, main = "Correlogram of Real GDP Growth")
dev.off()

png("out_data/NGDP_D_GR_corr.png", width = 800, height = 600)
acf(NGDP_D_GR, main = "Correlogram of GDP deflator growth")
dev.off()

png("out_data/PCPI_GR_corr.png", width = 800, height = 600)
acf(PCPI_GR, main = "Correlogram of inflation rate growth")
dev.off()

png("out_data/LUR_GR_corr.png", width = 800, height = 600)
acf(LUR_GR, main = "Correlogram of unemployment rate growth")
dev.off()

##### c) Are these time series stationary? Run an Augmented Dickey-Fuller
#        test, both on the logged levels, and on the first differences of the
#        log.

# Augmented Dickey-Fuller test for NGDP_R_log
summary(ur.df(NGDP_R_log, 
              type = "trend", 
              lags = 2, 
              selectlags = "Fixed"))

# COMMENT: the p-value, 0.42, is larger than alpha at 0.05, so H0:"data are non-stationary"
#isn't rejected.

# First differences of NGDP_R_log and Augmented Dickey-Fuller test for differences
NGDP_R_fd<- diff(NGDP_R_log, lag = 1, differences = 1)
summary(ur.df(NGDP_R_fd, 
              type = "trend", 
              lags = 2, 
              selectlags = "Fixed"))

# COMMENT: the p-value, 1.442e-06, is extremely lower than alpha at 0.05, so the
# null hypothesis can be strongly rejected with 95% confidence.

# Augmented Dickey-Fuller test for NGDP_D_log
summary(ur.df(NGDP_D_log, 
              type = "trend", 
              lags = 2, 
              selectlags = "Fixed"))

# COMMENT: the p-value, 3.625e-06, is extremely lower than alpha at 0.05, so the
# null hypothesis can be strongly rejected with 95% confidence.

# First differences of NGD_D_log and Augmented Dickey-Fuller test for 
# differences
NGDP_D_fd <- diff(NGDP_D_log, lag = 1, differences = 1)
summary(ur.df(NGDP_D_fd, 
              type = "trend", 
              lags = 2, 
              selectlags = "Fixed"))

# COMMENT: the p-value, 6.367e-05, is extremely lower than alpha at 0.05, so H0 can 
#be rejected with 95% confidence .

# Augmented Dickey-Fuller test for PCPI_log
summary(ur.df(PCPI_log, 
              type = "trend", 
              lags = 2, 
              selectlags = "Fixed"))

# COMMENT: the p-value, 3.544e-13, is extremely lower than alpha at 0.05, so the
# null hypothesis can be strongly rejected with 95% confidence.

# First differences of PCPI_log and Augmented Dickey-Fuller test for differences
PCPI_fd <- diff(PCPI_log, lag = 1, differences = 1)
summary(ur.df(PCPI_fd, 
              type = "trend", 
              lags = 2, 
              selectlags = "Fixed"))

# COMMENT: the p-value, 1.066e-06, is extremely lower than alpha at 0.05, so H0
#can be strongly rejected with 95% confidence.

# Augmented Dickey-Fuller test for LUR_log
summary(ur.df(LUR_log, 
              type = "trend", 
              lags = 2, 
              selectlags = "Fixed"))

# COMMENT: the p-value, 0.03152, is lower than alpha at 0.05, so H0 is rejected with
#95% confidence.

# First differences of LUR_log and Augmented Dickey-Fuller test for differences
LUR_fd <- diff(LUR_log, lag = 1, differences = 1)
summary(ur.df(LUR_fd, 
              type = "trend", 
              lags = 2, 
              selectlags = "Fixed"))

# COMMENT: the p-value, 3.887e-06, is extremely lower than alpha at 0.05, so the
# null hypothesis can be strongly rejected.

##### d) You may compute inflation using NGDP_D, or PCPI. Compare results in a) 
#        and b) using both measures of inflation. For the following question, 
#        use as a measure of inflation the one which is based on PCPI.

# Transform logged data to numeric
NGDP_D_log <- as.numeric(NGDP_D_log)
PCPI_log <- as.numeric(PCPI_log)


# Create a dataframe for the measures of inflation (logged)
Inflation <- data.frame(x = 1980:2023, y1 = PCPI_log,
                        y2 = NGDP_D_log)

# Plot the correlogram for inflation logs for the two measures of inflation
#Plot name PCPI_defl_log
p5 <- ggplot(imfweoCRI) +
  geom_line(aes(x = year, y = PCPI_log, colour = "CPI index")) +
  geom_line(aes(x = year, y = NGDP_D_log, colour = "GDP deflator")) +
  labs(title = "Log of inflation,Costa Rica (1980-2023)", y = "Log", x = "Year") +
  scale_colour_manual(values = c("CPI index" = "darkgreen", "GDP deflator" = "darkred"))

ggsave("out_data/PCPI_defl_log.png", plot = p5, width = 8, height = 6, dpi = 300)



# Plotting the correlogram for inflation growth rates for the two measures of inflation

Inflation_GR <- data.frame(x = 1981:2023, y1 = PCPI_GR,
                           y2 = NGDP_D_GR)
##Plot name PCPI_defl_GR
p6 <- ggplot(Inflation_GR, aes(x)) +
  geom_line(aes(y = y1, color = "deepgreen3"), lwd = 1) +
  geom_line(aes(y = y2, color= "steelblue"), lwd = 1) + 
  xlab("year") + ylab("Growth (yearly %)") + 
  ggtitle("Inflation growth rate, Costa Rica (1980-2023)") +     
  scale_color_manual(name= "Legend",labels = c("GDP deflator", "CPI index"), 
                     values = c("red", "blue")) 
ggsave("out_data/PCPI_defl_GR.png", plot = p6, width = 8, height = 6, dpi = 300)

## e) Estimate a "Phillips curve", where first differences of logged 
#     inflation are explained by their past history, and by the past levels 
#     of the unemployment rate. Find a suitable number of lags and comment 
#     your results.

# Y: First differences of logged inflation, PCPI_fd, is already known. Has to be transformed
#in time series
PCPI_fd <-ts(PCPI_fd, start = c(1980), end = c(2023),
             frequency = 1)
# X1: Past history of logged inflation: PCPI_log has to be transformed in a time series
PCPI_log <- ts(PCPI_log, start = c(1980), end = c(2023),
               frequency = 1 )
#X2: Past levels of the Unemployment rate
LUR_log <- ts(LUR_log, start = c(1980), end =c(2023), 
              frequency= 1)

### Find a suitable number of lags
# Create function BIC
BIC <- function(model) 
{  ssr <- sum(model$residuals^2)
t <- length(model$residuals)
npar <- length(model$coef)

return(round(c("p" = npar - 1,
               "BIC" = log(ssr/t) + npar * log(t)/t,
               "R2" = summary(model)$r.squared), 4))
}

# Loop the function BIC over multiple models 
order <- 1:6

# Compute BIC for different models

BICt <- sapply(order, function(x) 
  BIC(dynlm(PCPI_fd ~ L((PCPI_log), 1:x) + L((LUR_log), 1:x))))

BICt

# Select the model with the smallest BIC, which, as we can see, is the one with p=6

BICt[, which.min(BICt[2, ])]

# Estimate the relationship of the Phillips curve

phillips_t <- dynlm(PCPI_fd ~ PCPI_log + LUR_log )

# Create a table
stargazer(phillips_t, type = "text")


###########################THE END#######################################################

