knitr::opts_chunk$set(comment = NA, tidy = TRUE)

rm(list=ls())
load("../data/CO-temp-data.RData")
ls()

library(fields)
library(MBA)
library(geoR)
library(sf)
library(leaflet)

blue.red <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")

base.map <- leaflet(width = "100%") %>%
    addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
    addLayersControl(baseGroup = c("Satellite"), options = layersControlOptions(collapsed = FALSE))

pal <- colorNumeric(blue.red, domain = temp)

base.map %>%
    addCircleMarkers(lng = coords[, 1], lat = coords[, 2], col = pal(temp), stroke = FALSE,
        radius = 5, fillOpacity = 0.9, popup = paste("Mean min temp:", round(temp,
            1))) %>%
    addLegend("bottomright", pal = pal, values = temp, opacity = 0.9, title = "Temperature C")

lon2UTM <- function(longitude) {
    (floor((longitude + 180)/6) %% 60) + 1
} ##Function to find the UTM zone of a geographic coordinate in the western hemisphere (only longitude is needed).

Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }

  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
} ##Function to calculate the mode of the above data to find which zone contains most of the points
utmZone <- Mode(lon2UTM(coords$lon))
utmZone

coordsGeog <- st_as_sf(coords, coords=c("lon","lat"), crs="+proj=longlat +datum=WGS84 +no_defs")
coordsUTM <- st_transform(coordsGeog, crs="+proj=utm +zone=13 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
coordsUTM <- st_coordinates(coordsUTM)/1000
class(coordsUTM)

elev <- elev/1000
df <- data.frame(temp, elev, coordsUTM)
names(df) <- c("temp", "elev", "xUTM", "yUTM")

library("nimble")

lmCode <- nimbleCode({
        for (i in 1:N) {
                Y[i] ~ dnorm(mu[i], tausq)
                mu[i] <- inprod(X[i,1:p], beta[1:p])
            yFit[i] ~ dnorm(mu[i], tausq)
	    yResid[i] <- Y[i] - mu[i]
    }

    for (i in 1:p) { beta[i] ~ dflat()}
    tausq ~ dgamma(1.0E-5, 1.0E-5)
    sigma <- 1/sqrt(tausq)
}
)   

modelParameters <- c("beta", "tausq", "sigma", "yFit", "yResid")

nimbleData <- list(Y = df$temp)
N <- nrow(df)
X <- cbind(rep(1, times=N), df$xUTM, df$yUTM)
p <- ncol(X)
nimbleConstants <- list(N=N, X=X, p=p)
nimbleInits <- list(beta=rep(0, times=p), tausq=1.0, yFit=rep(0, times=N))

rModel <- nimbleModel(code=lmCode, constants=nimbleConstants, data=nimbleData, inits=nimbleInits)

n_iter = 11000 ## Number of iterations
n_chains = 1 ## Number of different chains
n_burn = 1000 ## Number of initial iterations ("burn" period) for MCMC to converge. These samples will be discarded and not used in posterior inference. 

mcmc.out <- nimbleMCMC(model = rModel, nchains=n_chains, niter = n_iter, nburnin=n_burn, monitor=modelParameters, WAIC=TRUE)

mcmc.out$WAIC

samps <- rbind(mcmc.out$samples)
dim(samps)

credibleIntervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))

CI_beta <- credibleIntervals[grep("beta", rownames(credibleIntervals)),]
CI_beta

CI_sigma <- credibleIntervals[grep("sigma", rownames(credibleIntervals)),]
CI_sigma

CI_yResid <- credibleIntervals[grep("yResid", rownames(credibleIntervals)),]
yResidPosteriorMedians <- CI_yResid[,grep("50%", colnames(CI_yResid))]

dMax <- max(rdist(coordsUTM))
dMax

vTemp <- variog(coords=coordsUTM, data=df$temp, uvec=(seq(0, 0.75*dMax, length=20)))

vResid <- variog(coords=coordsUTM, data=yResidPosteriorMedians, uvec=(seq(0, 0.75*dMax, length=20)))


par(mfrow=c(1,2))
plot(vTemp, xlab="Distance (km)")
plot(vResid, xlab="Distance (km)")

resid.surf <- mba.surf(cbind(coords, yResidPosteriorMedians), no.X = 200, no.Y = 200, extend = TRUE, sp = TRUE)$xyz.est

library(sp)
library(raster)

proj4string(resid.surf) <- "+proj=longlat +datum=WGS84"

resid.surf <- raster(resid.surf)

pal <- colorNumeric(blue.red, values(resid.surf), na.color = "transparent")

base.map %>%
    addRasterImage(resid.surf, colors = pal, opacity = 0.75, group = "Regression residuals") %>%
    addLegend("bottomright", pal = pal, values = values(resid.surf), opacity = 0.75,
        title = "<center>Regression<br> residuals</center>") %>%
    addLayersControl(baseGroup = c("Satellite"), overlayGroups = c("Regression residuals"),
        options = layersControlOptions(collapsed = FALSE))

