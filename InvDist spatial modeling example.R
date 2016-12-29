# Geog503-2015-SpatialModelling-InvDist.r
#
# Predicts mean annual precipitation in the GVRD 
# using inverse distance weighting interpolation
# and uses leave-one-out cross validation to assess
# predictive accuracy
#
# 2015-Mar-2 RDM 
#
# update 2016-Mar-4 RDM
#############################################################


#########################################################################
# clear work space and read in data
#########################################################################

rm(list = ls())

library(sp)
library(rgdal)
library(raster)

# specify working directory
setwd("c:/Teaching/Geog503/Data")


#########################################################################
# process climate data - generate Spatial* object
#########################################################################

# precipitation and temperature normals for gvrd. Data is in x,y,(var1),(var2) format
pt = read.csv("gvrdpt.csv")    

# change coordinate names in advance of coordinate transformation to UTM later in script
names(pt)[names(pt) == "LONG"] = "x"
names(pt)[names(pt) == "LAT"] = "y"

# convert to Spatial* object
pt$x = -pt$x  # negative longitude for west of Greenwich
coordinates(pt) = c("x", "y")
proj4string(pt) =  CRS("+proj=longlat +ellps=WGS84")

# convert to UTM coordinates
pt.utm = spTransform(pt, CRS("+proj=utm +zone=10"))
x = pt.utm$x
y = pt.utm$y
P = pt.utm$PANN
ELEV = pt.utm$ELEV

# add random noise to coordinates to avoid duplicate locations 
#
# ____ why you ask? there are actually stations which are so close 
#______ together that you can't really use them as separate points. 
#_______ if would be possible to have a nearest neighbor point of 0
# _______ add random noise to coordinates so that they don't fall on 
# ________exactly on the same location
n = length(x)
x = x + runif(n, -250, 250)
y = y + runif(n, -250, 250)

# create new SpatialPointsDataFrame with modified coordinates
pt2 = data.frame(x, y, pt.utm@data)
coordinates(pt2) = ~x + y
proj4string(pt2) = CRS("+proj=utm +zone=10")


#########################################################################
# get shoreline vectors from Geogratis and generate Spatial* object
#########################################################################

# specify directory (folder) containing shape files
shpdir = "c:/Teaching/Geog503/Data/Geogratis-GVRD/canvec_160222_122645_shp"

# specify layer to extract - shorelines
used.layer = c("hd_1440009_1")

# read layer into Spatial* objects and convert to UTM
sf1 = readOGR(dsn = shpdir, layer = used.layer)
sf1.utm = spTransform(sf1, CRS("+proj=utm +zone=10 +ellps=WGS84"))
class(sf1.utm)


#########################################################################
# define function for inverse distance interpolation
#########################################################################

#_________this is a home-made function - there is a better built in one
invdist = function(cp, pp, max.rad = 10000, k = 2) {

  # inverse distance interpolation
  # cp = matrix of (x, y, z) values for control points
  # pp = matrix of (x, y) values for prediction points
  # the function returns a matrix of (x, y, z, np) for prediction points
  #   where np is the number of control points used in the interpolation
  # if no control points lie within max.rad of the prediction point, a
  #   value of NA is returned
  
  xc = cp[, 1]
  yc = cp[, 2]
  zc = cp[, 3]
  xp = pp[, 1]
  yp = pp[, 2]
  nc = length(xc)
  np = length(xp)
  zp = numeric(np)
  nip = numeric(np)
  for (i in 1:np) {
    # compute distances between prediction point and control points
    dist = sqrt((xp[i] - xc)^2 + (yp[i] - yc)^2)
    # extract subset based on max.rad
    ip = which(dist < max.rad)
    nip[i] = length(ip)
    # error trap: check that there is at least one control point with dist < max.rad
    if (nip[i] == 0) { # this just returns an NA value if you have no control pts in your radius
      zp[i] = NA
    } else {
      vip = 1/(dist[ip]^k)
      wip = vip/sum(vip)
      zip = zc[ip]
      zp[i] = sum(wip*zip)
	}  
  }
  return(cbind(xp, yp, zp, nip))
}
 


#########################################################################
# test invdist
#########################################################################

# south-north transect
cpts = cbind(x, y, P)
xp = rep(495000, 5)
yp = seq(5430000, 5470000, 10000)
xyp = cbind(xp, yp)
xyzp = invdist(cpts, xyp)
xyzp

# west-east transect
cpts = cbind(x, y, P)
yp = rep(5455000, 5)
xp = seq(485000, 525000, 10000)
xyp = cbind(xp, yp)
xyzp = invdist(cpts, xyp)
xyzp

#########################################################################
# leave-one-out cross validation
#########################################################################

cpts = cbind(x, y, P)   # matrix of station data
nc = length(cpts[, 1])
invdist.loo = numeric()         # empty object to hold output     

mr = 10000  # value for max.rad (m)
kk = 2      # exponent in weighting function

# for assignment: in the loop, fit regression, calculate residuals, calculate prediction based on residuals? compare to thing that is there? 
# the point is to interpolate residuals
for (i in 1:nc) {
  # withhold point i from control points and make a prediction for it
  cp = cpts[-i, ]
  pp = rbind(cpts[i, 1:2])
  mod.1 = lm(cp[,3] ~ cp[,1] + cp[,2])
  res.mod = resid(mod.1)
  new.row = invdist(cp, pp, max.rad = mr, k = kk)
  # add results for point i to invdist.loo
  invdist.loo = rbind(invdist.loo, new.row)
  cpts = rbind(cp[,2:3], res.mod) # this gives you 
}

# error statistics
errors = P - invdist.loo[, 3]
MBE = mean(errors)
RMSE = sqrt(mean(errors^2))
MAE = mean(abs(errors))
MBE
RMSE
MAE


#########################################################################
# generate map of prediction error
#########################################################################\

par(mfrow = c(1, 1))
# plot spatial pattern of errors using blue + for positive and red - for negative
psym = ifelse(errors > 0, "+", "-")
pcol = ifelse(errors > 0, "blue", "red")
psize = sqrt(abs(errors)/100)
plot(pt.utm, 
  pch = psym, col = pcol, cex = psize  
)
lines(sf1.utm)

# add a legend
legsymsize = c(100, 300, 500, -100, -300, -500)
legpch = c("+", "+", "+", "-", "-", "-")
legend("bottom", bty = "n", bg = "white",
  ncol = 6,
  pt.cex = sqrt(abs(legsymsize)/100),
  pch = legpch,
  col = c("blue", "blue", "blue", "red", "red", "red"),
  legend = as.character(legsymsize)
)

#########################################################################
# plot a set of diagnostic plots
#########################################################################

par(mfrow = c(2, 2), mar = c(5, 5, 1, 1))

# plot spatial pattern of errors using blue + for positive and red - for negative
psym = ifelse(errors > 0, "+", "-")
pcol = ifelse(errors > 0, "blue", "red")
psize = sqrt(abs(errors)/100)
plot(pt.utm, 
  pch = psym, col = pcol, cex = psize  
)
lines(sf1.utm)

# add a legend
legsymsize = c(100, 300, 500, -100, -300, -500)

legpch = c("+", "+", "+", "-", "-", "-")
legend("bottom", bty = "n", bg = "white",
  ncol = 6,
  pt.cex = sqrt(abs(legsymsize)/100),
  pch = legpch,
  col = c("blue", "blue", "blue", "red", "red", "red"),
  legend = as.character(legsymsize)
)

# scatterplot of predicted and observed
plot(P, invdist.loo[, 3], 
  pch = 21, bg = "red",
  ylab = "Predicted MAP (mm)",
  xlab = "Observed MAP (mm)"  
)
abline(0, 1)
legend("bottomright", bty = "n", lty = 1, legend = "1:1 line")


# scatterplot of error vs elevation
plot(ELEV, errors, 
  pch = 21, bg = "red",
  ylab = "Prediction error (mm)",
  xlab = "Elevation (masl)"  
)
abline(lm((errors ~ ELEV)))
abline(h = 0, lty = 2)
legend("bottomright", bty = "n", lty = 1, legend = "regression line")


# add listing of interpolation options and statistics to fourth graph panel
#  - type = "n" suppresses plotting of points
#  - by = "n" suppresses plotting of box around graph
#  - xaxt = "n" and yaxt = "n" suppress plotting of axes
#  - ylab = "" and xlab = "" suppress axis labels
plot(ELEV, errors, type = "n", bty = "n",
  xaxt = "n", yaxt = "n",
  ylab = "",
  xlab = ""  
)

# create character strings to be added as a legend
mbetxt = paste("MBE = ", formatC(MBE, digits = 2), "mm")
rmsetxt = paste("RMSE = ", formatC(RMSE), "mm")
maetxt = paste("MAE = ", formatC(MAE), "mm")
mrtxt = paste("max.rad = ", formatC(mr), "m")
ktxt = paste("k = ", formatC(kk))
legend("left", bty = "n", 
  legend = c(mrtxt, ktxt, mbetxt, rmsetxt, maetxt)
)


#########################################################################
# generate map of mean annual precipitation
#########################################################################

par(mfrow = c(1, 1))

# define grid topology
xmin = bbox(pt.utm)[1, 1]
xmax = bbox(pt.utm)[1, 2]
ymin = bbox(pt.utm)[2, 1]
ymax = bbox(pt.utm)[2, 2]
dx = 100
nx = ceiling((xmax - xmin)/dx) + 1
ny = ceiling((ymax - ymin)/dx) + 1
gt = GridTopology(c(xmin, ymin), c(dx, dx), c(nx, ny))

# create SpatialGrid and generate coordinates
xygrid = SpatialGrid(grid = gt, proj4string = CRS("+proj=utm +zone=10"))
xy = coordinates(xygrid)

# interpolate mean annual precipitation for each grid point
Pid = invdist(cp = cpts, pp = xy, max.rad = 15000, k = 2)
head(Pid)
summary(Pid[, 4])
summary(Pid[, 3])

# create a raster object - first create a SpatialPointsDataFrame
Pid = as.data.frame(Pid)
coordinates(Pid) = ~xp + yp

# convert to SpatialPixelsDataFrame and then to raster
Pid.pix = as(Pid, "SpatialPixelsDataFrame")
Pid.ras = raster(Pid.pix, layer = 1)

# generate map
plot(Pid.ras,
  main = "Mean annual precipitation by inverse-distance interpolation"
)
lines(sf1.utm)
points(pt.utm, pch = 21, bg = "black")


