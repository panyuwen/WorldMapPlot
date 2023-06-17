###############################################################################
# World map using ggplot
# R 3.5.1
# The central/prime meridian can be shifted
# Enhanced aspect with graticules
# accept projection of known PROJ.4 string
###############################################################################

## by panyuwen
## 2020/06

## Reference:
## To find projections of known PROJ.4 string from, https://proj.org/usage/index.html, 
##                                                  https://proj.org/operations/projections/index.html
##                                               OR https://pro.arcgis.com/zh-cn/pro-app/help/mapping/properties/list-of-supported-map-projections.htm
## This script was initially borrowed from https://github.com/valentinitnelav/valentinitnelav.github.io/blob/master/gallery/Pacific%20centered%20world%20map%20with%20ggplot.R
## The code for the "Nine-dash line" was borrowed from https://blog.csdn.net/nikang3148/article/details/82871006
##

library(data.table)
library(ggplot2)
library(rgdal)
library(rgeos)
library(maps)
library(maptools)
library(raster)
library(metR)
library(reshape2)
library(akima)
library(sp)
library(scatterpie)
library(openxlsx)


# =============================================================================
# main functions
# =============================================================================

## coordinates projection, input shp data, and output data table
## using readOGR() to read the shp data
## centerlong stands for the central longitude in the final plot
## projection is a project of known PROJ.4 string, such as eqc, hammer, ortho, wintri, and so on.

## two functions provided here, for the same purpose
## shp2finaldf2() is more stable, but cost time and memory
## shp2finaldf() is much faster, but there remains some bugs, which should be obvious when u get the fig

## centerlong and centerlat indicate the geographic coordinates of the center in the final plot

shp2finaldf <- function(shp_data, centerlong, centerlat, projection){
    PROJ <- paste("+proj=",projection," +lon_0=",as.character(centerlong)," +lat_0=",as.character(centerlat)," +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep='')
    #PROJ <- "+proj=hammer +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
    # transform split country polygons in a data table that ggplot can use
    DT <- data.table(map_data(as(shp_data, "SpatialPolygonsDataFrame")))
    # project coordinates 
    DT[, c("X","Y") := data.table(project(cbind(long, lat), proj=PROJ))]

    return(DT)
}

shp2finaldf2 <- function(shp_data, centerlong, projection){
    # inspired from: https://stat.ethz.ch/pipermail/r-sig-geo/2015-July/023168.html

    # shift central/prime meridian towards west - positive values only
    shift <- 360 - centerlong
    # create "split line" to split country polygons
    WGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    split.line = SpatialLines(list(Lines(list(Line(cbind(180-shift,c(-90,90)))), ID="line")), proj4string=WGS84)
    # NOTE - in case of TopologyException' errors when intersecting line with country polygons,
    # apply the gBuffer solution suggested at:
    # http://gis.stackexchange.com/questions/163445/r-solution-for-topologyexception-input-geom-1-is-invalid-self-intersection-er
    # shp_data <- gBuffer(shp_data, byid=TRUE, width=0)

    # intersecting line with country polygons
    line.gInt <- gIntersection(split.line, shp_data)
    # create a very thin polygon (buffer) out of the intersecting "split line"
    bf <- gBuffer(line.gInt, byid=TRUE, width=0.000001)  
    # split country polygons using intersecting thin polygon (buffer)
    shp_data.split <- gDifference(shp_data, bf, byid=TRUE)
    # plot(shp_data.split) # check map
    class(shp_data.split) # is a SpatialPolygons object

    PROJ <- paste("+proj=",projection," +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep='')
    # transform split country polygons in a data table that ggplot can use
    Country.DT <- data.table(map_data(as(shp_data.split, "SpatialPolygonsDataFrame")))
    # Shift coordinates
    Country.DT[, long.new := long + shift]
    Country.DT[, long.new := ifelse(long.new > 180, long.new-360, long.new)]
    # project coordinates 
    Country.DT[, c("X","Y") := data.table(project(cbind(long.new, lat), proj=PROJ))]

    return(Country.DT)
}

## create graticules
## dis stands for the interval among E-W grid lines
## projection is a project of known PROJ.4 string
makeGraticules <- function(centerlong, centerlat, projection, longdis=45, latdis=30, minlat=-90, maxlat=90){
    # =============================================================================
    # Create graticules
    # =============================================================================
    # create a bounding box - world extent, Spatial-class
    b.box <- as(raster::extent(-180, 180, minlat, maxlat), "SpatialPolygons")
    # assign CRS to box
    WGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    proj4string(b.box) <- WGS84
    # create graticules/grid lines from box, Create N-S and E-W grid lines over a geographic region; create and plot corresponding labels
    longbreaks <- seq(-360, 360, by=longdis)
    longbreaks <- c(longbreaks[(longbreaks>=(centerlong-180)) & (longbreaks<=(centerlong+180))]) - centerlong
    longbreaks <- longbreaks[order(longbreaks)]
    longbreaks <- longbreaks[!duplicated(longbreaks)]
    
    latbreaks <- seq(-90, 90, by=latdis)
    latbreaks <- c(minlat, latbreaks[(latbreaks>minlat) & (latbreaks<maxlat)], maxlat)

    grid <- gridlines(b.box, easts  = longbreaks, norths = latbreaks)

    # give the PORJ.4 string for Eckert IV projection
    PROJ <- paste("+proj=",projection," +lon_0=0 +lat_0=",centerlat," +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep='')
    # transform graticules from SpatialLines to a data table that ggplot can use
    grid.DT <- data.table(map_data(SpatialLinesDataFrame(sl=grid, 
                                                       data=data.frame(1:length(grid)), 
                                                       match.ID = FALSE)))
    # project coordinates
    # assign matrix of projected coordinates as two columns in data table
    grid.DT[, c("X","Y") := data.table(project(cbind(long, lat), proj=PROJ))]
    grid.DT = grid.DT[is.finite(grid.DT$X),]

    return(grid.DT)
}

## convert "Nine-dash line" shp data
## works in the same way as shp2finaldf()
convertl9 <- function(shp_data,centerlong, centerlat, projection){
    l9 <- fortify(shp_data)
    #shift <- 360 - centerlong
    #l9$long.new <- l9$long + shift
    #l9$long.new <- ifelse(l9$long.new > 180, l9$long.new-360, l9$long.new)

    # project coordinates
    PROJ <- paste("+proj=",projection," +lon_0=",as.character(centerlong)," +lat_0=",as.character(centerlat)," +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep='')
    l9[, c("X","Y")] <- data.table(project(cbind(l9$long, l9$lat), proj=PROJ))

    return(l9)
}

## project a given dataframe
## must have <long> and <lat> in the columns
convertData <- function(df, centerlong, centerlat, projection){
    PROJ <- paste("+proj=",projection," +lon_0=",as.character(centerlong)," +lat_0=",as.character(centerlat)," +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep='')
    # project coordinates
    df[, c("X","Y")] <- data.table(project(cbind(df$long, df$lat), proj=PROJ))

    return(df)
}


# =============================================================================
# examples
# =============================================================================
# read the shapefile with readOGR

## world map
shp <- readOGR(dsn = "中国标准世界国家区划/ESRI世界地图中文版.shp")


## example 1
centerlong = 80
centerlat = 30
dis = 45
projection = 'ortho'

mapWorld <- shp2finaldf(shp, centerlong, centerlat, projection)
mapWorld <- mapWorld[which(is.finite(mapWorld$X)),]

china_map = readShapePoly('CHN_adm 中国行政区域空间文件，无台 湾/CHN_adm2.shp')
xj = china_map[which(china_map$NAME_1=='Xinjiang Uygur'),]
xj <- shp2finaldf(xj, centerlong, centerlat, projection)

# Nine Line
southChinaSea <- readOGR(dsn = "中国标准世界国家区划/南海九段线.shp")
southChinaSea <- convertl9(southChinaSea, centerlong, 0, projection)

# Wallace Line
hu = read.csv('中国标准世界国家区划/wh_line.txt',sep='\t',header=T)
hu <- convertData(hu, centerlong, 0, projection)
linelabel <- data.frame('label'=c("Wallace-Huxley's\nLine"),'lat'=c(-14.5),'long'=c(110))
linelabel <- convertData(linelabel, centerlong, 0, projection)

mapgrid <- makeGraticules(centerlong, centerlat, projection, 45)
outbound <- makeGraticules(0, 0, projection, 45)

## pie
data = read.csv('2_109513601.freq.txt',sep='\t',header=T)
data = convertData(data, centerlong, centerlat, projection)
data = data[which(is.finite(data$X)),]
data$radius <- ifelse(data$pop=='UYG', 550000, 350000)
data$REFfreq <- 1.0 - data$ALTfreq

# plot
p <- ggplot() +
  # background color
  geom_polygon(data = outbound[(long %in% c(-90,90) & region == "NS") | (long %in% c(-90,90) & lat %in% c(-90,90) & region == "EW")], aes(x = X, y = Y, group = group), colour = NA, fill = "#FEFEFE") +  ## change the color of ocean, rather than setting the color of panel background then mask the outer part 
  # world map
  geom_polygon(data = mapWorld,aes(x = X, y = Y, group = group), colour = "gray97",  fill = "gray80",size = 0.12) +
  # xinjiang map
  geom_polygon(data = xj,aes(x = X, y = Y, group = group), colour = NA,  fill = "#CE89E9",size = 0) +  
  # wallace line
  geom_path(data = hu, aes(x = X, y = Y, group = group), colour = "red", size = 0.25, linetype='longdash') +
  # Nine Line
  geom_path(data = southChinaSea,aes(x = X, y = Y, group = group), colour = "gray65",size = 0.2) +
  # pie plot
  geom_scatterpie(data=data,aes(x=X, y=Y, r=radius),cols=c('ALTfreq','REFfreq'),color='white', alpha=0.8, size=0.25, sorted_by_radius=F) +
  scale_fill_manual(values=c('#fb8a2e','#00B2EE'),breaks=c('ALTfreq','REFfreq'),labels=c('DAF','AAF')) +
  # add graticules
  geom_path(data = mapgrid,aes(x = X, y = Y, group = group),linetype = "solid", colour = "gray65", size = 0.15) +
  # add a bounding box (select graticules at edges)
  geom_path(data = outbound[(long %in% c(-90,90) & region == "NS") | (long %in% c(-90,90) & lat %in% c(-90,90) & region == "EW")], aes(x = X, y = Y, group = group), linetype = "solid", colour = "black", size = .3) +
  # ensures that one unit on the x-axis is the same length as one unit on the y-axis
  coord_equal(expand=T) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  #theme(panel.grid = element_blank(), panel.background = element_blank()) + 
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position=c(0.85,0.18),legend.background=element_rect(size=0.45,fill='white',color='black'),legend.title=element_blank(),legend.text=element_text(size=20),legend.key.size=unit(1.2,'cm')) +
  guides(fill=guide_legend(ncol=1))
  #theme_void() +
ggsave(file="example_1.202006.pdf",p,width=18, height=18, units="cm")



## example 2
centerlong = 60
centerlat = 0
projection = 'wintri'
minlat = -45
maxlat = 75

# background
cutting <- data.frame(long=c(-180,-180,+180,+180),lat=c(maxlat,minlat,minlat,maxlat))
cutting <- convertData(cutting, 0, 0, projection)
cutting$group <- c(1,2,2,1)

region <- data.frame(long=c(-180,-180,+180,+180),lat=c(maxlat,minlat,minlat,maxlat))
region <- convertData(region, 0, 0, projection)

outbound <- makeGraticules(0, 0, projection, minlat=minlat, maxlat=maxlat)

## map
mapWorld <- shp2finaldf2(shp, centerlong, projection)
mapWorld <- mapWorld[which((mapWorld$lat<maxlat) & (mapWorld$lat>minlat)),]

# Nine Line
southChinaSea <- readOGR(dsn = "中国标准世界国家区划/南海九段线.shp")
southChinaSea <- convertl9(southChinaSea, centerlong, centerlat, projection)

# Wallace Line
hu = read.csv('中国标准世界国家区划/wh_line.txt',sep='\t',header=T)
hu <- convertData(hu, centerlong, 0, projection)
linelabel <- data.frame('label'=c("Wallace-Huxley's\nLine"),'lat'=c(-14.5),'long'=c(110))
linelabel <- convertData(linelabel, centerlong, 0, projection)

## point
prop <- read.csv('proportion.txt',sep='\t',header=T)
prop <- convertData(prop, centerlong, centerlat, projection)
prop$arc <- ifelse(prop$arc > 0.00201, 0.00201, prop$arc)
prop <- prop[order(prop$arc,decreasing =FALSE),]

p <- ggplot() + 
  # background of the bounding box 
  geom_polygon(data = region, aes(x=X,y=Y), fill='#DFF2FC') +
  geom_polygon(data = outbound[(long %in% c(-180,180) & region == "NS") |(long %in% c(-180,180) & lat %in% c(minlat,maxlat) & region == "EW")], aes(x = X, y = Y, group = group), color=NA, fill='#DFF2FC') +
  # map
  geom_polygon(data = mapWorld,aes(x = X, y = Y, group = group), colour = "#FFF3E2",  fill = "#FFF3E2",size = 0.25) +
  # bounding box 
  geom_path(data = outbound[(long %in% c(-180,180) & region == "NS") |(long %in% c(-180,180) & lat %in% c(minlat,maxlat) & region == "EW")], aes(x = X, y = Y, group = group), linetype = "solid", colour = "black", size = .3) +
  geom_path(data = cutting,aes(x = X, y = Y, group = group), colour = "black",size = 0.3) +
  # Nine Line
  geom_path(data = southChinaSea,aes(x = X, y = Y, group = group), colour = "#FFF3E2",size = 0.65) +
  # wallace line
  geom_path(data = hu, aes(x = X, y = Y, group = group), colour = "red", size = 0.75, linetype='longdash') +
  geom_text(data=linelabel,aes(x = X, y = Y, label=label), size=4.2,color="gray55",family="Helvetica",fontface='italic') +
  # point
  geom_point(data = prop, aes(x = X, y = Y, fill = arc*100), size = 3.6, color = '#B2B2B2', shape = 21) +
  scale_fill_gradient2(high = "red", mid = "yellow", low = "blue", midpoint = 0.1,breaks = c(0.2,0.1,0.01), labels = c('>0.2%','0.1%','0.01%'),name='Ancestry') +
  # ensures that one unit on the x-axis is the same length as one unit on the y-axis
  coord_equal(ylim = c(min(cutting$Y), max(cutting$Y)),expand=T) + # same as coord_fixed(ratio = 1)
  theme(panel.grid = element_blank(), panel.background = element_blank()) + 
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=16), legend.position=c(0.835,0.4),legend.background=element_blank()) +
  theme(plot.title=element_text(size=32,hjust=0.5))
  #theme_void()
ggsave(file="example_2.202006.pdf",p,width=40.3, height=19.4, units="cm")



## example 3
centerlong = 110
centerlat = 0
projection = 'ortho'

mapWorld <- shp2finaldf(shp, centerlong, centerlat, projection)
mapWorld <- mapWorld[is.finite(mapWorld$X),]

# Nine Line
southChinaSea <- readOGR(dsn = "中国标准世界国家区划/南海九段线.shp")
southChinaSea <- convertl9(southChinaSea, centerlong, 0, projection)

# Wallace Line
hu = read.csv('中国标准世界国家区划/wh_line.txt',sep='\t',header=T)
hu <- convertData(hu, centerlong, 0, projection)
linelabel <- data.frame('label'=c("Wallace-Huxley's\nLine"),'lat'=c(-14.5),'long'=c(110))
linelabel <- convertData(linelabel, centerlong, 0, projection)

# grids, 1800 * 3600
xbreaks= seq(-180,180-0.1,by=0.1)
ybreaks = seq(90-0.1,-90,by=-0.1)
supergrid = data.frame(matrix(ncol = length(ybreaks), nrow = length(xbreaks)))
names(supergrid) = ybreaks
row.names(supergrid) = xbreaks

supergrid$long = as.numeric(row.names(supergrid))
supergrid = melt(supergrid,id.vars=c('long'))
names(supergrid) = c('long','lat','aasn')
supergrid$long = as.numeric(as.character(supergrid$long))
supergrid$lat = as.numeric(as.character(supergrid$lat))

# data
compdata = read.xlsx('popmean.xlsx')
compdata = compdata[!is.na(compdata$Latitude),]
compdata$aASN = compdata$aASN *100

prop = compdata[,c('Latitude','Longitude','aASN','Region')]
names(prop) = c('lat','long','aasn','Region')
prop$other = 1.0 - prop$aasn
prop <- convertData(prop, centerlong, centerlat, projection)
prop <- prop[is.finite(prop$X),]

mapgrid <- makeGraticules(centerlong, centerlat, projection, 45)
mapgrid = mapgrid[which((mapgrid$long>-90) & (grid.DT$long<90)),]

outbound <- makeGraticules(0, 0, projection, 45)

## make grid & filter data
griddata <- convertData(supergrid, centerlong, centerlat, projection)
griddata <- griddata[is.finite(griddata$X),]
inout = over(
    SpatialPoints(griddata[,c("long","lat")],proj4string=CRS(projection(shp))),
    as(shp,"SpatialPolygons")
)
griddata <- griddata[!is.na(inout),]

# new grid of 800*800, based on the converted coordinate
bound = data.frame('long'=c(-90,90,0,0),'lat'=c(0,0,90,-90))
bound <- convertData(bound, 0, 0, projection)
maxX = max(bound$X); maxY = max(bound$Y); minX = min(bound$X); minY = min(bound$Y)
xdis = as.integer((maxX - minX) / 800)
ydis = as.integer((maxY - minY) / 800)

# group into grids, ori data
prop$newX = as.integer(prop$X / xdis) * xdis
prop$newY = as.integer(prop$Y / ydis) * ydis
prop = prop[,c('long','lat','aasn','newX','newY')]
names(prop) = c('long','lat','aasn','X','Y')
prop = prop[order(prop$aasn,decreasing=T),]
# group into grids, empty grids
regrid = griddata
regrid$X = as.integer(regrid$X / xdis) * xdis
regrid$Y = as.integer(regrid$Y / ydis) * ydis

regrid = rbind(prop, regrid)
regrid = regrid[!duplicated(regrid[,c('X','Y')]),]

gridempty = regrid[is.na(regrid$aasn),]
gridvalue = regrid[!is.na(regrid$aasn),]

# geographic distance, for scaling
totlen = dist(matrix(c(minX,maxX,0,0),nrow=2))[1]
# imputation
annot <- function(xpos,ypos){
  tmp = gridvalue
  # distance to all the points in the ori data
  tmp$dis = apply(tmp,1,function(x){dist(matrix(c(as.numeric(x['X']),as.numeric(xpos),as.numeric(x['Y']),as.numeric(ypos)),nrow=2))[1]})
  # contrbutions from all the points in the ori data
  tmp$contri = tmp$aasn * (1.0 - (tmp$dis / totlen))
  # weighted average contribution
  value = sum(tmp$contri / (tmp$dis)^2.5) / sum(1.0 / (tmp$dis)^2.5)
}
gridempty$aasn = apply(gridempty,1,function(x){annot(x['X'],x['Y'])})

gridmerge = rbind(gridvalue, gridempty)
# better to save gridmerge
save(gridmerge, file='example_3.202006.RData')

## plot
p <- ggplot() +
  # add graticules
  geom_path(data = mapgrid,aes(x = X, y = Y, group = group),linetype = "solid", colour = "gray65", size = 0.15) +
  geom_tile(data=gridmerge, aes(x=X, y=Y, fill = aasn)) +
  #geom_polygon(data = mapWorld,aes(x = X, y = Y, group = group), colour = "gray85",  fill = NA,size = 0.08) +
  geom_path(data = hu, aes(x = X, y = Y, group = group), colour = "red", size = 0.25, linetype='longdash') +
  # add a bounding box (select graticules at edges)
  geom_path(data = outbound[(long %in% c(-90,90) & region == "NS") | (long %in% c(-90,90) & lat %in% c(-90,90) & region == "EW")], aes(x = X, y = Y, group = group), linetype = "solid", colour = "black", size = .3) +
  scale_fill_gradient2(high = "purple", mid = "#E066FF", low = "yellow", midpoint = 45,limits = c(0,100),breaks = c(90,50,10), name='aASN (%)') +
  coord_equal(expand=T) +
  theme(panel.grid = element_blank(), panel.background = element_blank()) + 
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  theme(legend.text=element_text(size=10),legend.position=c(0.35,0.3))
  #theme_void() +
ggsave(file="example_3.202006.pdf",p,width=18, height=12, units="cm")

