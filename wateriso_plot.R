map = function(froot){

  library(sp)
  library(rgdal)
  library(maptools)
  library(RColorBrewer)
  library(classInt)
  library(Grid2Polygons)
  library(raster)
    
  dbox = switch(Sys.info()["nodename"], "GJB-ZEN"="D:/Dropbox/", 
                "HYDROGEN"="C:/Users/gjbowen/Dropbox/")
  filedir = paste0(dbox, "Archived/Utilities/wateriso_plots/")
  setwd(paste0(dbox, "Archived/Utilities/IsotopeMaps/"))

  #setup - shapefile names
  wld.shp = "continents.shp"
  af.shp = "Africa.shp"
  as.shp = "Asia.shp"
  au.shp = "Australia.shp"
  eu.shp = "Europe.shp"
  na.shp = "NAmerica.shp"
  sa.shp = "SAmerica.shp"
  shps = c(wld.shp, af.shp, as.shp, au.shp, eu.shp, na.shp, sa.shp)

  #continent names
  conts = c("Global","Africa","Asia","Australia","Europe", "NAmer", "SAmer")
  
  #projections
  wld.proj = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"
  af.proj = "+proj=aea +lat_1=-32 +lat_2=12 +lat_0=-10 +lon_0=20 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
  as.proj = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=100 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
  au.proj = "+proj=aea +lat_1=-45 +lat_2=-15 +lat_0=-30 +lon_0=142 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
  eu.proj = "+proj=aea +lat_1=35 +lat_2=70 +lat_0=52.5 +lon_0=25 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
  na.proj = "+proj=aea +lat_1=20 +lat_2=70 +lat_0=45 +lon_0=-105 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
  sa.proj = "+proj=aea +lat_1=-40 +lat_2=0 +lat_0=-20 +lon_0=-65 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
  projs = c(wld.proj, af.proj, as.proj, au.proj, eu.proj, na.proj, sa.proj)

  #legend position
  wld.x = -15000000
  af.x = -4200000
  as.x = 3700000
  au.x = -4200000
  eu.x = -4200000
  na.x = -4600000
  sa.x = 1700000
  lpos.x = c(wld.x, af.x, as.x, au.x, eu.x, na.x, sa.x)
  wld.y = 900000
  af.y = 900000
  as.y = 2300000
  au.y = -100000
  eu.y = 3200000
  na.y = 0
  sa.y = -1000000
  lpos.y = c(wld.y, af.y, as.y, au.y, eu.y, na.y, sa.y)
   
  #legend titles
  if(substr(froot, 1, 1) == "H") {
    leg = c(expression(paste(delta^2, "H (\u2030)")), expression(paste(delta^2, "H 95% CI (\u2030)")))
  } else if(substr(froot, 1, 1) == "O") {
    leg = c(expression(paste(delta^18, "O (\u2030)")), expression(paste(delta^18, "O 95% CI (\u2030)")))
  } else {
    leg = c(rep("Isotope Value",2))
  }
  
  #plot dimensions
  w = c(7, 4.75, 5.5, 7, 5.5, 5, 4)
  h = c(5, 5, 5, 5, 5, 5, 5)
  
  #I/O file names
  froots = c(froot, paste0(substr(froot, 1, 3),"CI"))
  fnames = paste0(froots, ".asc")
  outroot = paste0(filedir, froot)
  mnames = paste0(filedir, conts, "/", froot, "_", conts, ".jpg")
  mnames = rbind(mnames,paste0(filedir, conts, "/", froot, "_", conts, "_CI.jpg"))
  
  for(i in 1:2){
    #read data
    grid = readAsciiGrid(fnames[i])
    plotvar = paste0("grid.sub$", fnames[i])
    
    #set projections
    grid=setproj(grid)
    
    for(j in 1:length(shps)){
      #get current polygon layer
      shp = readShapeSpatial(paste0(filedir, shps[j]))
      shp=setproj(shp)
      
      #make and clip raster
      rast = raster(grid)  #this works
      rast = crop(rast, extent(shp))
      rast = mask(rast, shp)
      grid.sub = as(rast, "SpatialGridDataFrame")
      rm(rast)
      
      #define classes and colorspace
      nclr = 9
      class = classIntervals(eval(parse(text=plotvar)), nclr, style = "equal", dataPrecision = 0.1)
      plotclr = rev(brewer.pal(nclr, "YlGnBu"))
      colcode = findColours(class, plotclr, digits=3)
      breaks = class$brks
      pal = attr(colcode, "palette")
      rm(class)
      rm(colcode)
      
      #convert to spatial polygons - THIS IS A MEMORY HOG!!!
      poly = Grid2Polygons(grid.sub, level=TRUE, at=breaks)
      rm(grid.sub)
      
      #transform projection
      poly.trans = spTransform(poly, CRS=CRS(projs[j]))
      shp.trans = spTransform(shp, CRS=CRS(projs[j]))
      rm(poly)
      
      #plot
      jpeg(mnames[i,j], width=w[j], height=h[j], units="in", pointsize=10, res=1200) #need to parameterize output size
      plot(poly.trans, border=rgb(0,0,0,max=255,alpha=20), col=plotclr)
      plot(shp.trans, lwd=0.5, add=TRUE)
      legent = paste(rev(round(breaks[1:nclr], digits=1)), "to", rev(round(breaks[1:nclr+1], digits=1)))
      legend(lpos.x[j], lpos.y[j], legend=legent, fill=rev(pal), box.col="white", 
             cex=0.75, title=leg[i])
      mtext("http://waterisotopes.org", side=1, cex=0.75, col="darkgrey")
      dev.off()
      rm(poly.trans)
      gc()
    }
    rm(grid)
    gc()
  }
}

setproj = function(spdf){
  proj4string(spdf) = CRS("+proj=longlat +ellps=WGS84")
  return(spdf)
}



#image(grid, "jandmap.asc", col=plotclr, axes=FALSE, breaks=class$brks)
#plot(countries, add=TRUE)
#legend("bottomleft", legend=names(rev(attr(colcode,"table"))), fill=rev(attr(colcode,"palette")), box.col="white", 
#       cex=0.8, title=expression(paste(delta^{2}, "H (\u2030)")))
