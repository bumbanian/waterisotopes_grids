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
  na.shp = "NAmerica.shp"
  shps = na.shp

  #projections
  na.proj = "+proj=aea +lat_1=20 +lat_2=70 +lat_0=45 +lon_0=-105 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
  projs = na.proj

  #legend position
  na.x = -4600000
  lpos.x = na.x
  na.y = 0
  lpos.y = na.y
   
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
  froots=c(froot, paste0(substr(froot, 1, 3),"CI"))
  fnames= paste0(froots, ".asc")
  outroot = paste0(filedir, froot)
  mnames= c(paste0(outroot, "_NAmer.jpg"))
  mnames= rbind(mnames,c(paste0(outroot, "_NAmer_CI.jpg")))
  
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
