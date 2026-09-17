map = function(froot){

  library(terra)
  library(RColorBrewer)
  library(classInt)
    
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
    grid = rast(fnames[i])
    
    #set projections
    grid=setproj(grid)
    
    for(j in 1:length(shps)){
      #get current polygon layer
      shp = vect(paste0(filedir, shps[j]))
      shp=setproj(shp)
      
      #clip raster
      rast.sub = crop(grid, ext(shp))
      rast.sub = mask(rast.sub, shp)
      
      #define classes and colorspace
      nclr = 9
      class = classIntervals(as.vector(values(rast.sub, na.rm=TRUE)), nclr, style = "equal", dataPrecision = 0.1)
      plotclr = rev(brewer.pal(nclr, "YlGnBu"))
      breaks = class$brks
      rm(class)
      
      #classify raster into color bins, then convert to smoothed contour polygons
      rcl = cbind(breaks[1:nclr], breaks[2:(nclr+1)], 1:nclr)
      classed = classify(rast.sub, rcl, include.lowest=TRUE)
      names(classed) = "class"
      poly = as.polygons(classed, dissolve=TRUE)
      rm(rast.sub, classed)
      
      #transform projection
      poly.trans = project(poly, projs[j])
      shp.trans = project(shp, projs[j])
      rm(poly)
      
      #map each polygon's class id back to its color
      polyclr = plotclr[poly.trans$class]
      
      #plot
      jpeg(mnames[i,j], width=w[j], height=h[j], units="in", pointsize=10, res=1200) #need to parameterize output size
      plot(poly.trans, border=rgb(0,0,0,max=255,alpha=20), col=polyclr)
      plot(shp.trans, lwd=0.5, add=TRUE)
      legent = paste(rev(round(breaks[1:nclr], digits=1)), "to", rev(round(breaks[1:nclr+1], digits=1)))
      legend(lpos.x[j], lpos.y[j], legend=legent, fill=rev(plotclr), box.col="white", 
             cex=0.75, title=leg[i])
      dev.off()
      rm(poly.trans)
      gc()
    }
    rm(grid)
    gc()
  }
}

setproj = function(x){
  crs(x) = "+proj=longlat +ellps=WGS84"
  return(x)
}



#image(grid, "jandmap.asc", col=plotclr, axes=FALSE, breaks=class$brks)
#plot(countries, add=TRUE)
#legend("bottomleft", legend=names(rev(attr(colcode,"table"))), fill=rev(attr(colcode,"palette")), box.col="white", 
#       cex=0.8, title=expression(paste(delta^{2}, "H (\u2030)")))
