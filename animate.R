#####
#Run after global d2H and d18O grids have been created to generate animations of
#seasonal cycle. Also outputs a set of clipped and masked .asc grids containing
#land area values for each continent except Antarctica. These are saved in the
#corresponding subdirectories of /IsotopeMaps/ and could be reused for other
#plotting and distribution purposes.
#
#current code removes jpg images created in directory 'clean-up' step (lines 
#207-209) but these could be retained and used for distribution, supplanting
#wateriso_plot.R.
#####

#set WGS84 projection, used throughout
setproj = function(spdf){
  proj4string(spdf) = CRS("+proj=longlat +ellps=WGS84")
  return(spdf)
}

#####

library(sp)
library(rgdal)
library(maptools)
library(RColorBrewer)
library(classInt)
library(Grid2Polygons)
library(raster)

dbox = switch(Sys.info()["nodename"], "GJB-ZEN"="D:/Dropbox/", 
              "HYDROGEN"="C:/Users/gjbowen/Dropbox/")

#for H
froots=c("H01c", "H02c", "H03c", "H04c", "H05c", "H06c", "H07c", "H08c", "H09c", "H10c", "H11c", "H12c")
#for O
froots=c("O01c", "O02c", "O03c", "O04c", "O05c", "O06c", "O07c", "O08c", "O09c", "O10c", "O11c", "O12c")

#number of color intervals
nclr = 9

#folder names for intermediate grids and animations
folder = c("Global/", "Africa/", "Asia/", "Australia/", "Europe/", "NAmer/", "SAmer/")

#shapefile names
wld.shp = "continents.shp"
af.shp = "Africa.shp"
as.shp = "Asia.shp"
au.shp = "Australia.shp"
eu.shp = "Europe.shp"
na.shp = "NAmerica.shp"
sa.shp = "SAmerica.shp"
shps = c(wld.shp, af.shp, as.shp, au.shp, eu.shp, na.shp, sa.shp)

#projections
wld.proj = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"
af.proj = "+proj=aea +lat_1=-32 +lat_2=12 +lat_0=-10 +lon_0=20 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
as.proj = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=100 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
au.proj = "+proj=aea +lat_1=-45 +lat_2=-15 +lat_0=-30 +lon_0=142 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
eu.proj = "+proj=aea +lat_1=35 +lat_2=70 +lat_0=52.5 +lon_0=25 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
na.proj = "+proj=aea +lat_1=20 +lat_2=70 +lat_0=45 +lon_0=-105 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
sa.proj = "+proj=aea +lat_1=-40 +lat_2=0 +lat_0=-20 +lon_0=-65 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
projs = c(wld.proj, af.proj, as.proj, au.proj, eu.proj, na.proj, sa.proj)

#legend position x and y coords
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
if(substr(froots[1], 1, 1) == "H") {
  leg = expression(paste(delta^2, "H (\u2030)"))
} else if(substr(froots[1], 1, 1) == "O") {
  leg = expression(paste(delta^18, "O (\u2030)"))
} else {
  leg = "Isotope Value"
}

#plot titles
tits = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

#plot dimensions
w = c(7, 4.75, 5.5, 7, 5.5, 5, 4)
h = c(5, 5, 5, 5, 5, 5, 5)

for(j in 1:7){
  setwd(paste0(dbox, "Archived/Utilities/IsotopeMaps/"))
  
  #I/O file names
  fnames= paste0(froots, ".asc")
  filedir = paste0(dbox, "Archived/Utilities/wateriso_plots/")
  outdir = paste0(filedir, folder[j], "anim/")
  if(dir.exists(outdir)){} else{
    dir.create(outdir)
  }
  if(dir.exists(folder[j])){} else{
    dir.create(folder[j])
  }
  gnames = paste0(folder[j], froots, ".asc")
  mnames = paste0(outdir, froots, ".jpg")
  
  #placeholder values get replaced below
  low = 100   #minimum for colorscale
  high = -500   #maximum for colorscale
  
  #get current polygon layer
  shp = readShapeSpatial(paste0(filedir, shps[j]))
  shp=setproj(shp)

  for(i in 1:length(fnames)){
    #read data
    grid = readAsciiGrid(fnames[i])
    plotvar = paste0("grid.sub$", fnames[i])
    
    #set projections
    grid=setproj(grid)
    
    #make and clip raster
    rast = raster(grid)  #this works
    rast = crop(rast, extent(shp))
    rast = mask(rast, shp)
    grid.sub = as(rast, "SpatialGridDataFrame")
    rm(rast)
    grid.sub@grid@cellsize[1] = grid.sub@grid@cellsize[2] #force square cells if numeric fuzz
    write.asciigrid(grid.sub, gnames[i])
    low = min(low, min(grid.sub@data, na.rm=TRUE))
    high = max(high, max(grid.sub@data, na.rm=TRUE))
    rm(grid.sub)
    
    rm(grid)
    gc()
  }
  
  #output low and high values for storage
  write.table(c(low, high), file=paste0(folder[j], substr(froots[1], 1, 1), "minmax.txt"), row.names=FALSE, col.names=FALSE)
  
  #define classes and colorspace - run from here if re-generating plots only
  lowhigh = read.table(paste0(folder[j], substr(froots[1], 1, 1), "minmax.txt"))
  low = lowhigh[1,1]
  high = lowhigh[2,1]
  pal = rev(brewer.pal(nclr, "YlGnBu"))
  breaks = seq(low, high, (high-low)/nclr)
  
  for(i in 1:length(fnames)){
    #read in grid and project
    grid.sub = read.asciigrid(gnames[i])
    grid.sub = setproj(grid.sub)
    
    #find color breaks represented on map
    gmin = min(grid.sub@data, na.rm=TRUE)
    gmax = max(grid.sub@data, na.rm=TRUE)
    imin = (gmin-low)/(high-low)*nclr
    imax = (gmax-low)/(high-low)*nclr
    imin = max(1, ceiling(imin))
    imax = min(nclr, ceiling(imax))
        
    #convert to spatial polygons - THIS IS A MEMORY HOG!!!
    poly = Grid2Polygons(grid.sub, level=TRUE, at=breaks)
    rm(grid.sub)
    
    #transform projection
    poly.trans = spTransform(poly, CRS=CRS(projs[j]))
    shp.trans = spTransform(shp, CRS=CRS(projs[j]))
    rm(poly)
        
    #plot
    jpeg(mnames[i], width=w[j], height=h[j], units="in", pointsize=10, res=600) #need to parameterize output size
    plot(poly.trans, border=rgb(0,0,0,max=255,alpha=20), col=pal[imin:imax])
    plot(shp.trans, lwd=0.5, add=TRUE)
    legent = paste(rev(round(breaks[1:nclr], digits=1)), "to", rev(round(breaks[1:nclr+1], digits=1)))
    legend(lpos.x[j], lpos.y[j], legend=legent, fill=rev(pal), box.col="white", 
           cex=0.75, title=leg)
    if(j==1){ #global vs regional
      text(par('usr')[1]+1000000, par('usr')[4]-400000, tits[i])} else{
        text(par('usr')[1]+350000, par('usr')[4]-200000, tits[i])
    }  
    mtext("http://waterisotopes.org", side=1, cex=0.75, col="darkgrey")
    dev.off()
    rm(poly.trans)
    gc()
  }
  
  setwd(outdir)
  magickroot = "\"C:/Program Files/ImageMagick/magick.exe\""
  
  for(i in 1:12){
    #prep images for animation
    cmd = paste0(magickroot, " ", froots[i], ".jpg -resize 1500x1500 -shave 100x100 +repage ", froots[i], ".gif")
    system(cmd)
  }
  #animate
  animstring = paste0(magickroot, " -delay 50 *c.gif ", substring(folder[j], 1, nchar(folder[j])-1), substr(froots[1], 1, 1), ".gif")
  system(animstring) 

  #clean up
  for(i in 1:12){
    file.remove(paste0(froots[i], ".jpg"))
    file.remove(paste0(froots[i], ".gif"))
  }
}
