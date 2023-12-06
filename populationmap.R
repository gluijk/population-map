# Function to build meaningful and minimalistic population like maps
# www.overfitting.net
# https://www.overfitting.net/


populatiomap=function(map,
          shape='circle', shapestyle='solid', mapstyle='solid', grid='none',
          inwidth=100, outwidth=100, overlap=1, gamma=1) {
    # map must be a matrix where...
    #   <0 values will be previosly clipped to 0
    #   >=0 values will be summed
    #   NA values will define the geographical limits of the map
    # shape='circle', 'square', 'none'
    # shapestyle='solid', 'outline', 'none'
    # mapstyle='solid', 'outline', 'none'
    # grid='none', 'centre', 'wrap'
    # inwidth: input grid size in pixels
    # outwidth: output grid size in pixels (inwidth=outwidth avoids resampling)
    # overlap: how much a square/circle can overlap its neighbours
    # gamma: output gamma lift curve
    
    require(raster)  # resample()
    require(tiff)  # save 16-bit TIFF's
    require(png)  # save 8-bit PNG's
    
    SOLIDVALUE=0.5  # background map and grid colours
    GRIDVALUE=0.25
    STROKE=3  # width of circles/squares
   
    # Clip negative values to 0 (will become part of the map)
    negatives=which(map<0 & !is.na(map))
    if (length(negatives)) {
        print(paste0("Warning: ",length(negatives)," negative values clipped to 0"))
        map[negatives]=0      
    }
    
    # Extend matrix to fit in output size (integer number of cells)
    DIMY=nrow(map)
    DIMX=ncol(map)
    NGRIDY=ceiling(DIMY/inwidth)  # wrap input map to preserve all data
    NGRIDX=ceiling(DIMX/inwidth)
    maptmp=matrix(NA, nrow=NGRIDY*inwidth, ncol=NGRIDX*inwidth)
    maptmp[1:DIMY, 1:DIMX]=map
    map=maptmp
    rm(maptmp)

    # Calculate solid map (0/1)
    solid=map
    solid[!is.na(solid)]=1  # set land areas to 1
    solid[is.na(solid)]=0  # set sea areas from NA to 0
    # Resample solid map only if needed
    if (outwidth!=inwidth &
        ((mapstyle=='solid'|mapstyle=='outline') | (grid=='centre'|grid=='wrap')) ) {
        print("Input and output grid sizes differ, resampling map...")
        # raster(solid) is a raster created with extent: 0, 1, 0, 1
        solidrs=raster(nrow=NGRIDY*outwidth, ncol=NGRIDX*outwidth, 
                        xmn=0, xmx=1, ymn=0, ymx=1)  # extent: 0, 1, 0, 1
        solid=resample(raster(solid), solidrs,
                        method='bilinear')  # bilinear is faster than ngb!
        rm(solidrs)
        solid=as.matrix(solid)
        solid[solid>=0.5]=1
        solid[solid<0.5]=0
        DIMY=nrow(solid)  # output new working dimensions
        DIMX=ncol(solid)
    }
    
    # Set NA values to 0
    map[is.na(map)]=0

    # Reduced map which will provide all output discrete values
    mapavg=matrix(0, nrow=NGRIDY, ncol=NGRIDX)
    for (i in 1:NGRIDX) {
        for (j in 1:NGRIDY) {
            mapavg[j,i]=sum(map[((j-1)*inwidth+1):(j*inwidth),
                                ((i-1)*inwidth+1):(i*inwidth)])
        }
    }
    
    # Write calculated maps
    writeTIFF((map/max(map))^(1/gamma), "map.tif",
              compression='LZW', bits.per.sample=16)
    writeTIFF((mapavg/max(mapavg))^(1/gamma), "mapavg.tif",
              compression='LZW', bits.per.sample=16)
    rm(map)
    mapout=matrix(0, nrow=NGRIDY*outwidth, ncol=NGRIDX*outwidth)
    
    # Draw output map shapes
    # shape='circle', 'square', 'none'
    # shapestyle='solid', 'outline', 'none'
    if ((shape=='circle'|shape=='square') & (shapestyle=='solid'|shapestyle=='outline')) {
        MAXD=outwidth*overlap  # max diameter of circles/squares into cells
        MAXPOP=max(mapavg)  # people on most populated cell
        THICK=round(STROKE/2)-1  # variable thickness not yet implemented
        for (i in 1:NGRIDX) {
            x0=outwidth*(i-1)+outwidth/2  # decimal figure
            for (j in 1:NGRIDY) {
                R=(mapavg[j,i]/MAXPOP)^0.5 * MAXD/2  # square area proportional to population
                if (R) {
                    y0=outwidth*(j-1)+outwidth/2  # decimal figure
                    if (shape=='circle' & shapestyle=='solid') {
                        for (x in round(x0-R):round(x0+R)) {
                            for (y in round(y0-R):round(y0+R)) {
                                if ( ((x-x0)^2 + (y-y0)^2 )^0.5 < R) mapout[y,x]=mapout[y,x]+1  # differentiated overlap
                            }
                        }
                    } else if (shape=='circle' & shapestyle=='outline') {
                        for (x in round(x0-R-1):round(x0+R+1)) {
                            for (y in round(y0-R-1):round(y0+R+1)) {
                                if ( ((x-x0)^2 + (y-y0)^2 )^0.5 < R+1 &
                                     ((x-x0)^2 + (y-y0)^2 )^0.5 > R-1) mapout[y,x]=1
                            }
                        }
                    } else if (shape=='square' & shapestyle=='solid') {
                        mapout[round(y0-R):round(y0+R), round(x0-R):round(x0+R)]=1
                    } else if (shape=='square' & shapestyle=='outline') {
                        mapout[round(y0-R-1):round(y0+R+1), round(x0-R-1):round(x0-R+1)]=1
                        mapout[round(y0-R-1):round(y0+R+1), round(x0+R-1):round(x0+R+1)]=1
                        mapout[round(y0-R-1):round(y0-R+1), round(x0-R-1):round(x0+R+1)]=1
                        mapout[round(y0+R-1):round(y0+R+1), round(x0-R-1):round(x0+R+1)]=1
                    }
                }
            }
        }
    }
    if (shape=='square' & overlap>1) print("Warning: 'square' plot with overlap>1 not recommended")
    mapplot=mapout  # preserve only shapes (needed for grid)
    
    # Draw solid/outline map
    # mapstyle='solid', 'outline', 'none'
    if (mapstyle=='solid') {
        writePNG(solid, "mapsolid.png")
        
        indices=which(solid==1 & mapout==0)  # plot solid
        mapout[indices]=SOLIDVALUE
    } else if (mapstyle=='outline') {
        # Calculate outline map (0/1) from solid map
        outline=solid*0
        # 1 pixel thickness outline
        outline[2:(DIMY-1), 2:(DIMX-1)]=
            abs(solid[1:(DIMY-2), 2:(DIMX-1)] -
                solid[2:(DIMY-1), 2:(DIMX-1)]) +
            abs(solid[2:(DIMY-1), 1:(DIMX-2)] -
                solid[2:(DIMY-1), 2:(DIMX-1)])
        # increase to 3 pixel thickness outline
        outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
            outline[1:(DIMY-2), 2:(DIMX-1)]+outline[3:(DIMY-0), 2:(DIMX-1)]+
            outline[2:(DIMY-1), 1:(DIMX-2)]+outline[2:(DIMY-1), 3:(DIMX-0)]
        outline[outline!=0]=1
        writePNG(outline, "mapoutline.png")
        
        indices=which(outline==1 & mapout==0)  # plot outline
        indices2=which(outline==1 & mapout!=0)  # invert outline on overlapping areas
        mapout[indices]=SOLIDVALUE
        mapout[indices2]=1-mapout[indices2]
    }
    
    # Draw grid
    # grid='none', 'centre', 'wrap'
    if (grid=='centre') {
        indices=which(  # plot grid
            ( (!(row(mapout)+round(outwidth/2))%%outwidth & col(mapout)%%2)
            | (!(col(mapout)+round(outwidth/2))%%outwidth & row(mapout)%%2))
            & solid==1 & mapplot==0)
        mapout[indices]=GRIDVALUE   
    } else if (grid=='wrap') {
        indices=which(  # plot grid
            ( (!row(mapout)%%outwidth & col(mapout)%%2)
            | (!col(mapout)%%outwidth & row(mapout)%%2))
            & solid==1 & mapplot==0)
        mapout[indices]=GRIDVALUE
    }
    
    return ((mapout/max(mapout))^(1/gamma))
}



# Usage example. Output map is non antialiased grayscale, colours must be added afterwards

library(raster)  # https://cran.r-project.org/web/packages/raster/raster.pdf

# https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/population-distribution-demography/geostat
europe=raster("ESTAT_OBS-VALUE-T_2021_V1-0.tiff")  # read GeoTIFF file
europe=populatiomap(as.matrix(europe), inwidth=100, outwidth=100,
                    shape='circle', shapestyle='solid',
                    mapstyle='solid', grid='none', overlap=1.2)
writeTIFF(europe, "europe.tif", compression='LZW', bits.per.sample=16)

pal=colorRampPalette(c(rgb(0,0,0.1), rgb(0.5,.5,0.5), rgb(1,1,0.7)))
image(t(europe[nrow(europe):1,]), col=pal(3), useRaster=TRUE,
      asp=nrow(europe)/ncol(europe), axes=FALSE)



# https://efi.int/knowledge/maps/woodproduction
wood=raster("woodprod_average.tif")  # read GeoTIFF file
wood=populatiomap(as.matrix(wood), inwidth=60, outwidth=60,
                    shape='circle', shapestyle='solid', gamma=2.2,
                    mapstyle='solid', grid='none', overlap=1)
writeTIFF(wood, "wood.tif", compression='LZW', bits.per.sample=16)

