# Function to build meaningful and minimalistic population like maps
# www.overfitting.net
# https://www.overfitting.net/


populatiomap=function(map, shape='circle', mapstyle='filled', grid='none',
                      inwidth=100, outwidth=100, overlap=1, gamma=1) {
    # shape='circle', 'square', 'none'
    # mapstyle='filled', 'contour', 'none'
    # grid='none', 'centre', 'wrap'
    # inwidth: input grid size in pixels
    # outwidth: output grid size in pixels (inwidth=outwidth avoids resampling)
    # overlap: how much a square or circle can overlap its neighbours
    # gamma: output gamma lift curve
    
    require(raster)  # resample
    require(tiff)  # save 16-bit TIFF's
    require(png)  # save 8-bit PNG's

    FILLEDVALUE=0.5  # background map and grid colours
    GRIDVALUE=0.25

    # Extend matrix to fit in output size (integer number of cells)
    DIMY=nrow(map)
    DIMX=ncol(map)
    NGRIDY=ceiling(DIMY/inwidth)  # wrap input map to preserve all data
    NGRIDX=ceiling(DIMX/inwidth)
    maptmp=matrix(NA, nrow=NGRIDY*inwidth, ncol=NGRIDX*inwidth)
    maptmp[1:DIMY, 1:DIMX]=map
    map=maptmp
    rm(maptmp)
    
    
    # Calculate map filled (0/1 filled)
    filled=map
    filled[!is.na(filled)]=1  # set land areas to 1
    filled[is.na(filled)]=0  # set sea areas from NA to 0
    # Resample filled map only if needed
    if (outwidth!=inwidth) {
        print("Warning: input and output grid sizes differ, resampling map...")
        # raster(filled) is a raster created with extent: 0, 1, 0, 1
        filledrs=raster(nrow=NGRIDY*outwidth, ncol=NGRIDX*outwidth, 
                        xmn=0, xmx=1, ymn=0, ymx=1)  # extent: 0, 1, 0, 1
        filled=resample(raster(filled), filledrs,
                        method='bilinear')  # bilinear is faster than ngb!
        rm(filledrs)
        filled=as.matrix(filled)
        filled[filled>=0.5]=1
        filled[filled<0.5]=0
        DIMY=nrow(filled)  # output new working dimensions
        DIMX=ncol(filled) 
    }
    
    
    # Process matrix
    map[is.na(map)]=0  # set sea areas from NA to 0
    # population=sum(map)  # 458 million people in EU
    
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
    
    
    # Draw output map shapes
    mapout=matrix(0, nrow=NGRIDY*outwidth, ncol=NGRIDX*outwidth)
    MAXD=outwidth*overlap  # max diameter of circles/squares into cells
    MAXPOP=max(mapavg)  # people on most populated cell
    # shape='circle', 'square', 'none'
    for (i in 1:NGRIDX) {
        x0=outwidth*(i-1)+outwidth/2  # decimal figure
        for (j in 1:NGRIDY) {
            R=(mapavg[j,i]/MAXPOP)^0.5 * MAXD/2  # square area proportional to population
            if (R) {
                y0=outwidth*(j-1)+outwidth/2  # decimal figure
                if (shape=='circle') {
                    for (x in round(x0-R):round(x0+R)) {
                        for (y in round(y0-R):round(y0+R)) {
                            if ( ((x-x0)^2 + (y-y0)^2 )^0.5 < R) mapout[y,x]=mapout[y,x]+1
                        }
                    }
                } else if (shape=='square') {
                    mapout[round(y0-R):round(y0+R), round(x0-R):round(x0+R)]=1
                }
            }
        }
    }
    if (shape=='square' & overlap>1)
        print("Warning: 'square' plot with overlap>1 not recommended")
    mapplot=mapout  # preserve only shapes (needed for grid)
    
    
    # Draw contour/filled map
    # mapstyle='filled', 'contour', 'none'
    if (mapstyle=='filled') {
        writePNG(filled, "mapfilled.png")
        
        indices=which(filled==1 & mapout==0)  # plot filled
        mapout[indices]=FILLEDVALUE
    } else if (mapstyle=='contour') {
        # Calculate map contour (0/1 contour) from filled
        contour=filled*0
        # 1 pixel thickness contour
        contour[2:(DIMY-1), 2:(DIMX-1)]=
            abs(filled[1:(DIMY-2), 2:(DIMX-1)] -
                filled[2:(DIMY-1), 2:(DIMX-1)]) +
            abs(filled[2:(DIMY-1), 1:(DIMX-2)] -
                filled[2:(DIMY-1), 2:(DIMX-1)])
        # increase to 3 pixel thickness contour
        contour[2:(DIMY-1), 2:(DIMX-1)]=contour[2:(DIMY-1), 2:(DIMX-1)]+
            contour[1:(DIMY-2), 2:(DIMX-1)]+contour[3:(DIMY-0), 2:(DIMX-1)]+
            contour[2:(DIMY-1), 1:(DIMX-2)]+contour[2:(DIMY-1), 3:(DIMX-0)]
        contour[contour!=0]=1
        writePNG(contour, "mapcontour.png")
        
        indices=which(contour==1 & mapout==0)  # plot contour
        indices2=which(contour==1 & mapout!=0)  # invert contour on overlapping areas
        mapout[indices]=FILLEDVALUE
        mapout[indices2]=1-mapout[indices2]
    }
    
    
    # Draw grid
    # grid='none', 'centre', 'wrap'
    if (grid=='centre') {
        indices=which(  # plot grid
            ( (!(row(mapout)+round(outwidth/2))%%outwidth & col(mapout)%%2)
            | (!(col(mapout)+round(outwidth/2))%%outwidth & row(mapout)%%2))
            & filled==1 & mapplot==0)
        mapout[indices]=GRIDVALUE   
    } else if (grid=='wrap') {
        indices=which(  # plot grid
            ( (!row(mapout)%%outwidth & col(mapout)%%2)
            | (!col(mapout)%%outwidth & row(mapout)%%2))
            & filled==1 & mapplot==0)
        mapout[indices]=GRIDVALUE
    }
    
    return ((mapout/max(mapout))^(1/gamma))
}


# Usage example. Map is built monochrome, colours must be added afterwards

library(raster)  # https://cran.r-project.org/web/packages/raster/raster.pdf

# https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/population-distribution-demography/geostat
europe=raster("ESTAT_OBS-VALUE-T_2021_V1-0.tiff")  # read GeoTIFF file
europe

# Output map
europe=populatiomap(as.matrix(europe))
writeTIFF(europe, "europeout.tif", compression='LZW', bits.per.sample=16)
