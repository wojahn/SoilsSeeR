#' SoilsSeeR Analysis Function
#' @param DSMW character: Path to FAO-UNESCO DSMW
#' @param ptz vector: list of lat/long points, must be Deg. Min. Sec.!
#' @export
# Written by John M. A. Wojahn October 2022
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn
# Provided under the Creative Commons Public License Attribution Non-Commercial
# Funded by Boise State University

SoilsInfoGrabbeR <- function(DSMW, ptz)
{
  oldw <- getOption("warn") #original opt
  options(warn = -1) # suppress warnings
  library("sp","rgdal") #load reqs
  message("Reading in FAO-UNESCO Digital Soils Map of the World")
  if(class(DSMW) != "character")
  {
    stop("Error: DSMW must be character!")
  }
  DSMW_sh <- rgdal::readOGR(DSMW, verbose = F) # read in DSMW
  message("Extracting and cleaning sample geospatial data")
  taxa <- ptz[,1] #IDs
  GeoCoords <- ptz[,2] #lat/long
  if(T %in% grepl("'",GeoCoords))
  {
    stop("Error: coordinates must be in Decimal, not Sexadecimal!")
  }
  # here we are converting the lat/long from DMS to DD
  GeoCoords_Expd <- as.data.frame(matrix(nrow=length(GeoCoords),ncol=2))
  for(i in 1:nrow(GeoCoords_Expd))
  {
    splitz <- as.character(unlist(strsplit(GeoCoords[i],split=" ")))
    latz <- splitz[1]

    #cd_lat <- char2dms(latz,chd="°",chm="'",chs="\"")
    #cd_lat <- as.numeric(cd_lat)
    longz <-splitz[2]
    #cd_long <- char2dms(longz,chd="°",chm="'",chs="\"")
    #cd_long <- as.numeric(cd_long)
    #if(grepl("S", latz))
    #{
    #  cd_long <- cd_long
    #}
    #if(grepl("W", latz))
    #{
    #  cd_long <- cd_lat
    #}
    #GeoCoords_Expd[i,1] <- cd_long
    #GeoCoords_Expd[i,2] <- cd_lat
    GeoCoords_Expd[i,1] <- longz
    GeoCoords_Expd[i,2] <- latz
  }
  colnames(GeoCoords_Expd) <- c("longitude","latitude")
  message("Creating spatial points object")
  require("sf")
  # Convert data frame to sf object
  # NA is OK because both use same projection
  # +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0
  # there is some strange error if I coerce DSMW_sh to projection
  # since it is empty in file, but visual inspection confirms
  # correctness
  pandan.sf.point <- st_as_sf(x = GeoCoords_Expd,
                          coords = c("longitude","latitude"),
                          crs = NA)
  # convert to sp object if needed
  message("Overlapping samples with FAO-UNESCO Digital Soils Map of the World")
  library(spatialEco)
  library(sp)
  pts.poly <- point.in.poly(pandan.sf.point, DSMW_sh, sp = T) # extract info
  pts_df <- as.data.frame(pts.poly) # make legible
  pts_df_old <- pts_df
  message("Checking answers for quality")
  pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                     max = nrow(pts_df), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     char = "=")   # Character used to create the bar
  for(i in 1:nrow(pts_df))
  {
    setTxtProgressBar(pb, i)
    if(is.na(pts_df[i,2]))
    {
      DSMW_f <- as(DSMW_sh,"sf")
      near <- sf::st_nearest_feature(pandan.sf.point[i,],DSMW_f,longlat=T)
      dist <- sf::st_distance(pandan.sf.point[i,],DSMW_f[near,],by_element = T)
      #assuming 110 k to degree on average...
      if(dist < 0.0091)
      {
        pts_df[i,2:ncol(pts_df)] <- as.data.frame(DSMW_f[near,])
      }else{
        pts_df[i,] <- "No_Soil_Data"
      }
    }
    close(pb)
  }
  taxa_pts_df <- as.data.frame(cbind(taxa,pts_df)) #bind to IDs
  taxa_pts_df <- taxa_pts_df[,-c(ncol(taxa_pts_df)-1,ncol(taxa_pts_df))]
  #message("Making interactive HTML map") BETA
  # interactive map:
  #library(mapview)
  #pandan.sf.point_pg <- sf::as_Spatial(pandan.sf.point)
  #mapOUT <- mapview(pandan.sf.point_pg)
  #mapshot(mapOUT,url = paste0(getwd(), "/IDs_map.html"))
  options(warn = oldw)
  return(taxa_pts_df)
}
