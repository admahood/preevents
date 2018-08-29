# this is to make figure 2 of the preevents proposal

library(sf)
library(raster)
library(ncdf4) # data is in netcdfs
library(rgdal) # target of the function is a point shapefile

# pdsi -------------------------------------------------------------------------
years = 1984:2017
dir.create("~/data/gridmet")
for(y in years){
  destfile <- path.expand(paste0("~/data/gridmet/pdsi_",y,".nc"))
  sourcefile <- paste0("http://www.northwestknowledge.net/metdata/data/pdsi_",y,".nc")
  if(!file.exists(destfile)){
    download.file(sourcefile, destfile)}
}


# fire -------------------------------------------------------------------------

mtbs_file <- "~/data/fire/mtbs/mtbs_perimeter_data/mtbs_perims_DD.shp"
destfile = "~/Downloads/mtbs.zip"
if(!file.exists(mtbs_file)){
  download.file("https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/burned_area_extent_shapefile/mtbs_perimeter_data.zip",
                destfile =destfile)
  unzip(destfile, exdir = "~/data/fire/mtbs/")
}

mtbs <- st_read(mtbs_file)

# earthquakes ------------------------------------------------------------------



# landslides -------------------------------------------------------------------