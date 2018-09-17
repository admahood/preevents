# this is to make figure 2 of the preevents proposal

library(sf)
library(raster)
library(ncdf4) # data is in netcdfs
library(rgdal) # target of the function is a point shapefile
library(dplyr)
library(foreach)
library(doParallel)
library(broom)
library(ggplot2) # needs to be from install_github("hadley/ggplot2")

# huc --------------------------------------------------------------------------
huc_file <- "~/data/background/hydrologic_units/huc250k_shp/huc250k.shp"
destfile <- "~/Downloads/huc.zip"

if(!file.exists(huc_file)){
  download.file("https://water.usgs.gov/GIS/dsdl/huc250k_shp.zip",
                destfile =destfile)
  unzip(destfile, exdir = "~/data/background/hydrologic_units")
  unlink(destfile)
}

huc <- st_read(huc_file)
huc1 <-readOGR(huc_file)

# pdsi -------------------------------------------------------------------------
years <- 1984:2016 # something'is messed up with 2017...its projection is upside down and backwards

# downloading
rawpath <- "~/data/gridmet/pdsi"
dir.create(rawpath, recursive = TRUE)
for(y in years){
  destfile <- path.expand(paste0(rawpath,"/pdsi_",y,".nc"))
  sourcefile <- paste0("http://www.northwestknowledge.net/metdata/data/pdsi_",y,".nc")
  if(!file.exists(destfile)){
    download.file(sourcefile, destfile)}
}

# calculating mean PDSI and applying a threshold of -2
respath<-"~/data/gridmet/pdsi_mean"
dir.create(respath)
threshold <- -2
rclmat <- matrix(c(-99, threshold, 1, threshold, 99, 0), nrow = 2, ncol=3, byrow=T)

c <- detectCores()-1
registerDoParallel(c)
foreach(y = years)%dopar%{
  stack(file.path(rawpath, paste0("pdsi_", y, ".nc"))) %>%
    calc(mean) %>%
    projectRaster(crs=crs(huc1)) %>%
    reclassify(rclmat, filename = paste0(respath, "/mean_pdsi_u-", threshold,"_", y, ".tif"))
  system(paste("echo",y))
}

rclmat10 <- matrix(c(-99, 9, 0, 9.1, 999, 1), nrow = 2, ncol=3, byrow=T)
drought_10yrs <- list.files(respath, full.names = TRUE) %>% # lapply(raster) #%>% lapply(extent)
  stack() %>%
  calc(sum) %>%
  reclassify(rclmat10) %>%
  aggregate(fun=max, fact=12) #%>% plot()
  
drought_10yrs[is.na(drought_10yrs[])] <- 0 


# fire -------------------------------------------------------------------------
mtbs_file <- "~/data/fire/mtbs/points/mtbs_fod_pts_DD.shp"
destfile <- "~/Downloads/mtbs.zip"
if(!file.exists(mtbs_file)){
  download.file("https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/fod_pt_shapefile/mtbs_fod_pts_data.zip",
                destfile =destfile)
  unzip(destfile, exdir = "~/data/fire/mtbs/points")
  unlink(destfile)
}
mtbs <- readOGR(mtbs_file) %>%
  spTransform(crs(drought_10yrs))
  
mtbs_10 <- rasterize(x=mtbs, y=drought_10yrs, fun= 'count', field=1) %>%
  reclassify(rclmat10) %>%
  calc(fun = function(x)x*10)#%>% plot(mtbs_10)
mtbs_10[is.na(mtbs_10[])] <- 0 


# floods -----------------------------------------------------------------------

flood_file <- "~/data/sheldus/flood_yi/SHELDUS/National_data/us_counties.shp"
floods10 <- st_read(flood_file)%>%
  st_transform(crs=st_crs(huc)) %>%
  mutate(flood10= ifelse(Floo_Incid>10, 100,0)) %>%
  dplyr::filter(flood10 > 0) %>%
  dplyr::select(flood10) %>%
  as("Spatial") %>%
  rasterize(x=., y=drought_10yrs, field=1, fun="max")%>%
  calc(fun = function(x)x*100) #%>% plot(floods10)

floods10[is.na(floods10[])] <- 0 

# states -----------------------------------------------------------------------
state_file <- "~/data/background/states/tl_2017_us_state.shp"
destfile <- "~/Downloads/states.zip"
if(!file.exists(state_file)){
  download.file("ftp://ftp2.census.gov/geo/tiger/TIGER2017/STATE/tl_2017_us_state.zip",
                destfile =destfile)
  unzip(destfile, exdir = "~/data/background/states")
  unlink(destfile)
}
states <- st_read(state_file) %>%
  dplyr::filter(REGION != 9,
                STUSPS != "AK",
                STUSPS != "HI") %>%
  st_transform(crs = crs(drought_10yrs, asText=TRUE))

# prepping overlay -------------------------------------------------------------

int <- stack(mtbs_10, floods10, drought_10yrs) %>%
  calc(sum) %>%
  mask(huc) %>%
  asFactor()

levs <- levels(int)[[1]]
levs$interactions = c("None", "Single", "Single", "Drought & Fire", "Single", 
                      "Flood & Fire", "Flood & Drought","Drought, Flood & Fire")
levs$number = c(0, 1, 1, 2, 1,2, 2,3)

colnames(levs) = c("ID", "interactions", "number")
levels(int) = levs
int_1 = deratify(int, 'interactions')
int_2 <- deratify(int, 'number')
int_2[is.na(int_2[])] <- 0 

watersheds <- raster::extract(int_2, huc, fun=max) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::rename(Watersheds_Affected = ".")
watersheds$num <- c("none", "Single", "Double", "Triple")
watersheds$Watersheds_Affected <- NULL

int_tab <- freq(int_1)[1:6,] %>% 
  as_tibble() %>%
  mutate(percent = count/sum(count)*100) %>%
  left_join(levels(int_1)[[1]], by = c("value" = "ID")) %>%
  dplyr::select(interactions,percent)

intp = rasterToPoints(int_1) %>%
  as_tibble() %>%
  rename(ID = interactions) %>%
  left_join(levels(int_1)[[1]]) %>%
  filter(interactions != "None") %>%
  left_join(int_tab) %>%
  mutate(interactions = as.factor(paste0(interactions, " ", round(percent), "%")))

intp$interactions <- factor(intp$interactions,levels(intp$interactions)[c(5,1,3,4,2)])

# plot -------------------------------------------------------------------------

wlabel = paste0("Watersheds Affected\nSingle: ",
                watersheds[watersheds$num=="Single",]$Freq,
                "\nDouble: ",watersheds[watersheds$num=="Double",]$Freq,
                "\nTriple: ",watersheds[watersheds$num=="Triple",]$Freq)
my_colors = c("lightgrey", "turquoise4", "orange", "khaki", "firebrick")

p <- ggplot() +
        geom_sf(data=huc, alpha = 0.5, lwd=0.25, color = "turquoise3", fill="white") +
        geom_sf(data=states, color = "grey40", fill="transparent") +
        geom_raster(data = intp, aes(x=x, y=y,fill = interactions), alpha=0.85) +
        scale_fill_manual(values = my_colors) +
        theme_void() +
        guides(fill=guide_legend(title=NULL)) +
        theme(legend.justification=c(0,0), legend.position=c(0,0),
              panel.grid.major = element_line(colour = 'transparent'),
              plot.title = element_text(hjust = 0.5)) +
        annotate("text", x=800000, y=380000,label = wlabel, size=3) +
        annotate("text", x=100000, y=3200000, label = "Disturbance Co-occurrence", size=8)

ggsave("co-occurrence_aea.png",limitsize = FALSE)

