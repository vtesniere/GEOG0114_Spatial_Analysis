
# in case there is a problem when running the code, we have attached the final version of the dataset for inner London
# this was used to create most of the maps and run the models

# load('GEOG0114 - code and datasets/Assignement/full_transformed_data.RData')
# full_3 <- full_transformed_data
# loading cultural venues
# load('GEOG0114 - code and datasets/Assignement/cultural_venues_inner_london.RData')


# packages to install/load

library(tidyverse)
library(sf)
library(tmap)
library(osmdata)
library(geodist)
library(dodgr)
library(spdep)
library(spatialreg)
library(rgeos)
library(terra)
library(raster)
library(dplyr)
library(Rfast)
library(jtools)
library(tidyverse)
library(timechange)
library(lubridate)
library(chron)
library(spgwr)

#inner london boroughs
inner_london <- c("Camden", "Greenwich", "Hackney", "Hammersmith and Fulham", "Islington", "Kensington and Chelsea",
                  "Lambeth", "Lewisham", "Southwark", "Tower Hamlets", "Wandsworth","Westminster", "City of London",
                  "Wandsworth and Westminster", "City of Westminster")

# inputs shapefiles and cleans names so they fit with the traditional borough names
London_lsoa <- read_sf('Assignement/Datasets 2/London LSOA Areas.shp')
London_lsoa$Cleaned_Names = substr(London_lsoa$LSOANAME,1,nchar(London_lsoa$LSOANAME)-5)

# adding borough shapefile to make the maps look better

boroughs <- read_sf('Assignement/Datasets 2/London Borough Areas.shp')
BOROUGHshp <- st_transform((boroughs%>% filter(BOROUGHN %in% inner_london)),4326)

#wards data
london_wards <- st_read('Assignement/london_wards.shp')
WARDSshp <- st_transform((london_wards%>% filter(DISTRICT %in% inner_london)),4326)

# we now proceed to get coordinates and extract road network from london
l_bbox <- c(-0.5103751, 51.2867602, 0.3340155, 51.6918741)
# as this bbox is too large and too much data to pull from OSM servers, we had to seperate into 6 smaller areas

l1_bbox_new <- c(-0.5103751, 51.2867602, 0.3340155, 51.4000000)
l2_bbox_new <- c(-0.5103751, 51.4000000, 0.3340155, 51.4800000)
l3_bbox_new <- c(-0.5103751, 51.4800000, 0.3340155, 51.5500000)
l4_bbox_new <- c(-0.5103751, 51.5500000, 0.3340155, 51.5900000)
l5_bbox_new <- c(-0.5103751, 51.5900000, 0.3340155, 51.6300000)
l6_bbox_new <- c(-0.5103751, 51.6300000, 0.3340155, 51.6918741)

# the below line of code was run for each bbox_new object, creating a new osmdata each time
osmdata <- opq(bbox = l6_bbox_new) %>%
  add_osm_feature(key = "highway", value = c("primary", "secondary", "tertiary",
                                             "residential", "path", "footway", "unclassified", "living_street", "pedestrian")) %>%
  osmdata_sf()

# extract the points, with their osm_id.
l1_roads_nodes <- osmdata$osm_points[, "osm_id"]

# extract the lines, with their osm_id, name, type of highway, max speed and
# oneway attributes
l1_roads_edges <- osmdata$osm_lines[, c("osm_id", "name", "highway")]

# this was repeated for all 6 boxes, making the final eadge a road datafile as follows

# putting all the london dataset together
london_roads_edges <- rbind(l1_roads_edges, l2_roads_edges)
london_roads_nodes <- rbind(l1_roads_nodes, l2_roads_nodes)
london_roads_edges <- rbind(london_roads_edges, l3_roads_edges)
london_roads_nodes <- rbind(london_roads_nodes, l3_roads_nodes)
london_roads_edges <- rbind(london_roads_edges, l4_roads_edges)
london_roads_nodes <- rbind(london_roads_nodes, l4_roads_nodes)
london_roads_edges <- rbind(london_roads_edges, l5_roads_edges)
london_roads_nodes <- rbind(london_roads_nodes, l5_roads_nodes)
london_roads_edges <- rbind(london_roads_edges, l6_roads_edges)
london_roads_nodes <- rbind(london_roads_nodes, l6_roads_nodes)

# such a big file that we saved it locally in case the data was lost 
# and/or if the data could not be pulled again if the servers were not responding
load('GEOG0114 - code and datasets/Assignement/london_edges.RData')
load('GEOG0114 - code and datasets/Assignement/london_nodes.RData')

# loading curltural amenities and putting them togather with their polygons

museums <- read.csv('Assignement/Museums_and_public_galleries.csv')
libraries <- read.csv('Assignement/Libraries.csv')
theatres <- read.csv('Assignement/Theatres.csv')
GMV <- read.csv('Assignement/Music_venues_grassroots.csv')

cultural_venues <- rbind(theatres, GMV)
cultural_venues <- rbind(cultural_venues, libraries)
cultural_venues <- rbind(cultural_venues, museums)

# checking and removing possible duplicates

n_occur <- data.frame(table(cultural_venues$name))
n_occur[n_occur$Freq > 1,]
cultural_venues <- unique(cultural_venues[ , 1:15 ] )
# remove line 200, 40 and 7
cultural_venues <- cultural_venues %>%  filter(!row_number() %in% c(7, 40, 200))

# creating geometry points and assigning the British grid coordinates
cultural_venues$polygon <- cultural_venues %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POINT")
cultural_venues$polygon <- st_transform(cultural_venues$polygon, 27700)


London_lsoa$centroid <- st_centroid(London_lsoa)

# finding the biggest and smallest lsoas in area size
London_lsoa$area <- st_area(London_lsoa) #Take care of units

hist(London_lsoa$area_sq_m)
London_lsoa <- rename(London_lsoa, area_sq_m=area)
options(scipen = 999) # this is to remove scientific notation, to put back to default make 0
clean_units <- function(z){
  attr(z,"units") <- NULL
  class(z) <- setdiff(class(z),"units")
  z
}
London_lsoa$area <- clean_units(London_lsoa$area)
mean(London_lsoa$area)

# adding the weights for walking
graph <- weight_streetnet(london_roads_edges, wt_profile = "foot")
# extract the largest connected graph component
graph_connected <- graph[graph$component == 1, ]

# create a distance matrix between area centroids and cultural venues
cent_to_cultural_calc <- dodgr_distances(graph_connected, from = st_coordinates(London_lsoa$geometry),
                                         to = st_coordinates(cultural_venues$polygon), shortest = TRUE, pairwise = FALSE, quiet = FALSE)

# getting list of lsoa's with no direct access to a museum
# getting the NA's is for all
no_direct_access <- prout %>% filter(is.na(rowSums(cent_to_cultural_calc <= 400)) %>% select(LSOANAME)
no_direct_access <- unlist(no_direct_access$LSOANAME)
length(no_direct_access) # 1094 lsoa's in greater london without any direct access

# mapping these out
tm_shape(London_lsoa$geometry) +
  tm_borders()+
  tm_shape(London_lsoa$centroid[London_lsoa$LSOANAME %in% no_direct_access,])+
  tm_dots(col = 'red')

# for this reason we decide to focus our ressearch on inner london, also to make it more computationally rapid
London_lsoa <- London_lsoa[London_lsoa$Cleaned_Names %in% inner_london,]

# outline of inner london and subsetting the points
outline <- London_lsoa %>% summarise(area = sum(area_sqkm))
outline <- st_transform(outline, 4326)
tm_shape(outline) +
  tm_polygons()
cultural_venues_prj <- st_transform(cultural_venues$polygon, 4326)
cultural_venues_ss <- cultural_venues_prj[outline,]
smth <- unlist(cultural_venues_ss$ward_code)
cultural_venues_out <- cultural_venues_prj %>% filter(!(ward_code %in% smth))

cent_to_cultural_calc_new <- dodgr_distances(graph_connected, from = st_coordinates(st_transform(London_lsoa$centroids, 4326)),
                                             to = st_coordinates(cultural_venues_ss), shortest = TRUE, pairwise = FALSE, quiet = FALSE)

distances_df <- as.data.frame(cent_to_cultural_calc_new)

# sorting out centroid distances to cultural amenities
London_lsoa$minimum <- apply(distances_df, 1, FUN = function(distances_df) (sort(distances_df))[1]) 
London_lsoa$second_lowest <- apply(distances_df, 1, FUN = function(distances_df) (sort(distances_df))[2])
London_lsoa$third_lowest <- apply(distances_df, 1, FUN = function(distances_df) (sort(distances_df))[3])

# time estimation in seconds when walking speed estimated a 5k per hour, this can be changed according to what predictions we decide to take
# walking distance from cultural amenities will here be used as a proxy for acces to culture

London_lsoa$min_walking_time <- full$minimum/(5000)*60
London_lsoa$min_walking_time <- times(full$min_walking_time / (24))
London_lsoa$avg_walking_time <- (full$minimum/5000*60)+(full$second_lowest/5000*60)+(full$third_lowest/5000*60)
London_lsoa$avg_walking_time <- full$avg_walking_time/3
London_lsoa$avg_walking_time <- times(full$avg_walking_time/24)


# now to get the indicators of deprivation that we will be using in our analysis
indicators_2 <- read_csv('Assignement/Indicators_2.csv')
indicators_2 <- indicators_2 %>% filter(`Local Authority District name (2019)` %in% inner_london)
indicators_2 <- rename(indicators_2, LSOANAME=`LSOA name (2011)`)

full_3 <- inner_join(test, indicators_2, by.x = "LSOANAME", by.y = "LSOANAME")


# map for education
plot1 <- tm_shape(full_3) +
  tm_fill("avg_walking_time", style = "quantile", n = 7, palette = "Oranges", title = "Average Walking Time (Seconds)") +
  tm_shape(st_transform((boroughs%>% filter(BOROUGHN %in% inner_london)),4326))+
  tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
  tm_text("BOROUGHN", size = "AREA") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

# map for education
plot2 <- tm_shape(full_3) +
  tm_fill("Education, Skills and Training Score", style = "quantile", n = 7, palette = "Blues", title = "Educational Deprivation Score") +
  tm_shape(st_transform((boroughs%>% filter(BOROUGHN %in% inner_london)),4326))+
  tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
  tm_text("BOROUGHN", size = "AREA") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)


# potentially remove if outliers from those residuals are too strong

###### full_3<-full_3[!(full_3$RESIDUALS>15),]

# Z score calculation attempt but not used in further research
full_3$z_score <- (full_3$avg_walking_time - mean(full_3$avg_walking_time)) / sd(full_3$avg_walking_time)

# putting values that go under 0 to positive by adding a constant and allow sqrt transformation
full_3$`Health Deprivation and Disability Score` <- full_3$`Health Deprivation and Disability Score`+3.2150
full_3$`Crime Score` <- full_3$`Crime Score`+2.343
full_3$`Children and Young People Sub-domain Score` <- full_3$`Children and Young People Sub-domain Score`+2.221
full_3$z_score <- full_3$z_score + 1.9682
full_3$`Income Score (rate)` <- full_3$`Income Score (rate)`*10

# sqrt so as to normalise the distributionsand then applying the first model with 4 variables with an R squared of 0.15
first_3 <- lm(formula = sqrt(avg_walking_time)~
                # sqrt(`Adult skills and English language proficiency indicator`)+
                sqrt(`Income Score (rate)`)+
                sqrt(`Health Deprivation and Disability Score`)+
                #sqrt(`Employment Score (rate)`) +
                sqrt(`Education, Skills and Training Score`)+
                sqrt(`Living Environment Score`),
              # sqrt(`Children and Young People Sub-domain Score`),
              data = full_3)

# Extract residuals from "modelLMR" object and dump into "spatialdatafile" and call the column "RESIDUALS"
full_3$RESIDUALS <- first_3$residuals

# Reporting basic summary measures to have an idea of its distribution before plotting them on map
summary(full_3$RESIDUALS)



# plotting the residuals to illustrate if spatial-autocorrellation is present
tm_shape(full_3) + tm_fill("RESIDUALS", style = "cont", midpoint = 0, palette = "-RdBu") +
  tm_shape(BOROUGHshp)+
  tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
  tm_text("BOROUGHN", size = "AREA") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)
# upderprediction of park and more exterior areas
# over_prediction central areas as most museums are concentrated in that area

#generate unique number for each row
full_3$ROWNUM <- 1:nrow(full_3)
# We need to coerce the sf spatialdatafile object into a new sp object
full_3.0 <- as(full_3, "Spatial")
# Create spatial weights matrix for areas
Weights <- poly2nb(full_3.0, row.names = full_3.0$ROWNUM)
WeightsMatrix <- nb2mat(Weights, style='B')
Residual_WeightMatrix <- mat2listw(WeightsMatrix , style='W')
# Run the test on the regression model output object "modelMLR" using lm.morantest()
lm.morantest(first_3, Residual_WeightMatrix, alternative="two.sided")

# there is evidence of high spatial autocorrelation therefore the use of appropriate spatial modelling is necessary
# moran's I at 0.6 (0.594) and p value under 0.01 threshold for significance

modelSLY <- lagsarlm(formula = sqrt(avg_walking_time)~
                       sqrt(Employment.Score..rate.)+
                       sqrt(Education..Skills.and.Training.Score)+
                       sqrt(Living.Environment.Score)+
                       #sqrt(Children.and.Young.People.Sub.domain.Score)
                       sqrt(Health.Deprivation.and.Disability.Score),
                     data = full_3.0, Residual_WeightMatrix)
summary(modelSLY)

# manually calculating the R-squared values for the spatial lag model
# this consider the residuals sum of squares (rss) divided by the total sum of squares (tss)
# the 1 represents predictions that are identical to the observed values
# we subtract rss/tss to account for all the errors present in the model and therefore give us a rsquared value
1-sum((modelSLY$residuals-mean(modelSLY$residuals))^2)/sum(((modelSLY$fitted.values)-mean(modelSLY$fitted.values))^2)


# extract the residuals for modelSLY object and dump back to original sf spatialdatafile object
full_3$RESID_SLY <- modelSLY$residuals
# use Moran's I test using moran.mc() function
moran.mc(full_3$RESID_SLY, Residual_WeightMatrix, 1000, zero.policy = T)

tm_shape(full_3) + tm_fill("RESID_SLY", style = "cont", midpoint = 0, palette = "-RdBu") +
  tm_shape(st_transform((boroughs%>% filter(BOROUGHN %in% inner_london)),4326))+
  tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
  tm_text("BOROUGHN", size = "AREA") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

# Interpretation of results using impacts
# impacts
Weights_2.0 <- as(Residual_WeightMatrix, "CsparseMatrix")
trMC <- trW(Weights_2.0, type="MC")
summary(impacts(modelSLY, tr = trMC, R=100), zstats=TRUE)

modelSER <- errorsarlm(formula =sqrt(avg_walking_time)~
                         sqrt(Employment.Score..rate.)+
                         sqrt(Education..Skills.and.Training.Score)+
                         sqrt(Living.Environment.Score)+
                         #sqrt(Children.and.Young.People.Sub.domain.Score)
                         sqrt(Health.Deprivation.and.Disability.Score), data = full_3.0,
                       Residual_WeightMatrix)
summary(modelSER)

# manually calculatin the r-squared value for the spatial error model
1-sum((modelSER$residuals-mean(modelSER$residuals))^2)/sum(((modelSER$fitted.values)-mean(modelSER$fitted.values))^2)

# extract the residuals for modelSLY object and dump back to original sf spatialdatafile object
full_3$RESID_SER <- modelSER$residuals
# use Moran's I test using moran.mc() function
moran.mc(full_3$RESID_SER, Residual_WeightMatrix, 1000, zero.policy = T)


tm_shape(full_3) + tm_fill("RESID_SER", style = "cont", midpoint = 0, palette = "-RdBu") +
  tm_shape(BOROUGHshp) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
  tm_text("BOROUGHN", size = "AREA") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)


# now on to the GWR to show how these variables can be illustrated

# finding bandwidth
bwG <- gwr.sel(sqrt(avg_walking_time) ~ sqrt(Employment.Score..rate.) + 
                 sqrt(Education..Skills.and.Training.Score) + sqrt(Living.Environment.Score) + 
                 sqrt(Health.Deprivation.and.Disability.Score), data = full_3.0, gweight = gwr.Gauss, verbose = FALSE)

# applying GWR model
gwrG <- gwr(sqrt(avg_walking_time) ~ sqrt(Employment.Score..rate.) + 
              sqrt(Education..Skills.and.Training.Score) + sqrt(Living.Environment.Score) + 
              sqrt(Health.Deprivation.and.Disability.Score), data = full_3.0, bandwidth = bwG, gweight = gwr.Gauss, hatmatrix = TRUE,
            se.fit = TRUE)

# results as data frame
gwr.data <- as.data.frame(gwrG$SDF)

full_3$GWR_education <- gwrG$SDF$sqrt.Education..Skills.and.Training.Score.
full_3$GWR_employment <- gwrG$SDF$sqrt.Employment.Score..rate..
full_3$R2 <- gwrG$SDF$localR2

full_3$tstat_education <- gwr.data$sqrt.Education..Skills.and.Training.Score./gwr.data$sqrt.Education..Skills.and.Training.Score._se
full_3$significant_education <- cut(full_3$tstat_education,
                                    breaks=c(min(full_3$tstat_education), -2, 2, max(full_3$tstat_education)),
                                    labels=c("Reduction: Significant","Not Significant", "Increase: Significant"))

full_3$tstat_employment <- gwr.data$sqrt.Employment.Score..rate../gwr.data$sqrt.Employment.Score..rate.._se
full_3$significant_employment <- cut(full_3$tstat_employment,
                                     breaks=c(min(full_3$tstat_employment), -2, 2, max(full_3$tstat_employment)),
                                     labels=c("Reduction: Significant","Not Significant", "Increase: Significant"))

tm_shape(full_3) + tm_fill("significant_education", style = "cont", midpoint = 0, palette = "-RdBu") +
  tm_shape(BOROUGHshp)+
  tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
  tm_text("BOROUGHN", size = "AREA") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("right")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

tm_shape(full_3) + tm_fill("significant_employment", style = "cont", midpoint = 0, palette = "-RdBu") +
  tm_shape(BOROUGHshp)+
  tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
  tm_text("BOROUGHN", size = "AREA") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

tm_shape(full_3) + tm_fill("R2", title = "Adaptive: Local R2", style = "cont", midpoint = 0.5, palette = "Spectral") +
  tm_shape(BOROUGHshp) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
  tm_text("BOROUGHN", size = "AREA") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("right")) +
  tm_layout(frame = FALSE, legend.title.size = 1, legend.text.size = 1)

summary(gwrG$SDF$localR2)


# figure 1
tm_shape(london_wards) +
  tm_polygons(alpha = 0, border.alpha = 0)+
  #tm_shape(cent$geometry) +
  #tm_dots()+
  tm_shape(boroughs)+
  tm_polygons(alpha = 1, border.alpha = 1, border.col = "black") +
  tm_shape(cultural_venues_ss$geometry)+
  tm_dots("red")+
  tm_shape(cultural_venues_out$geometry)+
  tm_dots("black")+
  # tm_text("BOROUGHN", size = "AREA", col = "blue", fontface = 'bold') +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

