# spatialML_package
This repository stores the spatial ml code 

## Packages 
```r
library("sf")            # Spatial data objects and methods
library("mapview")       # Interactive Map Viewing
library("ggmap")         # ggplot2 addon for base maps
library("cowplot")       # arrange ggplots
library("spatstat")      # KDE and other spatial functions
library("raster")        # cell-based spatial operations
library("tidyverse")     # data manipulation framework
library("Hmisc")         # using cut2() functions for ggplot legends
library("fitdistrplus")  # Distribution fitting functions
library("lubridate")     # Power tools for handling dates
library("tidycensus")    # download census data in tidy format
library("lwgeom")
library("Hmisc")
library("hrbrthemes")
library("gridExtra")
library("patchwork")
library("spdep")         # KNN functions
library("foreach")       # loops in parallel
library("doParallel")    # parallel backend
library("corrplot")      # correlation plot
library("ranger")        # randomforest implimentation      
library("glmnet")        # for Ridge and Lasso Regression
library("knitr")         # for kable table
library("kableExtra")    # nake nice tables in Rmarkdown
library("FNN")           # KNN for CPS vs. NN plots
library("groupdata2")
library("htmltools")
library("viridis")       # color palette
library("viridisLite")   # color palette
```

  ## Appendix 1: Data wrangling
  This appendix documents some key elements of loading, cleaning, and the initial manipulation of this project's data. As with each of the code appendices, the code blocks here do not show the entire source code, but do show most of the key functions and routines used to achieve these tasks. The code illustrated here follows the basic steps of:

1.  detect all `xls` and `csv` files in a directory and read them into a list
2.  Import spatial neighborhood data with `read_sf()` and `get_decennial()`
3.  create spatial fishnet grid with `st_make_grid()` and calculate spatial weights with `poly2nb()` and `nb2listw()`
4.intersect fishnet and census blocks to create populations estimates and weights per fishnet cell


This Code example of data import function. The inputs are `.xls`, `.xlsx`, and `.csv` files in a folder. The output is a `list` data type where each element of the list if one of the input data sets in the `sf` spatial data format.

```r
# Read data filenames
files <-list.files(file.path(base_dir,"/data"), pattern = "*\\.xls$|*\\.csv$")
var_list <- vector(mode = "list")
var_names <- NULL
# Loop over files and load
for(i in seq_along(files)){
  filename <- str_sub(files[i], start = 1, end = -5)
  sf_i <- tryCatch({
    if(tools::file_ext(files[i]) == "xls"){
      dat <- readxl::read_xls(file.path(base_dir,"data",files[i])) 
    } else if(tools::file_ext(files[i]) == "csv"){
      dat <- read.csv(file.path(base_dir,"data",files[i])) 
    }
    dat %>%
      filter(!is.na(X) | !is.na(Y)) %>%
      # Make files spatial `sf` object and assign geographic projection
      st_as_sf(., coords = c("X", "Y"), crs = 102747)
  }, error = function(e){
    cat(filename, "error = ",e$message,"\n")
    return(e)
  }
  )
  # add dataframes to list
  if(!inherits(sf_i, "error")){
    var_list[[length(var_list)+1]] <- sf_i
    var_names[length(var_list)] <- filename
  }
}
names(var_list) <- var_names
```


Code chunk demonstrating the reading of neighborhood and census block data from the web. The `read_sf()` function reads a `geojson` format file from the web and makes it into an `sf` files. This file is converted into a `raster` file type. the census block attribute for total population `P0010001` is read from the U.S. Census API using the `get_decennial` function. It is made into an `sf` object and calculations are made.
```r
# Get Richmond neighborhoods
nbr <- read_sf("https://data.richmondgov.com/resource/7juf-nwis.geojson")  %>%
  st_sf() %>%
  st_transform(102747)
# Create neighborhood outline
nbr_diss <- nbr %>%
  mutate(dissolve = 1) %>%
  # get rid of slivers
  st_buffer(., dist = 0.1) %>%
  group_by(dissolve) %>%
  summarise()

# Create neighborhood raster
nbr_rast_SP <- raster(as(nbr_diss, "Spatial"), nrows = 2000, ncol = 2000)

# Query Population data
vars10 <- c("P0010001") # total population
# Get total 2010 census population for blocks & calculate area
richmond_block <- get_decennial(geography = "block", variables = vars10, year = 2010,
                                summary_var = "P0010001", state = 51, county = 760, geometry = TRUE) %>%
  st_transform(crs = 102747)

# Calculate area
richmond_block <- richmond_block %>%
  mutate(acre = as.numeric(st_area(richmond_block)*2.29568e-5),
         # acre = units::set_units(acre, acre), 
         pop_acre_rate = value / acre) 
```


This block demonstrates the creation of the the spatial fishnet grid used as the unit of analysis for this project. The `sf_make_grid` function from the `sf` package makes this quite easy. A geographical intersect returns filters the fishnet to the cells that are with or touching the boundaries of Richmond City. Finally, the spatial neighborhood list is constructed. This is used by the spatial regression model to calculate spatial weights
```r
# set global parameter for fishnet grid dimensions
fishnet_grid_dim = 1200
# Fishnet creation
net <- st_make_grid(nbr, cellsize = fishnet_grid_dim) 

# Count CPS incidents per net cell
net_agg <- aggregate(cps_dissolve, net, sum) %>%
  tibble::rowid_to_column(.,"net_id")

# List of net cells IDs that intersect with Richmond Neighborhoods
net_intersect <- st_intersects(nbr, net_agg) 

# Extract Richmonds net cells based on intersect ID
net_Richmond <- net_agg[unique(unlist(net_intersect)),]

# Join neighborhood attributes to fishnet grid cells
net_hood <- st_join(net_Richmond, nbr, largest = TRUE)

# Calculate spatial neighborhood matrix
listw <- nb2listw(poly2nb(as(net_Richmond, "Spatial"), queen = TRUE))
```


In this final code example, the area weights population of each fishnet cell is calculated. The first step is to intersect the census blocks with the fishnet to create polygons of every intersecting region. After calculating the proportion of each Census block within each fishnet cell intersection and the matching proportion of that intersection, the `intersect_pop` populations are summed for each fishnet cell. The result is an estimate of the population of each fishnet cell based on the area and population of each census block that overlaps with it. Cells with zero population are dropped. Finally, the aggregated counts of maltreatment incidents are joined to the fishnet and a rate of incidents per 100 population is calculated.
```r
# Compute intersection of census blocks and fishnet
net_blocks_intersect <- st_intersection(richmond_block, net_Richmond)

# ... Calculate percent of each blocks area within each fishnet cell and divide block pop by that percent

# Summerise population by fishnet cell
fishnet_pop <- net_blocks_intersect %>% # xcc
  group_by(net_id) %>%
  summarise(net_pop = sum(intersect_pop)) %>%
  filter(net_pop > 0)   # <-  zeros or no zeros!!!!

# ... make cps_agg which is aggregate of each variable by cell

# Spatial join of fishnet_pop and fishnet_cps to then calculate rate for all CPS features
fishnet_pop_cps <- st_join(fishnet_pop, CPS_agg, join = st_equals) %>%
  mutate_at(vars(paste0("net_",CPS_vars)), funs(rate = ./(net_pop/100)))  %>% # cps per 100 person
  rename_at(vars( contains( "_rate")), funs(paste("rate", gsub("net_|_rate", "", .), sep = "_"))) %>% 
  replace(is.na(.), 0) # replace NA with zero

```


## Appendix 2: Feature engineering
This appendix documents some key elements of creating the variables, features, and data frames used in the train and test the machine learning models. As with each of the code appendices, the code blocks here do not show the entire source code, but do show most of the key functions and routines used to achieve these tasks. The code illustrated here follows the basic steps of:

1.  Compile relevant spatial point features from `var_list` into Protective and Risk specific spatial objects
2.  Create Nearest Neighbor, Aggregate Count, and Euclidean Distance features for Risk and Protective variables. Compose into single data frame of all features per fishnet cell
3.  Computer correlation between all features and dependent variable
4.  Select feature type with highest absolute correlation for each variable
5.  Combine selected features types for both Risk and Protective features and standardize independent variables


This code block illustrates the processes of selecting the pre-determined Risk and Protective variable from the `var_list` list that contained all for he candidate variables evaluated in this project. These steps are much more of a data processing step than an analytical step.
```r
# Create list of Protective variables
protective_class <- c("CommunityCenters","FireStations",
                      "HomelessShelters","Libraries","Parks","PointsOfInterest",
                      "PoliceStations","PublicSchools","ResourceOASIS","SNAP_WIC",
                      "VotingStations")
# Add column for variable set name and assign to new list
protective_vars <- list()
for(i in seq_along(protective_class)){
  dat <- var_list[[protective_class[i]]] %>%
    mutate(feature_name = protective_class[i],
           class = "protective") %>%
    dplyr::select(feature_name, class)
  protective_vars[[i]] <- dat
}

# Pull protective features from business features
Businesses_protective <- var_list[["BusinessProject"]] %>%
  filter(Classification == "PROTECTIVE") %>%
  mutate(feature_name = "BusinessProject",
         class = "protective") %>%
  dplyr::select(feature_name, class)
# Add to protective feature list
protective_vars[[length(protective_vars)+1]] <- Businesses_protective
# Turn into a spatial sf object with columsn for feature name, protective vs. risk class, and point coordinates
protective_vars <- do.call(rbind, protective_vars)
# Add dataframe to varaibles list
var_list[["Protective"]] <- protective_vars

## Simialr to above, but for Risk variables
risk_class <- c("BusStops")
risk_vars <- list()
for(i in seq_along(risk_class)){
  dat <- var_list[[risk_class[i]]] %>%
    mutate(feature_name = risk_class[i],
           class = "risk") %>%
    dplyr::select(feature_name, class)
  risk_vars[[i]] <- dat
}
# Business risks
Business_risk <- var_list[["BusinessProject"]] %>%
  filter(Classification == "RISK") %>%
  mutate(feature_name = "BusinessProject",
         class = "risk") %>%
  dplyr::select(feature_name, class)
risk_vars[[length(risk_vars)+1]] <- Business_risk
# Crime realted risks
CrimeData_risk <- var_list[["CrimeData"]] %>%
  filter(OFFENSE %in% c('DRUG/NARCOTIC VIOLATION','SIMPLE ASSAULT, DOMESTIC',
                        'Runaway','AGGRAVATED ASSAULT DOMESTIC')) %>%
  mutate(feature_name = "CrimeData",
         class = "risk") %>%
  dplyr::select(feature_name, class)
risk_vars[[length(risk_vars)+1]] <- CrimeData_risk
# Code violation risks
Violations_III_ks_risk <- var_list[["Violations_III_ks"]] %>%
  filter(CodeDsrp %in% c('General Violations','Unsafe Structure',
                         'Unfit Structure')) %>%
  mutate(feature_name = "Violations_III_ks_risk",
         class = "risk") %>%
  dplyr::select(feature_name, class)
risk_vars[[length(risk_vars)+1]] <- Violations_III_ks_risk
# Bind into dataframe
risk_vars <- do.call(rbind, risk_vars)
# Add to var_list
var_list[["Risk"]] <- risk_vars
```

This code block contains the three functions for turning the geographic point representation of each variable into a feature representation usable in regression or classification. This featureization is done in three different ways with the intent to capture different aspects of the spatial relationship between feature points (e.g. crime, code violations, certain businesses, etc...) and maltreatment events. The three ways of measuring this are 1) the mean Nearest Neighbor with `NN_point_features`, 2) the mean Euclidean distance to all points with `Euclidean_point_features`, 3) and the simple count of all points within each fishnet grid cell with `Aggregate_points_Features`. These function are written to use the `foreach()` and `%dopar` functions of the `foreach` package in order to execute the code using multiple processor cores. This helps to improve the efficiency of the computations. At the end of this code block, all of the features are joined together into a single data frame.
```r
# Functions to make thee feature types:
# Mean Nearest Neighbor distance
NN_point_features <- function(var_list, fishnet, k){
  NN_results <- foreach(i = seq_along(var_list),
                        .export=c('nn_function'),
                        .packages=c('raster', 'sf', 'dplyr', "FNN", "tibble", "tidyr")) %dopar% { 
                          feature <- names(var_list)[i]
                          # Create centroid of fishnet cell
                          fishnet_centroid_XY <- st_coordinates(st_centroid(fishnet))
                          dat <- var_list[[i]] 
                          # If more features the `k` nearest neighbors...
                          if(nrow(dat) >= k){
                            # get mean distance and id for `k` nearest neighbors ...
                            # from center of each fishnet cell
                            net_NN <- nn_function(fishnet_centroid_XY,
                                                  st_coordinates(dat)[,1:2], k) %>%
                              mutate(feature_name = paste0("NN_",feature),
                                     net_id = fishnet$net_id) %>%
                              left_join(., fishnet, by = "net_id") %>%
                              rename("value" = value.x) %>%
                              dplyr::select(-value.y) %>%
                              st_as_sf()
                          } else {
                            # If there are not enough features to make NN, then it is NA
                            net_NN <- data.frame(value = rep(NA, nrow(fishnet))) %>%
                              mutate(feature_name =  paste0("NN_",feature),
                                     net_id = fishnet$net_id) %>%
                              left_join(., fishnet, by = "net_id") %>%
                              rename("value" = value.x) %>%
                              dplyr::select(-value.y) %>%
                              st_as_sf()                       
                          }
                        }
  names(NN_results) <- paste0("NN_",names(var_list))
  return(NN_results)
}

# Nearest Euclidean Distance
Euclidean_point_features <- function(var_list, dist_raster, raster_mask, fishnet){
  ED_results <- foreach::foreach(i = seq_along(var_list), 
                                 .combine='comb', .multicombine=TRUE,
                                 .init=list(list(), list()),
                                 .export=c('distanceFromPoints', 'raster_to_fishnet'),
                                 .packages=c('raster', 'sf', 'dplyr')) %dopar% { 
                                   feature <- names(var_list)[i]
                                   # Calculate distance raster from all point features
                                   bs_dist <- distanceFromPoints(dist_raster, 
                                                                 sf::st_coordinates(var_list[[feature]]))
                                   # Clip raster to Richmond City
                                   bs_clip <- raster::mask(bs_dist, mask = as(raster_mask, "Spatial"))
                                   # Extract raster distance values within each fishnet cell and take mean
                                   fea_mean_dist <- raster_to_fishnet(bs_clip,fishnet,paste0("ed_",feature))
                                   list(fea_mean_dist, bs_clip)
                                 }
  # Add results to list
  dist_results <- ED_results[[1]]
  dist_rasters <- ED_results[[2]]
  names(dist_results) <- paste0("ED_",names(var_list))
  names(dist_rasters) <- paste0("ED_",names(var_list))
  return(list(dist_results, dist_raster))
}

# Aggregare Point Count
Aggregate_points_Features <- function(var_list, fishnet){
  agg_results <- foreach(i = seq_along(var_list),
                         .packages=c('raster', 'sf', 'dplyr')) %dopar% { 
                           feature <- names(var_list)[i]
                           dat <- var_list[[i]] %>%
                             mutate(value = 1) %>%
                             dplyr::select(value)
                           # Count feature points within each fishnet grid cell
                           net_agg <- aggregate(dat, fishnet, sum) %>%
                             mutate(feature_name = paste0("agg_",feature),
                                    net_id = fishnet$net_id)
                         }
  names(agg_results) <- paste0("agg_",names(var_list))
  return(agg_results)
}

#Make `ALL_FEATUES` dataframe by joining all Mean Nearest Neighbor, Mean Euclidean Distance, and Aggregare Count features
ALL_FEATURES <- full_join(NN_features, agg_features, by = "net_id") %>%
  full_join(.,ED_features, by = "net_id") %>%
  # Join in census data features
  full_join(.,sf1_features, by = "net_id")

# Clean up joined fields and replace `NA` with zeros
ALL_FEATURES <- ALL_FEATURES %>%
  dplyr::select(-cps_rate.y, -cps_rate.x.x, -cps_rate.y.y, 
                -cps_net.y, -cps_net.x.x, -cps_net.y.y,
                -net_pop.y, -net_pop.x.x, -net_pop.y.y) %>%
  dplyr::select(-contains("_CPS_")) %>%
  dplyr::rename(cps_net  = cps_net.x,
                cps_rate = cps_rate.x,
                net_pop  = net_pop.x) %>%
   mutate_all(funs(replace(., is.na(.), 0)))  %>%
  dplyr::rename_all(funs(make.names(.)))

```


Following the featurization demonstrated above, the pairwise Pearson's correlation of each feature and the dependent variable is computed using the `cor()` function. In order to perform feature selection on the three different ways of measuring the spatial relationships (e.g. mean Nearest Neighbor, mean Euclidean distance, and Aggregate count), the feature type with the highest absolute correlation is selected. This is achieved using the combination of `slice()` with `which.max()` functions.
```r
# Compute correlation between all pariwise features
cps_cor_ALL <- cor(ALL_FEATURES)
All_cors <- cps_cor_ALL[,"cps_net"]
# Compute p-values
p.mat_ALL <- cor.mtest(ALL_FEATURES)$p
p.mat_ALL <- p.mat_ALL[,which(colnames(cps_cor_ALL)=="cps_net")]
# Prepare data for plotting
cor_ALL_plot <- data.frame(feature = names(All_cors), 
                           cor = as.numeric(All_cors),
                           p_value   = p.mat_ALL) %>%
  filter(!(feature %in% c("cps_rate","cps_net","net_pop","net_cps","net_id"))) %>%
  filter(!(feature %in% grep("CPS", names(All_cors),value=T))) %>%
  arrange(desc(cor)) %>% 
  mutate(p_value = ifelse(p_value >= 0.1, "Not Significant", "Significant"))
cor_ALL_plot$feature <- factor(cor_ALL_plot$feature,
                               levels=cor_ALL_plot[order(cor_ALL_plot$cor,
                                                         decreasing=F),]$feature)

# For wach Protective variable, extract the feature type (mean NN, mean distance, aggregate) that as the highest absolute correlation
features_strong_protective_names <- cor_ALL_plot %>% 
  filter(feature %in% names(features_protective_all)) %>%
  mutate(prefix = str_extract(feature, "^[^_]+(?=_)"),
         suffix = str_extract(feature, "(?<=_)[^_].*"),
         feature = as.character(feature)) %>%
  group_by(suffix) %>%
  # Highest absolute correlation
  slice(which.max(abs(cor)))
features_protective_strong <- features_protective_all %>%
  dplyr::select(features_strong_protective_names$feature,
                NN_CPS_Accepted,
                cps_net, cps_rate, net_pop, net_id)

# Do the same for Protective
```

Finally, once the feature representation for each Risk and Protective variable is selected based on correlation, the features are joined into a single data frame. With some additional minor manipulation, this data frame will be the basis of the regression modeling that follows.
```r
# Join Risk and Protective features to selected census feaures to form data for regression models
og_dat <- full_join(features_risk_strong, features_census_select, by = "net_id") %>%
  full_join(., features_protective_strong, by = "net_id") %>% 
  dplyr::select(-net_pop.y, -cps_net.y, -cps_rate.y,
                -net_pop.x, -cps_net.x, -cps_rate.x,
                -NN_CPS_Accepted.y) %>% 
  rename("NN_CPS_Accepted" = NN_CPS_Accepted.x)

# Remove unneeded columns and center & scale all variables except the dependent variable
dat    <- og_dat %>% dplyr::select(-cps_rate, -net_pop, -net_id) %>%
  mutate_at(vars(-cps_net), scale_this)
# Add neighborhood name
net_hood <- st_join(net_Richmond, nbr, largest = TRUE)
og_dat$.block_id <- net_hood$name
```


## Appendix 3: Exploratory Analysis
This appendix documents some key elements of exploring and visualizing the features of this data set. The process of viewing, aggregating, and testing the data set in various undirected ways is referred to as Exploratory Data Analysis (EDA). As with each of the code appendices, the code blocks here do not show the entire source code, but do show most of the key functions and routines used to achieve these tasks. The code illustrated here follows the basic steps of:
  
1.  Mapping maltreatment event counts and rates by fishnet cell
2.  Plotting the redecoration of maltreatment events across a range of time scales
3.  Visualizing the spatial trends of all Risk and Protective features
4.  Calculating Global and Local Moran's I statistics


The first exploration of code is to map the count and rate of maltreatment incidents by fishnet cells. This done at a range of fishnet cells dimensions and for a variety of types of maltreatment events to better understand the spatial trends. The `make_cuts()` function uses the `Hmisc::cut2()` function to partition the data. In the `ggplot2` plot, the `ggmap::ggmap()` function allows for data to be overlain on base maps of various styles.
```r
# bin the countof maltreatment events
fishnet_pop_cps_cut <- fishnet_pop_cps %>%
  mutate(net_CPS_Accepted = ifelse(is.na(net_CPS_Accepted), 0, net_CPS_Accepted)) %>% 
  make_cuts(., "net_CPS_Accepted", cuts = "breaks", n_breaks = 10)
# plot maltreatment event counts on top of Richmond City basemap
CPS_COUNT_BY_FISHNET_PLOT <- ggmap(cps_base_map) +
  geom_sf(data = ll(fishnet_pop_cps_cut), aes(fill = cut_val), inherit.aes = FALSE, color = NA) +
  labs(title = "CPS Count per Fishnet Cell") +
  scale_fill_viridis_d(na.value = NA, option = "D", direction = 1, name = "CPS Count") +
  theme_bw()
```


This code example demonstrates a visualization of maltreatment events by normalized by month for each year in the data set. This shows each month of each years deviation from the average number of maltreatment events for that month across all years. This visualization also shows that for the years 2013 and 2017 only one-half of the years data is available. A number of examples of similar plots were viewed to understand the temporal patterns in maltreatment event redecoration. As with the vast majority of plots, the `ggplot2` library is used here.
```r
# Normalizes values by month
CPS_normalized_by_month <- st_drop_geometry(var_list[["CPS_Accepted"]]) %>%
    mutate(month = lubridate::month(RDate),
           year  = lubridate::year(RDate)) %>%
    group_by(year, month) %>%
    summarise(m_total = n()) %>%
    arrange(month, year) %>%
    dplyr::select(month, year, m_total) %>%
    ungroup() %>%
    group_by(month) %>%
    # Normalize
    mutate(m_mean = mean(m_total),
           m_sd   = sd(m_total),
           m_z    = (m_total - m_mean) / m_sd)
# Plot
CPS_LINE_NORMALIZED_plot <- ggplot(CPS_normalized_by_month, aes(x = as.factor(month), 
                                          y = m_z, group = year, 
                                          color = as.factor(year))) +
    geom_line() +
    geom_hline(yintercept = 0, color = "gray20", linetype = "dashed") +
    scale_y_continuous(limits = c(-2,2)) +
    theme_bw()
```


Additional spatial trends are describable from plotting Kernel Density Estimates (KDE) of each of the Risk and Protective features. A KDE is a spatial interpolation technique that can smooth out the spatial variation on counts to create a more continuous surface.  This surface can indicate larger scale trends that are not as evident with fishnet cell based features. Below is the code example for the Risk plot; the same is done for Protective features. The `geom_tile()` function in the `ggplot2` library allows from plotting raster data. The `dplyr` package `ntile()` function bins the values.
```r
# Risk KDE facet plot
RISK_KDE_FACET_PLOT <- ggmap(cps_base_map) +
  geom_tile(data = risk_plot_dat, 
            aes(x,y,fill = as.factor(ntile(value,brks)), 
                group = variable), alpha=0.5) +
  scale_fill_viridis_d(name = variable) +
  facet_wrap(~variable) +
  theme(
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    legend.key = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none"
  )
```


This code block shows the general method for simulating possible Global Moran's I values given the spatial structure of the fishnet grid. The `knearneigh()`, `knn2nb`, and `nb2listw` from the `spdep` package are the tools needed to convert from the fishnet to the spatial weights list. The `moran.mc()` function from the same package does the work of simulating permutations of the Global  Moran's statistic and returning `nsim` number of results (see note in code comments). The observed Global  Moran's I value is also returned. This visualization uses `geom_vline()` to plot the observed value in reference to the distribution of simulated values to demonstrate how likely or unlikely the observed value is.
```r
# Calculate fishent grid nearest neighbors using the "Queen" case (all directions)
# convert to a neighborhood graph
fishnet_knn <- knn2nb(knearneigh(fishnet_coords, k_direction))
# List of neighborhood weights
fishnet_Weights <- nb2listw(fishnet_knn, style="W")
# Simulate Global Moran's I with `999` simulations
globalMorans <- moran.mc(fishnet_pop_cps$net_CPS_Accepted, fishnet_Weights, nsim=999)

# Plot Global Moran's I simulated distribution and observed value
# note: the `globalMorans` I object will have 1000 rows, but only the first 999 are simulated. 
# The 1000th row is the observed value
GLOBAL_MORANS_PERMUTATION_plot <- ggplot(data.frame(res = globalMorans$res)[1:999,,0], aes(res)) + 
  geom_histogram(binwidth = 0.01) +
  geom_vline(aes(xintercept = globalMorans$statistic), colour = "red",size=1) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title="Observed and permuted Moran's I", x = "Simulated Moran's I Value") +
  theme_bw()
```


The final code example from the EDA phase of this project is the calculation and visualization of the Local Moran's I values. Using the same spatial neighborhood weights as above, the Local Moran's I is calculated with the `localmoran()` function of the `spdep` package. The data are joined together and plotted as three separate fishnet plots; the count of maltreatment events, binned Local Moran's I values, and the Local Moran's I p-value classified at alpha = `0.9` into significant and not significant clusters.
```r
# calculate Local Moran's I for each fishnet cell
localMorans  <- as.data.frame(localmoran(fishnet_pop_cps$net_CPS_Accepted, fishnet_Weights))
# Moran's I join
fishnet_pop_cps_morans        <- fishnet_pop_cps
fishnet_pop_cps_morans$Ii     <- localMorans$Ii
fishnet_pop_cps_morans$pvalue <- localMorans$`Pr(z > 0)`
fishnet_pop_cps_morans        <- cbind(fishnet_coords, fishnet_pop_cps_morans)

# Bin maltreatment event counts
fishnet_pop_cps_morans_cut <- make_cuts(fishnet_pop_cps_morans, "net_CPS_Accepted",
                                        cuts = "breaks", n_breaks = 10)
# Plot maltreatment events (assign to a variable)
plot_cps <- ggmap(cps_base_map) +
  geom_sf(data = ll(fishnet_pop_cps_morans_cut), aes(fill = cut_val),
          color = NA, inherit.aes = FALSE, alpha = 0.9) +
  scale_fill_viridis_d(na.value=NA, name = "Maltreatment Events", option = "D" ) +
  theme_void() +
  theme(
    legend.position = "right"
  )

# Bin Local Moran's I statistic
Ii_cut <- fishnet_pop_cps_morans %>%
  mutate(Ii_cut_val = as.character(Hmisc::cut2(.$Ii, 
                                               cuts = as.numeric(quantile(round(fishnet_pop_cps_morans$Ii,2), 
                                                                          na.rm=T, p = seq(0,1,0.25))))))
# plot binned Local Moran's I statistic (assign to a variable)
plot_Ii <- ggmap(cps_base_map) +
  geom_sf(data = ll(Ii_cut), aes(fill = Ii_cut_val),
          color = NA, inherit.aes = FALSE, alpha = 0.9) +
  scale_fill_viridis_d(na.value=NA, name = "Local Moran's I", option = "D")+
  theme_void()+
  theme(
    legend.position = "right"
  )

# Bin Local Moran's I p-value
p_cut <- fishnet_pop_cps_morans %>%
  mutate(pval_cut = ifelse(pvalue > 0.05, "Not Significant", "Significant"))
# Plot binned p-value (assign to a variable)
plot_p <- ggmap(cps_base_map) +
  geom_sf(data = ll(p_cut), aes(fill = pval_cut),
          color = NA, inherit.aes = FALSE, alpha = 0.9) +
  scale_fill_viridis_d(na.value=NA, name = "p-value", option = "D")    +
  theme_void()+
  theme(
    legend.position = "right"
  )
# use `cowplot` to put plots together
MORANS_I_P_plot <- cowplot::plot_grid(plot_cps, plot_Ii, plot_p, 
                                      rel_widths = c(0.9,0.9,0.9),
                                      ncol = 1, align = "v")
```


## Appendix 4: Modeling and Evaluation

This appendix documents some key elements of fitting machine learning models and evaluating their results. In these steps, the data manipulated and explored above is modeled to find the patterns and project those patterns back to the project area. From these projections, the goddness of fit is calculated and biases are uncovered. The purpose of this is to better understand the areas where factors contributing to maltreatment are present, but unreported, as well as understand conditions where the model may work better or worse. As with each of the code appendices, the code blocks here do not show the entire source code, but do show most of the key functions and routines used to achieve these tasks. The code illustrated here follows the basic steps of:
  
1.  Refine data sets and establish the cross-validation folds (one fold for each neighborhood)
2.  Fit three machine learning models (Poisson GLM, Random Forest, and Spatial Durbin) to each neighborhood
3.  Extract model predictions and join to the fishnet and neighborhood spatial data
4.  Use the prediction from three models to create a meta-model via model stacking with Random Forest
5.  Calculate goodness of fit for both each Quintilian of the data distribution and overall
5.  Download neighborhood statistical area shapes, divide into high/low classes, and compute error
6.  Compute Kernel Density Estimate (KDE) of data distribution and compare to meta-model predictions

The approach used in this report to facilitate cross-validation is to use a combination of `tibble` data frames and `purrr::map()` functions within `dplyr::mutate()`. This approach allows for a single table to include all CV folds, models, and derivative errors or metrics. By using this method, the data are compiled into a single object where they can be easily subset, viewed, and analyzed. The cons of this approach are that there is a bit more computational overhead (but not significant in this case) and the syntax may not be familiar to sum R users. In the code below, the first four lines create the neighborhood fixed effects and adds them to the regression data. Following this, the `groupdata2::fold()` function is used to create the CV fold index stratified by the neighborhood `.block_id`. Finally, the `tibble`, `purrr`, `dplyr` routine creates the CV fold table used in the model fitting section.

```r
# neighborhood fixed effect model
hood_matrix <- model.matrix(cps_net~.block_id,og_dat)
hood_model <- lm(sqrt(og_dat$cps_net) ~ hood_matrix)
dat$hood_fixed <- predict(hood_model, type = "response")^2
og_dat$hood_fixed <- predict(hood_model, type = "response")^2
# Create CV neighborhood fold index
all_hoods <- length(unique(net_hood$name))
n_folds = ifelse(n_folds == "LOOCV", all_hoods, n_folds)
folds_index <- groupdata2::fold(og_dat, k = n_folds, id_col = '.block_id')$.folds
# create tibble with all CV folds and assocaited data
cv_tbl <- tibble(folds = seq_len(n_folds),
                 train = NA, train_y = NA, train_index = NA, train_net_id = NA,
                 test  = NA, test_y  = NA, test_index  = NA, test_net_id  = NA)
for(k in seq_len(n_folds)){
  fold_i  <- which(folds_index == k)
  cv_tbl[k,]$train         <- list(dat[-fold_i,])
  cv_tbl[k,]$test          <- list(dat[ fold_i,])
  cv_tbl[k,]$train_y       <- list(og_dat[-fold_i,target_var])
  cv_tbl[k,]$test_y        <- list(og_dat[ fold_i,target_var])
  cv_tbl[k,]$train_index   <- list(setdiff(seq(1:nrow(dat)),fold_i))
  cv_tbl[k,]$test_index    <- list(fold_i)
  cv_tbl[k,]$train_net_id  <- list(og_dat[-fold_i,"net_id"])
  cv_tbl[k,]$test_net_id   <- list(og_dat[ fold_i,"net_id"])
  ```  
  
  The model fitting used helper functions to make the use of `purrr::map()` more streamlines. The `glm_fit`, `rf_fit`, and `score_model()` functions make a consistent interface for the `map()` to interact with the `glm()`, `ranger()`, and goodness of fit functions functions. The spatial model is fit directly with `errorsarlm()` because it doe not use the CV folds due to the neighborhood weight matrix. For each model, after it is fit, predictions are made, and the model is scored.
```r
  # Helper function for GLM fit
  glm_fit <- function(dat, formula, family){
    glm_model <- glm(formula, data = dat, family = family)
    return(glm_model)
  }
  # Helper funtion for Random Forest fit
  rf_fit <- function(dat, formula, mtry_add = 0, importance = "none"){
    mtry <- floor(sqrt(ncol(dat)-1))+mtry_add
    rf_model <- ranger(formula, data = dat, 
                       mtry = mtry,
                       splitrule = "variance",
                       importance = importance,
                       num.trees = 500,
                       min.node.size = 10)
    return(rf_model)
  }
  # Helper function for model scoring
  score_model <- function(dat){
    dat <- dat %>%
      mutate(R2     = map2_dbl(pred, test_y, r_squared),
             MAE    = map2_dbl(pred, test_y, mae),
             MAAPE  = map2_dbl(pred, test_y, maape),
             RMSE   = map2_dbl(pred, test_y, rmse),
             logdev = map2_dbl(pred, test_y, logdev_p))
    return(dat)
  }
  # Fit the Poisson GLM model
  po_cv_tbl <- cv_tbl %>%
    mutate(fit   = map(train, glm_fit, 
                       formula =  paste("cps_net ~ ."), 
                       family = "poisson"),
           pred  = map2(fit, test, lm_predict, sqrt = FALSE),
           mdl_nam = "GLm - Poisson") %>% 
    score_model()
  # Fit the Random Forest model
  rf_cv_tbl <- cv_tbl %>%
    mutate(fit   = map(train, rf_fit, formula = "cps_net ~ .", mtry_add = 2, importance = "impurity"),
           pred  = map2(fit, test, lm_predict),
           mdl_nam = "Random Forest") %>% 
    score_model()
  # Fit the Spatial Durbin model
  spat_durbin <- errorsarlm(sqrt(cps_net) ~ ., data = dat, listw, etype ="emixed")
  spat_durbin_tbl <- tibble(
    fit   = list(spat_durbin),
    pred  = map(fit, sar_pred),
    test_y= list(dat$cps_net),
    test_net_id = list(og_dat$net_id),
    mdl_nam = "Spatial Durbin - sqrt") %>% 
    score_model()
  ```
  
  This code section demonstrates the extraction of the out-of-fold predictions from the model tibbles and joining of these to each other. Only the example for the Poisson regression is show, but the code for the other two models is nearly identical. The `unnest()` function allows for the data bound in the tibble to be extracted. The sequence of `left_join()` creates a data set of each models predictions and the associated  neighborhood characteristics. These are the data that are used to fit the stacked meta-model.

```r
  # Extract predictions from Possion model (same is done for other two models)
  po_pred_dat <- po_cv_tbl %>%
    unnest(pred) %>%
    mutate(test_y = po_cv_tbl %>% unnest(test_y) %>% pull(test_y),
           test_net_id = po_cv_tbl %>% unnest(test_net_id) %>% pull(test_net_id))
  # Plot prediction map
  po_pred_geoplot <- model_pred_geoplot(po_pred_dat$pred,
                                        po_pred_dat$test_y,
                                        po_pred_dat$test_net_id,
                                        net_Richmond, cps_base_map, "po")
  # Join all predictions from three models into single data frame
  cps_preds <- og_dat %>% 
    dplyr::select(net_id, cps_net) %>% 
    left_join(., dplyr::select(po_pred_dat,
                               net_id = test_net_id,
                               pred_lm = pred), by = "net_id") %>%
    left_join(., dplyr::select(rf_pred_dat, 
                               net_id = test_net_id,
                               pred_rf = pred), by = "net_id") %>% 
    left_join(., dplyr::select(sarlm_pred_dat, 
                               net_id = test_net_id,
                               pred_sarlm = pred), by = "net_id") %>% 
    mutate_if(is.double, round, 2)
```
  
  As was done above, a `tibble` is created to contain the CV folds; one for each neighborhood. The fold index `fold_index` is recycled here to make sure that the out-of-fold predictions are used here. This assures that the meta-model is fit to data that are independently predicted from the model. The final block of code reuses the `rf_fit()` helper function to fit the candidate model predictions with the `Random Forest` algorithm using `ranger()`.

```r
  # Create CV neighborhood folds for meta-model
  cps_preds_cv_dat <- dplyr::select(cps_preds, -net_id)
  ens_cv_tbl <- tibble(folds = seq_len(n_folds),
                       train = NA, train_y = NA, train_index = NA, train_net_id = NA,
                       test  = NA, test_y  = NA, test_index  = NA, test_net_id  = NA)
  for(k in seq_len(n_folds)){
    fold_i  <- which(folds_index == k)
    ens_cv_tbl[k,]$train         <- list(cps_preds_cv_dat[-fold_i,])
    ens_cv_tbl[k,]$test          <- list(cps_preds_cv_dat[ fold_i,])
    ens_cv_tbl[k,]$train_y       <- list(cps_preds_cv_dat[-fold_i,target_var])
    ens_cv_tbl[k,]$test_y        <- list(cps_preds_cv_dat[ fold_i,target_var])
    ens_cv_tbl[k,]$train_index   <- list(setdiff(seq(1:nrow(cps_preds_cv_dat)),fold_i))
    ens_cv_tbl[k,]$test_index    <- list(fold_i)
    ens_cv_tbl[k,]$train_net_id  <- list(cps_preds[-fold_i,"net_id"])
    ens_cv_tbl[k,]$test_net_id   <- list(cps_preds[ fold_i,"net_id"])
  }
  # Fit meta-model with Random Forest
  ens_cv_tbl <- ens_cv_tbl %>%
    mutate(fit   = map(train, rf_fit, formula = "cps_net ~ pred_rf + pred_sarlm"),
           pred  = map2(fit, test, lm_predict),
           # pred  = map(pred, round),
           mdl_nam = "Meta-Model") %>% 
    score_model()
  ```
  
  These code examples show the compilation of the goodness of fit measures for each model. The errors are aggregated in two different ways; first over the quantiles of maltreatment event counts and secondly by mean and standard deviation across all neighborhoods. The quantile errors are calculated by grouping the predictions with the `quantile_error()` function. The last block of code uses the base `mean` and `sd` functions to calculate metrics grouped by each model. The `knitr::kable()` and `kableExtra` packages are used to present the findings.
  
```r  
  # Helper function for quantile error
  quantile_error <- function(pred,obs,quant){
    preds <- data.frame(pred = pred, obs = obs) %>%
      filter(quantile(seq(0,max(obs)), quant)>obs)
    return(preds)
  }
  # Join/bind model prediction tables
  models <- bind_rows(rf_cv_tbl, spat_durbin_tbl, ens_cv_tbl, po_cv_tbl)
  # Unnest predictions by model
  CV_preds_long <- models %>%
    group_by(mdl_nam) %>%
    unnest(pred, test_y) 
  ## Map over all quantiles to get error metrics
  quantile_errors <- CV_preds_long %>%
    nest(-mdl_nam) %>%
    mutate(q      = list(seq(0,1,0.01)),
           pred   = map(data, "pred"),
           test_y = map(data, "test_y")) %>%
    dplyr::select(-data) %>%
    unnest(q, .preserve = c(pred, test_y)) %>%
    filter(q != 0) %>% 
    mutate(q_dat  = pmap(list(pred, test_y, q), quantile_error),
           q_pred = map(q_dat, "pred"),
           q_obs  = map(q_dat, "obs"),
           q_RMSE = map2_dbl(q_pred, q_obs, rmse),
           q_MAE  = map2_dbl(q_pred, q_obs, mae),
           q_logdev  = map2_dbl(q_pred, q_obs, logdev_p),
           y_max  = quantile(seq(0,max(dat$cps_net)), q),
           q_cnt  = nrow(og_dat) - map_int(q_dat, nrow))
  
  # Map over all predictions grouped by model to calculate mean and sd for error metrics
  model_results <- models %>%
    dplyr::select("Model Name" = mdl_nam, R2, RMSE, MAE, logdev) %>%
    group_by(`Model Name`) %>%
    arrange(`Model Name`) %>%
    summarise(R2_mean      = mean(R2, na.rm=TRUE),
              R2_sd        = sd(R2, na.rm=TRUE),
              MAE_mean     = mean(MAE, na.rm=TRUE),
              MAE_sd       = sd(MAE, na.rm=TRUE),
              RMSE_mean    = mean(RMSE, na.rm=TRUE),
              RMSE_sd      = sd(RMSE, na.rm=TRUE),
              logdev_mean  = mean(logdev, na.rm=TRUE),
              logdev_sd    = sd(logdev, na.rm=TRUE)) 
  Model_Error_Results_table <- model_results %>%
    kable(., format = "html", digits = 3) %>%
    kable_styling()
  ```
  
  An important part of the evaluation of the meta-model is to compare it against important demographic factors that are not explicitly used in the model formulas. The `read_sf()` function retrieves spatial data for the neighborhood statistical groups (NSG) and `get_acs()` retrieves census data. The code that follows manipulates the census and NSG data by calculating percentages of population. After classifying the population percentage for both percent non-white and percent poverty, they are classified into high and low by diving at the 60the percentile. The errors are then aggregated to high and low classes using `aggregate()` and `make_cuts()`. Likewise the count of maltreatment events are aggregated in a similar manner. The last code chunk shows how the metrics are aggregated for poverty, but that same is done for percent non-white. 
  ```r
  # Download statarea
  nbr_statAreas <- read_sf("https://data.richmondgov.com/resource/8kyq-v9j2.geojson") %>%
    st_transform(crs = 102747) %>% 
    mutate(stat_area_id = id)
  
  # Download poverty and population data
  tract10 <- get_acs(geography = "tract", variables = c("B02001_001","B02001_002E","B17001_002"), 
                     year = 2010, state=51, county=760, geometry=T)
  # Aggregate data and make percentages
  tract10 <- tract10 %>%
    dplyr::select(variable,estimate) %>%
    as.data.frame() %>%
    spread(variable,estimate) %>%
    rename(TotalPop=B02001_001,
           NumberWhites=B02001_002,
           TotalPoverty=B17001_002) %>%
    mutate(percentNonWhite = ifelse(TotalPop > 0, ((TotalPop - NumberWhites) / TotalPop),0),
           percentPoverty  = ifelse(TotalPop > 0, TotalPoverty / TotalPop, 0),
           tract_id        = row_number()) %>%
    st_sf() %>%
    st_transform(102747) 
  tract10$tract_area <- st_area(tract10)
  
  # Aggreate mean errors to statareas
  stat_area_metric_logdev <- error_points %>%
    aggregate(., nbr_statAreas.spJoin, mean) %>%
    dplyr::select(logdev) %>% 
    mutate(logdev = round(logdev, 3)) %>% 
    make_cuts(., "logdev")
  stat_area_metric_MAE<- error_points %>%
    aggregate(., nbr_statAreas.spJoin, mean) %>%
    dplyr::select(MAE) %>% 
    mutate(MAE = round(MAE, 3)) %>% 
    make_cuts(., "MAE")
  
  # Aggregate sum of CPS incidents to statarea
  stat_area_cps <- error_points %>%
    aggregate(., nbr_statAreas.spJoin, sum) %>%
    dplyr::select(test_y)
  stat_area_errors <- stat_area_metric_logdev %>% 
    st_join(., stat_area_metric_MAE, join = st_equals) %>% 
    st_join(., stat_area_cps, join = st_equals) %>% 
    st_join(., nbr_statAreas.spJoin, join = st_equals)
  
  # Group by poverty and get median of statarea aggregate errors
  poverty_aggregate <- stat_area_errors %>% 
    group_by(poverty.percentile) %>% 
    summarise(med_dev = round(median(logdev),3),
              med_MAE = round(median(MAE),3),
              med_CPS = sum(test_y)) %>% 
    st_drop_geometry() %>% 
    dplyr::select(poverty.percentile, med_dev, med_MAE, med_CPS)
  ```
  
  The final code example from the meta-model evaluation phase is the comparison of prediction density between a KDE of maltreatment events and the meta-model predictions. The main parts of this section are the creation of the KDE using `spatstat::density.ppp`, the normalization of the predictions and KDE with the `bin_class()` function, and finally the plotting of spatial data and a bar plot for comparison. The `bin_class()` function takes in a vector of maltreatment event predictions or densities, bins them into quantiles with `.bincode()`, and then breaks the quantiles into five sensitivity classes with `cut()`. The code below shows the maltreatment event predictions aggregated by sensitivity class into the `p.summ` object. The approach is shown for aggregating the events by KDE sensitivity `kde.summ`. Finally, the percent of maltreatment events per sensitivity zone are counted for each of the meta-model predictions and KDE density sensitivity classes and plotted as a bar chart.
  ```r
  # Function to create sensitivity classes by binning and cutting quantities
  bin_class <- function(dat, bin_col = "pred", 
                        quantile_labels = 100, break_vec = c(-1, 30, 50, 70, 90, 100)){
    if(is(dat, "sf")){
      dat <- st_drop_geometry(dat)
    }
    pred_bin <- as.numeric(.bincode(dat[,bin_col]+1e-8, # wiggle factor to get above zero
                                    breaks = quantile(dat[,bin_col],
                                                      seq(0,1, by=(1/quantile_labels)), 
                                                      na.rm = TRUE,
                                                      labels = seq(1,quantile_labels,1))))
    pred_bin_class <- as.numeric(cut(pred_bin, 
                                     breaks = break_vec, 
                                     na.rm  = TRUE,
                                     labels = seq(1,length(break_vec)-1,1)))
    pred_bin_class <- ifelse(is.na(pred_bin_class), length(break_vec)-1, pred_bin_class)
  }
  
  error_geoplot$pred_bin_class <- bin_class(error_geoplot, "pred")
  # Count recorded incidents by meta-model prediction sensitivity classes
  p.summ <- error_geoplot %>%
    group_by(pred_bin_class) %>%
    dplyr::summarize(obs.total = sum(test_y),
                     obs.cnt = n()) %>% 
    rename(sens_group = pred_bin_class) %>%
    filter(!is.na(sens_group)) %>%
    identity()
  
  # Compute KDE
  cps_ppp <- as.ppp(st_coordinates(cps_dissolve), W = st_bbox(net_Richmond))
  cps_KDE <- spatstat::density.ppp(cps_ppp)
  # COnvert KDE to fishnet grid
  cps_KDE_tbl <- as.data.frame(cps_KDE) %>%
    st_as_sf(coords = c("x", "y"), crs = 102747) %>%
    aggregate(., net_Richmond, mean) %>%
    mutate(net_id = net_Richmond$net_id)
  
  error_geoplot$kde_bin_class  <- bin_class(cps_KDE_tbl, "value")
  # Count incident counts by KDE sensitvity classes
  kde.summ <- error_geoplot %>%
    group_by(kde_bin_class) %>%
    dplyr::summarize(kde.total = sum(test_y),
                     kde.cnt = n()) %>% 
    rename(sens_group = kde_bin_class) %>%
    filter(!is.na(sens_group)) %>%
    identity()
  
  # Compile counts of incidents from KDE and predictions then create percentages
  countComparisons <- merge(st_drop_geometry(p.summ), st_drop_geometry(kde.summ)) %>%
    mutate_if(is.double, round, 3) %>% 
    mutate(Category = rev(c("90% - 100%", "70% - 89%", "50% - 69%", 
                            "30% - 49%", "1% - 29%"))) %>%
    dplyr::mutate(kernelPct = round(kde.total / sum(kde.total),4),
                  fittedPct = round(obs.total / sum(obs.total), 4))
  
  countComparisonsLong <- countComparisons %>% 
    gather(Variable, Value, kernelPct:fittedPct)
  # Bat plot of results
  REALTIVE_RISK_BARPLOT_COMPARE_plot <- ggplot(data=countComparisonsLong, aes(Category,Value)) +
    geom_bar(aes(fill = Variable), position = "dodge", stat="identity", color = NA) +
    scale_fill_viridis_d(name = " ",
                         labels=c("Risk & Protective\nFeatures", "Kernel Density")) +
    labs(x= "Predicted Risk Levels",
         y="Percent of Test Set Cases",
         title= "Goodness of Fit: Spatial Risk Model vs. Kernel Density hotspot",
         caption = "Figure 1.6 - Kernel density comparison") +
    plotTheme() +
    theme(axis.line = element_blank())
  ```
  
  
  
  
  
