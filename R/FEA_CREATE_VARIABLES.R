##################################################################
#### Create Variables
##################################################################

library("sf")            # Spatial data objects and methods
library("spatstat")      # KDE and other spatial functions
library("raster")        # cell-based spatial operations
library("tidyverse")     # data manipulation framework
library("lubridate")     # Power tools for handling dates
library("tidycensus")    # Get census data

# requires all data in *.csv or *.xls files containing coordinate field names "X" and "Y"
# `crs` in the call to `st_as_sf()` needs to be set to the ESPG code of your data projection
# `base_dir` file path and many feature names are specified for the current project.

##1.1 Global Variables
mapviewOptions(basemaps = c("Stamen.TonerLite", "OpenStreetMap.DE"))
base_dir = "C:/projects/PAP_Virginia"


##2.1 Load Data
files <-list.files(file.path(base_dir,"/data"), pattern = "*\\.xls$|*\\.csv$")
var_list <- vector(mode = "list")
var_names <- NULL
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
      st_as_sf(., coords = c("X", "Y"), crs = 102747)
  }, error = function(e){
    cat(filename, "error = ",e$message,"\n")
    return(e)
  }
  )
  if(!inherits(sf_i, "error")){
    var_list[[length(var_list)+1]] <- sf_i
    var_names[length(var_list)] <- filename
  }
}
names(var_list) <- var_names


## 2.2 Protective Variables
## These are specified to the current project
protective_class <- c("CommunityCenters","FireStations",
                      "HomelessShelters","Libraries","Parks","PointsOfInterest",
                      "PoliceStations","PublicSchools","ResourceOASIS","SNAP_WIC",
                      "VotingStations")

protective_vars <- list()
for(i in seq_along(protective_class)){
  dat <- var_list[[protective_class[i]]] %>%
    mutate(feature_name = protective_class[i],
           class = "protective") %>%
    dplyr::select(feature_name, class)
  protective_vars[[i]] <- dat
}

Businesses_protective <- var_list[["BusinessProject"]] %>%
  filter(Classification == "PROTECTIVE") %>%
  mutate(feature_name = "BusinessProject",
         class = "protective") %>%
  dplyr::select(feature_name, class)

protective_vars[[length(protective_vars)+1]] <- Businesses_protective
protective_vars <- do.call(rbind, protective_vars)

var_list[["Protective"]] <- protective_vars


## 2.3 Risk variables
risk_class <- c("BusStops")
risk_vars <- list()
for(i in seq_along(risk_class)){
  dat <- var_list[[risk_class[i]]] %>%
    mutate(feature_name = risk_class[i],
           class = "risk") %>%
    dplyr::select(feature_name, class)
  risk_vars[[i]] <- dat
}

Business_risk <- var_list[["BusinessProject"]] %>%
  filter(Classification == "RISK") %>%
  mutate(feature_name = "BusinessProject",
         class = "risk") %>%
  dplyr::select(feature_name, class)
risk_vars[[length(risk_vars)+1]] <- Business_risk

CrimeData_risk <- var_list[["CrimeData"]] %>%
  filter(OFFENSE %in% c('DRUG/NARCOTIC VIOLATION','SIMPLE ASSAULT, DOMESTIC',
                        'Runaway','AGGRAVATED ASSAULT DOMESTIC')) %>%
  mutate(feature_name = "CrimeData",
         class = "risk") %>%
  dplyr::select(feature_name, class)
risk_vars[[length(risk_vars)+1]] <- CrimeData_risk

Violations_III_ks_risk <- var_list[["Violations_III_ks"]] %>%
  filter(CodeDsrp %in% c('General Violations','Unsafe Structure',
                         'Unfit Structure')) %>%
  mutate(feature_name = "Violations_III_ks_risk",
         class = "risk") %>%
  dplyr::select(feature_name, class)
risk_vars[[length(risk_vars)+1]] <- Violations_III_ks_risk

risk_vars <- do.call(rbind, risk_vars)

var_list[["Risk"]] <- risk_vars

# Violations Individual (min n = 100)
violations_list <- get_individual_features(var_list[["Violations_III_ks"]], "CodeNbr", "VIO", 100)
var_list <- c(var_list, violations_list)

# CrimeData Individual (min n = 100)
cimedata_list <- get_individual_features(var_list[["CrimeData"]], "OFFENSE", "CRIME", 100)
var_list <- c(var_list, cimedata_list)

# BusinessProject Individual (min n = 10)
businessproject_list <- get_individual_features(var_list[["BusinessProject"]], "BUSTYP", "BUSI", 10)
var_list <- c(var_list, businessproject_list)

## FILTER CPS_ACCEPTED TO REMOVE SUPES AND "unconfirmed"
CPS_filtered<- var_list[["CPS_Accepted"]] %>%
  filter(Disposition != "Unfounded - lack of evidence") %>% 
  group_by(Match_addr, RDate) %>% 
  mutate(n = row_number()) %>% 
  filter(n == 1) %>% 
  dplyr::select(-n) %>% 
  ungroup()
var_list[["CPS_Accepted"]] <- CPS_filtered

# CPS_ACCEPTED Individual (min n = 10)
CPS_Fatality_list <- get_individual_features(var_list[["CPS_Accepted"]], "ChildFatality", "CPS_Fatality", 10)
var_list <- c(var_list, CPS_Fatality_list)

CPS_Near_Fatality_list <- get_individual_features(var_list[["CPS_Accepted"]], "NearFatality", "CPS_NearFatality", 10)
var_list <- c(var_list, CPS_Near_Fatality_list)

CPS_SexualA_list <- get_individual_features(var_list[["CPS_Accepted"]], "SexualA", "CPS_SexualA", 10)
var_list <- c(var_list, CPS_SexualA_list)








