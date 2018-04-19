utils::globalVariables(c("out",
                         "coordinates"))
#' Download and clean daily climate records from the Global Historical Climate Database.
#'
#' @param region A spatial object defining the region
#' @param label A label for the region
#' @param calibration.years The years for defining climatology
#' @param google_maps_elevation_api_key the users Google Maps Elevation API key
#' @param raw_dir a directory in which to store the raw GHCN downloads
#' @param derived_dir a directory in which to store the processed GHCN output
#' @param force.redo Should the GHCN data be re-processed?
#'
#' @return A clean set of GHCN data
#' @export
prepare_ghcn <- function(region,
                         label,
                         calibration.years,
                         google_maps_elevation_api_key,
                         raw_dir = "./data/raw_data/ghcn/",
                         derived_dir = "./data/derived_data/ghcn/",
                         force.redo = FALSE) {
  # Keep only the clean stations
  if (!force.redo & file.exists(paste0(derived_dir,"/ghcn_data_final.Rds"))) {
    GHCN.data.final <- readr::read_rds(paste0(derived_dir,"/ghcn_data_final.Rds"))

    return(GHCN.data.final)
  }

  ## GHCN DATA ##
  # These are the GHCN stations that will be used to calibrate the interpolation model.
  GHCN.data <- FedData::get_ghcn_daily(
    template = region,
    label = label,
    elements = c("tmin", "tmax", "prcp"),
    years = calibration.years,
    raw.dir = raw_dir,
    extraction.dir = derived_dir,
    standardize = TRUE,
    force.redo = force.redo
  )

  GHCN.stations <- GHCN.data[[1]] # The spatial station data
  GHCN.data <- GHCN.data[[2]] # The actual temperature data
  GHCN.stations <- GHCN.stations[!duplicated(GHCN.stations$ID), ] # Remove duplicated station location data
  GHCN.data <- GHCN.data[as.character(GHCN.stations$ID)] # Ensure that the order of the two datasets are the same

  # Some stations have the same "location", but not the same data. Remove these points.
  nonduplicates <- !duplicated(sp::coordinates(GHCN.stations)) * !duplicated(sp::coordinates(GHCN.stations), fromLast = T)
  GHCN.stations <- GHCN.stations[nonduplicates, ]
  GHCN.data <- GHCN.data[nonduplicates]

  # Remove empty records
  GHCN.data <- GHCN.data[!(GHCN.data %>% sapply(length) < 2)]

  ## Clean the GHCN data
  if (force.redo | !file.exists(paste0(derived_dir,"/ghcn_data_clean.Rds"))) {

    # Get run length encoding of missing data
    all.rles <- do.call(c, lapply(GHCN.data, function(test) {
      tryCatch(getMissingRLE(test), error = function(e) NULL)
    }))

    # Calculate the cutoff length of data gaps.
    # Years with gaps longer than "cutoff" will be dropped
    cutoff <- calcGapCutoff(
      rleVector = all.rles,
      pLevel = 0.95
    )

    GHCN.data.clean <- lapply(GHCN.data, function(station.data) {
      # cat(station.data$TMAX$STATION[[1]],"\n")
      test <- ghcnCleaner(
        data.list = station.data,
        min.years = 10,
        year.range = calibration.years,
        na.cutoff = cutoff
      )
    })
    names(GHCN.data.clean) <- names(GHCN.data)
    GHCN.data.clean <- GHCN.data.clean[!sapply(GHCN.data.clean, is.null)]
    readr::write_rds(GHCN.data.clean,
      path = paste0(derived_dir,"/ghcn_data_clean.Rds"),
      compress = "gz"
    )
  }
  GHCN.data.clean <- readr::read_rds(paste0(derived_dir,"/ghcn_data_clean.Rds"))


  # Keep only the clean stations
  if (force.redo | !file.exists(paste0(derived_dir,"/ghcn_data_final.Rds"))) {
    GHCN.stations <- GHCN.stations[GHCN.stations$ID %in% names(GHCN.data.clean), ]

    # Get the station elevations using google maps API
    GHCN.stations$elevation <- rgbif::elevation(latitude = GHCN.stations@coords[, 2],
                                                longitude = GHCN.stations@coords[, 1],
                                                key = google_maps_elevation_api_key)[, "elevation"]

    # Get all stations averages over the calibration period
    GHCN.data.averages <- lapply(GHCN.data.clean, function(station) {
      return(lapply(station, calcDailyMeanSD))
    })

    # Create a final dataset
    GHCN.data.final <- list(
      spatial = GHCN.stations,
      weather = GHCN.data.clean,
      climatology = GHCN.data.averages
    )

    readr::write_rds(GHCN.data.final,
      path = paste0(derived_dir,"/ghcn_data_final.Rds"),
      compress = "gz"
    )
  }
  GHCN.data.final <- readr::read_rds(paste0(derived_dir,"/ghcn_data_final.Rds"))

  return(GHCN.data.final)
}

# climate_plotter(data = GHCN.data.final, station = "AM000037789", element = "TMIN")
