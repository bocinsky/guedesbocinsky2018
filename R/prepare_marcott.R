prepare_marcott <- function(calibration.years){
  # Get Mann et al. 2008 infilled instrumental temperature data, 
  # and extract the 1961--1990 period.
  message("Downloading Mann et al. 2008 infilled instrumental temperature data.")
  dir.create(out("DATA/mann2008/"), showWarnings = FALSE, recursive = TRUE)
  FedData::download_data(url="ftp://ftp.ncdc.noaa.gov/pub/data/paleo/contributions_by_author/mann2008/instrument.zip",destdir=out("DATA/mann2008/"))
  unzip(out("DATA/mann2008/instrument.zip"), exdir=out("DATA/mann2008/"))
  # Get just the Northern Hemisphere HAD CRU V3 data.
  iHAD_NH_reform <- readr::read_table(out("DATA/mann2008/iHAD_NH_reform"), col_names = c("Year","Temperature"))
  iHAD_NH_reform <- iHAD_NH_reform[iHAD_NH_reform$Year %in% calibration.years,]
  message("Calculating standard deviation for 1961--1990 calibration period.")
  calib.sd <- sd(iHAD_NH_reform$Temperature)
  
  # Get the Marcott et al. 2013 Northern Hemisphere temperature reconstruction, 
  # which is referenced to the Mann et al. 2008 infilled instrumental average from 1961--1990
  message("Downloading Marcott et al. 2013 Northern Hemisphere temperature reconstruction.")
  dir.create(out("DATA/marcott2013/"), showWarnings = FALSE, recursive = TRUE)
  FedData::download_data(url="http://www.sciencemag.org/content/suppl/2013/03/07/339.6124.1198.DC1/Marcott.SM.database.S1.xlsx",destdir=out("DATA/marcott2013/"))
  marcott2013 <- readxl::read_excel(out("DATA/marcott2013/Marcott.SM.database.S1.xlsx"),
                                    sheet="TEMPERATURE STACKS",
                                    range = "U3:W571",
                                    col_types = "numeric")
  names(marcott2013) <- c("YearBP","Temperature","Uncertainty")
  message("Transforming temperature deviations to standard scores.")
  marcott2013 <- marcott2013[(marcott2013$YearBP <= 5510)&(!is.na(marcott2013$Temperature)),]
  marcott2013$Z_Lower <- (marcott2013$Temperature-marcott2013$Uncertainty)/calib.sd
  marcott2013$Z <- marcott2013$Temperature/calib.sd
  marcott2013$Z_Upper <- (marcott2013$Temperature+marcott2013$Uncertainty)/calib.sd
  
  # Write the final Marcott dataset for the reconstruction
  message("Exporting the standard score data for later use.")
  readr::write_csv(marcott2013,out("DATA/marcott2013/MARCOTT2013_Z.csv"))
  return(marcott2013)
}

