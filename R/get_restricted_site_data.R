#' Get the restricted archaeological site and chronometric data from tDAR
#'
#' @param tdar_user Your tDAR user name
#' @param tdar_password Your tDAR password
#'
#' @return A data.frame of site and chronometric data
#' @export
get_restricted_site_data <- function(tdar_user = Sys.getenv('tdar_un'),
                                     tdar_password = Sys.getenv('tdar_pw')){

  # Log in to tDAR API
  tdar::tdar_login(tdar_user = tdar_user,
                   tdar_password = tdar_password)

  on.exit(
    # Log out of tDAR
    tdar::tdar_logout()
  )

  data_path <- tdar::tdar_download_resource(428089,
                                            out_dir = tempdir(),
                                            overwrite = TRUE)

  # httr::GET("https://core.tdar.org/filestore/download/428089",
  #           httr::write_disk(paste0(tempdir(),"/sites.xlsx"),
  #                            overwrite = TRUE))

  # Read in the site/chronometric data
  out <- data_path %>% # This gets the file path for the resource ID
    readxl::read_xlsx(sheet = "chronometric_data",
                      col_types = c("text",
                                    "numeric",
                                    "numeric",
                                    "numeric",
                                    "logical",
                                    "text",
                                    "text",
                                    "text",
                                    "logical",
                                    "numeric",
                                    "numeric",
                                    "text",
                                    "numeric",
                                    "numeric",
                                    "text",
                                    rep("logical",18))
    ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(`Exclude?` = ifelse(is.na(`Exclude?`), FALSE, `Exclude?`)) %>%
    dplyr::filter(!is.na(Site),
                  !is.na(Longitude),
                  !is.na(Latitude),
                  !(is.na(`14C age BP`) & is.na(`Age range lower (BP)`)),
                  !is.na(`14C date on cereal?`),
                  !`Exclude?`) %>%
    dplyr::select(-`Exclude?`, -Notes) %>%
    dplyr::group_by(Site,Period)

labs <- data_path %>% # This gets the file path for the resource ID
  readxl::read_xlsx(sheet = "14C Labs")

  refs <- data_path %>% # This gets the file path for the resource ID
  readxl::read_xlsx(sheet = "References")



  unlink(data_path)

  return(list(data = out,
              labs = labs,
              refs = refs))

}
