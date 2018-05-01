utils::globalVariables(c("VEGE_ID",
                         "Vegetation_ID",
                         "Vegetation_formation_and_sub-formation",
                         "Vegetation_type",
                         "crop"))
#' Get the restricted vegetation_atlas data from tDAR
#'
#' @param tdar_user Your tDAR user name
#' @param tdar_password Your tDAR password
#'
#' @return A Simple Features Data Frame of vegetation atlas data
#' @export
get_restricted_vegetation_atlas_data <- function(tdar_user = Sys.getenv('tdar_un'),
                                     tdar_password = Sys.getenv('tdar_pw')){

  # Log in to tDAR API
  tdar::tdar_login(tdar_user = tdar_user,
                   tdar_password = tdar_password)

  on.exit(
    # Log out of tDAR
    tdar::tdar_logout()
  )

  httr::GET("https://core.tdar.org/filestore/download/442484",
            httr::write_disk(paste0(tempdir(),"/vegetation.zip"),
                                                         overwrite = TRUE))
  paste0(tempdir(),"/vegetation.zip") %>%
    utils::unzip(exdir = paste0(tempdir(),"/vegetation/"))

  out <- sf::st_read(paste0(tempdir(),"/vegetation/vegetation.shp")) %>%
    dplyr::select(VEGE_ID) %>%
    dplyr::rename(`Vegetation_ID` = VEGE_ID) %>%
    dplyr::mutate(`Vegetation_ID` = `Vegetation_ID` %>%
                    as.character() %>%
                    stringr::str_replace_all("[0-9]{1,3}",
                                             function(x){
                                               stringr::str_pad(x,
                                                                width = 3,
                                                                side = "left",
                                                                pad = "0")
                                             })) %>%
    dplyr::group_by(`Vegetation_ID`) %>%
    dplyr::summarise() %>%
    sf::st_cast()

codes <- readr::read_csv("../inst/vegetation_codes.csv") %>%
  dplyr::mutate(`Vegetation_ID` = `Vegetation_ID` %>%
                  stringr::str_replace_all("[0-9]{1,3}",
                                           function(x){
                                             stringr::str_pad(x,
                                                              width = 3,
                                                              side = "left",
                                                              pad = "0")
                                           }))

crops <- list(tibble::tibble(crop = "wheat",
                             `Vegetation_ID` = c("556a","556b","557", "558", "558a",
                                                 "558e", "559", "560", "560a", "560e",
                                                 "560d", "560f", "560g", "560i", "561",
                                                 "562", "564", "564a", "564b", "565",
                                                 "566", "568", "568g", "569", "570",
                                                 "570a", "419a+556", "550+83", "550+069a", "568+568c")),
              tibble::tibble(crop = "barley",
                             `Vegetation_ID` = c("548")),
              tibble::tibble(crop = "rice",
                             `Vegetation_ID` = c("568+568c", "549a", "550a", "555", "558c",
                                                 "560i", "564", "564a", "564c", "565",
                                                 "566", "567", "568", "568g", "568h",
                                                 "568i", "569", "570", "570a", "571",
                                                 "571g", "571i", "571j", "572", "572a",
                                                 "572l", "573", "568+568c")),
              tibble::tibble(crop = "millet",
                             `Vegetation_ID` = c("560","419a+556", "556"))#,
              # tibble::tibble(crop = "millet and sorghum",
              #                Vegetation_ID = c("560","419a+556", "556", "550+069","550",
              #                                  "563", "570", "550+083","550+069a")),
              # tibble::tibble(crop = "potato",
              #                Vegetation_ID = c("547","551"))
              ) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(codes)


  out %<>%
    dplyr::right_join(crops) %>%
    dplyr::mutate(`Vegetation_formation_and_sub-formation` = factor(`Vegetation_formation_and_sub-formation`),
                  Vegetation_type = factor(Vegetation_type)) %>%
    dplyr::filter(!is.na(crop)) %>%
   sf::st_transform(4326) %>%
    dplyr::arrange(Vegetation_ID)

   # out %>%
   #   sf::st_transform(4326) %>%
   #   sf::st_write(dsn = "/Users/bocinsky/Desktop/china_vegetation_atlas.gpkg",
   #                delete_dsn = TRUE)

  unlink(paste0(tempdir(),"/vegetation.zip"))
  unlink(paste0(tempdir(),"/vegetation"),
         recursive = TRUE)

  return(out)

}
