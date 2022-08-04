#' @title Correct Degenerate NWM files
#' @description Degenerate = having lost the physical, mental, or moral qualities considered normal and desirable; showing evidence of decline. NWM files do not store extents or projections...
#' @param path path to NWM file
#' @param subds subdatasets to extract
#' @return SpatRast object
#' @export
#' @importFrom terra ext crs rast

correct_nwm_spatial = function(path, ...){

  template = list(
    ext = ext(
      -2303999.62876143,
      2304000.37123857,
      -1920000.70008381,
      1919999.29991619
    ),
    crs = 'PROJCS["Sphere_Lambert_Conformal_Conic",GEOGCS["GCS_Sphere",DATUM["D_Sphere",SPHEROID["Sphere",6370000.0,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-97.0],PARAMETER["standard_parallel_1",30.0],PARAMETER["standard_parallel_2",60.0],PARAMETER["latitude_of_origin",40.000008],UNIT["Meter",1.0]];-35691800 -29075200 126180232.640845;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision'
  )


  data = suppressWarnings({ rast(path, ...) })
  terra::ext(data) <-  template$ext
  terra::crs(data) <-  template$crs

  data
}

#' Aggregated CFE variables frome NWM
#' @param gpkg a geopackage with aggregation units
#' @param catchment_name the layer name of the aggregation units
#' @param flowline_name the layer name of the flowpath units
#' @param single_layer should only the top layer of a multilayer parameter be processed?
#' @param dir directory of expected nwm_grid_file (see hyAggregate::get_nwm_grids())
#' @param precision the precision of the computations
#' @return NULL
#' @export
#' @importFrom sf read_sf st_drop_geometry st_is_empty
#' @importFrom zonal weight_grid execute_zonal
#' @importFrom dplyr select mutate filter bind_cols rename inner_join group_by summarize across everything right_join
#' @importFrom tidyr unnest_longer
#' @importFrom RNetCDF open.nc var.get.nc close.nc
#' @importFrom stats weighted.mean complete.cases setNames
#' @importFrom logger log_info log_success
#' @importFrom terra sources


aggregate_cfe_noahowp = function(gpkg = NULL,
                                 dir = NULL,
                                 catchment_name = "aggregate_divides",
                                 flowline_name  = "aggregate_flowpaths",
                                 ID = "id",
                                 precision = 9,
                                 add_weights_to_gpkg = TRUE,
                                 add_to_gpkg = TRUE,
                                 overwrite = FALSE,
                                 log = TRUE) {

  if(!is.logical(log)){
    log_appender(appender_file(log))
    verbose = TRUE
  } else {
    log_appender(appender_console)
    verbose = log
  }

  lyr = "cfe_noahowp_attributes"

  .SD <- . <- .data <-  NULL

  if(is.null(gpkg)){
    stop('gpkg cannot be missing', call. = FALSE)
  }

  if(!overwrite){
   if(lyr %in% st_layers(gpkg)$name){
     return(gpkg)
   }
  }

  cats = read_sf(gpkg, catchment_name)

  # cats = filter(cats, st_is_empty(cats))

  if(is.null(dir)){
    stop("dir cannot be NULL")
  }

  hyaggregate_log("INFO", "Getting NWM grid data", verbose)

  data = get_nwm_grids(dir = dir, spatial = TRUE)

  hyaggregate_log("INFO", glue("Building weighting grid from: {terra::sources(data)[1]}"),  verbose)

  nwm_w_1000m = weight_grid(data, cats,  ID = ID, progress = FALSE)

  if(add_weights_to_gpkg){
    write_sf(nwm_w_1000m, gpkg, "nwm1km_weights")
  }

  soils_exe = list()

  ### soil_properties
  soil_mode_var = c("bexp", "IVGTYP", "ISLTYP")
  soil_gm_var   = c("dksat", "psisat")
  soil_mean_var = c("slope", "smcmax", "smcwlt", "refkdt", 'cwpvt', 'vcmx25', 'mp', 'mfsno')

  soils_exe[[1]] = zonal::execute_zonal(
    data = data,
    w = nwm_w_1000m,
    ID = ID,
    drop = ID,
    subds = grepl(paste0(soil_mode_var, collapse = "|"), names(data)),
    fun = "mode"
  )

  hyaggregate_log("INFO", glue('Getting mode: {paste(soil_mode_var, collapse = ", ")}'),  verbose)

  soils_exe[[2]] = execute_zonal(
    data = data,
    w = nwm_w_1000m,
    ID = ID,
    drop = ID,
    subds = grepl(paste0(soil_gm_var, collapse = "|"), names(data)),
    fun = zonal:::geometric_mean
  )

  hyaggregate_log("INFO", glue('Getting geometric_mean: {paste(soil_gm_var, collapse = ", ")}'),  verbose)

  soils_exe[[3]] = execute_zonal(
    data = data,
    w = nwm_w_1000m,
    ID = ID,
    drop = ID,
    subds = grepl(paste0(soil_mean_var, collapse = "|"), names(data)),
    fun = "mean"
  )

  hyaggregate_log("INFO", glue('Getting mean: {paste(soil_mean_var, collapse = ", ")}'),  verbose)

  exe <- cbind(soils_exe[[1]], soils_exe[[2]], soils_exe[[3]])

  exe[[ID]] = cats[[ID]]

  ####

  if (!is.null(flowline_name)) {

    crosswalk <- st_drop_geometry(read_sf(gpkg, flowline_name))

    crosswalk = select(crosswalk, .data$id, .data$member_comid) %>%
      mutate(comid = strsplit(.data$member_comid, ",")) %>%
      unnest_longer(col = c("comid")) %>%
      mutate(comid = as.integer(.data$comid)) %>%
      filter(!duplicated(.))

    gwnc = open.nc(get_nwm_grids(dir, spatial = FALSE))

    on.exit(close.nc(gwnc))

    vars      = c("Area_sqkm", "ComID", "Coeff",  "Zmax")
    vars_mode = c("Area_sqkm", "ComID", "Expon")


    gwparams_means = suppressMessages({
      lapply(vars, function(x)
        var.get.nc(gwnc, x)) %>%
        bind_cols() %>%
        setNames(vars) %>%
        rename(comid = .data$ComID) %>%
        mutate(comid = as.integer(.data$comid)) %>%
        inner_join(select(crosswalk, .data$id, .data$comid), by = 'comid') %>%
        filter(complete.cases(.)) %>%
        filter(!duplicated(.)) %>%
        group_by(.data$id) %>%
        summarize(across(everything(), ~ round(
          weighted.mean(.x, w = .data$Area_sqkm, na.rm = TRUE),
          precision
        ))) %>%
        select(-.data$comid, -.data$Area_sqkm)
    })

    hyaggregate_log("INFO", glue('Getting weighted mean: {paste(vars, collapse = ", ")}'),  verbose)

    getmode = function(x){
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }

    gwparams_mode = suppressMessages({
      lapply(vars_mode, function(x)
        var.get.nc(gwnc, x)) %>%
        bind_cols() %>%
        setNames(vars_mode) %>%
        rename(comid = .data$ComID) %>%
        inner_join(select(crosswalk, .data$id, .data$comid), by = 'comid') %>%
        filter(complete.cases(.)) %>%
        filter(!duplicated(.)) %>%
        group_by(.data$id) %>%
        summarize(Expon = getmode(floor(.data$Expon)))
    })

    hyaggregate_log("INFO", glue('Getting mode: {paste(vars_mode, collapse = ", ")}'),  verbose)

    traits = left_join(gwparams_means, gwparams_mode, by = ID) %>%
      setNames(c('id', paste0('gw_', names(.)[-1]))) %>%
      right_join(exe, by = ID)
  }

  names(traits) = gsub("_Time=1", "", names(traits))

  if(add_to_gpkg){
    write_sf(traits, gpkg, lyr)
    return(gpkg)
  } else {
    return(traits)
  }

  log_appender(appender_console)

}


