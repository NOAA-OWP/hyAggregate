#' @title Build Headwater Collapse Table
#' @description  Identifies small (pathlength or area) headwater catchments and returns a data.frame
#' with the current ID and the feature ID it should collapse into (becomes).
#' Headwaters are segments in which there are no inflows (!ID %in% toID).
#' @param network_list  a list containing flowpath and catchment `sf` objects
#' @param min_area_sqkm The minimum allowable size of the output hydrofabric catchments
#' @param min_length_km The minimum allowable length of the output hydrofabric flowlines
#' @return a 2 column data.frame with {id, becomes}
#' @export
#' @importFrom dplyr mutate filter group_by mutate ungroup distinct
#' @importFrom sf st_set_geometry st_geometry st_buffer st_intersects
#' @importFrom nhdplusTools get_node

build_collapse_table = function(network_list,
                                min_area_sqkm  = 3,
                                min_length_km  = 1) {
  bad =  mutate(
    network_list$flowpaths,
    hw = ifelse(!id %in% toid, TRUE, FALSE),
    small = areasqkm < min_area_sqkm | lengthkm < min_length_km
  ) |>
    filter(hw, small) |>
    # needed for VPU02 - DO NOT REMOVE!!!!
    st_cast("MULTILINESTRING")

  outlets =  st_buffer(st_set_geometry(bad, st_geometry(get_node(bad, "end"))), 1)

  emap = st_intersects(outlets, network_list$flowpaths)

  df = data.frame(
    id       = rep(outlets$id, times = lengths(emap)),
    toid     = rep(outlets$toid, times = lengths(emap)),
    touches  = network_list$flowpaths$id[unlist(emap)],
    id_type  = rep(outlets$type, times = lengths(emap))
  ) %>%
    filter(!.data$id == .data$touches) %>%
    filter(is.na(id_type)) %>%
    group_by(id) |>
    mutate(becomes = ifelse(any(toid == touches), toid, touches)) |>
    ungroup()   |>
    distinct(id, becomes) %>%
    filter(!id %in% becomes)

  df
}

#' Collapse Headwaters
#' @description  This function identifies small (pathlength or area) headwater catchments and
#' collapses them into the existing network until none remain.
#' Headwaters are those segments in which there are no inflows (!ID %in% toID).
#' @param network_list  a list containing flowpath and catchment `sf` objects
#' @param min_area_sqkm The minimum allowable size of the output hydrofabric catchments
#' @param min_length_km The minimum allowable length of the output hydrofabric flowlines
#' @param verbose should messages be emitted?
#' @param cache_file If not NULL results will be written to a provide path (.gpkg)
#' @return a list containing flowpath and catchment `sf` objects
#' @importFrom dplyr filter left_join mutate bind_rows
#' @importFrom nhdplusTools get_node

collapse_headwaters = function(network_list,
                               min_area_sqkm  = 3,
                               min_length_km  = 1,
                               verbose = TRUE,
                               cache_file = NULL) {

  start <- nrow(network_list$flowpaths)

  mapping_table <- build_collapse_table(network_list, min_area_sqkm, min_length_km)

  count = 0

  while (nrow(mapping_table) > 0) {

    count = count + 1

    hyaggregate_log("INFO", glue("Collapsing: {nrow(mapping_table)} features (round {count})"), verbose)

    fl = filter(network_list$flowpaths, !id %in% mapping_table$id)

    cat = filter(network_list$catchments,
                 id %in% c(mapping_table$id, mapping_table$becomes)) |>
      left_join(mapping_table, by = "id") |>
      mutate(id = ifelse(is.na(becomes), id, becomes)) |>
      geos_union_polygon_hyaggregate("id") |>
      bind_rows(filter(
        network_list$catchments,
        !id %in% c(mapping_table$id, mapping_table$becomes)
      ))

    network_list = prepare_network(list(flowpaths = fl, catchments = cat))

    mapping_table = build_collapse_table(network_list, min_area_sqkm, min_length_km)
  }

  hyaggregate_log("SUCCESS", glue("Collapsed {start - nrow(network_list$flowpaths)} features."), verbose)

  if (!is.null(cache_file)) {
    write_hydrofabric_package(network_list, cache_file, "hw_agg_cat", "hw_agg_fp", verbose)
  }

  return(network_list)
}


#' @title Prepare Hydrologic Network
#' @details This function adds an area, length, hydrosequence, streamorder and contributing drainage area
#' metric to the flowpath list element of network_list.
#' @details tot_drainage_areasqkm can only be added when there are no NA areas
#' @param network_list a list with flowpath and catchment data
#' @return a list containing flowpath and catchment `sf` objects
#' @export
#' @importFrom nhdplusTools get_streamorder calculate_total_drainage_area
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @noRd

prepare_network = function(network_list) {

  names(network_list$flowpaths)  = tolower(names(network_list$flowpaths))
  names(network_list$catchments) = tolower(names(network_list$catchments))

  # Add a hydrosequence to the flowpaths
  network_list$flowpaths = add_hydroseq(flowpaths = network_list$flowpaths)

  # Add area and length measures to the network list
  network_list = add_measures(network_list$flowpaths, network_list$catchments)

  if (!any(is.na(network_list$flowpaths$areasqkm))) {
    network_list$flowpaths$tot_drainage_areasqkm = calculate_total_drainage_area(st_drop_geometry(
      select(
        network_list$flowpaths,
        ID = id,
        toID = toid,
        area = areasqkm
      )
    ))
  }

  network_list$flowpaths$order = network_list$flowpaths %>%
    st_drop_geometry() %>%
    flush_prefix(c("id", "toid")) %>%
    select(ID = id, toID = toid) %>%
    get_streamorder()

  check_network_validity(network_list)

}

#' @title Aggregate Network to Distribution
#' @description This function aggregated a network to a desired size distribution while
#' enforcing minimum flowpath legnths and catchment area. Addtionally a set of explicit nexus
#' locations can be provided over which the netwrok cannot be aggregated.
#' @param gf A path to a hydrofabric file (see US reference fabric)
#' @param ideal_size_sqkm The ideal size of catchments (default = 10 sqkm)
#' @param min_length_km The minimum allowable length of flowpath features (default = 1 km)
#' @param min_area_sqkm The minimum allowable area of catchment features (default = 3 sqkm)
#' @param outfile of not NULL, where to write the output files
#' @param verbose print status updates. Default = TRUE
#' @param overwrite overwrite existing gf file. Default is FALSE
#' @param nexus_topology should a hy-features cat-->nex topology be enforced? default = TRUE
#' @param nexus_locations a data.frame with columns specifiying the ID, and the nexus type.
#' @param log a filepath to write to or booleen (TRUE = print to console; FALSE = no messages)
#' @return if outfile = TRUE, a file path, else a list object
#' @export
#' @importFrom sf st_transform read_sf st_set_crs write_sf st_layers
#' @importFrom dplyr left_join filter
#' @importFrom nhdplusTools get_sorted calculate_total_drainage_area get_streamorder
#' @importFrom logger log_success log_info log_appender appender_file appender_console

aggregate_network_to_distribution = function(gpkg = NULL,
                                             catchment_name = "refactored_divides",
                                             flowpath_name = "refactored_flowpaths",
                                             ideal_size_sqkm = 10,
                                             min_length_km = 1,
                                             min_area_sqkm  = 3,
                                             nexus_topology = TRUE,
                                             routelink_path = NULL,
                                             nexus_locations = NULL,
                                             outfile = NULL,
                                             log = TRUE,
                                             overwrite = FALSE,
                                             cache_file = NULL) {

  if(!is.logical(log)){
    log_appender(appender_file(log))
    verbose = TRUE
  } else {
    log_appender(appender_console)
    verbose = log
  }

  hyaggregate_log("INFO", glue("ideal_size_sqkm --> {ideal_size_sqkm}"), verbose)
  hyaggregate_log("INFO", glue("min_length_km --> {min_length_km}"), verbose)
  hyaggregate_log("INFO", glue("min_area_sqkm --> {min_area_sqkm}"), verbose)

  if(!is.null(outfile)){
    hyaggregate_log("INFO", glue("outfile --> {outfile}\n"), verbose)
  }

  hyaggregate_log("INFO", glue("\n--- Read in data from {gpkg} ---\n"), verbose)

  if (!is.null(gpkg)) {

    if (file.exists(outfile) & overwrite) {
      unlink(outfile)
    } else if (file.exists(outfile)) {
      hyaggregate_log("WARN", glue("{outfile} already exists and overwrite is FALSE"), verbose)
      return(outfile)
    }

    network_list = read_hydrofabric_package(gpkg,
                                            catchment_name = catchment_name,
                                            flowpath_name = flowpath_name,
                                            crs = 5070)

    if (!is.null(nexus_locations)) {
      network_list$flowpaths  = left_join(network_list$flowpaths,  nexus_locations, by = "ID")
      network_list$catchments = left_join(network_list$catchments, nexus_locations, by = "ID")
    } else {
      network_list$flowpaths$type  = NA
      network_list$flowpaths$value  = NA
    }
  }

  # Create a network list and check if DAG and connected

  network_list$flowpaths  =   network_list$flowpaths[!duplicated(network_list$flowpaths), ]
  network_list$catchments =   network_list$catchments[!duplicated(network_list$catchments), ]

  network_list <- drop_extra_features(prepare_network(network_list), verbose)

  if (!is.null(cache_file)) {
    write_hydrofabric_package(network_list, cache_file, "base_cat", "base_fp", verbose)
  }

  # Perform first aggregation step!
  hyaggregate_log("INFO", "\n---  Aggregate Along Mainstem ---\n", verbose)

  network_list = aggregate_along_mainstems(network_list,
                                           ideal_size_sqkm,
                                           min_area_sqkm,
                                           min_length_km,
                                           verbose = verbose,
                                           cache_file = cache_file)

  hyaggregate_log("INFO", "\n--- Collapse Network Inward ---\n", verbose)

  network_list  = collapse_headwaters(network_list,
                                      min_area_sqkm,
                                      min_length_km,
                                      verbose = verbose,
                                      cache_file = cache_file)

  if (nexus_topology) { network_list = apply_nexus_topology(network_list, verbose = verbose) }


  if (!is.null(routelink_path)) {
    hyaggregate_log("INFO", glue("\n--- Creating Routing Table based on: {routelink_path} ---\n "), verbose)

    network_list$waterbody_attributes <-
      length_average_routelink(flowpaths = network_list$flowpaths,
                               rl_path = routelink_path)
  }


  if (is.null(outfile)) {
    return(network_list)
  } else {

    write_hydrofabric_package(network_list,
                              outfile,
                              catchment_name  = "aggregate_divides",
                              flowpath_name   ="aggregate_flowpaths",
                              verbose)

    if (!is.null(network_list$nex)) { write_sf(network_list$nex, outfile, "nexus") }

    if (!is.null(network_list$flowpath_edge_list)) {
      write_sf(network_list$flowpath_edge_list, outfile, "flowpath_edge_list")
    }

    if (!is.null(network_list$waterbody_attributes)) {
      write_sf(network_list$waterbody_attributes, outfile, "flowpath_attributes")
    }

    hyaggregate_log("SUCCESS", "Done!", verbose)
  }

  log_appender(appender_console)
}


#' Aggregate along network mainstems
#' @description Given a set of ideal catchment sizes, plus the
#' minimum allowable catchment size and segment length, aggregate the network along mainstems.
#' @param network_list a list containing flowline and catchment `sf` objects
#' @param ideal_size The ideal size of output hydrofabric catchments
#' @param min_area_sqkm The minimum allowable size of the output hydrofabric catchments
#' @param min_length_km The minimum allowable length of the output hydrofabric flowlines
#' @param term_cut cutoff integer to define terminal IDs
#' @return a list containing aggregated and validated flowline and catchment `sf` objects
#' @export
#' @importFrom dplyr filter group_by arrange mutate ungroup select
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr %>% cur_group_id n
#' @importFrom logger log_info

aggregate_along_mainstems = function(network_list,
                                     ideal_size_sqkm,
                                     min_area_sqkm,
                                     min_length_km,
                                     verbose = TRUE,
                                     cache_file = NULL) {

  # Requires having a set, id, toid, levelpathid, hydroseq, member_comid, type, n
  index_table = network_list$flowpaths %>%
    st_drop_geometry() %>%
    group_by(.data$levelpathid) %>%
    arrange(.data$hydroseq) %>%
    mutate(
      ind = cs_group(
        .data$areasqkm,
        .data$lengthkm,
        .data$type,
        ideal_size_sqkm,
        min_area_sqkm,
        min_length_km
      )
    ) %>%
    ungroup()   %>%
    group_by(.data$levelpathid, .data$ind) %>%
    mutate(set = cur_group_id(), n = n()) %>%
    ungroup() %>%
    select(set, id, toid, levelpathid, hydroseq, member_comid, type, value, n)

  v = aggregate_sets(
    flowpaths  = filter(network_list$flowpaths,  id %in% index_table$id),
    cat        = filter(network_list$catchments, id %in% index_table$id),
    index_table
  ) %>%
    prepare_network()

  hyaggregate_log("SUCCESS",
                  glue("Merged to idealized catchment size of {ideal_size_sqkm} sqkm: {nrow(network_list$flowpaths) - nrow(v$flowpaths)} features removed"),
                  verbose)

  if (!is.null(cache_file)) {
    write_hydrofabric_package(v, cache_file, "mainstem_agg_cat", "mainstem_agg_fp", verbose)
  }

  return(v)

}
