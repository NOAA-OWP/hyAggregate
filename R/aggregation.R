#' Build Collapse Table
#' @description  This function identifies small (pathlength or area) headwater catchments and returns a data.frame
#' with the current ID and the the feature ID it should collapse into (becomes).
#' Headwaters are those segments in which there are no inflows (!ID %in% toID).
#' @param network_list  a list containing flowpath and catchment `sf` objects
#' @param min_area_sqkm The minimum allowable size of the output hydrofabric catchments
#' @param min_length_km The minimum allowable length of the output hydrofabric flowlines
#' @return a 2 column data.frame with {id,becomes}
#' @export
#' @importFrom dplyr mutate filter group_by mutate ungroup distinct
#' @importFrom sf st_set_geometry st_geometry st_buffer st_intersects
#' @importFrom nhdplusTools get_node

build_collapse_table = function(network_list,
                                min_area_sqkm = 3,
                                min_length_km  = 1){

  bad =  mutate(network_list$flowpaths, hw = ifelse(!id %in% toid, TRUE, FALSE),
                small = areasqkm < min_area_sqkm | lengthkm < min_length_km) |>
    filter(hw, small)

  outlets =  st_buffer(st_set_geometry(bad, st_geometry(get_node(bad, "end"))), 1)

  emap = st_intersects(outlets, fl)

  data.frame(
    id       = rep(outlets$id, times = lengths(emap)),
    toid     = rep(outlets$toid, times = lengths(emap)),
    touches  = fl$id[unlist(emap)],
    id_type  = rep(outlets$type, times = lengths(emap))
  ) %>%
    filter(!.data$id == .data$touches) %>%
    filter(is.na(id_type)) %>%
    group_by(id) |>
    mutate(becomes = ifelse(any(toid == touches), toid, touches)) |>
    ungroup()   |>
    distinct(id, becomes) %>%
    filter(!id %in% becomes)
}

#' Collapse Headwaters
#' @description  This function identifies small (pathlength or area) headwater catchments and
#' collapses them into the existing network uptil none remain.
#' Headwaters are those segments in which there are no inflows (!ID %in% toID).
#' @param network_list  a list containing flowpath and catchment `sf` objects
#' @param min_area_sqkm The minimum allowable size of the output hydrofabric catchments
#' @param min_length_km The minimum allowable length of the output hydrofabric flowlines
#' @param verbose should messages be emitted?
#' @return a list containing flowpath and catchment `sf` objects
#' @export
#' @importFrom logger log_info log_success
#' @importFrom dplyr filter left_join mutate bind_rows
#' @importFrom nhdplusTools get_node

collapse_hw = function(network_list,
                       min_area_sqkm = 3,
                       min_length_km  = 1,
                       verbose = TRUE){

  start = nrow(network_list$flowpaths)

  mapping_table = build_collapse_table(network_list, min_area_sqkm, min_length_km)

  count = 0

  while(nrow(mapping_table) > 0){
    count = count + 1

    if(verbose){
      log_info("Collapsing: {nrow(mapping_table)} features (round {count})")
    }

    fl = filter(network_list$flowpaths, !id %in% mapping_table$id)

    cat = filter(network_list$catchments, id %in% c(mapping_table$id, mapping_table$becomes)) |>
      left_join(mapping_table, by = "id") |>
      mutate(id = ifelse(is.na(becomes), id, becomes)) |>
      geos_union_polygon_hyaggregate("id") |>
      bind_rows(filter(network_list$catchments, !id %in% c(mapping_table$id, mapping_table$becomes)))

    network_list = check_network_validity(flowpaths = fl, cat = cat) |>
      prep_for_ngen()

    mapping_table = build_collapse_table(network_list, min_area_sqkm, min_length_km)
  }

  if(verbose){
    log_success("Collapsed {start - nrow(network_list$flowpaths)} features.")
  }

  network_list
}


#' Prepare Network
#' This function adds a hydrosequence to the flowpath list element and adds area and length calculations.
#' @param network_list a list with flowpath and catchment data
#' @return a list containing flowpath and catchment `sf` objects
#' @export
#' @noRd

prep_for_ngen = function(network_list){
  names(network_list$flowpaths)  = tolower(names(network_list$flowpaths))
  names(network_list$catchments) = tolower(names(network_list$catchments))

  # Add a hydrosequence to the flowpaths
  network_list$flowpaths = add_hydroseq(flowpaths = network_list$flowpaths)

  # Add area and length measures to the network list
  network_list = add_measures(network_list$flowpaths, network_list$catchments)

  network_list
}

#' @title Aggregate Network
#' @param gf A path to a refactored geofarbric file (see US reference fabric) <-
#' @param flowpaths an sf object
#' @param catchments an sf object
#' @param ideal_size_sqkm The ideal size of catchments (default = 10 sqkm)
#' @param min_length_km The minimum allowable length of flowpaths (default = 1 km)
#' @param min_area_sqkm The minimum allowable size of catchments (default = 3 sqkm)
#' @param outfile of not NULL, where to write the output files
#' @param verbose print status updates. Default = TRUE
#' @param overwrite overwrite existing gf file. Default is FALSE
#' @param nexus_topology should a hy-features cat-->nex topology be enforced? default = TRUE
#' @param nexus_locations a data.frame with columns specifiying the ID, and the nexus type.
#' @return if outfile = TRUE, a file path, else a list object
#' @export
#' @importFrom sf st_transform read_sf st_set_crs write_sf st_layers
#' @importFrom dplyr left_join filter
#' @importFrom nhdplusTools get_sorted calculate_total_drainage_area get_streamorder
#' @importFrom logger log_success log_info

aggregate_network_to_distribution = function(gf = NULL,
                                             ideal_size_sqkm = 10,
                                             min_length_km = 1,
                                             min_area_sqkm  = 3,
                                             outfile = NULL,
                                             verbose = TRUE,
                                             overwrite = FALSE,
                                             nexus_topology = TRUE,
                                             routelink_path = NULL,
                                             nexus_locations = NULL) {

  if (!is.null(gf)) {

    if(file.exists(outfile) & overwrite){
      unlink(outfile)
    } else if(file.exists(outfile)){
      warning(outfile, " already exists and overwrite is FALSE", call. = FALSE)
      return(outfile)
    }

    # Make sure needed layers exist
    if (all(!all(c("reconciled", "divides") %in% sf::st_layers(gf)$name),
        !all(c("refactored_flowpaths", "refactored_divides") %in% sf::st_layers(gf)$name))) {
      stop("Make sure you are using a refactored product!")
    }

    # Read Data into R
    flowpaths  <- tryCatch({
      st_transform(read_sf(gf, "refactored_flowpaths"), 5070)
    }, error = function(e){
      st_transform(read_sf(gf, "reconciled"), 5070)
    })

    catchments  <- tryCatch({
      st_transform(read_sf(gf, "refactored_divides"), 5070)
    }, error = function(e){
      st_transform(read_sf(gf, "divides"), 5070)
    })

    if(!is.null(nexus_locations)){
      flowpaths  = left_join(flowpaths,  nexus_locations, by = "ID")
      catchments = left_join(catchments, nexus_locations, by = "ID")
    } else {
        flowpaths$type  = NA
        catchments$type = NA
    }
  }

  # Create a network list and check if DAG and connected
  network_list <- check_network_validity(flowpaths = flowpaths, cat = catchments) %>%
    prep_for_ngen()

  # Perform first aggregation step!
  network_list = aggregate_along_mainstems(network_list,
                                           ideal_size_sqkm,
                                           min_area_sqkm,
                                           min_length_km,
                                           verbose = verbose)

  network_list = collapse_hw(network_list,
                             min_area_sqkm,
                             min_length_km,
                             verbose = verbose)

  if (verbose) {
    log_success("Network is valid.")
  }

  network_list$flowpaths$tot_drainage_areasqkm = calculate_total_drainage_area(st_drop_geometry(
    select(
      network_list$flowpaths,
      ID = id,
      toID = toid,
      area = areasqkm
    )
  ))

  if (verbose) {
    log_success("Total Upstream Drainage Computed and Added.")
  }

  network_list$flowpaths$order = get_streamorder(st_drop_geometry(
    select(
      network_list$flowpaths,
      ID = id,
      toID = toid
    )
  ))

  if (verbose) {
    log_success("Stream Order Computed and Added.")
  }

  if (nexus_topology) {

    if (verbose) {
      log_info("Applying HY_feature topology...")
    }

    network_list$flowpaths =  assign_nex_ids(network_list$flowpaths)

    network_list$catchment_edge_list <- get_catchment_edges_terms(network_list$flowpaths)

    network_list$flowpath_edge_list  <- get_catchment_edges_terms(network_list$flowpaths, catchment_prefix = 'wb-')

    #network_list$waterbody_edge_list <- get_waterbody_edges_terms(network_list$flowpaths)

    network_list$nex =  left_join(
      get_nexus(fp = network_list$flowpaths),
      network_list$catchment_edge_list,
      by = "id"
    )

    network_list$catchments = get_catchment_data(network_list$catchments, network_list$catchment_edge_list)

    network_list$flowpaths  = get_flowpath_data(fline = network_list$flowpaths, catchment_edge_list = network_list$catchment_edge_list)

    if(!is.null(routelink_path)){
      if (verbose) {
        log_info("Creating Routeing Table based on: {routelink_path}")
      }

      network_list$waterbody_attributes =  length_average_routelink(flowpaths = network_list$flowpaths,
                                                                    rl_path = get_routelink_path())

      if (verbose) {
        log_success("Done!")
      }
    }
  }


  if (is.null(outfile)) {
    return(network_list)
  } else {
    if (verbose) {
      log_info("Writing data to: {outfile}")
    }

    write_sf(network_list$flowpaths,  outfile, "aggregate_flowpaths")
    write_sf(network_list$catchments, outfile, "aggregate_divides")

    if (!is.null(network_list$nex)) {
      write_sf(network_list$nex, outfile, "nexus")
    }

    # if (!is.null(network_list$catchment_edge_list)) {
    #   write_sf(network_list$catchment_edge_list, outfile, "catchment_edge_list")
    # }

    if (!is.null(network_list$flowpath_edge_list)) {
      write_sf(network_list$flowpath_edge_list, outfile, "flowpath_edge_list")
    }

    # if (!is.null(network_list$waterbody_edge_list)) {
    #   write_sf(network_list$waterbody_edge_list, outfile, "waterbody_edge_list")
    # }

    if (!is.null(network_list$waterbody_attributes)) {
      write_sf(network_list$waterbody_attributes, outfile, "flowpath_attributes")
    }

    if (verbose) {
      log_success("Done!")
    }

    return(outfile)
  }
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
                                     verbose = TRUE) {

  # Requires having a set, id, toid, levelpathid, hydroseq, member_comid, type, n
  index_table = network_list$flowpaths %>%
    st_drop_geometry() %>%
    group_by(.data$levelpathid) %>%
    arrange(.data$hydroseq) %>%
    mutate(ind = cs_group(.data$areasqkm, .data$lengthkm, .data$type, ideal_size_sqkm, min_area_sqkm, min_length_km)) %>%
    ungroup()   %>%
    group_by(.data$levelpathid, .data$ind) %>%
    mutate(set = cur_group_id(), n = n()) %>%
    ungroup() %>%
    select(set, id, toid, levelpathid, hydroseq, member_comid, type, n)

  v = aggregate_sets(
    flowpaths  = filter(network_list$flowpaths,  id %in% index_table$id),
    cat        = filter(network_list$catchments, id %in% index_table$id),
    index_table
  )

  v = check_network_validity(flowpaths = v$flowpaths, cat = v$catchments) %>%
    prep_for_ngen()

  if (verbose) {
    log_info("Merged to idealized catchment size of {ideal_size_sqkm} sqkm: {nrow(network_list$flowpath) - nrow(v$flowpath)} features removed")
  }

  v

}
