#' Apply Nexus Topology
#' This function enforces the nexus-->flowpath topology and adds a catchment and flowpath
#' edgelist to the network_list object. Additonally, nexus locations are identified and
#' added as well.
#' @param network_list  a list containing flowpath and catchment `sf` objects
#' @param nexus_prefix  character prefix for nexus IDs
#' @param terminal_nexus_prefix character prefix for terminal nexus IDs
#' @param catchment_prefix character prefix for catchment IDs
#' @param waterbody_prefix character prefix for catchment IDs
#' @param term_cut terminal IDs begin above a defined threshold
#' @return list
#' @export


apply_nexus_topology = function(network_list,
                                nexus_prefix = "nex-",
                                terminal_nexus_prefix = "tnex-",
                                catchment_prefix = "cat-",
                                waterbody_prefix = "wb-",
                                term_cut = 1e9,
                                verbose = TRUE){

  hyaggregate_log("INFO", "\n--- Applying HY_feature topology ---\n", verbose)

  nl = list()

  nl$flowpaths     <- assign_nex_ids(fline = network_list$flowpaths, term_cut = term_cut)

  nl$catchment_edge_list <- get_catchment_edges_terms(flowpaths = nl$flowpaths,
                                                      nexus_prefix = nexus_prefix,
                                                      terminal_nexus_prefix = terminal_nexus_prefix,
                                                      catchment_prefix = catchment_prefix,
                                                      term_cut = term_cut )

  nl$flowpath_edge_list  <- get_catchment_edges_terms(nl$flowpaths,
                                                      nexus_prefix = nexus_prefix,
                                                      terminal_nexus_prefix = terminal_nexus_prefix,
                                                      catchment_prefix = waterbody_prefix,
                                                      term_cut = term_cut)

  nex =  get_nexus(fline = nl$flowpaths,
                   term_cut = term_cut,
                   nexus_prefix = nexus_prefix,
                   terminal_nexus_prefix = terminal_nexus_prefix)

  nl$nex =  left_join(nex,  nl$catchment_edge_list, by = "id")

  hyaggregate_log("INFO", glue("Created {nrow(network_list$nex)} nexus locations"), verbose)

  nl$catchments <- get_catchment_data(network_list$catchments,
                                      nl$catchment_edge_list,
                                      catchment_prefix = catchment_prefix)

  nl$flowpaths  <- get_flowpath_data( fline = nl$flowpaths,
                                      catchment_edge_list = nl$catchment_edge_list,
                                      waterbody_prefix = waterbody_prefix,
                                      catchment_prefix = catchment_prefix)

  nl

}



#' Get catchment edge list
#' @description get a edge list for catchments
#' @param flowpaths  sf data.frame containing hyRefactor or hyAggregate output.
#' @param nexus_prefix  character prefix for nexus IDs
#' @param terminal_nexus_prefix character prefix for terminal nexus IDs
#' @param catchment_prefix character prefix for catchment IDs
#' @param term_cut terminal IDs begin above a defined threshold
#' @return data.frame
#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select mutate left_join bind_rows `%>%`

get_catchment_edges_terms = function(flowpaths,
                                     nexus_prefix = "nex-",
                                     terminal_nexus_prefix = "tnex-",
                                     catchment_prefix = "cat-",
                                     term_cut = 1e9) {

  fline = select(st_drop_geometry(flowpaths), id, toid)

  fline = flush_prefix(fline, c("id", "toid"))

  obj1 = fline %>%
    mutate(id = paste0(catchment_prefix, .data$id),
           toid = paste0(
             ifelse(.data$toid > term_cut, terminal_nexus_prefix, nexus_prefix), .data$toid
           ))

  obj2 =  data.frame(id = unique(fline$toid)) %>%
    left_join(mutate(select(fline, id), toid = id), by = "id") %>%
    mutate(toid = ifelse(is.na(.data$toid), 0, .data$toid)) %>%
    mutate(id =  paste0(
      ifelse(.data$id > term_cut, terminal_nexus_prefix, nexus_prefix),
      .data$id
    ),
    toid = paste0(catchment_prefix, .data$toid))

  bind_rows(obj1, obj2)
}


assign_nex_ids = function(fline, term_cut = 1e9) {

  term_node = filter(fline, toid == 0 | is.na(toid)) %>%
    mutate(toid = term_cut + 1:n())


  no_term = filter(fline, !id %in% term_node$id)

  bind_rows(term_node, no_term) %>%
    rename_geometry("geometry")

}


#' @title get nexuses
#' @title get nexuses
#' @param fline sf data.frame NHDPlus Flowlines or hyRefactor output.
#' @param nexus_prefix character prefix for nexus IDs
#' @importFrom sf st_coordinates st_as_sf st_crs
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by filter ungroup select n row_number rename
#' @export
#'

get_nexus <- function(fline, term_cut = 1e9,
                      nexus_prefix = "nex-",
                      terminal_nexus_prefix = "tnx-") {

  nexus <- fline %>%
    st_cast("MULTILINESTRING") %>%
    st_coordinates() %>%
    as.data.frame()

  if("L2" %in% names(nexus)) {
    nexus <- rename(nexus, GG = .data$L2)
  } else {
    nexus <- rename(nexus, GG = .data$L1)
  }

  fline <- check_nexus(fline)

  nexus <- nexus %>%
    group_by(.data$GG) %>%
    filter(row_number() == n()) %>%
    ungroup() %>%
    select(.data$X, .data$Y) %>%
    st_as_sf(coords = c("X", "Y"), crs = st_crs(fline))

  nexus$id <- fline$to_nID
  nexus$type <- ifelse(is.na(fline$type), "infered", fline$type)
  nexus$type_id <- fline$value

  nexus$id = ifelse(nexus$id >= term_cut, paste0(terminal_nexus_prefix, nexus$id), paste0(nexus_prefix, nexus$id))

  if(length(unique(nexus$id)) < nrow(nexus)) {
    nexus <- group_by(nexus, .data$id) %>%
      filter(row_number() == 1) %>%
      ungroup()
  }

  return(nexus)
}

check_nexus <- function(fline) {

  fline$from_nID <- fline$id

  fline <- left_join(fline,
                     select(st_drop_geometry(fline), .data$id, to_nID = .data$from_nID),
                     by = c("toid" = "id"))

  fline$to_nID[is.na(fline$to_nID)] <- fline$toid[is.na(fline$to_nID)]

  fline

}


#' Get waterbody edge list
#' @description get a edge list for waterbodies
#' @param flowpaths  sf data.frame containing hyRefactor or hyAggregate output.
#' @param wb_prefix  character prefix for waterbody IDs
#' @param terminal_wb_prefix character prefix for terminal waterbody IDs
#' @param term_cut terminal IDs begin above a defined threshold
#' @return data.frame
#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select mutate

get_waterbody_edges_terms = function(flowpaths,
                                     wb_prefix = "wb-",
                                     terminal_wb_prefix = "twb-",
                                     term_cut = 1e9) {

  fline = select(st_drop_geometry(flowpaths), id, toid)

  fline = flush_prefix(fline, "id")
  fline = flush_prefix(fline, "toid")

  fline %>% select(.data$id, .data$toid) %>%
    mutate(id = paste0(
      ifelse(.data$id > term_cut, terminal_wb_prefix, wb_prefix),
      .data$id
    ),
    toid = paste0(
      ifelse(.data$toid > term_cut, terminal_wb_prefix, wb_prefix),
      .data$toid
    ))
}


#' Get Catchment Data
#' @description get a edge list for waterbodies
#' @param catchment  sf data.frame containing hyRefactor or hyAggregate output.
#' @param wb_prefix  character prefix for waterbody IDs
#' @param terminal_wb_prefix character prefix for terminal waterbody IDs
#' @param cutoff terminal IDs begin above a defined threshold
#' @return data.frame
#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select mutate

get_catchment_data = function(catchment,
                              catchment_edge_list,
                              catchment_prefix = "cat-") {
  catchment %>%
    mutate(id = paste0(catchment_prefix, .data$id)) %>%
    left_join(catchment_edge_list,  by = "id")
}


get_flowpath_data = function(fline,
                             catchment_edge_list,
                             waterbody_prefix = "wb-",
                             catchment_prefix = "cat-") {

  if ("main_id" %in% names(fline)) {
    fline = rename(fline, levelpathid = main_id)
  }

  if (!"slope" %in% names(fline)) {
    fline = add_slope(flowpaths = fline)
  }


  select(
      fline,
      id = .data$id,
      lengthkm = .data$lengthkm,
      slope_percent = .data$slope,
      main_id = .data$levelpathid,
      member_comid = .data$member_comid,
      tot_drainage_areasqkm = tot_drainage_areasqkm,
      order = order
    ) %>%
    mutate(id = paste0(waterbody_prefix, .data$id)) %>%
    mutate(realized_catchment = gsub(waterbody_prefix, catchment_prefix, id)) %>%
    left_join(catchment_edge_list, by = c("realized_catchment" = "id"))

}
