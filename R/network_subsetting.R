#origin = 'wb-16354'
#gpkg = '/Users/mjohnson/github/hyAggregate/data/ngen_01.gpkg'


#' Find ID from location
#' @param gpkg path to a hydrofabric
#' @param pt a spatial point (sf)
#' @return a waterbody ID (character)
#' @export
#' @importFrom sf read_sf st_transform

find_origin = function(gpkg, pt) {
  tmp = read_sf(gpkg,  'aggregate_catchment')[st_transform(pt, 5070), ]
  gsub("cat-", "wb-", tmp$id)
}

#' Subset the upstream protion of a network
#' @param gpkg path to a hydrofabric
#' @param origin the ID to begin navigation
#' @param flowpath_edgelist layer name of flowpath edgelist in gpkg
#' @param flowpath_name layer name of flowpaths in gpkg
#' @param catchment_name layer name of cathcments in gpkg
#' @param  mainstem should only the mainstem flowpath be retruned (default = FALSE)
#' @param export_gpkg a path to write the data to. If NULL a list is returned
#' @export
#' @importFrom nhdplusTools get_sorted
#' @importFrom sf read_sf
#' @importFrom dplyr filter

subset_network = function(gpkg,
                          origin,
                          flowpath_edgelist = 'flowpath_edge_list',
                          flowpath_name     = 'aggregate_flowpaths',
                          catchment_name    = 'aggregate_divides',
                          mainstem = FALSE,
                          export_gpkg = NULL) {
  trace = get_sorted(read_sf(gpkg, flowpath_edgelist),
                     split = TRUE,
                     outlets = origin)

  flowpaths = filter(read_sf(gpkg,  flowpath_name),  id %in% trace$id)
  divides   = filter(read_sf(gpkg,  catchment_name),
                     id %in% flowpaths$realized_catchment)

  if ("nexus" %in% st_layers(gpkg)$name) {
    nexus     = filter(read_sf(gpkg,  "nexus"), id %in% divides$toid)
  }
  nexus     = filter(read_sf(gpkg,  "nexus"), id %in% divides$toid)

  if (mainstem) {
    tmp = filter(flowpaths, id == origin)
    flowpaths = filter(flowpaths, main_id == tmp$main_id)
  }

  ll = list(flowpaths = flowpaths,
            divides = divides,
            nexus = nexus)

  if (!is.null(export_gpkg)) {
    if (length(ll) > 0) {
      names = names(ll)

      for (i in 1:length(ll)) {
        logger::log_info("Writing ", names[i])
        write_sf(ll[[i]], export_gpkg, names[i])
      }
    }

    return(gpkg)
  } else {
    return(ll)
  }
}

#pt = data.frame(x = 2141136, y = 2824888) |>
#   st_as_sf(coords = c("x", "y"), crs = 5070)
#
# set = subset_network(gpkg, find_origin(gpkg, pt))
# mapview::mapview(set)  +pt
