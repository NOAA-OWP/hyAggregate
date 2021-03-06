#origin = 'wb-16354'
#gpkg = '/Users/mjohnson/github/hyAggregate/data/ngen_01.gpkg'

#' Check if geopackage layer exists
#' @param gpkg path to geopackage
#' @param name layer name
#' @return boolean
#' @export

layer_exists = function(gpkg, name){

  if(!file.exists(gpkg)){ return(FALSE) }

  n = sf::st_layers(gpkg)$name

  if(name %in% n){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#' Find ID from location
#' @param gpkg path to a hydrofabric
#' @param pt a spatial point (sf)
#' @return a waterbody ID (character)
#' @export
#' @importFrom sf read_sf st_transform

find_origin = function(gpkg, pt, catchment_name = "aggregate_divides") {
  tmp = read_sf(gpkg,  catchment_name)[st_transform(pt, 5070), ]
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
                          attribute_layers = NULL,
                          export_gpkg = NULL) {

  trace = get_sorted(read_sf(gpkg, flowpath_edgelist),
                     split = TRUE,
                     outlets = origin)

  ll = list()

  ll[['flowpaths']] = filter(read_sf(gpkg,  flowpath_name),  id %in% trace$id)

  ll[['divides']]   = filter(read_sf(gpkg,  catchment_name),
                     id %in% ll$flowpaths$realized_catchment)

  if ("nexus" %in% st_layers(gpkg)$name) {
    ll[['nexus']]     = filter(read_sf(gpkg,  "nexus"), id %in% ll$divides$toid)
  }

  if (mainstem) {
    tmp = filter(ll$flowpaths, id == origin)
    ll[['flowpaths']] = filter(ll[['flowpaths']], main_id == tmp$main_id)
  }

  ll$flowpath_edge_list =  get_catchment_edges_terms(ll$flowpaths, catchment_prefix = 'wb-')

  if(!is.null(attribute_layers)){
    for(i in 1:length(attribute_layers)){
      if(!layer_exists(gpkg, attribute_layers[i])){
        tmp = read_sf(gpkg, attribute_layers[i])
        ll[[attribute_layers[i]]] = filter(tmp, id %in% divides$id )
      }
    }
  }

  if(layer_exists(gpkg, "flowpath_attributes")){
    tmp = read_sf(gpkg, "flowpath_attributes")
    ll[["flowpath_attributes"]] = filter(tmp, id %in% ll$flowpaths$id )
  }

  if (!is.null(export_gpkg)) {
    if (length(ll) > 0) {
      names = names(ll)

      for (i in 1:length(ll)) {
        logger::log_info("Writing ", names[i])
        write_sf(ll[[i]], export_gpkg, names[i])
      }
    }

    return(export_gpkg)
  } else {
    return(ll)
  }
}


