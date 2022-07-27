#' Flush existing ID prefixes
#' Given a data object and column, remove a prefix and adjoining "-"
#' @param input input data object
#' @param col column to remove prefix from
#' @return data object with updated column
#' @export
flush_prefix = function(input, col) {
  for (i in col) {
    input[[i]] = as.numeric(gsub(".*-", "", input[[i]]))
  }
  input
}

#' Add hydrosequence
#' @param flowpaths sf object (LINESTRING)
#' @return sf object
#' @export
#' @importFrom  nhdplusTools get_sorted rename_geometry

add_hydroseq = function(flowpaths) {

  flowpaths$terminalID = NULL
  flowpaths$terminalid = NULL
  flowpaths$hydroseq   = NULL
  flowpaths$toid = ifelse(is.na(flowpaths$toid), 0, flowpaths$toid)
  topo = get_sorted(st_drop_geometry(select(flowpaths, id, toid)), split = FALSE)
  gc()
  topo['hydroseq'] = 1:nrow(topo)

  left_join(flowpaths, select(topo, id, hydroseq), by = "id")
}
#' Add Measures to Flowlines and Catchments
#' @param flowpaths sf object (LINESTRING)
#' @param cat sf object (POLYGON)
#' @return list
#' @export
#' @importFrom dplyr select left_join
#' @importFrom sf st_drop_geometry
#' @importFrom nhdplusTools rename_geometry

add_measures = function(flowpaths, cat) {
  flowpaths$lengthkm  = add_lengthkm(flowpaths)
  cat$areasqkm = add_areasqkm(cat)
  flowpaths$areasqkm = NULL
  flowpaths = left_join(flowpaths,
                 select(st_drop_geometry(cat), id, areasqkm),
                 by = "id")
  list(flowpaths  = rename_geometry(flowpaths, "geometry"),
       catchments = rename_geometry(cat, "geometry"))
}


#' @title Flowpaths to linestrings
#' @description Takes an input list of flowpaths object and converts
#' all possible features to LINESTRING
#' @param fl a LINESTRING/MULTILINESTRING `sf` flowlines object
#' @return a LINESTRING `sf` flowlines object
#' @export
#' @importFrom sf st_geometry_type st_geometry st_line_merge
#' @importFrom dplyr bind_rows

flowpaths_to_linestrings = function(fl){
  bool = (st_geometry_type(sf::st_geometry(fl)) == "MULTILINESTRING")
  multis = fl[bool, ]
  if(nrow(multis) > 0){
    sf::st_geometry(multis) = st_line_merge(sf::st_geometry(multis))
  }
  singles = fl[!bool, ]

  bind_rows(multis, singles)
}

#' Compute length in kilometers
#' @param x LINESTRING sf object
#' @return numeric vector
#' @export
#' @importFrom units set_units drop_units
#' @importFrom sf st_length

add_lengthkm = function (x) { drop_units(units::set_units(st_length(x), "km")) }

#' Compute area in square kilometers
#' @param x POLYGON sf object
#' @return numeric vector
#' @export
#' @importFrom units set_units drop_units
#' @importFrom sf st_area

add_areasqkm = function (x) { drop_units(set_units(st_area(x), "km2")) }

#' Aggregate Sets by Index Table
#' @param flowpaths LINESTRING flowpaths
#' @param cat POLYGON catchments
#' @param index_table index table to aggregate with
#' @return a list of catchments and flowpaths that have been validated
#' @export
#' @importFrom dplyr group_by mutate slice_max ungroup select left_join everything filter bind_rows rename `%>%`
#' @importFrom sf st_as_sf
#' @importFrom nhdplusTools rename_geometry get_sorted
#'
aggregate_sets = function(flowpaths, cat, index_table) {

  set_topo = index_table %>%
    group_by(set) %>%
    mutate(member_comid  = paste(.data$member_comid, collapse = ","),
           type  = paste(.data$type[!is.na(.data$type)], collapse = ","),
           type = ifelse(.data$type == "", NA, .data$type)) %>%
    slice_max(hydroseq) %>%
    ungroup() %>%
    select(set, id = toid, levelpathid, member_comid, type) %>%
    left_join(select(index_table, toset = set, id), by = "id") %>%
    select(-id)

  ####

  single_flowpaths = filter(index_table, n == 1) %>%
    left_join(flowpaths, by = "id") %>%
    st_as_sf() %>%
    select(set) %>%
    rename_geometry("geometry")

  flowpaths_out  = filter(index_table, n > 1) %>%
    left_join(flowpaths, by = "id") %>%
    st_as_sf() %>%
    select(set) %>%
    filter(!sf::st_is_empty(.)) %>%
    geos_union_linestring_hyaggregate('set') %>%
    rename_geometry("geometry") %>%
    bind_rows(single_flowpaths) %>%
    select(set) %>%
    left_join(set_topo, by = "set") %>%
    rename(id = set, toid = toset)

  ####

  single_catchments = filter(index_table, n == 1) %>%
    left_join(cat, by = "id") %>%
    st_as_sf() %>%
    select(set) %>%
    rename_geometry("geometry")

  catchments_out  = filter(index_table, n != 1) %>%
    left_join(cat, by = "id") %>%
    st_as_sf() %>%
    select(set) %>%
    filter(!sf::st_is_empty(.)) %>%
    geos_union_polygon_hyaggregate('set') %>%
    rename_geometry("geometry") %>%
    bind_rows(single_catchments) %>%
    select(set) %>%
    left_join(set_topo, by = "set") %>%
    rename(id = set, toid = toset)

  catchments_out$toid = ifelse(is.na(catchments_out$toid), 0, catchments_out$toid)

  check_network_validity(flowpaths = flowpaths_out, cat = catchments_out) %>%
    prep_for_ngen()
}

##' @title Fast LINESTRING union
#' @description Wayyyy faster then either data.table, or sf based line merging
#' @param lines lines to merge
#' @param ID ID to merge over
#' @return LINESTRING sf object
#' @export
#' @importFrom terra aggregate vect
#' @importFrom dplyr select
#' @importFrom sf st_as_sf

geos_union_linestring_hyaggregate = function (lines, ID)  {
  aggregate(vect(lines), by = eval(ID)) %>%
    st_as_sf() %>%
    select(!!ID) %>%
    flowpaths_to_linestrings()
}

#' @title Fast POLYGON Union
#' @description This is significantly faster then sf::st_union or summarize
#' @param poly sf POLYGON object
#' @param ID the column name over which to union geometries
#' @return sf object
#' @export
#' @importFrom terra aggregate vect makeValid
#' @importFrom dplyr select
#' @importFrom sf st_as_sf st_collection_extract st_geometry_type st_make_valid

geos_union_polygon_hyaggregate = function(poly, ID) {

  poly = makeValid(vect(poly)) %>%
    aggregate(by = eval(ID)) %>%
    st_as_sf() %>%
    select(!!ID)

  if (any(grepl("COLLECTION",  st_geometry_type(poly)))) {
    poly = st_collection_extract(poly, "POLYGON")
  }
  return(poly)
}

#' Enforces area and length grouping
#' @description This function takes a vector of area's and length's and returns a
#' grouping vector that enforces the grouping of lengths and areas less then defined thresholds
#' @param l a vector of lengths
#' @param a a vector of areas
#' @param lthres a minimum length that must be achieved
#' @param athres a minimum length that must be achieved
#' @return a vector of length(a) containing grouping indexes
#' @export

agg_length_area   <- function(l, a, lthres, athres) {

  ids = 1:length(l)

  if(length(ids) != 1){

    if(!is.null(lthres)){
      for (i in 1:(length(l)-1)) {
        if (l[i] < lthres) {
          ids[(i+1):length(l)] = ids[(i+1):length(l)] - 1
          l[i+1] = l[i] + l[i+1]
          l[i]   = l[i+1]
          a[i+1] = a[i] + a[i+1]
          a[i] =   a[i+1]
        }
      }
    }

    if(!is.null(athres)){
      for (i in 1:(length(a)-1)) {
        if (a[i] < athres) {
          ids[(i+1):length(a)] = ids[(i+1):length(a)] - 1
          a[i+1] = a[i] + a[i+1]
          a[i] =   a[i+1]
        }
      }
    }

    if(is.null(athres)){ athres = 0 }
    if(is.null(lthres)){ lthres = 0 }

    if(a[length(a)] < athres | l[length(l)] < lthres){
      ids[length(ids)] = pmax(1, ids[length(ids)] - 1)
    }
  }

  return (ids)
}

splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% (pos +1)) ))

#' Index a Vector by cumulative Sum
#'
#' @param a a vector of values
#' @param athres the ideal size of each index. Cummulative sums will get as close to this value without exceeding it
#' @return a vector of length(a)
#' @export

assign_id = function(a, athres
                     #, amin
                     ){

  cumsum <- 0
  group  <- 1
  result <- numeric()

  for (i in 1:length(a)) {
    cumsum <- cumsum + a[i]
    if (cumsum > athres) {
      group <- group + 1
      cumsum <- a[i]
    }
    result = c(result, group)
  }

  # if(a[1] < amin & length(a) > 1){
  #   result[1] = result[2]
  # }
  #
  # if(a[length(a)] < amin & length(a) > 1){
  #   result[length(result)] = result[length(result) - 1]
  # }

  return (result)
}


#' Cumulative sum area grouping
#' @description This function takes a vector of areas and lengths and returns a
#' index vector that combines them towards an ideal aggregate area (athres). While enforcing a minimum area (amin) and length (lmin).
#' Additionally, this function can take a set of indexes to exclude over which the network cannot be aggregated.
#' @param areas a vector of areas
#' @param lengths a vector of lengths
#' @param exclude a vector of equal length to areas and lengths. Any non NA value will be used to enforce an aggregation break
#' @param ideal_size_sqkm a vector of areas
#' @param amin a threshold, or target, cumulative size
#' @param lmin a threshold, or target, cumulative size
#' @return a vector of length(areas) containing grouping indexes
#' @export

cs_group <- function(areas, lengths, exclude, ideal_size_sqkm, amin, lmin) {

  areas[is.na(areas)] = 0
  lengths[is.na(lengths)] = 0

  if(length(areas) == 1){ return(1) }

  break_index = which(!is.na(exclude))

  if(length(break_index) != 0){
    sub_areas = splitAt(areas, break_index)
    sub_lengths = splitAt(lengths, break_index)
    #splitAt(index_table$id, break_index)
  } else {
    sub_areas = list(areas)
    sub_lengths = list(lengths)
  }

  if(all(lengths(sub_areas) != lengths(sub_areas))){
    stop("Yuck~")
  }

   o1 = lapply(sub_areas, assign_id, athres = ideal_size_sqkm)
   #lengths(o1) == lengths(sub_areas)
   o2 = lapply(1:length(sub_areas),   function(i) { pinch_sides(   x = sub_areas[[i]],   ind = o1[[i]], thres = amin) })
   #lengths(o2) == lengths(sub_areas)
   o3 = lapply(1:length(sub_lengths), function(i) { pinch_sides(   x = sub_lengths[[i]], ind = o2[[i]], thres = lmin) })
   #lengths(o3) == lengths(sub_areas)
   o4 = lapply(1:length(sub_areas),   function(i) { middle_massage(x = sub_areas[[i]],   ind = o3[[i]], thres = amin) })
   #lengths(o4) == lengths(sub_areas)
   o5 = lapply(1:length(sub_lengths), function(i) { middle_massage(x = sub_lengths[[i]], ind = o4[[i]], thres = lmin) })
   #lengths(o5) == lengths(sub_areas)

    for(i in 1:length(o5)){ o5[[i]] = o5[[i]] + 1e9*i }

    unlist(o5)

}

#' Re-index the edges of vector by threshold
#' Merge the outside edges of a vector if they are less then the provides threshold.
#' @param x vector of values
#' @param ind current index values
#' @param thres threshold to evaluate x
#' @return a vector of length(x) containing grouping indexes
#' @export

pinch_sides = function(x, ind, thres){
  # i = 2
  # x = sub_areas[[i]]; ind = o[[i]]
  tmp_areas = unlist(lapply(split(x, ind), sum))

  if(length(tmp_areas) == 1){ return(ind) }
  #
  n = as.numeric(names(tmp_areas))

  if(tmp_areas[1] < thres){
    names(tmp_areas)[1] = names(tmp_areas[2])
  }

  if(tmp_areas[length(tmp_areas)] < thres){
    names(tmp_areas)[length(tmp_areas)] = names(tmp_areas[length(tmp_areas) - 1])
  }

  n2 = as.numeric(names(tmp_areas))

  n2[match(ind, n)]
}


#' Re-index the interior of vector by threshold
#' Merge the interior values of a vector if they are less then the provided threshold.
#' Merging will look "up" and "down" the vector and merge into the smaller of the two.
#' @param x vector of values
#' @param ind current index values
#' @param thres threshold to evaluate x
#' @return a vector of length(x) containing grouping indexes
#' @export

middle_massage = function(x, ind, thres){

  # i = 5
  # x = sub_areas[[i]]; ind = o2[[i]]

  tmp_areas = unlist(lapply(split(x, ind), sum))

  if(length(tmp_areas) == 1){ return(ind) }

  n = as.numeric(names(tmp_areas))

  if(any(tmp_areas < thres)){
    tmp = which(tmp_areas < thres)

    for(j in 1:length(tmp)){
      base = as.numeric(tmp[j])
      edges = c(base - 1, base + 1)
      becomes = names(which.min(tmp_areas[edges]))
      names(tmp_areas)[base] = becomes
    }
  }

  n2 = as.numeric(names(tmp_areas))

  n2[match(ind, n)]
}


#' Check Network Validity
#' **INTERNAL** function that validates a flowpath and catchment network
#' @param flowpaths a LINESTRING `sf` flowpaths object
#' @param cat a POLYGON `sf` catchments object
#' @param term_cut cutoff integer to define terminal IDs
#' @return a list containing flowline and catchment `sf` objects
#' @export
#' @importFrom dplyr mutate select left_join
#' @importFrom sf st_drop_geometry

check_network_validity     <- function(flowpaths, cat, term_cut = 1e9, check = TRUE){

  names(flowpaths) = tolower(names(flowpaths))
  names(cat) = tolower(names(cat))

  flowpaths$toid    = ifelse(is.na(flowpaths$toid), 0, flowpaths$toid)
  DAG        = network_is_dag(flowpaths)
  CONNECTION = sum(!(flowpaths$toid %in% flowpaths$id | flowpaths$toid > term_cut | flowpaths$toid == 0)) == 0

  if(!check){ return(list(flowpaths = fl, catchments = cat))}

  if(all(DAG,  CONNECTION)){

    flowpaths  = mutate(flowpaths, lengthkm = add_lengthkm(flowpaths))

    if(!is.null(cat)){
      cat =  mutate(cat, areasqkm = add_areasqkm(cat))  %>%
        select(.data$id, .data$areasqkm)

      if('areasqkm' %in% names(flowpaths)){
        flowpaths = flowpaths %>%
          select(-.data$areasqkm) %>%
          left_join(st_drop_geometry(cat), by = "id")
      } else {
        flowpaths =  left_join(flowpaths, st_drop_geometry(cat), by = "id")
      }
    }

    return(list(flowpaths = flowpaths, catchments = cat))

  } else {
    if(!DAG){ stop("Network is not a graph.")}
    if(!CONNECTION){stop("All toIDs are not present in network")}

    return(NULL)
  }
}

#' Check if network is DAG
#' Checks in a `sf` flowline network is a DAG (Directed acyclic graph).
#' @param fl a LINESTRING `sf` flowlines object
#' @param ID the name of the ID column in `fl`
#' @param toID the name of the toID column in `fl`
#' @return boolean
#' @export
#' @importFrom igraph graph_from_data_frame is.dag
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select

network_is_dag = function(fl, ID = "id", toID = "toid"){
  st_drop_geometry(select(fl, !!ID, !!toID)) %>%
    graph_from_data_frame(directed = TRUE) %>%
    is.dag()
}

#' Get Nexus Locations
#'
#' @param fp a sf flowpath object
#' @param term_cut cutoff integer to define terminal IDs
#' @return POINT sf object
#' @export
#' @importFrom dplyr left_join select mutate filter select bind_rows group_by ungroup arrange rename slice_head slice_min
#' @importFrom nhdplusTools get_node rename_geometry
#' @importFrom sf st_intersects st_drop_geometry st_as_sf

get_nexus_locations = function(fp, term_cut =  1e9){

  fp$toid = ifelse(is.na(fp$toid), 0, fp$toid)

  term_node = filter(fp, .data$toid > term_cut | .data$toid == 0) %>%
    rename_geometry("geometry") %>%
    mutate(geometry = get_node(., "end")$geometry) %>%
    dplyr::slice_min(.data$hydroseq) %>%
    select(id = .data$toid)


  if(nrow(fp) <= 2){
    nex = term_node
  } else {
    nex = fp %>%
      left_join(st_drop_geometry(select(., toid = .data$id, ds_toID = .data$toid)), by = c("toid")) %>%
      filter(.data$id %in% unique(.data$toid)) %>%
      rename_geometry("geometry") %>%
      rmapshaper::ms_explode() %>%
      mutate(geometry = get_node(., "start")$geometry) %>%
      select(.data$id, .data$toid)
  }

  imap = st_intersects(nex, fp)

  df = data.frame(
    id       = rep(nex$id, times = lengths(imap)),
    touches  = fp$id[unlist(imap)]) %>%
    mutate(cond = ifelse(.data$id == .data$touches, "move","aaa")) %>%
    group_by(.data$id) %>%
    arrange(.data$cond) %>%
    slice_head(n = 1) %>%
    ungroup()

  to_move = filter(fp, .data$id %in% filter(df, .data$cond == "move")$id) %>%
    select(.data$id) %>%
    rename_geometry("geometry")

  to_keep = filter(nex, .data$id %in% filter(df, .data$cond != "move")$id) %>%
    select(.data$id) %>%
    rename_geometry("geometry")

  fp_ends = bind_rows(
    select(fp, .data$id, .data$hydroseq, .data$toid) %>%
      rename_geometry("geometry") %>%
      rmapshaper::ms_explode() %>%
      mutate(geometry = get_node(., "end")$geometry, pos = "end"),
    select(fp, .data$id, .data$hydroseq, .data$toid) %>%
      rename_geometry("geometry") %>%
      mutate(geometry = get_node(., "start")$geometry, pos = "start")
  )

  imap = st_intersects(to_move, fp_ends)

  df = data.frame(
    id       = rep(to_move$id, times = lengths(imap)),
    touches  = fp_ends$id[unlist(imap)],
    pos      = fp_ends$pos[unlist(imap)],
    hs       = fp_ends$hydroseq[unlist(imap)]) %>%
    left_join(rename(fp_ends, touches = .data$id), by = c('touches', 'pos')) %>%
    group_by(.data$id) %>%
    filter(.data$id == .data$toid) %>%
    dplyr::slice_max(.data$hydroseq) %>%
    ungroup() %>%
    st_as_sf() %>%
    select(.data$id) %>%
    bind_rows(to_keep) %>%
    bind_rows(term_node) %>%
    filter(!duplicated(.))

  df
}


add_grid_mapping = function(gpkg = NULL,
                                 catchment_name = "aggregate_divides",
                                 template = '/Users/mjohnson/Downloads/AORC-OWP_2012063021z.nc4',
                                 grid_name = NULL,
                                 add_to_gpkg = TRUE){

  out = weight_grid(rast(template),
                    geom = sf::read_sf(gpkg, catchment_name),
                    ID = "id")

  if(add_to_gpkg){
    if(is.null(grid_name)){ stop("To write this file to a gpkg, a `grid_name` must be provided ...") }
    write_sf(out, gpkg, grid_name)
  } else{
    return(out)
  }


}

