sbtools::authenticate_sb("jjohnson@lynker.com", "Mj7-franklin-109034")

ids = nhdplusTools::get_boundaries()$VPUID

#system('git push --mirror https://github.com/NOAA-OWP/hyAggregate.git')



library(zonal)

devtools::load_all()

#get_vaa()

for (i in 12:21) {
  VPU = ids[i]
  logger::log_info("Processing {VPU}")
  base = '/Volumes/Transcend/ngen/CONUS-hydrofabric/'
  g01 = get_reference_fabric(VPU = VPU, dir = glue::glue("{base}refactor"), overwrite = FALSE)
  g02 = glue::glue("{base}ngen/ngen_{VPU}.gpkg")
  logger::log_info(basename(g01), " --> ", basename(g02))

  aggregate_network_to_distribution(
    gf = g01,
    outfile        = g02,
    routelink_path = get_routelink_path(),
    overwrite      = FALSE
  )

  aggregate_cfe_noahowp(gpkg = g02,
                        dir = '/Volumes/Transcend/nwmCONUS-v216/',
                        add_to_gpkg = TRUE)

  if (!dir.exists(glue::glue('{base}/ngen/ngen_{VPU}'))) {
    write_ngen_dir(g02,  export_shapefiles = TRUE)
  }

  logger::log_success("Finished {VPU}")
   gc()
}


vpu = nhdplusTools::get_boundaries()[1:21,]


t= list.files('/Volumes/Transcend/ngen/CONUS-hydrofabric/ngen', pattern = "gpkg")
t=gsub("ngen_", "", gsub(".gpkg", "", t))

ifelse(vpu$VPUID %in% t, "red", "gray")
plot(vpu$geometry, col = ifelse(vpu$VPUID %in% t, "red", "gray"),
     main = glue::glue("Beta Testing Completed\n{Sys.Date()}"))



vpu = nhdplusTools::get_boundaries()[1:21,]


t= list.files('/Volumes/Transcend/ngen/CONUS-hydrofabric/ngen', pattern = "gpkg")
t=gsub("ngen_", "", gsub(".gpkg", "", t))

ifelse(vpu$VPUID %in% t, "red", "gray")
plot(vpu$geometry, col = ifelse(vpu$VPUID %in% t, "red", "gray"),
     main = glue::glue("Currently Completed VPU Regions\n{Sys.Date()}"))



aggregate_cfe_noahowp(gpkg = '/Users/mjohnson/Downloads/ngen_12.gpkg',
                      dir = '/Volumes/Transcend/nwmCONUS-v216/',
                      add_to_gpkg = TRUE)

sf::st_layers('/Users/mjohnson/Downloads/ngen_12.gpkg')

st_layers('/Users/mjohnson/Downloads/ngen_12.gpkg')

mypath = "..."
catchment_name = "subset_divides",
flowpath_name  = "subset_flowpaths"

write_sf(flowpaths, mypath, flowpath_name)
write_sf(catchments, mypath, catchment_name)
write_ngen_dir(g02,
               catchment_name = catchment_name,
               flowpath_name  = flowpath_name)
