library(hyAggregate)
library(sf)

gpkg = '/Volumes/Transcend/ngen/CONUS-hydrofabric/ngen/ngen_12.gpkg'

## Example 1: Subset based on ID:
set = subset_network(gpkg,
                     origin = 'wb-10977')
sapply(set, nrow)

# Write out subset based on ID to gpkg:

set = subset_network(gpkg,
                     origin = 'wb-10977',
                     export_gpkg = "data/test.gpkg")

st_layers("data/test.gpkg")

#Subset based on ID and find mainstem only:
set = subset_network(gpkg,
                     origin = 'wb-10977',
                     mainstem = TRUE)
sapply(set, nrow)


## Example 2: Subset based on location:
pt = data.frame(x = -275151, y = 845494) |>

set = subset_network(gpkg, find_origin(gpkg, pt))
sapply(set, nrow)

st_layers(gpkg)

### Subset attribtues sets:

set = subset_network(gpkg,
                     find_origin(gpkg, pt),
                     attribute_layers = "cfe_noahowp_attributes")
sapply(set, nrow)







-275151,845494
