#' Authenticated Upload to s3
#' You must authenticate s3 for this working session to use this function. Also,
#' the aws.s3 package is required on your machine and is NOT a hyAggregate dependency.
#' @param path a path to a file or directory
#' @param bucket the s3 bucket to upload to
#' @param prefix an optional file object prefix (think of sub directory)
#' @param verbose should messages be emitted?

upload_to_aws = function(path,
                         bucket = "formulations-dev",
                         prefix = "hf_1.0",
                         verbose = TRUE){

  check_pkg("aws.s3")

  prefix = gsub("/$", "", prefix)

  if(!file_test("-f", path)){
    f = list.files(path,
                   full.names = TRUE,
                   recursive = TRUE)

    o = file.path(prefix, basename(path), gsub("/", "", gsub(path, "", f)))

    fin = lapply(1:length(f), function(i) {
      aws.s3::put_object(file = f[i],
                         object = o[i],
                         bucket = bucket,
                         multipart = TRUE,
                         show_progress = verbose)
    })

  } else {

    o = file.path(prefix, basename(path))

    aws.s3::put_object(file = path,
                       object = o,
                       bucket = bucket,
                       multipart = TRUE,
                       show_progress = verbose)
  }

  hyaggregate_log("SUCCESS", glue("{length(o)} file(s) uplaoded to s3!"), verbose)

}


#' Find Processing Unit (Vector or Raster)
#' @param location sf object
#' @param pu either "vpu" or "rpu" (default is "vpu")
#' @return intersection processing units (sf object)
#' @export
#' @importFrom nhdplusTools get_boundaries
#' @importFrom sf st_transform
#' @importFrom dplyr filter slice_min
#' @importFrom sbtools item_file_download

find_pu = function(location, pu = "vpu"){ get_boundaries(pu)[st_transform(location, 4269),] }

#' Extract File Extension
#' @details returns file extension
#' @param x a file path
#' @param prefix character string to precede extracted extension.  default = "". (Usfull if you want to keep the ".")
#' @return character string
#' @export

.getExtension = function (x, prefix = "") {
  ext <- strsplit(basename(x), split = "\\.")[[1]]
  return(paste0(prefix, ext[length(ext)]))
}

#' Download Reference Fabric Data by VPU ID
#'
#' @param VPU a VPU ID
#' @param type either 'refactored' (default) or 'reference'
#' @param dir directory path to save data to
#' @param overwrite should existing files be overwritten? (default = FALSE)
#' @return file path
#' @export
#' @importFrom nhdplusTools get_boundaries
#' @importFrom dplyr filter slice_min
#' @importFrom jsonlite fromJSON
#' @importFrom sbtools item_file_download
#' @importFrom httr GET write_disk

get_reference_fabric = function(VPU = "01",
                                type = "refactored",
                                dir  = NULL,
                                overwrite = FALSE) {
  if (is.null(dir)) {
    stop("`dir` cannot be NULL", call. = FALSE)
  }

  if(!VPU %in% get_boundaries()$VPUID){
    stop(VPU, " is not a valid VPU ID", call. = FALSE)
  }

  path = ifelse(
    type == "refactored",
    'https://www.sciencebase.gov/catalog/item/61fbfdced34e622189cb1b0a',
    #reference
    'https://www.sciencebase.gov/catalog/item/61295190d34e40dd9c06bcd7'
  )

  xx = fromJSON(paste0(path, '?format=json'), simplifyDataFrame = TRUE)

  find = slice_max(filter(xx$files, grepl(VPU, xx$files$name)), dateUploaded)

  out = file.path(dir, find$name)

  out2 = gsub("zip", "gpkg", out)

  if (file.exists(out2) & !overwrite) {
    return(out2)
  } else {
    tryCatch({
      sbtools::item_file_download(
        sb_id = basename(path),
        names = find$name,
        destinations = out,
        overwrite_file = TRUE
      )}, error = function(e){
        httr::GET(find$url, httr::write_disk(out, overwrite = TRUE))
      })
  }


  if (.getExtension(out) == "zip") {
    unzip(out, exdir = dir)
    unlink(out)
    return(out2)
  } else {
    return(out)
  }
}

#' Authenticated Upload to ScienceBase
#' @param x A string vector of paths to files to be uploaded
#' @param item a character ScienceBase ID corresponding to the item (default = '629a4246d34ec53d276f446d')
#' @param username Sciencebase username, prompts user if not supplied
#' @param password Sciencebase password, prompts user if not supplied
#' @importFrom sbtools authenticate_sb item_append_files

upload_aggregate_gpkg = function(x, item = "629a4246d34ec53d276f446d",
                                 username, password) {

  authenticate_sb(username, password)
  item_append_files(item, x)
}
