library(terra)
library(data.table)

setwd(this.path::here())

# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

replace_nodata <- function(r, nodata = -9999) {
  for (x in nodata) {
    r <- terra::ifel(r == x, NA, r)
  }
  r
}

check_same_geometry <- function(...) {
  rs <- list(...)
  ref <- rs[[1]]
  
  for (i in seq_along(rs)[-1]) {
    if (!terra::compareGeom(ref, rs[[i]],
                            stopOnError = FALSE,
                            crs = TRUE,
                            ext = TRUE,
                            rowcol = TRUE,
                            res = TRUE,
                            orig = TRUE)) {
      stop("Input rasters do not have matching geometry: raster ", i)
    }
  }
  
  invisible(TRUE)
}

# ---------------------------------------------------------------------
# Read the zone rasters
# ---------------------------------------------------------------------

input_dir <- "./inputs/Chloe/OkaWen_LFTFC/LFTFC_MaxEffect/Input/Zone_01/FOA"

fvt  <- rast(file.path(input_dir,
                       "landfire2020_FOA_Z01_FVT_120m.tif"))

fvc  <- rast(file.path(input_dir,
                       "landfire2020_FOA_Z01_FVC_120m.tif"))

fvh  <- rast(file.path(input_dir,
                       "landfire2020_FOA_Z01_FVH_120m.tif"))

fdist <- rast(file.path(input_dir,
                        "landfire2020_FOA_Z01_FDist121_120m.tif"))

bps <- rast(file.path(input_dir,
                      "landfire2020_FOA_Z01_BPS_120m.tif"))

check_same_geometry(fvt, fvc, fvh, fdist, bps)

# ---------------------------------------------------------------------
# Minimal preprocessing
# ---------------------------------------------------------------------

# Convert standard nodata values to NA
fvt   <- replace_nodata(fvt,   -9999)
fvc   <- replace_nodata(fvc,   -9999)
fvh   <- replace_nodata(fvh,   -9999)
fdist <- replace_nodata(fdist, -9999)

# BPS convention used by the collaborator:
#   -9999 = missing
#   -1111 = valid special value recoded to 8888
bps <- replace_nodata(bps, -9999)
bps <- terra::ifel(bps == -1111, 8888, bps)

# Combine the rasters
stack <- c(fvt, fvc, fvh, fdist, bps)
names(stack) <- c("fvt", "fvc", "fvh", "fdist", "bps")

# Remove categorical raster levels without changing cell values
for (i in seq_len(nlyr(stack))) {
  levels(stack[[i]]) <- NULL
}