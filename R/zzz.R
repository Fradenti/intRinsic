# Clean after unload
.onUnload <- function(libpath) {
  library.dynam.unload("intRinsic", libpath)
}
