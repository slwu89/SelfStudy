

.onLoad <- function(libname, pkgname) {

  XRJulia::RJulia()
  XRJulia::juliaAddToPath("./inst/julia/")
  XRJulia::juliaUsing("MutationDriftSelection")

  assign(x = "TIMESTAMP",  as.numeric(Sys.time()),asNamespace('MutationDriftSelection'))


}
