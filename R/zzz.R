# Code copied from: https://github.com/luukvdmeer/sfnetworks/blob/master/R/zzz.R
.onLoad <- function(...) {
  has_ggplot2_3.0 <-
    requireNamespace("ggplot2", quietly = TRUE) &&
    utils::packageVersion("ggplot2") >= "3.0.0"

  if (has_ggplot2_3.0) {
    s3_register("ggplot2::autoplot", "twonn_linfit")
    s3_register("ggplot2::autoplot", "twonn_bayes")
    s3_register("ggplot2::autoplot", "gride_bayes")
    s3_register("ggplot2::autoplot", "gride_mle")
    s3_register("ggplot2::autoplot", "gride_evolution")
    s3_register("ggplot2::autoplot", "Hidalgo")
    s3_register("ggplot2::autoplot", "twonn_mle")
  }
  invisible()
}

# Register S3 helper
# As suggested here: https://vctrs.r-lib.org/reference/s3_register.html
# (To avoid taking a dependency on vctrs for this one function,
# please feel free to copy and paste the function source into your own
# package.)
s3_register <- function(generic, class, method = NULL) {
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)

  pieces <- strsplit(generic, "::")[[1]]
  stopifnot(length(pieces) == 2)
  package <- pieces[[1]]
  generic <- pieces[[2]]

  caller <- parent.frame()

  get_method_env <- function() {
    top <- topenv(caller)
    if (isNamespace(top)) {
      asNamespace(environmentName(top))
    } else {
      caller
    }
  }
  get_method <- function(method, env) {
    if (is.null(method)) {
      get(paste0(generic, ".", class), envir = get_method_env())
    } else {
      method
    }
  }

  method_fn <- get_method(method)
  stopifnot(is.function(method_fn))

  # Always register hook in case package is later unloaded & reloaded
  setHook(packageEvent(package, "onLoad"),
          function(...) {
            ns <- asNamespace(package)

    # Refresh the method, it might have been updated by `devtools::load_all()`
            method_fn <- get_method(method)

            registerS3method(generic, class, method_fn, envir = ns)
          })

  # Avoid registration failures during loading (pkgload or regular)
  if (!isNamespaceLoaded(package)) {
    return(invisible())
  }

  envir <- asNamespace(package)

  # Only register if generic can be accessed
  if (exists(generic, envir)) {
    registerS3method(generic, class, method_fn, envir = envir)
  }

  invisible()
}
# nocov end


# Clean after unload
.onUnload <- function(libpath) {
  library.dynam.unload("intRinsic", libpath)
}
