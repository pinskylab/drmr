.onAttach <- function(libname, pkgname) {
  tryCatch({
    cmdstan_ver <- instantiate::stan_cmdstan_version()
    if (!is.null(cmdstan_ver) && package_version(cmdstan_ver) < "2.36") {
      packageStartupMessage("Note: Your CmdStan version is outdated. Please update to 2.36 or higher.")
    }
  }, error = function(e) {
    packageStartupMessage("CmdStan not found.")
  })
}
