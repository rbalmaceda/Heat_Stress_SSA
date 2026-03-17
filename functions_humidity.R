
tdps2huss <- function(tdps, ps,formula ="bohren") {
  
  if (isMultigrid(tdps) | isMultigrid(ps)) 
    stop("Multigrids are not an allowed input")
  
  stopifnot(isGrid(tdps), isGrid(ps))
  
  tdps %<>% redim(member = TRUE)
  ps   %<>% redim(member = TRUE)
  
  suppressMessages(checkDim(tdps, ps))
  checkSeason(tdps, ps)
  
  # --- Units check ---
  # tdps must be degC
  u.td <- getGridUnits(tdps)
  if (u.td != "degC") {
    if (!ud.are.convertible(u.td, "degC")) {
      stop("Non compliant tdps units (should be convertible to degC)")
    }
    message("[", Sys.time(), "] Converting tdps units ...")
    tdps <- udConvertGrid(tdps, new.units = "degC") %>% redim(member = TRUE)
  }
  # ps must be Pa
  u.ps <- getGridUnits(ps)
  if (u.ps != "Pa") {
    if (!ud.are.convertible(u.ps, "Pa")) {
      stop("Non compliant ps units (should be convertible to Pa)")
    }
    message("[", Sys.time(), "] Converting ps units ...")
    ps <- udConvertGrid(ps, new.units = "Pa") %>% redim(member = TRUE)
  }
  
  coords <- getCoordinates(tdps)
  n.mem  <- getShape(tdps, "member")
  
  message("[", Sys.time(), "] Calculating specific humidity from dew point ...")
  
  epsilon <- 0.621981  # Rv/Rd
  
  huss <- tdps
  
  l <- lapply(1:n.mem, function(x) {
    
    td <- subsetGrid(tdps, members = x, drop = TRUE) %>% 
      redim(member = FALSE) %>% 
      extract2("Data") %>% 
      array3Dto2Dmat()
    
    p <- subsetGrid(ps, members = x, drop = TRUE) %>% 
      redim(member = FALSE) %>% 
      extract2("Data") %>% 
      array3Dto2Dmat()
    
    # Td in °C → K
    T <- td + 273.15  
    
    T0 <- 273.16  # triple point (K)
    
    iceMask   <- which(T < T0)
    waterMask <- which(T >= T0)
    
    e <- T  
    if (formula=='bohren'){
      # 
      es0 <- 611 # Pa
      e[iceMask] <- es0 * exp((6293/T0) - (6293/e[iceMask]) -
                                0.555 * log(abs(e[iceMask]/T0)))
      e[waterMask] <- es0 * exp((6808/T0) - (6808/e[waterMask]) -
                                  5.09 * log(abs(e[waterMask]/T0)))
    }
    # e is now vapor pressure in hPa 
    
    # ---- specific humidity formulation) ----
    w <- epsilon * e / (p - e)
    q <- w / (1 + w)
   # q[q < 0] <- 0
    huss$Data <- mat2Dto3Darray(q, x = coords$x, y = coords$y)
    return(huss)
  })
  
  huss <- suppressWarnings(bindGrid(l, dimension = "member"))
  
  # ---- metadata ----
  huss$Variable$varName <- "huss"
  huss$Variable$level <- NULL
  attr(huss$Variable, "units") <- "kg/kg"
  attr(huss$Variable, "longname") <- "specific_humidity"
  attr(huss$Variable, "description") <- "Specific humidity computed from dew-point temperature and surface pressure using WMO08 formulation"
  attr(huss, "origin") <- "climate4R-convertR adaptation"
  
  message("[", Sys.time(), "] Done.")
  invisible(huss)
}

###################################

huss2hurs=function (huss, ps, tas) 
{
  if (isMultigrid(huss) | isMultigrid(ps) | isMultigrid(tas)) 
    stop("Multigrids are not an allowed input")
  stopifnot(isGrid(huss) | isGrid(ps) | isGrid(tas))
  huss %<>% redim(member = TRUE)
  ps %<>% redim(member = TRUE)
  tas %<>% redim(member = TRUE)
  suppressMessages(checkDim(huss, ps, tas))
  checkSeason(huss, ps, tas)
  u.huss <- getGridUnits(huss)
  if (u.huss != "kg/kg") {
    if (!ud.are.convertible(u.huss, "kg/kg")) {
      stop("Non compliant huss units (should be convertible to 'kg/kg')")
    }
    message("[", Sys.time(), "] Converting units ...")
    huss <- udConvertGrid(huss, new.units = "kg/kg")
  }
  huss %<>% redim(member = TRUE)
  message("[", Sys.time(), "] Calculating Relative humidity ...")
  ws <- suppressMessages(tas2ws(tas, ps)) %>% redim(member = TRUE)
  tas <- NULL
  coords <- getCoordinates(huss)
  n.mem <- getShape(huss, "member")
  l <- lapply(1:n.mem, function(x) {
    a <- subsetGrid(huss, members = x, drop = TRUE) %>% redim(member = FALSE) %>% 
      extract2("Data") %>% array3Dto2Dmat()
    b <- subsetGrid(ws, members = x, drop = TRUE) %>% redim(member = FALSE) %>% 
      extract2("Data") %>% array3Dto2Dmat()
    w <- a/(1 - a)
    a <- 100 * w/b
    a[a > 100] <- 100
    a[a < 0] <- 0
    ps$Data <- mat2Dto3Darray(a, x = coords$x, y = coords$y)
    a <- b <- NULL
    return(ps)
  })
  huss <- ws <- ps <- NULL
  hurs <- suppressWarnings(bindGrid(l, dimension = "member"))
  hurs$Variable$varName <- "hurs"
  hurs$Variable$level <- NULL
  attr(hurs$Variable, "units") <- "%"
  attr(hurs$Variable, "longname") <- "Surface_air_relative_humidity"
  attr(hurs$Variable, "description") <- "Estimated relative humidity from saturation pressure and specific humidity"
  attr(hurs, "origin") <- paste0("Calculated with R package 'convertR' v", 
                                 packageVersion("convertR"))
  attr(hurs, "URL") <- "https://github.com/SantanderMetGroup/convertR"
  message("[", Sys.time(), "] Done.")
  invisible(hurs)
}

#################################
tas2ws=function (tas, ps) 
{
  if (isMultigrid(ps) | isMultigrid(tas)) 
    stop("Multigrids are not an allowed input")
  stopifnot(isGrid(ps) | isGrid(tas))
  ps %<>% redim(member = TRUE)
  tas %<>% redim(member = TRUE)
  suppressMessages(checkDim(ps, tas))
  checkSeason(ps, tas)
  u.ps <- getGridUnits(ps)
  u.tas <- getGridUnits(tas)
  if (u.ps != "Pa") {
    if (!ud.are.convertible(u.ps, "Pa")) {
      stop("Non compliant ps units (should be convertible to 'Pascals')")
    }
    message("[", Sys.time(), "] Converting units ...")
    ps <- udConvertGrid(ps, new.units = "Pa") %>% redim(member = TRUE)
  }
  if (u.tas != "K") {
    if (!ud.are.convertible(u.tas, "K")) {
      stop("Non compliant tas units (should be convertible to 'Kelvin')")
    }
    message("[", Sys.time(), "] Converting units ...")
    tas <- udConvertGrid(tas, new.units = "K") %>% redim(member = TRUE)
  }
  message("[", Sys.time(), "] Calculating saturation pressure ...")
  source(file.path(find.package(package = "convertR"), "constants.R"), 
         local = TRUE)
  coords <- getCoordinates(ps)
  n.mem <- getShape(tas, "member")
  wt <- ps
  l <- lapply(1:n.mem, function(x) {
    a <- subsetGrid(tas, members = x, drop = TRUE) %>% redim(member = FALSE) %>% 
      extract2("Data") %>% array3Dto2Dmat()
    b <- subsetGrid(ps, members = x, drop = TRUE) %>% redim(member = FALSE) %>% 
      extract2("Data") %>% array3Dto2Dmat()
    iceMask <- which(a < T0)
    waterMask <- which(a >= T0)
    a[iceMask] <- es0 * exp((6293/T0) - (6293/a[iceMask]) - 
                              0.555 * log(abs(a[iceMask]/T0)))
    a[waterMask] <- es0 * exp((6808/T0) - (6808/a[waterMask]) - 
                                5.09 * log(abs(a[waterMask]/T0)))
    a <- (Rd/Rv) * (a/(b - a))
    wt$Data <- mat2Dto3Darray(a, x = coords$x, y = coords$y)
    a <- b <- NULL
    return(wt)
  })
  tas <- ps <- NULL
  wt <- suppressWarnings(bindGrid(l, dimension = "member"))
  wt$Variable$varName <- "ws"
  wt$Variable$level <- NULL
  attr(wt$Variable, "units") <- "atm"
  attr(wt$Variable, "longname") <- "water_vapour_saturation_pressure"
  attr(wt$Variable, "description") <- "Estimated water vapour saturation pressure"
  attr(wt, "origin") <- paste0("Calculated with R package 'convertR' v", 
                               packageVersion("convertR"))
  attr(wt, "URL") <- "https://github.com/SantanderMetGroup/convertR"
  message("[", Sys.time(), "] Done.")
  invisible(wt)
}
