#' @title Run a one-dimensional burial and thermal history model
#'
#' @description
#' Forward models the burial, compaction and temperature history of a
#' sedimentary column. The basin is described by one row per stratigraphic unit
#' and one row per erosion event; deposition phases, hiatuses, layer stacks,
#' maximum-burial windows and grain thicknesses are all derived from that
#' description, so nothing has to be respecified per time interval.
#'
#' @details
#' Depths are measured from the sediment-water interface, so the water column is
#' not part of the burial depth and a changing water depth does not compact or
#' decompact the column; it only shifts the column relative to sea level.
#'
#' Porosity follows \eqn{\phi(z) = \phi_0 e^{-c z}} and compaction is
#' irreversible: the thickness of a layer is the smaller of the thickness that
#' solves the decompaction equation at its current burial depth and the
#' thickness it had at its deepest burial. An exhumed layer therefore keeps the
#' porosity it acquired at maximum burial. The decompaction equation is solved
#' by iteration rather than on a search grid, so layer boundaries are not quantised.
#'
#' During deposition grain material is added at a constant rate, which makes a
#' unit reach exactly its maximum-burial thickness at the end of its own
#' deposition phase, so no step appears at the phase boundary.
#'
#' Temperatures are the steady-state conductive solution for the layered column.
#' With \code{radiogenic = "cumulative"} the surface heat flow is the basal heat
#' flow plus the radiogenic production of the whole column and the heat flow
#' entering each layer is reduced by the production of the layers above it. With
#' \code{radiogenic = "legacy"} every layer instead receives the basal heat flow
#' plus its own production over the full column thickness.
#'
#' @param site a basin description: a list with the elements \code{units},
#'   \code{erosions}, \code{water_depth}, \code{heat_flow}, \code{surface_temp}
#'   and \code{dz}, and optionally \code{name} and \code{source_units}. Defaults to \code{ardennes}.
#' @param units data frame of stratigraphic units, oldest first. Columns:
#'   \code{name}, \code{start} and \code{end}, \code{thickness}, \code{phi0},
#'   \code{ck}, \code{K}, \code{A} and optionally \code{rho}. Defaults to
#'   \code{site$units}.
#' @param erosions data frame of erosion events with \code{start}, \code{end}
#'   and \code{amount}, or \code{NULL}. Defaults to \code{site$erosions}.
#' @param water_depth water depth in m as a two column matrix or data frame,
#'   a csv path, or a single number. Defaults to \code{site$water_depth}.
#' @param heat_flow heat flow into the base of the column in mW/m2, in the same
#'   formats. Defaults to \code{site$heat_flow}.
#' @param surface_temp surface temperature in degrees Celsius, as a two column
#'   matrix or data frame, the path of a csv file, or a single number.
#'   Defaults to \code{site$surface_temp}.
#' @param dz depth step of the returned temperature raster (m). Defaults to \code{site$dz}.
#' @param k_water thermal conductivity of pore water (W/m/K).
#' @param radiogenic how radiogenic heat is carried through the column, either
#'   \code{"cumulative"} or \code{"legacy"}.
#' @param verbose print a one line summary of the run.
#' @param plot whether to plot the results.
#'
#' @return
#' An object of class \code{"basin_model"}.
#'
#' @examples
#' condroz_example <- list(
#'   name = "Huy (Condroz inlier)",
#'   units = data.frame(
#'     name = c(
#'       "layer 1", "Mousty eq. (unit 2)", "layer 3",
#'       "Huy-Vitrival-Bruyere Fms. (unit 5)", "layer 5", "layer 6",
#'       "layer 7", "layer 8"
#'     ),
#'     start = c(520, 497, 485, 467, 450, 440, 433, 400),
#'     end = c(497, 485, 480, 450, 440, 433, 425, 300),
#'     thickness = c(1.50, 1.00, 0.05, 0.50, 0.05, 0.50, 0.50, 3.50),
#'     phi0 = c(0.490, 0.602, 0.595, 0.595, 0.602, 0.602, 0.602, 0.419),
#'     ck = c(0.270, 0.462, 0.450, 0.450, 0.462, 0.462, 0.462, 0.399),
#'     K = c(5.500, 2.700, 2.875, 2.875, 2.700, 2.700, 2.700, 3.720),
#'     A = c(1.200, 1.600, 1.575, 1.575, 1.600, 1.600, 1.600, 0.850),
#'     rho = c(2650, 2720, 2702.5, 2702.5, 2720, 2720, 2720, 2680),
#'     stringsAsFactors = FALSE
#'   ),
#'   erosions = data.frame(
#'     start = c(425, 300, 230, 45),
#'     end = c(400, 230, 45, 0),
#'     amount = c(0.5, 3.0, 0.5, 0.5)
#'   ),
#'   water_depth = cbind(
#'     age = c(
#'       520, 497, 485, 480, 467, 450, 440, 425,
#'       400, 325, 310, 300, 275, 230, 0
#'     ),
#'     m = c(
#'       0, 50, 250, 0, 0, 200, 100, 0,
#'       0, 250, 0, 0, 0, 0, 0
#'     )
#'   ),
#'   heat_flow = cbind(
#'     age = c(
#'       520, 497, 485, 480, 467, 450, 440, 425,
#'       400, 325, 300, 275, 230, 0
#'     ),
#'     mWm2 = c(
#'       80, 70, 60, 50, 50, 60, 60, 70,
#'       60, 70, 70, 60, 60, 63
#'     )
#'   ),
#'   surface_temp = data.frame(
#'     age = seq(0, 520, by = 10),
#'     temp = c(
#'       11.0411784, 0.1500778, 3.8685811, 6.5473048, 12.1965383,
#'       17.0861816, 17.1607056, 20.4808655, 21.1662140, 23.6852976,
#'       25.2335841, 23.0593262, 19.7055817, 22.0795466, 20.6099243,
#'       20.9470520, 20.6926142, 15.7622960, 16.8022995, 20.8179728,
#'       23.1166636, 26.5647456, 26.6191177, 38.6240616, 39.1549352,
#'       39.2777100, 33.0727590, 26.5889257, 27.4612198, 24.3741633,
#'       27.8295568, 27.3116150, 31.9937617, 33.6753845, 27.6599223,
#'       25.8550212, 28.2169037, 30.6980921, 33.2465744, 30.7135722,
#'       28.8724314, 29.5901311, 24.2159042, 28.4941711, 20.8871002,
#'       16.4026972, 11.2207235, 15.0190506, 13.3535004, 11.9391505,
#'       6.1531270, 7.4877523, 6.3858999
#'     )
#'   ),
#'   source_units = c(2L, 4L)
#' )
#'
#' res <- run_basin_model(condroz_example)
#'
#' @export

run_basin_model <- function(site,
                            units = site$units,
                            erosions = site$erosions,
                            water_depth = site$water_depth,
                            heat_flow = site$heat_flow,
                            surface_temp = site$surface_temp,
                            dz = 10,
                            k_water = 0.6,
                            radiogenic = c("cumulative", "legacy"),
                            verbose = FALSE,
                            plot = FALSE) {

  radiogenic <- match.arg(radiogenic)

  need <- c("name", "start", "end", "thickness", "phi0", "ck", "K", "A")
  miss <- setdiff(need, names(units))
  if (length(miss)) {
    stop("the unit table is missing the column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  if (any(units$start <= units$end)) {
    stop("every unit needs start larger than end", call. = FALSE)
  }
  if (any(units$thickness <= 0)) {
    stop("thickness must be positive", call. = FALSE)
  }
  if (any(units$phi0 <= 0 | units$phi0 >= 1)) {
    stop("phi0 must lie between 0 and 1", call. = FALSE)
  }
  if (any(units$ck <= 0)) {
    stop("ck must be positive", call. = FALSE)
  }
  if (any(units$K <= 0)) {
    stop("K must be positive", call. = FALSE)
  }

  units <- units[order(-units$start), , drop = FALSE]
  rownames(units) <- NULL
  n <- nrow(units)

  thk  <- units$thickness * 1000
  phi0 <- units$phi0
  cc   <- units$ck / 1000
  Kgr  <- units$K
  Arad <- units$A * 1e-6

  ph <- list()
  for (i in seq_len(nrow(units))) {
    ph[[length(ph) + 1L]] <- list(type = "deposition", start = units$start[i],
                                  end = units$end[i], unit = i, amount = 0)
  }
  if (!is.null(erosions) && nrow(erosions) > 0) {
    for (i in seq_len(nrow(erosions))) {
      ph[[length(ph) + 1L]] <- list(type = "erosion", start = erosions$start[i],
                                    end = erosions$end[i], unit = NA_integer_,
                                    amount = erosions$amount[i])
    }
  }
  ph <- ph[order(-vapply(ph, function(p) p$start, 0))]

  phases <- list()
  for (k in seq_along(ph)) {
    p <- ph[[k]]
    if (p$start <= p$end) {
      stop("phase start must be the older age", call. = FALSE)
    }
    if (k > 1) {
      prev <- phases[[length(phases)]]$end
      if (p$start > prev) {
        stop("phases overlap", call. = FALSE)
      }
      if (p$start < prev) {
        phases[[length(phases) + 1L]] <- list(type = "hiatus", start = prev,
                                              end = p$start, unit = NA_integer_,
                                              amount = 0)
      }
    }
    phases[[length(phases) + 1L]] <- p
  }

  t_base <- phases[[1]]$start
  t_top  <- phases[[length(phases)]]$end
  ages   <- seq(t_base - 1, t_top)
  nA     <- length(ages)

  idx_of <- function(a) t_base - a

  if (is.null(water_depth)) {
    stop("water_depth is missing", call. = FALSE)
  }
  if (is.character(water_depth)) {
    if (!file.exists(water_depth)) {
      wd <- rep(20, length(ages))
    } else {
      water_depth <- utils::read.csv(water_depth)
      water_depth <- as.matrix(water_depth)
      if (ncol(water_depth) > 2) {
        water_depth <- water_depth[, 1:2, drop = FALSE]
      }
      water_depth <- water_depth[order(water_depth[, 1]), , drop = FALSE]
      wd <- stats::approx(water_depth[, 1], water_depth[, 2],
                          xout = ages, rule = 2)$y
    }
  } else if (is.numeric(water_depth) && is.null(dim(water_depth)) && length(water_depth) == 1L) {
    wd <- rep(water_depth, length(ages))
  } else {
    water_depth <- as.matrix(water_depth)
    if (ncol(water_depth) > 2) {
      water_depth <- water_depth[, 1:2, drop = FALSE]
    }
    water_depth <- water_depth[order(water_depth[, 1]), , drop = FALSE]
    wd <- stats::approx(water_depth[, 1], water_depth[, 2],
                        xout = ages, rule = 2)$y
  }

  if (is.null(heat_flow)) {
    stop("heat_flow is missing", call. = FALSE)
  }
  if (is.character(heat_flow)) {
    if (!file.exists(heat_flow)) {
      qbase <- rep(20e-3, length(ages))
    } else {
      heat_flow <- utils::read.csv(heat_flow)
      heat_flow <- as.matrix(heat_flow)
      if (ncol(heat_flow) > 2) {
        heat_flow <- heat_flow[, 1:2, drop = FALSE]
      }
      heat_flow <- heat_flow[order(heat_flow[, 1]), , drop = FALSE]
      qbase <- stats::approx(heat_flow[, 1], heat_flow[, 2] * 1e-3,
                             xout = ages, rule = 2)$y
    }
  } else if (is.numeric(heat_flow) && is.null(dim(heat_flow)) && length(heat_flow) == 1L) {
    qbase <- rep(heat_flow * 1e-3, length(ages))
  } else {
    heat_flow <- as.matrix(heat_flow)
    if (ncol(heat_flow) > 2) {
      heat_flow <- heat_flow[, 1:2, drop = FALSE]
    }
    heat_flow <- heat_flow[order(heat_flow[, 1]), , drop = FALSE]
    qbase <- stats::approx(heat_flow[, 1], heat_flow[, 2] * 1e-3,
                           xout = ages, rule = 2)$y
  }

  if (is.null(surface_temp)) {
    stop("surface_temp is missing", call. = FALSE)
  }
  if (is.character(surface_temp)) {
    if (!file.exists(surface_temp)) {
      tsurf <- rep(20, length(ages))
    } else {
      surface_temp <- utils::read.csv(surface_temp)
      surface_temp <- as.matrix(surface_temp)
      if (ncol(surface_temp) > 2) {
        surface_temp <- surface_temp[, 1:2, drop = FALSE]
      }
      surface_temp <- surface_temp[order(surface_temp[, 1]), , drop = FALSE]
      tsurf <- stats::approx(surface_temp[, 1], surface_temp[, 2],
                             xout = ages, rule = 2)$y
    }
  } else if (is.numeric(surface_temp) && is.null(dim(surface_temp)) && length(surface_temp) == 1L) {
    tsurf <- rep(surface_temp, length(ages))
  } else {
    surface_temp <- as.matrix(surface_temp)
    if (ncol(surface_temp) > 2) {
      surface_temp <- surface_temp[, 1:2, drop = FALSE]
    }
    surface_temp <- surface_temp[order(surface_temp[, 1]), , drop = FALSE]
    tsurf <- stats::approx(surface_temp[, 1], surface_temp[, 2],
                           xout = ages, rule = 2)$y
  }

  wd[wd < 0] <- 0

  rem <- rep(0, n)
  stack <- integer(0)
  zt0 <- rep(NA_real_, n)
  zb0 <- rep(NA_real_, n)

  for (p in phases) {
    if (p$type == "deposition") {
      stack <- c(p$unit, stack)
      rem[p$unit] <- thk[p$unit]
    }

    if (p$type == "erosion") {
      left <- p$amount * 1000

      while (left > 1e-9 && length(stack) > 0) {
        u <- stack[1]
        take <- min(left, rem[u])
        rem[u] <- rem[u] - take
        left <- left - take

        if (rem[u] <= 1e-9) {
          stack <- stack[-1]
        }
      }
    }

    z <- 0
    for (u in stack) {
      z1 <- z
      z <- z + rem[u]

      if (is.na(zb0[u]) || z > zb0[u]) {
        zb0[u] <- z
        zt0[u] <- z1
      }
    }
  }

  pore_volume <- function(phi0, cc, z1, z2) {
    (phi0 / cc) * (exp(-cc * z1) - exp(-cc * z2))
  }

  solid_full <- (zb0 - zt0) - pore_volume(phi0, cc, zt0, zb0)

  if (any(is.na(solid_full)) || any(solid_full <= 0)) {
    stop("zero or negative grain thickness detected", call. = FALSE)
  }

  hmin   <- rep(Inf, n)
  zt     <- rep(NA_real_, n)
  zb     <- rep(NA_real_, n)
  solid  <- rep(0, n)
  active <- rep(FALSE, n)
  stack  <- integer(0)

  layers <- vector("list", nA)
  tot    <- numeric(nA)

  for (p in phases) {
    steps <- seq(p$start - 1, p$end)
    N <- length(steps)

    for (j in seq_len(N)) {
      age <- steps[j]
      k <- idx_of(age)

      if (p$type == "deposition") {
        u <- p$unit
        if (!active[u]) {
          active[u] <- TRUE
          stack <- c(u, stack)
        }
        solid[u] <- (j / N) * solid_full[u]
        hmin[u] <- Inf

      } else if (p$type == "erosion") {
        left <- (p$amount * 1000) / N

        while (left > 1e-9 && length(stack) > 0) {
          u <- stack[1]
          if (!is.finite(hmin[u])) {
            break
          }
          take <- min(left, hmin[u])
          zt[u] <- zt[u] + take
          hmin[u] <- zb[u] - zt[u]
          solid[u] <- hmin[u] - pore_volume(phi0[u], cc[u], zt[u], zb[u])
          left <- left - take

          if (hmin[u] <= 1e-6) {
            active[u] <- FALSE
            stack <- stack[-1]
          }
        }
      }

      ns <- length(stack)
      if (ns == 0) {
        tot[k] <- 0
        next
      }

      z_top <- numeric(ns)
      z_bas <- numeric(ns)
      hh <- numeric(ns)
      por <- numeric(ns)
      z <- 0

      for (i in seq_len(ns)) {
        u <- stack[i]
        a <- (phi0[u] / cc[u]) * exp(-cc[u] * z)
        h <- solid[u]

        for (iter in seq_len(500)) {
          hn <- solid[u] + a * (1 - exp(-cc[u] * h))
          if (abs(hn - h) < 1e-10) {
            h <- hn
            break
          }
          h <- hn
        }

        hc <- h
        if (hc < hmin[u]) {
          hmin[u] <- hc
          zt[u] <- z
          zb[u] <- z + hc
        }

        hh[i] <- hmin[u]
        z_top[i] <- z
        z <- z + hh[i]
        z_bas[i] <- z

        por[i] <- if (hh[i] > 0) {
          pore_volume(phi0[u], cc[u], zt[u], zb[u]) / hh[i]
        } else {
          0
        }
      }

      tot[k] <- z
      us <- stack
      Kb <- (1 - por) * Kgr[us] + por * k_water
      Ai <- Arad[us]

      Tt <- numeric(ns)
      Tb <- numeric(ns)
      Tm <- numeric(ns)
      qt <- numeric(ns)

      qs <- qbase[k] + sum(Ai * hh)
      Tp <- tsurf[k]

      for (i in seq_len(ns)) {
        qi <- if (radiogenic == "cumulative") qs else qbase[k] + Ai[i] * z
        qt[i] <- qi
        Tt[i] <- Tp
        Tb[i] <- Tp + qi * hh[i] / Kb[i] - Ai[i] * hh[i]^2 / (2 * Kb[i])
        Tm[i] <- Tp + qi * (hh[i] / 2) / Kb[i] - Ai[i] * (hh[i] / 2)^2 / (2 * Kb[i])

        if (radiogenic == "cumulative") {
          qs <- qs - Ai[i] * hh[i]
        }
        Tp <- Tb[i]
      }

      layers[[k]] <- data.frame(
        age = age,
        unit = us,
        name = units$name[us],
        z_top = z_top,
        z_base = z_bas,
        z_mid = z_top + hh / 2,
        h = hh,
        porosity = por,
        K = Kb,
        A = Ai,
        q_top = qt,
        T_top = Tt,
        T_base = Tb,
        T_mid = Tm,
        water_depth = wd[k],
        stringsAsFactors = FALSE
      )
    }
  }

  profile <- do.call(rbind, layers[!vapply(layers, is.null, TRUE)])
  zmax <- max(tot + wd)
  depth_grid <- seq(0, ceiling(zmax / dz) * dz, by = dz)
  temp_grid <- matrix(NA_real_, nrow = nA, ncol = length(depth_grid))

  for (k in seq_len(nA)) {
    L <- layers[[k]]
    v <- rep(NA_real_, length(depth_grid))
    v[depth_grid <= wd[k]] <- tsurf[k]

    if (!is.null(L)) {
      zz <- depth_grid - wd[k]
      for (i in seq_len(nrow(L))) {
        sel <- zz > L$z_top[i] & zz <= L$z_base[i]
        if (any(sel)) {
          d <- zz[sel] - L$z_top[i]
          v[sel] <- L$T_top[i] + L$q_top[i] * d / L$K[i] - L$A[i] * d^2 / (2 * L$K[i])
        }
      }
    }
    temp_grid[k, ] <- v
  }

  res <- list(
    units = units,
    phases = phases,
    ages = ages,
    water_depth = wd,
    heat_flow = qbase * 1000,
    surface_temp = tsurf,
    total_thickness = tot,
    profile = profile,
    depth_grid = depth_grid,
    temp_grid = temp_grid,
    radiogenic = radiogenic,
    site_name = if (is.null(site$name)) NA_character_ else site$name,
    source_units = site$source_units
  )

  class(res) <- "basin_model"

  if (verbose) {
    cat("1-D basin model", if (!is.na(res$site_name)) paste0(": ", res$site_name) else "", "\n", sep = "")
  }

  invisible(res)
}


