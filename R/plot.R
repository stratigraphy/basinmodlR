#' @title Composite burial and temperature figure
#'
#' @description
#' Draws the standard figure of a model run: a geological period bar, an event
#' bar, the temperature raster with the water column, the horizon lines and
#' labelled isotherms, a colour key, and one generation-window panel per
#' selected unit. The function also calculates and prints phase boundary
#' diagnostics internally, comparing the change in layer depth and temperature
#' over the single million years that straddles every phase boundary with the
#' background rate on either side of it.
#'
#' @param res an object of class \code{"basin_model"}, from
#'   \code{\link{run_basin_model}}.
#' @param window numeric of length two, the temperature window in degrees
#'   Celsius used for the generation panels.
#' @param window_units which units get a generation panel, given as row numbers
#'   of \code{res$units} or as unit names. \code{NULL} uses
#'   \code{res$source_units} when the basin description provides it, and draws
#'   no panel otherwise.
#' @param n_levels number of colour levels in the raster.
#' @param scale how the colours are spread over the temperature range, either
#'   \code{"quantile"} (equalised, the default) or \code{"linear"}.
#' @param contour_n approximate number of isotherms.
#' @param xmax oldest age shown (Ma). \code{NULL} uses the start of the model.
#' @param tick spacing of the age axis labels (Myr).
#' @param new_device open a new graphics window before plotting. Leave this
#'   \code{FALSE} when writing to a file device.
#' @param use the temperature used for the generation window test, one of
#'   \code{"T_mid"} (the middle of the layer, the default), \code{"T_top"} or
#'   \code{"T_base"}.
#' @param span half width, in Myr, of the interval on either side of a boundary
#'   used to measure the background rate for continuity checks.
#'
#' @return \code{res}, invisibly. Called for its side effect, the figure and diagnostic print.
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
#' plot_basin_model(res)
#'
#' @export
#' @importFrom graphics abline
#' @importFrom graphics arrows
#' @importFrom graphics axis
#' @importFrom graphics box
#' @importFrom graphics lines
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics rect

plot_basin_model <- function(res,
                             window = c(95, 135),
                             window_units = NULL,
                             n_levels = 250,
                             scale = c("quantile", "linear"),
                             contour_n = 20,
                             xmax = NULL,
                             tick = 25,
                             new_device = FALSE,
                             use = c("T_mid", "T_top", "T_base"),
                             span = 5) {

  scale <- match.arg(scale)
  use <- match.arg(use)

  turbo_pal <- grDevices::colorRampPalette(c(
    "#30123B", "#4145AB", "#4675ED", "#39A2FC", "#1BCFD4", "#24ECA6", "#61FC6C",
    "#A4FC3B", "#D1E834", "#F3C63A", "#FE9B2D", "#F36315", "#D93806", "#B11901",
    "#7A0403"
  ))

  ics_periods <- data.frame(
    name = c("Quaternary", "Neogene", "Paleogene", "Cretaceous", "Jurassic",
             "Triassic", "Permian", "Carboniferous", "Devonian", "Silurian",
             "Ordovician", "Cambrian"),
    abbr = c(NA, "Neo.", "Pale.", "Cret.", "Jura.", "Tria.", "Perm.", "Carb.",
             "Devo.", "Silu.", "Ordo.", "Camb."),
    base = c(2.58, 23.03, 66, 143.1, 201.4, 251.9, 298.9, 358.9, 419.2, 443.8,
             486.9, 538.8),
    top = c(0, 2.58, 23.03, 66, 143.1, 201.4, 251.9, 298.9, 358.9, 419.2, 443.8,
            486.9),
    col = c("#F9F97F", "#FFE619", "#FD9A52", "#7FC64E", "#34B2C9", "#812B92",
            "#F04028", "#67A599", "#CB8C37", "#B3E1B6", "#009270", "#7FA056"),
    stringsAsFactors = FALSE
  )

  d_cc <- do.call(rbind, lapply(split(res$profile, res$profile$unit), function(x) {
    x <- x[order(-x$age), ]
    if (nrow(x) < 2) return(NULL)
    data.frame(age = x$age[-1], unit = x$unit[-1],
               dz = abs(diff(x$z_base)), dT = abs(diff(x$T_base)))
  }))

  bnd <- vapply(res$phases[-length(res$phases)], function(p) p$end, 0)

  out_cc <- do.call(rbind, lapply(bnd, function(b) {
    at <- d_cc[d_cc$age == b - 1, ]
    nb <- d_cc[abs(d_cc$age - (b - 1)) <= span & d_cc$age != b - 1, ]
    data.frame(boundary_Ma = b,
               jump_z = if (nrow(at)) max(at$dz) else NA_real_,
               bg_z = if (nrow(nb)) stats::median(nb$dz) else NA_real_,
               jump_T = if (nrow(at)) max(at$dT) else NA_real_,
               bg_T = if (nrow(nb)) stats::median(nb$dT) else NA_real_)
  }))

  if (!is.null(out_cc) && nrow(out_cc) > 0) {
    out_cc$ratio_z <- round(out_cc$jump_z / pmax(out_cc$bg_z, 1e-9), 2)
    out_cc$ratio_T <- round(out_cc$jump_T / pmax(out_cc$bg_T, 1e-9), 2)
    cat("depth and temperature step at each phase boundary versus the background rate\n")
    print(out_cc, row.names = FALSE, digits = 3)
  }

  if (is.null(xmax)) xmax <- max(res$ages) + 1

  xf <- function(a) xmax - a
  o <- order(-res$ages)
  ages <- res$ages[o]
  xx <- xf(ages)
  z <- res$temp_grid[o, , drop = FALSE]
  wd <- res$water_depth[o]
  dg <- res$depth_grid
  xtick <- seq(0, floor(xmax / tick) * tick, by = tick)

  if (is.null(window_units)) window_units <- res$source_units
  if (is.character(window_units)) {
    window_units <- match(window_units, res$units$name)
  }
  window_units <- window_units[!is.na(window_units)]
  nw <- length(window_units)

  if (scale == "quantile") {
    br <- unique(stats::quantile(
      z,
      probs = seq(0, 1, length.out = n_levels + 1),
      na.rm = TRUE
    ))
  } else {
    br <- seq(
      min(z, na.rm = TRUE),
      max(z, na.rm = TRUE),
      length.out = n_levels + 1
    )
  }

  cols <- turbo_pal(length(br) - 1)

  if (new_device) grDevices::dev.new(width = 10, height = 9)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  m <- rbind(c(1, 0), c(2, 0), c(3, 4))
  if (nw > 0) {
    for (i in seq_len(nw)) {
      m <- rbind(m, c(4 + i, 0))
    }
  }

  graphics::layout(
    m,
    heights = c(1, 1, 10, rep(1, nw)),
    widths = c(5, 1)
  )

  graphics::par(mar = c(0.5, 4, 0.5, 0.5))
  plot(
    NA,
    xlim = c(0, xmax),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n"
  )

  graphics::mtext("Period", side = 2, line = 1, las = 2, cex = 0.75)

  for (i in seq_len(nrow(ics_periods))) {
    b <- min(ics_periods$base[i], xmax)
    tp <- ics_periods$top[i]

    if (tp >= xmax) next

    graphics::rect(
      xf(b), 0,
      xf(tp), 1,
      col = ics_periods$col[i]
    )

    if (!is.na(ics_periods$abbr[i]) && (b - tp) > xmax / 40) {
      graphics::text(
        xf((b + tp) / 2),
        0.5,
        ics_periods$abbr[i],
        cex = 0.8
      )
    }
  }

  graphics::box()

  plot(
    NA,
    xlim = c(0, xmax),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n"
  )

  graphics::mtext("Event", side = 2, line = 1, las = 2, cex = 0.75)

  ecol <- c(
    deposition = "gold",
    erosion = "grey80",
    hiatus = "white"
  )

  for (i in seq_along(res$phases)) {
    p <- res$phases[[i]]

    graphics::rect(
      xf(p$start), 0,
      xf(p$end), 1,
      col = ecol[[p$type]]
    )

    graphics::text(
      xf((p$start + p$end) / 2),
      0.5,
      i,
      cex = 0.8
    )
  }

  graphics::box()

  graphics::par(mar = c(4, 4, 1, 0.5))

  max_plot_depth <- max(c(dg, res$profile$z_base + res$profile$water_depth), na.rm = TRUE)

  graphics::image(
    x = xx,
    y = dg,
    z = z,
    col = cols,
    breaks = br,
    useRaster = TRUE,
    xlim = c(0, xmax),
    ylim = c(max_plot_depth, 0),
    xaxt = "n",
    xlab = "Age (Ma)",
    ylab = "Depth below sea level (m)"
  )

  graphics::axis(
    1,
    at = xf(xtick),
    labels = xtick,
    las = 2
  )

  graphics::polygon(
    c(xx, rev(xx)),
    c(rep(0, length(xx)), rev(wd)),
    col = "lightblue2",
    border = NA
  )

  graphics::lines(xx, wd, lwd = 2)

  for (u in sort(unique(res$profile$unit))) {
    d <- res$profile[res$profile$unit == u, ]
    d <- d[order(-d$age), ]

    graphics::lines(
      xf(d$age),
      d$z_base + d$water_depth,
      lwd = 2
    )
  }

  graphics::contour(
    x = xx,
    y = dg,
    z = z,
    add = TRUE,
    drawlabels = TRUE,
    labcex = 0.6,
    method = "flattest",
    col = "grey15",
    levels = pretty(range(z, na.rm = TRUE), contour_n)
  )

  graphics::box()

  graphics::par(mar = c(4, 1, 1, 4))

  nb <- length(br) - 1

  graphics::image(
    x = 1,
    y = seq_len(nb),
    z = matrix(seq_len(nb), nrow = 1),
    col = cols,
    useRaster = TRUE,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = ""
  )

  at <- round(seq(1, nb, length.out = 6))

  graphics::axis(
    4,
    at = at,
    labels = round(br[at + 1]),
    las = 2
  )

  graphics::mtext(
    "Temperature (C)",
    side = 4,
    line = 2.5,
    cex = 0.75
  )

  graphics::box()

  if (nw > 0) {
    mw_ages <- res$ages
    present <- sort(unique(res$profile$unit))
    mw <- matrix(0, nrow = length(mw_ages), ncol = length(present),
                 dimnames = list(NULL, res$units$name[present]))

    for (j in seq_along(present)) {
      d_mw <- res$profile[res$profile$unit == present[j], ]
      d_mw <- d_mw[order(-d_mw$age), ]
      tt <- d_mw[[use]]
      peak <- cummax(tt)
      hot <- tt >= peak - 1e-9 & tt > window[1] & tt < window[2]
      mw[match(d_mw$age, mw_ages), j] <- as.numeric(hot)
    }

    graphics::par(mar = c(0.5, 4, 0.5, 0.5))

    for (u in window_units) {
      plot(
        NA,
        xlim = c(0, xmax),
        ylim = c(-0.1, 1.2),
        xaxs = "i",
        yaxs = "i",
        xlab = "",
        ylab = "",
        xaxt = "n",
        yaxt = "n"
      )

      j <- which(present == u)

      if (length(j)) {
        graphics::lines(
          xf(res$ages),
          mw[, j]
        )
      }

      graphics::axis(2, at = c(0, 1), labels = NA)

      graphics::mtext(
        c("no-gen", "gen"),
        side = 2,
        at = c(0, 1),
        line = 0.5,
        las = 2,
        cex = 0.7
      )

      graphics::text(
        xf(xmax / 3),
        0.5,
        res$units$name[u],
        cex = 0.9
      )

      graphics::box()
    }
  }

  invisible(res)
}

