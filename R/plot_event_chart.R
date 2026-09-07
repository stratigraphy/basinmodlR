#' Plot a petroleum-system event chart from a basin-model result
#'
#' @description
#' Draws a petroleum-system event chart from the output of
#' [run_basin_model()]. The chart shows the timing of source-rock
#' deposition, reservoir and seal deposition, overburden development,
#' uplift/erosion, optional trap formation, generation/migration/accumulation,
#' preservation, and the critical moment.
#'
#' Source rocks are displayed as separate bars and may be assigned
#' different shades of red. Generation/migration/accumulation is displayed
#' as one sub-row per source rock, allowing a single source rock to have
#' multiple separate generation events without creating additional rows.
#'
#' Overburden begins after the first specified source rock has finished being
#' deposited. From that point onward, every subsequent deposition phase is
#' represented as overburden, regardless of whether the deposited material is
#' source rock, reservoir rock, seal rock, or another unit. Hiatuses, uplift,
#' and erosion interrupt the overburden record and therefore appear as gaps.
#'
#' The generation/migration/accumulation interval is calculated from the
#' modelled temperature history of each specified source rock. Only periods
#' within \code{generation_window} during which temperature reaches or exceeds
#' all previous temperatures in the source-rock history are plotted.
#'
#' Preservation is plotted from the supplied critical moment to the present
#' only when source rock, reservoir rock, seal rock, and trap formation have
#' all been specified.
#'
#' @param res An object returned by [run_basin_model()], containing the model
#'   units, geological phases, ages, and temperature profiles.
#' @param source_units Source-rock units. May be supplied as unit indices or
#'   as unit names corresponding to \code{res$units$name}.
#' @param reservoir_units Reservoir-rock units. May be supplied as unit
#'   indices or as unit names corresponding to \code{res$units$name}.
#' @param seal_units Seal-rock units. May be supplied as unit indices or as
#'   unit names corresponding to \code{res$units$name}.
#' @param trap_formation Optional data frame defining trap-formation
#'   intervals. The data frame must contain columns named \code{start} and
#'   \code{end}. Defaults to \code{NULL}.
#' @param generation_window Numeric vector of length two defining the
#'   temperature window used to identify generation/migration/accumulation,
#'   in degrees C. Defaults to \code{c(95, 135)}.
#' @param critical_moment Numeric value giving the age of the petroleum-system
#'   critical moment in Ma.
#' @param temperature_column Character string giving the name of the
#'   temperature column in \code{res$profile}. Defaults to \code{"T_mid"}.
#' @param xmax Optional maximum age shown on the x-axis. If \code{NULL}, the
#'   maximum model age in \code{res$ages} is used.
#' @param tick Spacing between x-axis age ticks in Ma. Defaults to 25.
#' @param col Default colour used for event bars. Defaults to
#'   \code{"darkred"}.
#' @param critical_col Colour used for the critical-moment marker.
#'   Defaults to \code{"black"}.
#' @param main Plot title. Defaults to an empty string.
#' @param xlab Label for the x-axis. Defaults to \code{"Age (Ma)"}.
#'
#' @return
#' Invisibly returns a named list containing the event intervals used to
#' construct the chart.
#'
#' @details
#' Ages are plotted from old to young, with the present at the right-hand
#' side of the chart. Internally, the x-axis is transformed so that increasing
#' plotted x-position corresponds to decreasing geological age.
#'
#' Source-rock bars are kept separate even when their depositional intervals
#' overlap. If multiple source rocks overlap, they are stacked within the
#' Source Rock row.
#'
#' Generation events are handled differently: all events belonging to the
#' same source rock are drawn in the same generation sub-row. Consequently,
#' repeated generation or expulsion events from one source rock do not create
#' additional source-rock rows.
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
#' plot_event_chart(
#'   res,
#'   source_units = c(2, 4),
#'   reservoir_units = c(5),
#'   seal_units = c(5),
#'   trap_formation = data.frame(start = 320, end = 295),
#'   critical_moment = 310
#' )
#'
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
#' @importFrom grDevices colorRampPalette

plot_event_chart <- function(
    res,
    source_units,
    reservoir_units,
    seal_units,
    trap_formation = NULL,
    generation_window = c(95, 135),
    critical_moment,
    temperature_column = "T_mid",
    xmax = NULL,
    tick = 25,
    col = "darkred",
    critical_col = "black",
    main = "",
    xlab = "Age (Ma)"
) {



  # Resolve unit names or indices to indices in res$units



  resolve_units <- function(x, units) {

    if (is.character(x)) {

      idx <- match(x, units$name)

      if (anyNA(idx)) {
        stop("Unknown unit name(s).")
      }

      idx

    } else {

      x

    }

  }

  units <- res$units
  n <- nrow(units)

  src_idx <- resolve_units(source_units, units)
  res_idx <- resolve_units(reservoir_units, units)
  seal_idx <- resolve_units(seal_units, units)



  # Define colours for individual source rocks



  red_shades <- c(
    "firebrick",
    "darkred",
    "indianred",
    "brown",
    "tomato",
    "coral"
  )

  if (length(src_idx) == 1) {


    source_unit_colors <- col


  } else if (length(src_idx) <= length(red_shades)) {


    source_unit_colors <- red_shades[seq_along(src_idx)]


  } else {


    source_unit_colors <-
      colorRampPalette(c("firebrick", "darkred"))(length(src_idx))


  }



  # Construct event intervals



  event_data <- list()

  # Source Rock

  if (length(src_idx) > 0) {


    event_data[["Source Rock"]] <- lapply(src_idx, function(i) {

      data.frame(
        start = units$start[i],
        end = units$end[i],
        unit = i
      )
    })


  } else {


    event_data[["Source Rock"]] <- list()


  }

  # Reservoir Rock

  if (length(res_idx) > 0) {


    event_data[["Reservoir Rock"]] <- lapply(res_idx, function(i) {

      data.frame(
        start = units$start[i],
        end = units$end[i]
      )
    })


  } else {


    event_data[["Reservoir Rock"]] <- list()


  }

  # Seal Rock--

  if (length(seal_idx) > 0) {


    event_data[["Seal Rock"]] <- lapply(seal_idx, function(i) {

      data.frame(
        start = units$start[i],
        end = units$end[i]
      )
    })


  } else {


    event_data[["Seal Rock"]] <- list()


  }

  # Overburden-

  #

  # Overburden begins after the first specified source rock has finished

  # being deposited.

  #

# Every subsequent deposition phase is represented as overburden,
# irrespective of the identity of the deposited unit.
# Hiatuses, uplift, and erosion are not deposition phases and therefore
# remain as gaps in the overburden bar.



  if (length(src_idx) > 0) {


    overburden_onset <- max(units$end[src_idx])

    deposition_phases <- res$phases[
      vapply(
        res$phases,
        function(p) p$type == "deposition",
        logical(1)
      )
    ]

    over_intervals <- list()

    for (p in deposition_phases) {

      # Deposition entirely before overburden onset is ignored.
      if (p$end >= overburden_onset) {
        next
      }

      # Truncate a deposition phase crossing the onset.
      start <- min(p$start, overburden_onset)
      end <- p$end

      if (start > end) {

        over_intervals[[length(over_intervals) + 1]] <-
          data.frame(
            start = start,
            end = end
          )
      }
    }

    # Sort from old to young.
    if (length(over_intervals) > 1) {

      over_intervals <- over_intervals[
        order(
          vapply(
            over_intervals,
            function(x) -x$start,
            numeric(1)
          )
        )
      ]
    }

    # Merge directly adjacent deposition phases only.
    # Because non-depositional phases are absent from over_intervals,
    # hiatuses, uplift, and erosion remain visible as gaps.
    if (length(over_intervals) > 1) {

      merged <- list(over_intervals[[1]])

      for (j in 2:length(over_intervals)) {

        current <- over_intervals[[j]]
        last <- merged[[length(merged)]]

        if (isTRUE(all.equal(last$end, current$start))) {

          merged[[length(merged)]]$end <- current$end

        } else {

          merged[[length(merged) + 1]] <- current
        }
      }

      over_intervals <- merged
    }

    event_data[["Overburden"]] <- over_intervals


  } else {


    event_data[["Overburden"]] <- list()


  }

  # Uplift/Erosion

  erosions <- res$phases[
    vapply(
      res$phases,
      function(p) p$type == "erosion",
      logical(1)
    )
  ]

  if (length(erosions) > 0) {


    event_data[["Uplift/Erosion"]] <- lapply(erosions, function(p) {

      data.frame(
        start = p$start,
        end = p$end
      )
    })


  } else {


    event_data[["Uplift/Erosion"]] <- list()


  }

  # Trap Formation

  if (!is.null(trap_formation)) {


    tf <- as.data.frame(trap_formation)

    event_data[["Trap Formation"]] <- lapply(
      seq_len(nrow(tf)),
      function(i) {

        data.frame(
          start = tf$start[i],
          end = tf$end[i]
        )
      }
    )


  } else {


    event_data[["Trap Formation"]] <- list()


  }

  # Generation/Migration/Accumulation
  # A source rock may generate more than once. All generation intervals
  # belonging to the same source rock are retained and later plotted in
  # one common sub-row.

  gen_intervals <- list()

  if (length(src_idx) > 0) {


    lo <- min(generation_window)
    hi <- max(generation_window)

    for (u in src_idx) {

      prof_u <- res$profile[
        res$profile$unit == u,
        ,
        drop = FALSE
      ]

      prof_u <- prof_u[order(-prof_u$age), ]

      temp <- prof_u[[temperature_column]]
      ages <- prof_u$age

      peak <- cummax(temp)

      in_window <-
        (temp >= peak - 1e-9) &
        (temp > lo) &
        (temp < hi)

      if (!any(in_window)) {
        next
      }

      runs <- rle(in_window)

      end_idx <- cumsum(runs$lengths)

      start_idx <- c(
        1,
        end_idx[-length(end_idx)] + 1
      )

      for (k in which(runs$values)) {

        i_start <- start_idx[k]
        i_end <- end_idx[k]

        gen_intervals[[length(gen_intervals) + 1]] <- data.frame(
          start = ages[i_start],
          end = ages[i_end],
          unit = u
        )
      }
    }


  }

  event_data[["Generation/Migration/Accumulation"]] <-
    gen_intervals

  # Preservation
  # Preservation is only shown when all petroleum-system elements have
  # been supplied.



  all_units_specified <-
    !is.null(source_units) &&
    length(source_units) > 0 &&
    !is.null(reservoir_units) &&
    length(reservoir_units) > 0 &&
    !is.null(seal_units) &&
    length(seal_units) > 0 &&
    !is.null(trap_formation) &&
    nrow(as.data.frame(trap_formation)) > 0

  if (all_units_specified) {


    event_data[["Preservation"]] <-
      data.frame(
        start = critical_moment,
        end = 0
      )


  } else {


    event_data[["Preservation"]] <- list()


  }

  # Plot configuration

  category_labels <- c(
    "Source Rock",
    "Reservoir Rock",
    "Seal Rock",
    "Overburden",
    "Uplift/Erosion",
    "Trap Formation",
    "Generation/Migration/Accumulation",
    "Preservation",
    "Critical Moment"
  )

  n_cat <- length(category_labels)

  if (is.null(xmax)) {
    xmax <- max(res$ages)
  }

  xf <- function(age) xmax - age

  # Ensure all event categories use lists of data frames.

  normalize <- function(x) {


    if (is.data.frame(x)) {
      return(list(x))
    }

    if (is.list(x) && all(sapply(x, is.data.frame))) {
      return(x)
    }

    if (is.list(x) && length(x) == 0) {
      return(list())
    }

    stop(
      "Each category must be a data.frame or list of data.frames."
    )


  }

  event_data <- lapply(event_data, normalize)



  # Plot



  old_mar <- par("mar")

  par(
    mar = c(5, 2, 4, 8)
  )

  on.exit(
    par(mar = old_mar)
  )

  plot(
    NA,
    xlim = c(0, xmax),
    ylim = c(0.5, n_cat + 0.5),
    xaxs = "i",
    yaxs = "i",
    xlab = "",
    ylab = "",
    axes = FALSE,
    main = main)

  # X-axis

  xticks <- seq(
    0,
    floor(xmax / tick) * tick,
    by = tick
  )

  axis(
    1,
    at = xf(xticks),
    labels = xticks,
    lwd = 1
  )

  mtext(
    xlab,
    side = 1,
    line = 2.5,
    cex = 0.9
  )

  # Horizontal category boundaries

  for (i in seq_len(n_cat)) {


    abline(
      h = i + 0.5,
      lwd = 1,
      col = "grey70"
    )


  }

  abline(
    h = 0.5,
    lwd = 1,
    col = "grey70"
  )

  abline(
    h = n_cat + 0.5,
    lwd = 1,
    col = "grey70"
  )

  # Determine whether intervals overlap

  intervals_overlap <- function(df_list) {


    if (length(df_list) < 2) {
      return(FALSE)
    }

    for (i in seq_len(length(df_list) - 1)) {

      for (j in (i + 1):length(df_list)) {

        if (
          df_list[[i]]$start <= df_list[[j]]$end &&
          df_list[[j]]$start <= df_list[[i]]$end
        ) {
          return(TRUE)
        }
      }
    }

    FALSE


  }



  # Draw event categories



  for (i in seq_len(n_cat - 1)) {


    cat_name <- names(event_data)[i]
    intervals <- event_data[[cat_name]]

    n_bars <- length(intervals)

    if (n_bars == 0) {
      next
    }

    y_center <- n_cat - i + 1

    y_bottom <- y_center - 0.5
    y_top <- y_center + 0.5


    # Generation/Migration/Accumulation
    # One sub-row is allocated to each source rock. Multiple generation
    # intervals from the same source rock are drawn in that same sub-row.

    if (cat_name == "Generation/Migration/Accumulation") {

      if (
        !all(
          vapply(
            intervals,
            function(x) "unit" %in% names(x),
            logical(1)
          )
        )
      ) {
        stop(
          "Generation intervals must contain a 'unit' column."
        )
      }

      generation_units <- unique(
        vapply(
          intervals,
          function(x) x$unit[1],
          numeric(1)
        )
      )

      n_source_rows <- length(generation_units)

      if (n_source_rows == 0) {
        next
      }

      sub_height <- 1 / n_source_rows

      for (s in seq_along(generation_units)) {

        unit_id <- generation_units[s]

        unit_intervals <- intervals[
          vapply(
            intervals,
            function(x) identical(x$unit[1], unit_id),
            logical(1)
          )
        ]

        y0 <- y_bottom + (s - 1) * sub_height
        y1 <- y0 + sub_height

        col_idx <- which(src_idx == unit_id)

        if (length(col_idx) == 0) {
          col_use <- col
        } else {
          col_use <- source_unit_colors[col_idx[1]]
        }

        for (b in seq_along(unit_intervals)) {

          rect(
            xf(unit_intervals[[b]]$start),
            y0,
            xf(unit_intervals[[b]]$end),
            y1,
            col = col_use,
            border = NA
          )
        }
      }

      next
    }


    # All other categories

    should_stack <-
      n_bars > 1 &&
      intervals_overlap(intervals)

    if (!should_stack) {

      for (b in seq_len(n_bars)) {

        col_use <- col

        if (
          cat_name == "Source Rock" &&
          "unit" %in% names(intervals[[b]])
        ) {

          unit_id <- intervals[[b]]$unit

          col_idx <- which(src_idx == unit_id)

          if (length(col_idx) > 0) {
            col_use <- source_unit_colors[col_idx[1]]
          }
        }

        rect(
          xf(intervals[[b]]$start),
          y_bottom,
          xf(intervals[[b]]$end),
          y_top,
          col = col_use,
          border = NA
        )
      }

    } else {

      sub_height <- 1 / n_bars

      for (b in seq_len(n_bars)) {

        y0 <- y_bottom + (b - 1) * sub_height
        y1 <- y0 + sub_height

        col_use <- col

        if (
          cat_name == "Source Rock" &&
          "unit" %in% names(intervals[[b]])
        ) {

          unit_id <- intervals[[b]]$unit

          col_idx <- which(src_idx == unit_id)

          if (length(col_idx) > 0) {
            col_use <- source_unit_colors[col_idx[1]]
          }
        }

        rect(
          xf(intervals[[b]]$start),
          y0,
          xf(intervals[[b]]$end),
          y1,
          col = col_use,
          border = NA
        )
      }
    }


  }



  # Critical moment



  x_crit <- xf(critical_moment)

  lines(
    c(x_crit, x_crit),
    c(0.5, 1.5 - 0.1),
    lwd = 2,
    col = critical_col
  )

  arrows(
    x0 = x_crit,
    y0 = 1.5 - 0.15,
    x1 = x_crit,
    y1 = 1.5,
    length = 0.08,
    lwd = 2,
    col = critical_col
  )

  # Right-hand category labels

  y_positions <- seq(
    n_cat,
    1,
    by = -1
  )

  axis(
    4,
    at = y_positions,
    labels = category_labels,
    las = 1,
    tick = FALSE,
    line = -0.5
  )

  box()

  invisible(event_data)
}
