#' @title basinmodlR: one-dimensional burial and thermal history modelling
#'
#' @description
#' \pkg{basinmodlR} forward models the burial, compaction and temperature history
#' of a one-dimensional sedimentary column through geological time. The whole
#' basin history is described by two small tables, one row per stratigraphic
#' unit and one row per erosion event, from which everything else is derived.
#'
#' @section Describing a basin:
#' A basin is a list with four elements (see \code{\link{brabant}},
#' \code{\link{condroz}} and \code{\link{ardennes}} for worked examples):
#' \describe{
#'   \item{\code{units}}{one row per stratigraphic unit, oldest first, with the
#'     deposition interval (\code{start}, \code{end}, in Ma, \code{start} being
#'     the older age), the thickness at maximum burial (\code{thickness}, km),
#'     the depositional porosity (\code{phi0}), the porosity decay coefficient
#'     (\code{ck}, 1/km), the matrix thermal conductivity (\code{K}, W/m/K), the
#'     radiogenic heat production (\code{A}, uW/m3) and the grain density
#'     (\code{rho}, kg/m3).}
#'   \item{\code{erosions}}{one row per erosion event, with \code{start},
#'     \code{end} and the amount of section removed (\code{amount}, km). Erosion
#'     is stripped off the top of the column and may cut through several units.}
#'   \item{\code{water_depth}}{water depth (m) as (age, value) knots.}
#'   \item{\code{heat_flow}}{heat flow into the base of the column (mW/m2) as
#'     (age, value) knots.}
#' }
#' Deposition phases follow from the unit intervals, erosion phases from the
#' erosion table, and every remaining gap in the timeline automatically becomes a
#' hiatus. Nothing has to be repeated per phase:
#' layer stacks, layer identifiers, maximum-burial windows and the amount of
#' section previously removed are all worked out by \code{\link{run_basin_model}}.
#'
#' @section How the model works:
#' Porosity follows \eqn{\phi(z) = \phi_0 e^{-c z}}, with \eqn{z} measured from
#' the sediment-water interface: the water column is not part of the burial
#' depth. For each unit the model first determines the deepest position it ever
#' occupies, and from the thickness there it computes the invariant grain (solid)
#' thickness. At every time step the thickness of a layer whose top sits at depth
#' \eqn{z_1} then solves
#' \deqn{h = h_{solid} + (\phi_0/c)(e^{-c z_1} - e^{-c (z_1 + h)}),}
#' which is iterated to convergence. Compaction is irreversible, so the thickness
#' used is the smaller of that solution and the thickness the layer had at its
#' maximum burial; an exhumed layer simply keeps its porosity. During deposition
#' grain material is added at a constant rate, so a unit reaches exactly its
#' maximum-burial thickness at the end of its own deposition phase and hands over
#' to the next phase without a step.
#'
#' Temperatures are the steady-state conductive solution for the layered column.
#' Bulk conductivity of a layer is the porosity-weighted mean of the matrix and
#' pore-water conductivities. Surface heat flow is the basal heat flow plus the
#' radiogenic production of the whole column, and the heat flow entering each
#' layer is reduced by the production of the layers above it, so that within a
#' layer of thickness \eqn{h}
#' \deqn{T_{base} = T_{top} + q h / K - A h^2 / (2 K).}
#' Setting \code{radiogenic = "legacy"} reverts to the older formulation in which
#' every layer receives the basal heat flow plus its own production over the full
#' column thickness; use it to reproduce results calibrated against that scheme.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{run_basin_model}}}{run the model; returns a tidy layer
#'     table plus a temperature raster on a common age-depth grid.}
#'   \item{\code{\link{plot_basin_model}}}{composite figure: geological periods,
#'     event bar, temperature raster with horizon lines and isotherms, colour key
#'     and generation-window panels. Also available as
#'     \code{plot(<basin_model>)}.}
#'   \item{\code{\link{plot_event_chart}}}{figures that states the events ocurring
#'   in the time domain.}
#' }
#'
#' @section Units:
#' Input thicknesses and erosion amounts are in km, porosity coefficients in
#' 1/km, radiogenic heat production in uW/m3 and heat flow in mW/m2. Everything
#' the model returns is in metres, degrees Celsius, W/m/K and Ma.
#'
#' @examples
#' res <- run_basin_model(brabant)
#' res
#' head(res$profile)
#'
#' @docType package
#' @name basinmodlR-package
#' @aliases basinmodlR
#' @keywords internal
"_PACKAGE"
