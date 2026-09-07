# basinmodlR 0.1.0

* First release.
* `run_basin_model()` replaces the phase-by-phase scripts: a basin is described
  by one row per stratigraphic unit and one row per erosion event, and the
  phase list, hiatuses, layer stacks, maximum-burial windows, grain thicknesses
  and previously eroded thicknesses are all derived from that.
* Burial depths are measured from the sediment-water interface, so a changing
  water depth no longer compacts or decompacts the column.
* The decompaction equation is solved by iteration instead of on a 1 m search
  grid.
* Compaction is irreversible through a single rule, so the ad-hoc
  `previous_erosion` branches are gone and erosion may cut through several
  units.
* Radiogenic heat is carried down the column properly; `radiogenic = "legacy"`
  reproduces the earlier formulation.
* `plot_basin_model()` draws the raster, the horizon lines and the sea floor on
  one shared age grid, so they coincide.
* Input tables for the Brabant Massif, Huy (Condroz inlier) and the Ardennes
  are included as `brabant`, `condroz` and `ardennes`.
