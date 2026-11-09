Figure 1: platform schematic + AFM micrographs (image data, not computed). AFM = Supplementary Fig. S6.
Tile design: see ../Tile_Tile_Master.ipynb.

Fig 1C — Fitted melting curves (computed):
  make_melting_curves.py   regenerates "Melting Curves.svg" (the normalized two-state fits of every
                           construct, 30<=Tm<=50, dx>=1) and writes the raw plotted data to
                           Melting_Curves_data.csv (+ a copy at ../Figure_Source_Data/Fig1C_MeltingCurves.csv).
                           Wired into ../run_all.sh. Fit params come from the authoritative Filtered_DF_Tm;
                           unaffected by the tile-concentration (C_T) correction (melt-fit Tm/dx only).
