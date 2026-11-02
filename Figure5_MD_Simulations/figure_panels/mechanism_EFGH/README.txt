Fig 5 E-H binder-position summaries (HOOMD), from plots_robust/mechanism/.
Below: which COLUMNS are plotted for each manuscript Fig. 5 panel.

E -- tip-gap distributions  (manuscript Fig. 5e; "Panel C" in the source SVG)
  panelC_histograms_raw.csv
     x (histogram): tip_gap
     series/hue:    chain_bend   (linker stiffness 4, 10, 20, 40)
     fixed:         temperature = 2.0 (kT=2), rod_bend = 80
  panelC_histograms_binned.csv  (pre-binned version of the same)
     x: bin_center     y: density_chain4 / density_chain10 / density_chain20 / density_chain40
  figure: panelC_histograms.svg

F -- scaffold registration spread  (manuscript Fig. 5f)
  registration_summary.csv   (x axis = rod_bend, log; series/hue = temperature 1.25/1.5/1.75)
     Fig. 5f panel:   intra_sigma   +/- intra_sigma_reperr    (anchor-spacing spread sigma)
     companion panels in the SVG:
        intra_mean  +/- intra_mean_reperr    (mean anchor spacing vs rigidity)
        dS_reg_kB   +/- dS_reg_reperr         (registration entropy gain)
  registration_raw.csv       (per-frame swarm / distribution)
     plotted: value  where metric == "intra_spacing"   (x of the spacing distribution), grouped by rod_bend
     (metric == "transverse_misreg" is the confounded transverse component; not plotted in F)
  figure: registration.svg

G -- mean tip gap across scaffold x linker  (manuscript Fig. 5g)
  heatmap_grid.csv
     color: tipgap_mean
     y axis: chain_bend (linker 4,10,20,40)   x axis: rod_bend (scaffold 4,20,80,240)
     fixed:  temperature = 2.0

H -- tip-gap spread across scaffold x linker  (manuscript Fig. 5h)
  heatmap_grid.csv
     color: tipgap_sigma        same x/y axes as G,  temperature = 2.0
  (Same file's anchorspacing_mean / anchorspacing_sigma are the companion anchor-spacing
   heatmaps, plotted at temperature = 1.5. tipgap_n_unbound is the per-cell sample count, not plotted.)
  figure: heatmap_grid.svg  (4-up: tipgap sigma, tipgap mean, anchor sigma, anchor mean)

mechanism.csv -- per-bond enthalpy / unbound-gap summary vs stiffness
  columns: axis, stiffness, dH_bridge_perbond, dH_bend_perbond, unbound_gap_mean, unbound_gap_std
  supports the mechanism discussion (mechanism.svg); not one of panels E-H.
