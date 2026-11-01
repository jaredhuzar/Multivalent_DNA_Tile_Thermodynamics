# Phenomenological $T_m$ models --- equations

## Shared two-state core

$$T_m(^\circ\mathrm{C}) = \frac{1000\sum_k g_k\,\Delta H_{NN,k}}{\Delta S_{\mathrm{eff}}} - 273.15, \qquad g_k = e^{-\bar d_k/\lambda}$$

Here $\bar d_k$ is the mean distance of bond $k$ to the other bonds, $N$ is the number of bonds (2 per binding site), and every $\Delta S_{\mathrm{eff}}$ also includes the bimolecular term $+R\ln(C_T/2)\approx-39.6$ cal/mol·K.

## Family A --- occupancy-weighted

**A1 --- mean:** $\displaystyle \Delta S_{\mathrm{eff}} = \sum_k g_k\big(\Delta S_{NN,k}+\Delta S_{Bond}\big) + \mathrm{base}$

**A2 --- N-scaled $\Delta S_{Bond}$ (best):** $\displaystyle \Delta S_{\mathrm{eff}} = \sum_k g_k\,\Delta S_{NN,k} + N\,\Delta S_{Bond} + \mathrm{base}$

**A3 --- g on $\Delta H$ only:** $\displaystyle \Delta S_{\mathrm{eff}} = \sum_k \Delta S_{NN,k} + N\,\Delta S_{Bond} + \mathrm{base}$ (the weight $g_k$ is applied only in the $\Delta H$ sum)

**A4 --- nearest-neighbor g:** as A1 with $g_k = e^{-d_{\min,k}/\lambda}$ ($d_{\min,k}$ = distance to the single closest bond)

**A5 --- per-type $\Delta S_{Bond}$:** $\displaystyle \Delta S_{\mathrm{eff}} = \sum_k g_k\big(\Delta S_{NN,k}+\Delta S_{Bond}^{\,\mathrm{type}}\big) + \mathrm{base}$

**A6 --- split $\lambda_H/\lambda_S$:** $\Delta H$ uses $g(\lambda_H)$ and $\Delta S$ uses $g(\lambda_S)$

**A7 --- occupancy + rotational:** A2 with $\displaystyle -\,k_{rot}R\ln\!\Big(1+\sum_k (z_k-\bar z)^2\Big)$ added to $\Delta S_{\mathrm{eff}}$

## Family B --- elastic / mechanical-coupling

**B1 --- g = 1 (no weighting):** $\displaystyle \Delta S_{\mathrm{eff}} = \sum_k \Delta S_{NN,k} + N\,\Delta S_{Bond} + \mathrm{base}$ with $g_k\equiv1$

**B2 --- g = 1 + rotational:** B1 with $\displaystyle -\,k_{rot}R\ln\!\Big(1+\sum_k (z_k-\bar z)^2\Big)$

**B3 --- g-weighted + WLC tethering:** $\displaystyle \Delta S_{\mathrm{eff}} = \sum_k g_k\big(\Delta S_{NN,k}+\Delta S_{teth}^{\,\mathrm{type}}\big) + \mathrm{base}$

**B4 --- WLC tethering only:** $\displaystyle \Delta S_{\mathrm{eff}} = \sum_k \big(\Delta S_{NN,k}+\Delta S_{teth}^{\,\mathrm{type}}\big) + \mathrm{base}$

## Mechanical terms

Rotational (elastic Green's-function):

$$\Delta S_{rot} = -\,k_{rot}\,R\ln\!\Big(1+\sum_k (z_k-\bar z)^2\Big)$$

WLC tethering entropy (derived; $\ell_p=1$ nm, rise $0.63$ nm/nt):

$$\Delta S_{teth}^{\,\mathrm{type}} = R\ln\frac{V_2}{V_1^{2}\,C_0}, \qquad V = \tfrac{4}{3}\pi\langle r^2\rangle^{3/2}$$

The 6nt bond is asymmetric (8 nt on one tile, 7 nt on the other):

$$\Delta S_{teth}^{\,6} = R\ln\frac{V(1.5)}{V(8)\,V(7)\,C_0}$$
