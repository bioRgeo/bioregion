## Formulas Reference

### Rho (ρ)

$$\rho_{ij} = \frac{n_{ij} - \frac{n_i n_j}{n}}{\sqrt{\frac{n - n_j}{n-1} \cdot \left(1-\frac{n_j}{n}\right) \cdot \frac{n_i n_j}{n}}}$$

Where:
- $n$ = total number of sites
- $n_i$ = number of sites where species $i$ is present
- $n_j$ = number of sites in bioregion $j$
- $n_{ij}$ = occurrences of species $i$ in bioregion $j$

### Affinity (A)

$$A_i = \frac{R_i}{Z}$$

Where:
- $R_i$ = occurrence/range size of species $i$ in its bioregion
- $Z$ = total number of sites in the bioregion

### Fidelity (F)

$$F_i = \frac{R_i}{D_i}$$

Where:
- $R_i$ = occurrence/range size of species $i$ in its bioregion
- $D_i$ = total range size of species $i$ (across all bioregions)

### Indicator Value (IndVal)

$$IndVal_i = F_i \cdot A_i$$

### Participation Coefficient (C)

$$C_i = 1 - \sum_{s=1}^{N_M}{\left(\frac{k_{is}}{k_i}\right)^2}$$

Where:
- $k_{is}$ = number of links of node $i$ to nodes in bioregion $s$
- $k_i$ = total degree of node $i$

### Within-bioregion Degree Z-score (z)

$$z_i = \frac{k_i - \overline{k_{si}}}{\sigma_{k_{si}}}$$

Where:
- $k_i$ = links of node $i$ to nodes in its bioregion
- $\overline{k_{si}}$ = average degree in bioregion
- $\sigma_{k_{si}}$ = standard deviation of degrees in bioregion
