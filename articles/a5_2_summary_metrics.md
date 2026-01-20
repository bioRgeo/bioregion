# 5.2 Summary metrics

In this vignette, we describe two functions to compute summary metrics:

- metrics calculated for each species and/or site
  [`site_species_metrics()`](https://bioRgeo.github.io/bioregion/reference/site_species_metrics.md)
- metrics calculated for each bioregion
  [`bioregion_metrics()`](https://bioRgeo.github.io/bioregion/reference/bioregion_metrics.md)

## 1. Terminology clarification

The `bioregion` is focused on bioregionalization, i.e. clustering of
geographical areas on the basis of species data. However, there are
several cases where species can also become part of the clustering (for
example, in bipartite network clustering), which poses terminology
issues.

To be conceptually accurate, we have chosen to name species clusters as
‘chorotypes’:

- **Bioregion**: A group of sites with similar species composition,
  identified through clustering analysis. Bioregions are geographic
  units.

- **Chorotype**: A group of species with similar distributions within
  the study area. Chorotypes are biological units. This generally
  corresponds to the concept of “regional chorotype” sensu (Baroni
  Urbani *et al.*, 1978), as clarified by (Fattorini, 2015). Note that
  when clustering on worldwide ranges, the concept becomes “global
  chorotypes” (see (Fattorini, 2015) for further details).

### Possible cases of chorotypes

| Clustering scenario          | Site clusters |       Species clusters        | Conceptual basis                                           |
|:-----------------------------|:-------------:|:-----------------------------:|:-----------------------------------------------------------|
| Site-only clustering         |  Bioregions   |               —               | Sites grouped by compositional similarity                  |
| Bipartite network clustering |  Bioregions   | Chorotypes (same cluster IDs) | Sites and species grouped by shared network structure      |
| Species-only clustering      |       —       |          Chorotypes           | Species grouped by distributional similarity               |
| Post-hoc species assignment  |  Bioregions   |     Chorotypes (derived)      | Species assigned to bioregions based on specificity/IndVal |

#### Bipartite network clustering

In bipartite network clustering, both sites and species are assigned to
the **same clusters** (network modules). A species assigned to cluster 1
belongs to the same bioregion as sites assigned to cluster 1. We use the
term **chorotype** to refer to the set of species assigned to a given
bioregion, but it is important to understand that:

> **In bipartite clustering, bioregion ID = chorotype ID.** They are two
> perspectives on the same network partition: bioregion refers to the
> sites in a cluster, chorotype refers to the species in that same
> cluster.

#### Site-only clustering with post-hoc species assignment

Species can be secondarily assigned to bioregions based on metrics such
as maximum specificity or IndVal. Here, **chorotype** refers to the
group of species most strongly associated with a given bioregion. Unlike
bipartite clustering, this assignment is derived rather than intrinsic
to the clustering algorithm.

#### Species-only clustering

When clustering species directly (e.g., by distributional similarity),
the resulting groups are true **chorotypes** in the regional sense
(Fattorini, 2015): species with similar distributions within the study
area.

## 2. Example data

We use the vegetation dataset included in the `bioregion`.

``` r
data("vegedf")
data("vegemat")

# Calculation of (dis)similarity matrices
vegedissim <- dissimilarity(vegemat, metric = c("Simpson"))
vegesim <- dissimilarity_to_similarity(vegedissim)
```

## 3. Bioregionalization

We use the same three bioregionalization algorithms as in the
[visualization
vignette](https://biorgeo.github.io/bioregion/articles/a5_visualization.html),
i.e., non-hierarchical, hierarchical, and network bioregionalizations.
In addition, we include a network bioregionalization algorithm based on
a bipartite network, which assigns clusters to both sites and species.
We chose three bioregions for the non-hierarchical and hierarchical
bioregionalizations.  

``` r
# Non hierarchical bioregionalization
vege_nhclu <- nhclu_kmeans(vegedissim, 
                           n_clust = 3, 
                           index = "Simpson",
                           seed = 1)
vege_nhclu$cluster_info 
```

    ##     partition_name n_clust
    ## K_3            K_3       3

``` r
# Hierarchical bioregionalization
set.seed(1)
vege_hclu <- hclu_hierarclust(dissimilarity = vegedissim,
                              index = "Simpson",
                              method = "average", 
                              n_clust = 3,
                              optimal_tree_method = "best",
                              verbose = FALSE)
vege_hclu$cluster_info
```

    ##   partition_name n_clust requested_n_clust output_cut_height
    ## 1            K_3       3                 3            0.5625

``` r
# Network bioregionalization
set.seed(1)
vege_netclu <- netclu_walktrap(vegesim,
                               index = "Simpson")
vege_netclu$cluster_info 
```

    ##     partition_name n_clust
    ## K_3            K_3       3

``` r
# Bipartite network bioregionalization
install_binaries(verbose = FALSE)
vege_netclubip <- netclu_infomap(vegedf,
                                 seed = 1, 
                                 bipartite = TRUE)
vege_netclubip$cluster_info
```

    ##     partition_name n_clust
    ## K_8            K_8       8

## 4. Metric components

Before diving into specific metrics, we can understand the core terms
using a simple example. Consider a study area with **4 sites** and **4
species**, where sites have been assigned to **2 bioregions**.

### 4.1 Species-derived metrics

The following diagram shows the site-species matrix where sites are
grouped by bioregion. Marginal sums give us all the core terms needed to
compute metrics:

                              Species
                       sp1   sp2   sp3   sp4      n_b 
                     ┌─────┬─────┬─────┬─────┐
              Site A │  1  │  1  │  ·  │  ·  │
         B1   ───────┼─────┼─────┼─────┼─────┤     2
              Site B │  1  │  1  │  1  │  ·  │
     Bioregion ══════╪═════╪═════╪═════╪═════╪══════
              Site C │  ·  │  1  │  1  │  1  │
         B2   ───────┼─────┼─────┼─────┼─────┤     2
              Site D │  ·  │  ·  │  1  │  1  │
                     └─────┴─────┴─────┴─────┘
                      
         n_sb            sp1   sp2   sp3   sp4     n_b
       (per bioregion) ┌─────┬─────┬─────┬─────┐
                  B1   │  2  │  2  │  1  │  0  │   2
                       ├─────┼─────┼─────┼─────┤
                  B2   │  0  │  1  │  2  │  2  │   2
                       └─────┴─────┴─────┴─────┘
         n_s (total)      2     3     3     2      n = 4
         K_s (# bioreg)   1     2     2     1      K = 2

| Term        | Meaning                                         | Where to find it                 |
|:------------|:------------------------------------------------|:---------------------------------|
| \\n\\       | Total number of sites                           | Bottom-right corner (4)          |
| \\K\\       | Total number of bioregions                      | Bottom-right corner (2)          |
| \\n_b\\     | Sites in bioregion \\b\\                        | Right margin per bioregion row   |
| \\n_s\\     | Sites where species \\s\\ occurs                | Bottom margin per species column |
| \\K_s\\     | Number of bioregions where species \\s\\ occurs | Bottom margin \\n_s\\            |
| \\n\_{sb}\\ | Sites in bioregion \\b\\ with species \\s\\     | The \\n\_{sb}\\ summary table    |

#### Examples of calculations

From the \\n\_{sb}\\ table, all species-per-bioregion metrics follow
directly:

**Specificity** (fraction of species’ occurrences in a bioregion):
\\A\_{sp1,B1} = \frac{n\_{sp1,B1}}{n\_{sp1}} = \frac{2}{2} = 1.00 \quad
\text{(sp1 is exclusive to B1)}\\ \\A\_{sp2,B1} =
\frac{n\_{sp2,B1}}{n\_{sp2}} = \frac{2}{3} = 0.67 \quad \text{(sp2
mostly in B1)}\\

**Fidelity** (fraction of bioregion’s sites with the species):
\\B\_{sp2,B1} = \frac{n\_{sp2,B1}}{n\_{B1}} = \frac{2}{2} = 1.00 \quad
\text{(sp2 in all B1 sites)}\\ \\B\_{sp3,B1} =
\frac{n\_{sp3,B1}}{n\_{B1}} = \frac{1}{2} = 0.50 \quad \text{(sp3 in
half of B1)}\\

**IndVal** (indicator value = Specificity × Fidelity):
\\IndVal\_{sp1,B1} = 1.00 \times 1.00 = 1.00 \quad \text{(perfect
indicator of B1)}\\ \\IndVal\_{sp2,B1} = 0.67 \times 1.00 = 0.67\\

### 4.2 Site-derived metrics

The following diagram shows the same site-species matrix, but now
**species are grouped by cluster** (chorotype). We compute how many
species from each cluster occur in each site:

                                    Chorotypes 
                            ┌─── C1 ───┐ ┌─── C2 ───┐
                              sp1   sp2   sp3   sp4
                            ┌─────┬─────┬─────┬─────┐
                     Site A │  1  │  1  │  ·  │  ·  │  2
                            ├─────┼─────┼─────┼─────┤
       Sites         Site B │  1  │  1  │  1  │  ·  │  3
                            ├─────┼─────┼─────┼─────┤
                     Site C │  ·  │  1  │  1  │  1  │  3
                            ├─────┼─────┼─────┼─────┤
                     Site D │  ·  │  ·  │  1  │  1  │  2
                            └─────┴─────┴─────┴─────┘
                      n_c          2           2         n = 4
                      
           n_gc                C1      C2          n_g
         (per cluster)     ┌───────┬───────┐
                    Site A │   2   │   0   │        2
                           ├───────┼───────┤
                    Site B │   2   │   1   │        3
                           ├───────┼───────┤
                    Site C │   1   │   2   │        3
                           ├───────┼───────┤
                    Site D │   0   │   2   │        2
                           └───────┴───────┘
         n_c                   2       2           n = 4

| Term        | Meaning                                          | Where to find it              |
|:------------|:-------------------------------------------------|:------------------------------|
| \\n\\       | Total number of species                          | Bottom-right corner (4)       |
| \\n_c\\     | Species in cluster \\c\\                         | Bottom margin per cluster     |
| \\n_g\\     | Species present in site \\g\\                    | Right margin per site row     |
| \\n\_{gc}\\ | Species from cluster \\c\\ present in site \\g\\ | The \\n\_{gc}\\ summary table |

**NOTE:** in bipartite clustering, bioregion and chorotypes can be the
**exact same clusters.** Nevertheless, we use different terms here to
avoid confusion in the calculation of metrics.

#### Examples of calculations

**Specificity** of Site A for C1 (fraction of site’s species belonging
to C1): \\A\_{A,C1} = \frac{n\_{A,C1}}{n_A} = \frac{2}{2} = 1.00 \quad
\text{(Site A has only C1 species)}\\

**Specificity** of Site B for C1: \\A\_{B,C1} = \frac{n\_{B,C1}}{n_B} =
\frac{2}{3} = 0.67 \quad \text{(Site B mostly has C1 species)}\\

**Fidelity** of Site A for C1 (fraction of C1 species present in Site
A): \\B\_{A,C1} = \frac{n\_{A,C1}}{n_C1} = \frac{2}{2} = 1.00 \quad
\text{(Site A has all C1 species)}\\

**Fidelity** of Site C for C1: \\B\_{C,C1} = \frac{n\_{C,C1}}{n\_{C1}} =
\frac{1}{2} = 0.50 \quad \text{(Site C has half of C1 species)}\\

## 5. List of site/species metrics included in the package

### Metrics per cluster

#### When clusters are assigned to sites (`cluster_on = "site"` or `cluster_on = "both"`)

| Metric        | Entity  | Cluster type | Based on      | Occ | Ab  | Formula (occurrence)                                       | Interpretation                                                                        |
|:--------------|:--------|:-------------|:--------------|:---:|:---:|:-----------------------------------------------------------|:--------------------------------------------------------------------------------------|
| Specificity   | Species | Bioregion    | Co-occurrence |  ✓  |  ✓  | \\A\_{sb} = \frac{n\_{sb}}{n_s}\\                          | Fraction of species’ occurrences in bioregion                                         |
| NSpecificity  | Species | Bioregion    | Co-occurrence |  ✓  |  ✓  | \\\bar{A}\_{sb} = \frac{n\_{sb}/n_b}{\sum_k n\_{sk}/n_k}\\ | Size-normalized specificity                                                           |
| Fidelity      | Species | Bioregion    | Co-occurrence |  ✓  |  ✓  | \\B\_{sb} = \frac{n\_{sb}}{n_b}\\                          | Fraction of bioregion’s sites with species                                            |
| IndVal        | Species | Bioregion    | Co-occurrence |  ✓  |  ✓  | \\A\_{sb} \times B\_{sb}\\                                 | Indicator value (specificity × fidelity)                                              |
| NIndVal       | Species | Bioregion    | Co-occurrence |  ✓  |  ✓  | \\\bar{A}\_{sb} \times B\_{sb}\\                           | Size-normalized indicator value                                                       |
| Rho           | Species | Bioregion    | Co-occurrence |  ✓  |  ✓  | See section 6.1.1                                          | Standardized contribution index                                                       |
| CoreTerms     | Species | Bioregion    | Co-occurrence |  ✓  |  ✓  | \\n\\, \\n_b\\, \\n_s\\, \\n\_{sb}\\                       | Raw counts for custom calculations                                                    |
|               |         |              |               |     |     |                                                            |                                                                                       |
| Richness      | Site    | —            | Co-occurrence |  ✓  |  —  | \\S_g = n_g\\                                              | Number of species                                                                     |
| Rich_Endemics | Site    | Bioregion    | Co-occurrence |  ✓  |  —  | \\E_g = \sum{K_s}\\                                        | Number of endemic species in the site (i.e., species occurring in only one bioregion) |
| Prop_Endemics | Site    | Bioregion    | Co-occurrence |  ✓  |  —  | \\\bar{PctEnd}\_{g} = \frac{E_g}{S_g}\\                    | Proportion of endemic species in the site                                             |
|               |         |              |               |     |     |                                                            |                                                                                       |
| MeanSim       | Site    | Bioregion    | Similarity    |  —  |  —  | \\\frac{1}{n_b - \delta} \sum\_{g' \neq g} sim\_{gg'}\\    | Mean similarity to bioregion                                                          |
| SdSim         | Site    | Bioregion    | Similarity    |  —  |  —  | See section 6.2.1                                          | SD of similarity to bioregion                                                         |

#### When clusters are assigned to species (`cluster_on = "species"` or `cluster_on = "both"`)

| Metric       | Entity | Cluster type | Based on      | Occ | Ab  | Formula (occurrence)                                       | Interpretation                           |
|:-------------|:-------|:-------------|:--------------|:---:|:---:|:-----------------------------------------------------------|:-----------------------------------------|
| Specificity  | Site   | Chorotype    | Co-occurrence |  ✓  |  ✓  | \\A\_{gc} = \frac{n\_{gc}}{n_g}\\                          | Fraction of site’s species in cluster    |
| NSpecificity | Site   | Chorotype    | Co-occurrence |  ✓  |  ✓  | \\\bar{A}\_{gc} = \frac{n\_{gc}/n_c}{\sum_k n\_{gk}/n_k}\\ | Size-normalized specificity              |
| Fidelity     | Site   | Chorotype    | Co-occurrence |  ✓  |  ✓  | \\B\_{gc} = \frac{n\_{gc}}{n_c}\\                          | Fraction of cluster’s species in site    |
| IndVal       | Site   | Chorotype    | Co-occurrence |  ✓  |  ✓  | \\A\_{gc} \times B\_{gc}\\                                 | Indicator value (specificity × fidelity) |
| NIndVal      | Site   | Chorotype    | Co-occurrence |  ✓  |  ✓  | \\\bar{A}\_{gc} \times B\_{gc}\\                           | Size-normalized indicator value          |
| Rho          | Site   | Chorotype    | Co-occurrence |  ✓  |  ✓  | See section 6.2.2                                          | Standardized contribution index          |
| CoreTerms    | Site   | Chorotype    | Co-occurrence |  ✓  |  ✓  | \\n\\, \\n_c\\, \\n_g\\, \\n\_{gc}\\                       | Raw counts for custom calculations       |

### Metrics in bioregionalization/clustering

These metrics summarize how an entity is distributed across *all*
clusters, rather than in relation to each individual cluster.

#### When `cluster_on = "site"` (or `"both"`)

| Metric     | Entity  | Based on      | Occ | Ab  | Formula                                           | Interpretation                              |
|:-----------|:--------|:--------------|:---:|:---:|:--------------------------------------------------|:--------------------------------------------|
| P          | Species | Co-occurrence |  ✓  |  ✓  | \\1 - \sum_k \left(\frac{n\_{sk}}{n_s}\right)^2\\ | Evenness of species across bioregions (0–1) |
| Silhouette | Site    | Similarity    |  —  |  —  | \\\frac{a_g - b_g}{\max(a_g, b_g)}\\              | Fit to assigned vs. nearest bioregion       |

#### When `cluster_on = "species"` (or `"both"`)

| Metric | Entity | Based on      | Occ | Ab  | Formula                                           | Interpretation                           |
|:-------|:-------|:--------------|:---:|:---:|:--------------------------------------------------|:-----------------------------------------|
| P      | Site   | Co-occurrence |  ✓  |  ✓  | \\1 - \sum_k \left(\frac{n\_{gk}}{n_g}\right)^2\\ | Evenness of site across chorotypes (0–1) |

## 6. Metrics per cluster

### 6.1 Species-per-bioregion metrics

These metrics are computed when sites have clusters (i.e.,
`cluster_on = "site"` (or `"both"`)). In the following example, we
compute all metrics
(`bioregion_metrics = c("Specificity", "NSpecificity", "Fidelity", "IndVal", "NIndVal", "Rho", "CoreTerms")`).
To compute these metrics, we need to provide `comat`.

#### 6.1.1 Co-occurrence metrics: occurrence version

The occurrence metrics are computed when `data_type = "occurrence"`. By
default, the function will detect the type of data used for the
clustering. However, this parameter can be overriden by users, such that
occurrence metrics can be calculated for abundance clustering, and
vice-versa. Users can also specify `data_type = "both"` if they want to
obtain both versions of co-occurrence metrics.

``` r
nsb <- site_species_metrics(bioregionalization = vege_nhclu,
                            bioregion_metrics = c("Specificity", "NSpecificity",
                                                  "Fidelity", "IndVal", "NIndVal",
                                                  "Rho", 
                                                  "CoreTerms"),
                            bioregionalization_metrics = NULL,
                            data_type = "occurrence",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = NULL,
                            index = NULL, # Name of similarity column
                            verbose = FALSE)

nsb$species_bioregions[1:10,]
```

    ##    Species Bioregion n_sb n_s n_b Specificity_occ NSpecificity_occ Fidelity_occ
    ## 1    10001         1   27 254 358      0.10629921       0.05586158  0.075418994
    ## 2    10001         2   97 254 150      0.38188976       0.47897510  0.646666667
    ## 3    10001         3  130 254 207      0.51181102       0.46516332  0.628019324
    ## 4    10002         1    5  19 358      0.26315789       0.16739339  0.013966480
    ## 5    10002         2    1  19 150      0.05263158       0.07990244  0.006666667
    ## 6    10002         3   13  19 207      0.68421053       0.75270417  0.062801932
    ## 7    10003         1   82  86 358      0.95348837       0.92219928  0.229050279
    ## 8    10003         2    0  86 150      0.00000000       0.00000000  0.000000000
    ## 9    10003         3    4  86 207      0.04651163       0.07780072  0.019323671
    ## 10   10004         1  205 438 358      0.46803653       0.31989959  0.572625698
    ##      IndVal_occ NIndVal_occ    Rho_occ
    ## 1  0.0080169797 0.004213024 -15.645265
    ## 2  0.2469553806 0.309737230   8.383617
    ## 3  0.3214272129 0.292131557   9.721763
    ## 4  0.0036753896 0.002337896  -2.097447
    ## 5  0.0003508772 0.000532683  -1.704102
    ## 6  0.0429697432 0.047271276   3.842178
    ## 7  0.2183967780 0.211230003   8.947452
    ## 8  0.0000000000 0.000000000  -5.090898
    ## 9  0.0008987754 0.001503396  -5.293785
    ## 10 0.2680097446 0.183182724  -2.194976

##### Specificity (occurrence)

The specificity \\A\_{sb}\\ of species \\s\\ for bioregion \\b\\ (De
Cáceres & Legendre, 2009) is defined as

\\A\_{sb} = \frac{n\_{sb}}{n_s}\\

and measures the fraction of occurrences of species \\s\\ that belong to
bioregion \\b\\. It therefore reflects the uniqueness of a species to a
particular bioregion.

##### NSpecificity (occurrence)

A normalized version that accounts for the size of each bioregion is
also available, as defined in (De Cáceres & Legendre, 2009):

\\\bar{A}\_{sb} = \frac{n\_{sb}/n_b}{\sum\_{k=1}^K n\_{sk}/n_k}\\

It corresponds to a normalized specificity value that adjusts for
differences in bioregion size.

##### Fidelity (occurrence)

The fidelity \\B\_{sb}\\ of species \\s\\ for bioregion \\b\\ (De
Cáceres & Legendre, 2009) is defined as

\\B\_{sb} = \frac{n\_{sb}}{n_b}\\

and measures the fraction of sites in bioregion \\b\\ where species
\\s\\ is present. It therefore reflects the frequency of occurrence of a
species within a bioregion.

##### IndVal (occurrence)

The indicator value \\IndVal\_{sb}\\ of species \\s\\ for bioregion
\\b\\ can be defined as the product of specificity and fidelity (De
Cáceres & Legendre, 2009):

\\IndVal\_{sb} = A\_{sb} \times B\_{sb}\\

This index quantifies the strength of association between a species and
a bioregion by combining its specificity (uniqueness to that bioregion)
and fidelity (consistency of occurrence within that bioregion). High
IndVal values identify species that are both frequent and restricted to
a single bioregion, making them good indicators of that region.

##### NIndVal (occurrence)

A normalized version of the indicator value is also available:

\\\bar{IndVal}\_{sb} = \bar{A}\_{sb} \times B\_{sb}\\

This normalization adjusts for differences in bioregion size, allowing
more comparable indicator values across regions with unequal sampling
effort or extent.

##### Rho (occurrence)

The contribution index \\\rho\\ can also be calculated following
(Lenormand *et al.*, 2019):

\\\rho\_{sb} = \frac{n\_{sb} - n_s\frac{n_b}{n}}{\sqrt{\frac{n_b(n -
n_b)}{n - 1} \frac{n_s}{n}(1 - \frac{n_s}{n}) }}\\

This index measures the deviation between the observed number of
occurrences of species \\s\\ in bioregion \\b\\ and the expected value
under random association, providing a standardized measure of
contribution to the bioregional structure.

#### Co-occurrence metrics: abundance version

The occurrence metrics are computed when `data_type = "occurrence"`. By
default, the function will detect the type of data used for the
clustering. However, this parameter can be overriden by users, such that
occurrence metrics can be calculated for abundance clustering, and
vice-versa.

The abundance version of these metrics can also be computed when
`data_type = "abundance"` (or `data_type = "both"`). In this case the
core terms and associated metrics are:

- \\w\_{sb}\\ is the sum of abundances of species **s** in sites of
  bioregion **b**.
- \\w_s\\ is the total abundance of species **s**.  
- \\w_b\\ is the total abundance of all species present in sites of
  bioregion **b**.

``` r
wsb <- site_species_metrics(bioregionalization = vege_nhclu,
                            bioregion_metrics = c("Specificity", "NSpecificity",
                                                  "Fidelity",
                                                  "IndVal", "NIndVal",
                                                  "Rho",
                                                  "CoreTerms"),
                            bioregionalization_metrics = NULL,
                            data_type = "abundance",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = NULL, # Name of similarity column
                            index = NULL,
                            verbose = FALSE)

wsb$species_bioregions[1:10,]
```

    ##    Species Bioregion w_sb  w_s     w_b Specificity_abund NSpecificity_abund
    ## 1    10001         1   85 6255 1889243        0.01358913        0.006665761
    ## 2    10001         2 3037 6255 1081424        0.48553157        0.568417434
    ## 3    10001         3 3133 6255 1271474        0.50087930        0.424916804
    ## 4    10002         1   81  295 1889243        0.27457627        0.178777214
    ## 5    10002         2    3  295 1081424        0.01016949        0.015803023
    ## 6    10002         3  211  295 1271474        0.71525424        0.805419763
    ## 7    10003         1  439  444 1889243        0.98873874        0.980682689
    ## 8    10003         2    0  444 1081424        0.00000000        0.000000000
    ## 9    10003         3    5  444 1271474        0.01126126        0.019317311
    ## 10   10004         1 2227 7476 1889243        0.29788657        0.187592363
    ##    Fidelity_abund IndVal_abund NIndVal_abund Rho_abund
    ## 1    4.499157e-05 1.024878e-03  0.0005027250 -6.805816
    ## 2    2.808334e-03 3.139771e-01  0.3675766076  4.731084
    ## 3    2.464069e-03 3.145619e-01  0.2668559641  3.255761
    ## 4    4.287432e-05 3.834864e-03  0.0024968885 -1.446975
    ## 5    2.774120e-06 6.779661e-05  0.0001053535 -1.568655
    ## 6    1.659491e-04 4.491935e-02  0.0505819175  3.003468
    ## 7    2.323682e-04 2.264709e-01  0.2246256438  4.175535
    ## 8    0.000000e+00 0.000000e+00  0.0000000000 -2.204186
    ## 9    3.932444e-06 2.176089e-04  0.0003732814 -2.624519
    ## 10   1.178779e-03 1.705775e-01  0.1074202078 -5.624022

##### Specificity (abundance)

\\A\_{sb} = \frac{w\_{sb}}{w_s}\\

##### NSpecificity (abundance)

\\\bar{A}\_{sb} = \frac{w\_{sb}/n_b}{\sum\_{k=1}^K w\_{sk}/n_k}\\

##### Fidelity (abundance)

\\B\_{sb} = \frac{w\_{sb}}{w_b}\\

##### IndVal (abundance)

\\IndVal\_{sb} = A\_{sb} \times \frac{n\_{sb}}{n_b}\\ Note that the
fidelity based on occurrence is used here (De Cáceres & Legendre, 2009).

##### NIndVal (abundance)

\\\bar{IndVal}\_{sb} = \bar{A}\_{sb} \times \frac{n\_{sb}}{n_b}\\

Note that the fidelity based on occurrence is used here (De Cáceres &
Legendre, 2009).

##### Rho (abundance)

\\\rho\_{sb} = \frac{\mu\_{sb} - \mu_s}{\sqrt{\left(\frac{n -
n_b}{n-1}\right) \left(\frac{{\sigma_s}^2}{n_b}\right)}}\\ where

- \\\mu\_{sb} = \frac{w\_{sb}}{n_b}\\ the average abundance of species
  \\s\\ in bioregion \\b\\ (as in **NSpecificity** and **NIndVal**)
- \\\mu_s = \frac{w_s}{n}\\ the average abundance of species \\s\\
- \\\sigma_s\\ the associated standard deviation.

### 6.2 Site metrics

For sites, two types of metrics can be computed, depending on whether
the clustering is based on site or species:

- if the clustering is based on sites (`cluster_on = "site"` (or
  `"both"`)), then richness and similarity-based metrics can be computed
- if the clustering is based on species (`cluster_on = "species"` (or
  `"both"`)), then we can also compute metrics that are typically
  applied at the species level, such as affinity, fidelity, IndVal and
  other similar metrics. The conceptual interpretation differs in this
  case.

#### 6.2.1 Diversity & endemicity site metrics

When clusters are assigned to sites (bioregions), we can compute basic
diversity metrics:

- Richness = number of species in the site
- Rich_Endemics = number of species in the site that are endemic to a
  single region (i.e., occur in only one bioregion)
- Prop_Endemics = proportion of endemic species, i.e. ratio between
  Rich_Endemics and Richness

``` r
sim_metrics <- site_species_metrics(bioregionalization = vege_nhclu,
                            bioregion_metrics = c("Richness", "Rich_Endemics",
                                                  "Prop_Endemics"),
                            bioregionalization_metrics = NULL,
                            data_type = "occurrence",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = vegesim,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

sim_metrics$site_bioregions[1:10,]
```

    ##    Site Bioregion Richness Rich_Endemics Prop_Endemics
    ## 1    35         1      129             0    0.00000000
    ## 2    35         2      129             2    0.01550388
    ## 3    35         3      129             0    0.00000000
    ## 4    36         1      867             0    0.00000000
    ## 5    36         2      867            24    0.02768166
    ## 6    36         3      867             0    0.00000000
    ## 7    37         1      659             0    0.00000000
    ## 8    37         2      659            28    0.04248862
    ## 9    37         3      659             0    0.00000000
    ## 10   38         1      562             0    0.00000000

#### 6.2.2 Similarity-based site metrics

To compute similarity-based metrics for sites, we need to provide the
site similarity matrix (`vegesim`).

These metrics include the average similarity of each site to the sites
of  
each bioregion (\\MeanSim\\) and the associated standard deviation
(\\SdSim\\).  
When computing the average similarity, the focal site itself is not
included  
in the calculation for its own bioregion.

``` r
sim_metrics <- site_species_metrics(bioregionalization = vege_nhclu,
                            bioregion_metrics = c("MeanSim", "SdSim"),
                            bioregionalization_metrics = NULL,
                            data_type = "occurrence",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = vegesim,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

sim_metrics$site_bioregions[1:10,]
```

    ##    Site Bioregion    MeanSim      SdSim
    ## 1    35         1 0.05505001 0.03334749
    ## 2    35         2 0.55413137 0.34023380
    ## 3    35         3 0.25070767 0.19491881
    ## 4    36         1 0.21852554 0.06086231
    ## 5    36         2 0.63174150 0.19172348
    ## 6    36         3 0.43691643 0.11753321
    ## 7    37         1 0.12094558 0.05025788
    ## 8    37         2 0.57121940 0.23598380
    ## 9    37         3 0.31929608 0.13362203
    ## 10   38         1 0.08133525 0.03614405

##### MeanSim

Let \\g\\ be a site and \\b\\ a bioregion with sites \\g' \in b\\, then:

\\MeanSim\_{gb} = \frac{1}{n_b - \delta\_{g \in b}} \sum\_{g' \in b, g'
\neq g} sim\_{gg'}\\ where \\sim\_{gg'}\\ is the similarity between
sites \\g\\ and \\g'\\, \\n_b\\ is the number of sites in bioregion
\\b\\, and \\\delta\_{g \in b}\\ is 1 if site \\g\\ belongs to bioregion
\\b\\ (to exclude itself), 0 otherwise.

##### SdSim

The standard deviation of similarities of site \\g\\ to bioregion \\b\\
is:

\\SdSim\_{gb} = \sqrt{\frac{1}{n_b - 1 - \delta\_{g \in b}} \sum\_{g'
\in b, g' \neq g} \left( sim\_{gg'} - MeanSim\_{gb} \right)^2}\\ where
\\sim\_{gg'}\\ is the similarity between sites \\g\\ and \\g'\\, \\n_b\\
is the number of sites in bioregion \\b\\, and \\\delta\_{g \in b}\\ is
1 if site \\g\\ belongs to bioregion \\b\\ (to exclude itself), 0
otherwise.

#### 6.2.3 Chorotype/Cluster-based site metrics

In the following example we compute only metrics for sites, on the basis
of species clusters (`cluster_on = "species"`).

``` r
gc <- site_species_metrics(bioregionalization = vege_netclubip,
                            bioregion_metrics = c("Specificity", "NSpecificity",
                                                  "Fidelity",
                                                  "IndVal", "NIndVal",
                                                  "Rho",
                                                  "CoreTerms"),
                            bioregionalization_metrics = "P",
                            data_type = "both",
                            cluster_on = "species",
                            comat = vegemat,
                            similarity = NULL,
                            index = NULL,
                            verbose = FALSE)

gc$site_chorotypes[1:10,]
```

    ##    Site Chorotypes n_gc n_g  n_c Specificity_occ NSpecificity_occ Fidelity_occ
    ## 1    35          1    1 129 2219     0.007751938      0.002901223 0.0004506534
    ## 2    35          2  121 129  873     0.937984496      0.892297163 0.1386025200
    ## 3    35          3    7 129  430     0.054263566      0.104801614 0.0162790698
    ## 4    35          4    0 129  137     0.000000000      0.000000000 0.0000000000
    ## 5    35          5    0 129   16     0.000000000      0.000000000 0.0000000000
    ## 6    35          6    0 129   16     0.000000000      0.000000000 0.0000000000
    ## 7    35          7    0 129    3     0.000000000      0.000000000 0.0000000000
    ## 8    35          8    0 129    3     0.000000000      0.000000000 0.0000000000
    ## 9    36          1  241 867 2219     0.277970012      0.116252995 0.1086074808
    ## 10   36          2  585 867  873     0.674740484      0.717275560 0.6701030928
    ##      IndVal_occ  NIndVal_occ     Rho_occ  w_gc   w_g     w_c Specificity_abund
    ## 1  3.493438e-06 1.307446e-06 -13.9811951     1   423 2873216       0.002364066
    ## 2  1.300070e-01 1.236746e-01  19.1029702   411   423 1149003       0.971631206
    ## 3  8.833604e-04 1.706073e-03  -2.2372244    11   423  103264       0.026004728
    ## 4  0.000000e+00 0.000000e+00  -2.2676901     0   423   85285       0.000000000
    ## 5  0.000000e+00 0.000000e+00  -0.7621238     0   423   16026       0.000000000
    ## 6  0.000000e+00 0.000000e+00  -0.7621238     0   423    9049       0.000000000
    ## 7  0.000000e+00 0.000000e+00  -0.3294281     0   423    3356       0.000000000
    ## 8  0.000000e+00 0.000000e+00  -0.3294281     0   423    2942       0.000000000
    ## 9  3.018962e-02 1.262594e-02 -22.1362445  4139 43753 2873216       0.094599227
    ## 10 4.521457e-01 4.806486e-01  34.7507328 38926 43753 1149003       0.889676136
    ##    NSpecificity_abund Fidelity_abund IndVal_abund NIndVal_abund   Rho_abund
    ## 1        0.0009070715   3.480421e-07 1.065375e-06  4.087749e-07  -7.9728244
    ## 2        0.9476029112   3.577014e-04 1.346705e-01  1.313402e-01  11.3128933
    ## 3        0.0514900173   1.065231e-04 4.233328e-04  8.382096e-04  -1.8400980
    ## 4        0.0000000000   0.000000e+00 0.000000e+00  0.000000e+00  -1.2815040
    ## 5        0.0000000000   0.000000e+00 0.000000e+00  0.000000e+00  -0.4306870
    ## 6        0.0000000000   0.000000e+00 0.000000e+00  0.000000e+00  -0.4306870
    ## 7        0.0000000000   0.000000e+00 0.000000e+00  0.000000e+00  -0.1861645
    ## 8        0.0000000000   0.000000e+00 0.000000e+00  0.000000e+00  -0.1861645
    ## 9        0.0387188108   1.440546e-03 1.027418e-02  4.205153e-03 -14.0235465
    ## 10       0.9255703214   3.387807e-02 5.961747e-01  6.202275e-01  20.9066178

``` r
gc$site_chorological[1:10,]
```

    ##    Site     P_occ    P_abund
    ## 1    35 0.1171805 0.05525097
    ## 2    36 0.4653281 0.19928153
    ## 3    37 0.3498196 0.08246363
    ## 4    38 0.2925811 0.06936005
    ## 5    39 0.3870672 0.18829826
    ## 6    84 0.4950958 0.41403147
    ## 7    85 0.4611738 0.31729045
    ## 8    86 0.4555509 0.17223117
    ## 9    87 0.5147554 0.29310903
    ## 10   88 0.4965267 0.30487674

## 7. Metrics over the entire bioregionalization (i.e., over all clusters)

### 7.1 Site metrics

Based on \\MeanSim\\, it is possible to derive aggregated metrics that
assess  
how well a site fits within its assigned bioregion relative to others.

For now, only the Silhouette index (Rousseeuw, 1987) is proposed.

#### Silhouette

The Silhouette index for a site \\g\\ is defined as:

\\Silhouette_g = \frac{a_g - b_g}{\max(a_g, b_g)}\\

where:

- \\a_g\\ is the average similarity of site \\g\\ to all other sites in
  its own bioregion,  
- \\b_g\\ is the average similarity of site \\g\\ to all sites belonging
  to the nearest bioregion.

This index reflects how strongly a site is associated with its assigned
bioregion relative to the most similar alternative bioregion, ranging
from -1, when the site may be misassigned (i.e., more similar to another
bioregion than its own), to 1, when the site is well matched to its own
bioregion, and around 0 when the site lies near the boundary between
bioregions.

``` r
sil_metrics <- site_species_metrics(bioregionalization = vege_nhclu,
                            bioregion_metrics = NULL,
                            bioregionalization_metrics = "Silhouette",
                            data_type = "occurrence",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = vegesim,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

sil_metrics$site_bioregionalization[1:10,]
```

    ##    Site Silhouette
    ## 1    35 0.54756636
    ## 2    36 0.30839365
    ## 3    37 0.44102725
    ## 4    38 0.54025960
    ## 5    39 0.37117843
    ## 6    84 0.05207862
    ## 7    85 0.26743108
    ## 8    86 0.35331378
    ## 9    87 0.25830421
    ## 10   88 0.28035626

#### Site participation coefficient

We can compute the participation coefficient \\P_s\\ of a species \\s\\
to the bioregionalization as described in (Denelle *et al.*, 2020),
available in both its occurrence and abundance versions.

These metrics measure whether a site has species from a single region or
from multiple regions - useful when investigating transition zones
(Leroy *et al.*, 2019). There are ranging from 0 to 1. Values close to 0
indicate that the site only has species from a single chorotype (i.e.,
not a transition zone), whereas values close to 1 indicate that the site
has species evenely distributed across multiple chorotypes (i.e., likely
a transition zone).

#### P (occurrence)

\\ P_s = 1 - \sum\_{k=1}^K \left(\frac{n\_{sk}}{n_s}\right)^2 \\

``` r
p_occ_site <- site_species_metrics(bioregionalization = vege_netclubip,
                            bioregion_metrics = NULL,
                            bioregionalization_metrics = "P",
                            data_type = "occurrence",
                            cluster_on = "species",
                            comat = vegemat,
                            similarity = NULL,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

p_occ_site$site_chorological[1:10,]
```

    ##    Site     P_occ
    ## 1    35 0.1171805
    ## 2    36 0.4653281
    ## 3    37 0.3498196
    ## 4    38 0.2925811
    ## 5    39 0.3870672
    ## 6    84 0.4950958
    ## 7    85 0.4611738
    ## 8    86 0.4555509
    ## 9    87 0.5147554
    ## 10   88 0.4965267

#### P (abundance)

\\ P_s = 1 - \sum\_{k=1}^K \left(\frac{w\_{sk}}{w_s}\right)^2 \\

``` r
p_ab_site <- site_species_metrics(bioregionalization = vege_netclubip,
                            bioregion_metrics = NULL,
                            bioregionalization_metrics = "P",
                            data_type = "abundance",
                            cluster_on = "species",
                            comat = vegemat,
                            similarity = NULL,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

p_ab_site$site_chorological[1:10,]
```

    ##    Site    P_abund
    ## 1    35 0.05525097
    ## 2    36 0.19928153
    ## 3    37 0.08246363
    ## 4    38 0.06936005
    ## 5    39 0.18829826
    ## 6    84 0.41403147
    ## 7    85 0.31729045
    ## 8    86 0.17223117
    ## 9    87 0.29310903
    ## 10   88 0.30487674

### 7.2 Species metrics

We can compute the participation coefficient \\P_s\\ of a species \\s\\
to the bioregionalization for species as well.

#### P (occurrence)

\\ P_s = 1 - \sum\_{k=1}^K \left(\frac{n\_{sk}}{n_s}\right)^2 \\

``` r
p_occ_sp <- site_species_metrics(bioregionalization = vege_netclubip,
                            bioregion_metrics = NULL,
                            bioregionalization_metrics = "P",
                            data_type = "occurrence",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = NULL,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

p_occ_sp$species_bioregionalization[1:10,]
```

    ##    Species      P_occ
    ## 1    10001 0.59811520
    ## 2    10002 0.18836565
    ## 3    10003 0.06787453
    ## 4    10004 0.29062155
    ## 5    10005 0.16566163
    ## 6    10006 0.29637642
    ## 7    10007 0.44808648
    ## 8    10008 0.51907253
    ## 9    10009 0.00000000
    ## 10   10010 0.07249541

#### P (abundance)

\\ P_s = 1 - \sum\_{k=1}^K \left(\frac{w\_{sk}}{w_s}\right)^2 \\

``` r
p_ab_sp <- site_species_metrics(bioregionalization = vege_netclubip,
                            bioregion_metrics = NULL,
                            bioregionalization_metrics = "P",
                            data_type = "abundance",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = NULL,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

p_ab_sp$species_bioregionalization[1:10,]
```

    ##    Species    P_abund
    ## 1    10001 0.50904954
    ## 2    10002 0.09041080
    ## 3    10003 0.02230947
    ## 4    10004 0.31682118
    ## 5    10005 0.04592151
    ## 6    10006 0.19857965
    ## 7    10007 0.44297818
    ## 8    10008 0.48782839
    ## 9    10009 0.00000000
    ## 10   10010 0.03312698

These metrics measure how evenly a species is distributed among
bioregions. There are ranging from 0 to 1. Values close to 0 indicate
that the species is largely restricted to a single bioregion, while
values close to 1 indicate that the species is evenly distributed across
multiple bioregions.

Calculations on both occurrence & abundance at the same time:

``` r
ps <- site_species_metrics(bioregionalization = vege_nhclu,
                           bioregion_metrics = NULL,
                            bioregionalization_metrics = "P",
                            data_type = "both",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = NULL,
                            index = NULL,
                            verbose = FALSE)

ps$species_bioregionalization[1:10,]
```

    ##    Species      P_occ    P_abund
    ## 1    10001 0.58091016 0.51319436
    ## 2    10002 0.45983380 0.41291583
    ## 3    10003 0.08869659 0.02226889
    ## 4    10004 0.59334668 0.55362167
    ## 5    10005 0.53267897 0.50874080
    ## 6    10006 0.59944565 0.48784442
    ## 7    10007 0.60353798 0.59088266
    ## 8    10008 0.61624257 0.53496844
    ## 9    10009 0.44444444 0.44444444
    ## 10   10010 0.18447026 0.11565410

## 8. Bioregion metrics & spatial coherence

At the granularity of bioregions, we can calculate the number of sites
it contains and the number of species present in those sites. The number
and proportion of endemic species are also computed. Endemic species are
defined as those occurring only in sites assigned to a particular
bioregion (i.e., species that occur in only one bioregion).

``` r
bioregion_summary <- bioregion_metrics(bioregionalization = vege_nhclu,
                                       comat = vegemat)
bioregion_summary
```

    ##   Bioregion Site_number Species_number Endemics Percentage_Endemic
    ## 1         2         150           2688      133           4.947917
    ## 2         3         207           3090       58           1.877023
    ## 3         1         358           2821      407          14.427508

We use the metric of spatial coherence as in (Divíšek *et al.*, 2016),
except that we replace the number of pixels per bioregion with the area
of each coherent part.

The spatial coherence is expressed in percentage, and has the following
formula:

\\SC_j = 100 \times \frac{LargestPatch_j}{Area_j}\\

where \\j\\ is a bioregion.

Here is an example with the vegetation dataset.

``` r
# Spatial coherence
vegedissim <- dissimilarity(vegemat)
hclu <- nhclu_kmeans(dissimilarity = vegedissim, n_clust = 4)
vegemap <- map_bioregions(hclu, vegesf, write_clusters = TRUE, plot = FALSE)

bioregion_metrics(bioregionalization = hclu, comat = vegemat, map = vegemap,
col_bioregion = 2) 
```

    ##   Bioregion Site_number Species_number Endemics Percentage_Endemic Coherence
    ## 1         2         128           2527       90           3.561535  49.21875
    ## 2         1         169           2983       45           1.508548  56.21302
    ## 3         4         298           2936       56           1.907357  98.99329
    ## 4         3         120           2262       67           2.961981  79.16667

The bioregion 4 is almost constituted of one homogeneous block, which is
why the spatial coherence is very close to 100 %.

``` r
ggplot(vegemap) +
  geom_sf(aes(fill = as.factor(K_4))) +
  scale_fill_viridis_d("Bioregion") +
  theme_bw() +
  theme(legend.position = "bottom")
```

![](a5_2_summary_metrics_files/figure-html/unnamed-chunk-16-1.png)

## 9. References

Baroni Urbani C, Ruffo S & Vigna Taglianti A (1978) Materiali per una
biogeografia italiana fondata su alcuni generi di coleotteri
cicindelidi, carabidi e crisomelidi. *Memorie della Società Entomologica
Italiana* 56, 35–92.

De Cáceres M & Legendre P (2009) Associations between species and groups
of sites: Indices and statistical inference. *Ecology* 90, 3566–3574.

Denelle P, Violle C & Munoz F (2020) Generalist plants are more
competitive and more functionally similar to each other than specialist
plants: Insights from network analyses. *Journal of Biogeography* 47,
1922–1933.

Divíšek J, Storch D, Zelený D & Culek M (2016) Towards the spatial
coherence of biogeographical regionalizations at subcontinental and
landscape scales. *Journal of biogeography* 43, 2489–2501.

Fattorini S (2015) On the concept of chorotype. *Journal of
Biogeography* 42, 2246–2251.

Lenormand M, Papuga G, Argagnon O, Soubeyrand M, Alleaume S & Luque S
(2019) Biogeographical network analysis of plant species distribution in
the mediterranean region. *Ecology and evolution* 9, 237–250.

Leroy B, Dias MS, Giraud E *et al.* (2019) Global biogeographical
regions of freshwater fish species. *Journal of Biogeography* 2407–2419.

Rousseeuw PJ (1987) Silhouettes: A graphical aid to the interpretation
and validation of cluster analysis. *Journal of Computational and
Applied Mathematics* 20, 53–65.
