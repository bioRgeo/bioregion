# Package index

## Utils

- [`install_binaries()`](https://bioRgeo.github.io/bioregion/reference/install_binaries.md)
  : Download, unzip, check permissions, and test the bioregion's binary
  files

- [`mat_to_net()`](https://bioRgeo.github.io/bioregion/reference/mat_to_net.md)
  : Create a data.frame from a contingency table

- [`net_to_mat()`](https://bioRgeo.github.io/bioregion/reference/net_to_mat.md)
  : Create a contingency table from a data.frame

- [`site_species_subset()`](https://bioRgeo.github.io/bioregion/reference/site_species_subset.md)
  :

  Extract a subset of sites or species from a `bioregion.clusters`
  object

## Pairwise similarity and dissimilarity metrics

- [`similarity()`](https://bioRgeo.github.io/bioregion/reference/similarity.md)
  : Compute similarity metrics between sites based on species
  composition
- [`similarity_to_dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/similarity_to_dissimilarity.md)
  : Convert similarity metrics to dissimilarity metrics
- [`dissimilarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity.md)
  : Compute dissimilarity metrics (beta-diversity) between sites based
  on species composition
- [`dissimilarity_to_similarity()`](https://bioRgeo.github.io/bioregion/reference/dissimilarity_to_similarity.md)
  : Convert dissimilarity metrics to similarity metrics
- [`as_bioregion_pairwise()`](https://bioRgeo.github.io/bioregion/reference/as_bioregion_pairwise.md)
  : Convert a matrix or list of matrices to a bioregion (dis)similarity
  object
- [`bind_pairwise()`](https://bioRgeo.github.io/bioregion/reference/bind_pairwise.md)
  : Combine and enrich bioregion (dis)similarity object(s)
- [`betapart_to_bioregion()`](https://bioRgeo.github.io/bioregion/reference/betapart_to_bioregion.md)
  : Convert betapart dissimilarity to bioregion dissimilarity
  (DEPRECATED)

## Clustering

### Hierarchical clustering

- [`hclu_hierarclust()`](https://bioRgeo.github.io/bioregion/reference/hclu_hierarclust.md)
  : Hierarchical clustering based on dissimilarity or beta-diversity
- [`cut_tree()`](https://bioRgeo.github.io/bioregion/reference/cut_tree.md)
  : Cut a hierarchical tree
- [`hclu_diana()`](https://bioRgeo.github.io/bioregion/reference/hclu_diana.md)
  : Divisive hierarchical clustering based on dissimilarity or
  beta-diversity
- [`hclu_optics()`](https://bioRgeo.github.io/bioregion/reference/hclu_optics.md)
  : OPTICS hierarchical clustering algorithm

### Non-hierarchical clustering

- [`nhclu_clara()`](https://bioRgeo.github.io/bioregion/reference/nhclu_clara.md)
  : Non-hierarchical clustering: CLARA
- [`nhclu_clarans()`](https://bioRgeo.github.io/bioregion/reference/nhclu_clarans.md)
  : Non-hierarchical clustering: CLARANS
- [`nhclu_dbscan()`](https://bioRgeo.github.io/bioregion/reference/nhclu_dbscan.md)
  : Non-hierarchical clustering: DBSCAN
- [`nhclu_kmeans()`](https://bioRgeo.github.io/bioregion/reference/nhclu_kmeans.md)
  : Non-hierarchical clustering: K-means analysis
- [`nhclu_pam()`](https://bioRgeo.github.io/bioregion/reference/nhclu_pam.md)
  : Non-hierarchical clustering: Partitioning Around Medoids
- [`nhclu_affprop()`](https://bioRgeo.github.io/bioregion/reference/nhclu_affprop.md)
  : Non-hierarchical clustering: Affinity Propagation

### Network clustering

- [`netclu_beckett()`](https://bioRgeo.github.io/bioregion/reference/netclu_beckett.md)
  : Community structure detection in weighted bipartite networks via
  modularity optimization
- [`netclu_infomap()`](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.md)
  : Infomap community finding
- [`netclu_greedy()`](https://bioRgeo.github.io/bioregion/reference/netclu_greedy.md)
  : Community structure detection via greedy optimization of modularity
- [`netclu_labelprop()`](https://bioRgeo.github.io/bioregion/reference/netclu_labelprop.md)
  : Finding communities based on propagating labels
- [`netclu_leiden()`](https://bioRgeo.github.io/bioregion/reference/netclu_leiden.md)
  : Finding communities using the Leiden algorithm
- [`netclu_leadingeigen()`](https://bioRgeo.github.io/bioregion/reference/netclu_leadingeigen.md)
  : Finding communities based on the leading eigenvector of the
  community matrix
- [`netclu_louvain()`](https://bioRgeo.github.io/bioregion/reference/netclu_louvain.md)
  : Louvain community finding
- [`netclu_oslom()`](https://bioRgeo.github.io/bioregion/reference/netclu_oslom.md)
  : OSLOM community finding
- [`netclu_walktrap()`](https://bioRgeo.github.io/bioregion/reference/netclu_walktrap.md)
  : Community structure detection via short random walks

## Bioregionalization analysis

- [`site_species_metrics()`](https://bioRgeo.github.io/bioregion/reference/site_species_metrics.md)
  : Calculate metrics for sites and species relative to bioregions and
  chorotypes
- [`bioregion_metrics()`](https://bioRgeo.github.io/bioregion/reference/bioregion_metrics.md)
  : Calculate contribution metrics for bioregions
- [`bioregionalization_metrics()`](https://bioRgeo.github.io/bioregion/reference/bioregionalization_metrics.md)
  : Calculate metrics for one or several bioregionalizations
- [`find_optimal_n()`](https://bioRgeo.github.io/bioregion/reference/find_optimal_n.md)
  : Search for an optimal number of clusters in a list of
  bioregionalizations
- [`compare_bioregionalizations()`](https://bioRgeo.github.io/bioregion/reference/compare_bioregionalizations.md)
  : Compare cluster memberships among multiple bioregionalizations

## Visualization

- [`map_bioregions()`](https://bioRgeo.github.io/bioregion/reference/map_bioregions.md)
  : Create a map of bioregions
- [`bioregion_colors()`](https://bioRgeo.github.io/bioregion/reference/bioregion_colors.md)
  : Add color palettes to bioregion cluster objects
- [`exportGDF()`](https://bioRgeo.github.io/bioregion/reference/exportGDF.md)
  : Export a network to GDF format for Gephi visualization

## Data

- [`fishdf`](https://bioRgeo.github.io/bioregion/reference/fishdf.md) :
  Spatial distribution of fish in Europe (data.frame)
- [`fishmat`](https://bioRgeo.github.io/bioregion/reference/fishmat.md)
  : Spatial distribution of fish in Europe (co-occurrence matrix)
- [`fishsf`](https://bioRgeo.github.io/bioregion/reference/fishsf.md) :
  Spatial distribution of fish in Europe
- [`vegedf`](https://bioRgeo.github.io/bioregion/reference/vegedf.md) :
  Spatial distribution of Mediterranean vegetation (data.frame)
- [`vegemat`](https://bioRgeo.github.io/bioregion/reference/vegemat.md)
  : Spatial distribution of Mediterranean vegetation (co-occurrence
  matrix)
- [`vegesf`](https://bioRgeo.github.io/bioregion/reference/vegesf.md) :
  Spatial distribution of Mediterranean vegetation (spatial grid)
