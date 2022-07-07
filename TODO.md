# bioRgeo

**<span style="color:blue"><font size="4">TO DO list</span></font>**

## 1. Package V1

### 1.1 Améliorations générales

* glossaire

* see also

* bien se mettre d'accord sur les args

* faire la vignette

* écrire article blog sur binaires

* écrire article blog sur format matrix/network

* homogénéiser la terminologie

* Bien vérifier les packages importés (il y en aura à supprimer...)

* Bien vérifier les potentiels problèmes character, numeric, factor pour les sites/species des data.frame (net_to_mat et mat_to_net)  

* Bien vérifier les potentiels problèmes NA, integer, numeric pour les présences et abondances

* plant CBNmed ask permission 

* Bien vérifier utilisation similarité pour réseaux mais dissimilarité pour clustering classique (cas particulier euclidean)

* vérifier les NAs, négatif et 0 pour les réseaux

* créer une fonction générique pour les controles ?

* créer fonction is_bipartite ?

* créer une fonction minimaliste pour faire des cartes


### 1.2 Amélirations spécifiques dans des fonctions


#### clustering_nonhierarchical

* permettre de lancer plusieurs dbscan et optics avec différentes valeurs de paramètres ?

#### clustering_hierarchical

* Ajouter l'algorithme `diana` ?

* Vérifier arguments en entrée standard

* 

#### partition_metrics

* décider si MARS utile ou pas. Switcher vers segmented?

* Voir si rounding nécessaire avec la nouvelle version de elbow

* séparer en deux fonctions, l'une pour calculer les métriques, l'autre pour chercher un nombre optimal

* Implémenter la parallélisation de la recherche de clusters 

* Fusionner le hiérarchique et le non hiérarchique dans cette fonction ?

## 2. Package V2

* Faire une fonction d'evaluation des bioregionalisations: a quel point la partition qu'on obtient est significative ou non (exemple: metrique de modularite pour les reseaux, l'ideal serait d'avoir une metrique simple et commune a tous les algos)

## 3. Other

* rendre les vignettes le plus simple et explicite possible

