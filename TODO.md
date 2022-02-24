# bioRgeo

**<span style="color:blue"><font size="4">TO DO list</span></font>**

## 1. Package V1

### 1.1 Améliorations générales

* bien se mettre d'accord sur les args

* faire la vignette

* écrire article blog sur bin()

* homogénéiser la terminologie

* Bien vérifier les packages importés (il y en aura à supprimer...)

* Bien vérifier les potentiels problèmes character, numeric, factor pour les sites/species des data.frame (net_to_mat et mat_to_net)  

* Bien vérifier les potentiels problèmes NA, integer, numeric pour les présences et abondances

* executables (Mac)

* decide name and order column outputs algo community detection (+ description in Details)
  
* plant CBNmed ask permission 

* Bien vérifier utilisation similarité pour réseaux mais dissimilarité pour clustering classique (cas particulier euclidean)


### 1.2 Amélirations spécifiques dans des fonctions

#### mat_to_net
* Se mettre d'accord sur les noms des deux colonnes de df dans mat_to_net appelé "Object" pour l'instant

#### spproject
* Bien vérifier les formules de projection dans "spproject" 

* renommer spproject ?

* remove prodmat ?

* Ajouter transformation data.frame en matrix dans spproject plutôt qu'une erreur

#### clustering_nonhierarchical

* Ajouter une recherche le nombre optimal de clusters

#### clustering_hierarchical

* Ajouter l'algorithme `diana` ?

* Vérifier arguments en entrée standard

* 

#### find_nclust_tree

* décider si MARS utile ou pas. Switcher vers segmented?

* Voir si rounding nécessaire avec la nouvelle version de elbow

* Implémenter la parallélisation de la recherche de clusters 

* Fusionner le hiérarchique et le non hiérarchique dans cette fonction ?

## 2. Package V2

* Faire une fonction d'evaluation des bioregionalisations: a quel point la partition qu'on obtient est significative ou non (exemple: metrique de modularite pour les reseaux, l'ideal serait d'avoir une metrique simple et commune a tous les algos)

## 3. Other

* rendre les vignettes le plus simple et explicite possible

