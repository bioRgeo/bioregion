# bioRgeo

**<span style="color:blue"><font size="4">TO DO list</span></font>**

## 1. Package V1

### 1.1 Améliorations générales

* A harmoniser
    -> stop( , call. = FALSE)
    -> \code{} pour le code dans les docs et `...` dans les vignettes
    -> mettre un point "." à la fin des phrases dans la documentation (?)
    
* glossaire

* faire la vignette

* Bien vérifier les packages importés (il y en aura à supprimer...)

* Bien vérifier utilisation similarité pour réseaux mais dissimilarité pour clustering classique (cas particulier euclidean)

* créer une fonction minimaliste pour faire des cartes


### 1.2 Amélirations spécifiques dans des fonctions

#### optics

* permettre de lancer plusieurs optics avec différentes valeurs de paramètres ?

#### clustering_hierarchical

* Ajouter l'algorithme `diana` ?

#### partition_metrics

* décider si MARS utile ou pas. Switcher vers segmented?

* Voir si rounding nécessaire avec la nouvelle version de elbow

* Implémenter la parallélisation de la recherche de clusters 

* Ajouter d'autres métriques, basées par exemple sur l'abondance

### 1.3 Amélioration des classes

### 1.4 Fonctions génériques

* Envisager une fonction générique plot qui s'adapte au type de cluster pour la classe bioRgeo.clusters

## 2. Package V2

## 3. Other


