# bioRgeo

**<span style="color:blue"><font size="4">TO DO list</span></font>**

## 1. Package V1

* gérer l'intégration continue avec "GitHub Actions" au lieu de travis (switch to gitlab?)

* Bien vérifier les packages importés (il y en aura à supprimer...)

* Gérer exécutables + permission exécution pour le CRAN et autre, binman (https://cran.r-project.org/web/packages/binman/vignettes/binman-Basics.html) ?

* Bien vérifier les potentiels problèmes character, numeric, factor pour les sites/species des data.frame (df_to_contingency et contingency_to_df)  

* Bien vérifier les potentiels problèmes NA, integer, numeric pour les présences et abondances

* Se mettre d'accord sur les noms des deux colonnes de df dans contingency_to_df appelé "Object" pour l'instant

* Bien vérifier les formules de projection dans "spproject" 

* renommer spproject

* Bien vérifier utilisation similarité pour réseaux mais dissimilarité pour clustering classique (cas particulier euclidean)

* executables (Mac)

* decide name and order column outputs algo community detection (+ description in Details)
  
* plant CBNmed ask permission 

## 2. Package V2

* Faire une fonction d'evaluation des bioregionalisations: a quel point la partition qu'on obtient est significative ou non (exemple: metrique de modularite pour les reseaux, l'ideal serait d'avoir une metrique simple et commune a tous les algos)

## 3. Other
* rendre les vignettes le plus simple et explicite possible

