# Éléments Finis en C++

#### Objectif

Écrire [un code aux éléments finis en C++](course/code.md) pour résoudre des 
[problèmes de Poisson](course/poisson.md) 2D avec des éléments triangulaires
linéaires. 

#### Contact

Paul Cupillard : paul.cupillard@univ-lorraine.fr

#### Liens utiles

Cours sur [Arche](http://arche.univ-lorraine.fr/course/view.php?id=61482)

Vidéo de Gilbert Strang : [Finite element method](https://www.youtube.com/watch?v=WwgrAH-IMOk)

Cours de Grégoire Allaire : [Approximation numérique et optimisation](http://www.cmap.polytechnique.fr/~allaire/map411/polycopie-map411.pdf)

Générateurs de maillages triangulaires : [Gmsh](http://gmsh.info/),[Triangle](https://www.cs.cmu.edu/~quake/triangle.html)

# Utilisation du code fem2A_NOGUES
### Utilisation des commades du code
Liste des commandes :
-t (pour les tests)
-s (pour les simulations)
-h (pour obtenir de l'aide)

### Tests
Pour tester les fonctions du code :
- choisir les fonctions à tester dans main.cpp/run_tests()
- définir les booléens test (t_nomfonction) sur true pour lancer le test, false sinon

/!\ Attention
\Pour les tests de certaines fonctions, il est nécessaire de décommenter des lignes directement dans le fichier src\fem.cpp. C'est le cas pour :
- test du constructeur de ElementMapping
- test de la fonction local_to_global_matrix
- test de la fonction local_to_global_vector


### Simulations
Pour lancer les simulations du code :
- choisir les simulations à lancer dans main.cpp/run_simu()
- définir les booléens test (s_nom_simulation) sur true pour lancer la simulation, false sinon
Le choix du maillage se fait dans le main.cpp/run_simu() dans la boucle if associé au test.

Nos résultats de simulations sont stockées dans le dossier data/output.

Pour afficher les simulations, il est nécessaire d'utiliser Medit.
