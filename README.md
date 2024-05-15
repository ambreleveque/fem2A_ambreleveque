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

#### Etudiant 

Ambre Lévêque : ambre.leveque4@etu.univ-lorraine.fr

#### Date de rendu

15/05/2024

#### Problèmes résolus 

Probléme de Dirichlet pur et problème de Dirichlet avec terme source

#### Comment utiliser le projet ?

Étant donné que le projet a été réalisée sous Linux, les explications suivantes sont valables pour le terminal de commandes de Linux.
Dans le fichier ./Documents/fem2A_ambreleveque/main.cpp, vous avez le choix de compiler et d'afficher ou non les huits tests mis en place en mettant les booléens true (pour compiler et afficher) ou false (pour ne pas compiler ni afficher) dans la fonction run_tests( ) du fichier main.cpp.
Les neufs tests mis en place sont les suivants :
  * Le test t_lmesh permet de vérifier que les maillages sont bien téléchargés. Il affiche les coordonnées des vertex et des attributs du maillage (vertice, segment ou triangle)
  * Le test t_io permet de vérifier que le maillage est bien sauvegardé. vous pouvez modifier le maillage que vous souhaitez tester.
  * Le test t_quadrature permet de vérifier que la calcul de la quadrature est correct. Il affiche le nombre de points de la quadrature et la somme des poids. Deux tests sont réalisés, un avec un ordre égal à 0 et un autre avec un ordre égal à 2.
  * Le test t_mapping permet de vérifier que les méthodes de la classe ElementMapping sont focntionnelles et correctes
  * Le test t_ShapeFunction permet de vérifier que les méthodes de la classe ShapeFunctions sont fonctionnelles et correctes.
  * Le test t_Ke permet de calculer le Ke (local) et le K (global) en fonction de l'ordre de la quadrature 0 ou 2.
  * Le test t_src_term permet de calculer le Fe (local) et le F (global) en fonction de l'ordre de la quadrature 0 ou 2.
  
De la même façon que les tests, le fichier main.cpp vous permet aussi de choisir quelle simulation vous souhaitez lancer parmi les 2 réalisées. À la suite de la fonction run_tests( ), vous trouverez la fonction run_simu( ) qui permet, en changeant les booléens devant les différentes simulations, de compiler ou pas celles-ci.
Deux chois de simulations sont possibles. En rentrant true devant simu_pure_dirichlet, la simulation du problème de Dirichlet pure sera lancée, pour ne plus la visualiser rentrez false. En rentrant true devant simu_dirichlet_source_term, la simulation du problème de Dirichlet avec terme source constant égal à 1 sera lancée, entrez false pour l'annuler.
Sous Linux, la compilation du projet se fait par la commande "make" dans le terminal.
Pour lancer l'exécutable, plusieurs options sont possibles. Si vous souhaitez lancer les tests, entrez la commande suivante : ./build/fem2a -t  . 
Pour lancer les simulations utilisez la commande suivante : ./build/fem2a -s   .
L'affichage graphique des résultats se fait grâce à Medit. Une fois que vous avez compilé et exécuté, la commande : Medit/build/medit data/output/simu.mesh depuis le dossier fem2A_ambreleveque vous permet d'afficher graphiquement la solution souhaité (dirichlet_with_source_term ou bien pure_dirichlet). Une fenêtre s'ouvre, appuyer sur la touche m de votre clavier pour visualiser le maillage à l'aide du fichier .bb portant le même nom de la simulation.
