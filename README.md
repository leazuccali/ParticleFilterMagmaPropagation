README at [Voir le PDF](https://github.com/leazuccali/ParticleFilterMagmaPropagation/blob/PF_V4/README.pdf)
or download with https://raw.githubusercontent.com/leazuccali/ParticleFilterMagmaPropagation/PF_V4/README.pdf

Structure of one simulation

BEFORE SIMULATION

**0. Definition of simulation parameters.** global_data, grid_data, p_reference_data, resampling_data are 5 dictionnaries that defines the simulation. Each dictionnary is describes in the pdf.

**1. Generation of various elements**

* (R, $\gamma$) grids

-> To generate independently with creation_grilles_RG.py

* truth particle : creation_fichier_particule_reference.py
* a set of n particles : creation_fichier_particule_joblib.py

-> Can be generate with the independently of directly by running main.py

**2. Simulation of propagation**

By running main.py 

Outputs are stored step by step (step_X) in a general folder output_test. 

**3. Results of simulation**

Particle evolution (in output_test) is illustrated with plot_particles_positions.py, genere_animation_gif.py and plot_distribution_values.py.





# Fonctionnement du code et description des fonctions principales 

--------------------------------------------------------------------------------------------------
## Une fonction de création des trajectoires pour l'initialisaton:
-------------------------------------------
**Dossier trajectoire_particules.py** 

Sous-fonctions :  
*calcul_parametres  
fonction_watanabe  
calcul_sig1_sig3  
streamplot_trajectoire*  


Fonction principale   
*trajectoire_une_particule(start_point, R, G, mu, E, vol, temps_courant, pas_temps, xmin, xmax, zmin, zmax, pas_trajectoire, P_load, rayon_load, pas_vect, norme)*

Entrées :  
start_point : point de départ de la particule  
R, G, mu, E, vol : paramètres de la particules   
temps_courant : temps début trajectoire
pas_temps : pas de temps   
xmin, xmax : étendue axe abscisses (m)  
zmin, zmax : étendue axe ordonnées (m)  
pas_trajectoire : discrétisation de l'espace (m)   
P_load, rayon_load : données charge ou décharge  
pas_vect : discrétisation de représentation des vecteurs sigma1 (si besoin)  
norme : passage du mètres à une autre unité (si besoin). Pas défaut norme = 1  

Sorties :  
vec_Xt_eff : coordonnées X trajectoire complète  
vec_Zt : coordonnée Z trajectoire complète  
vec_temps_eff : temps depuis le début de l'expérience quand (X,Z) est atteint  
longueur : longueur maximale du dike  
ouverture : ouverture du dike  
vitesse : vitesse de propagation  


### Initialisation des n particules et de la particule vérité :

**Dossier creation-fichier-particule-reference.py**  
*creation_trajectoire_reference():  
Appel à trajectoire_une_particule*

Sortie : sauvegarde_trajectoire_ref.json dans dossier data_initialisation


**Dossier creation_fichier_particule_joblib.py**  
Sous-fonctions :   
*creation_particules(Nb_particules):  
traiter_particule(p, sauvegarde_param_physiques, temps_courant, pas_temps, xmin, xmax, zmin, zmax, pas_trajectoire, P_load, rayon_load, pas_vect, norme):  
Appel à trajectoire_une_particule*  


Fonction principale :   
*creation_trajectoires(Nb_particules):*

Sortie : dossier particule_n 
- fichier param_init_particule_0.json  
- fichier trajectoire_init_particule_0.json  
**dans dossier data_initialisation**
  
--------------------------------------------------------------------------------------------------


--------------------------------------------------------------------------------------------------
## Gestion de l'expérience d'assimilation
Dans **main.py**

Sous-fonctions de gestion des dossiers :  
*copier_fichiers_particule  
propager_particule  
copy_files_for_particle  
propag_deplac_particule  
traiter_deplacement  
donnes_particule_avant_selec  
donnees_particule_apres_selec  
propag_deplac_particule_apres*






## FONCTION PRINCIPALE DE L'EXPÉRIENCE

*fonction_principale(nb_particules, pas_temps, nb_propag_libre, nb_pas_max , int_assimilation, type_selection, type_observations, pas, nobs, dossier_output)*  
Entrées : 

nb_particules : nombre de particules (ex : 100)  
pas_temps : pas de temps (ex : 60)  
nb_propag_libre : Nombre de pas_temps en propagation libre initiale (pas_temps * nb_propag_libre = temps en secondes) (ex : 60)  
nb_pas_max : Nombre maximal de pas_temps de l'expérience (ex : 550)  
int_assimilation : Nombre de pas_temps entre 2 assimilations (ex : 5)  
type_selection : 'systematic','stratified','multinomial,'residual'  
type_observations : 'regulier','normal','random'  
Si regulier : Pas de discrétiation de l'espace (tous les pas mètres)  
pas : pas de discrétisation SI REGULIER  
nobs : nombre d'observations considérées SI 'normal ou 'random    
#Dossier de sortie :  
dossier_output : nom du dossier de sortie (ex: "output_sauvegarde_plots")  

Sorties :  
**Dossier : output_data_and_figures**  
Un sous-dossier par fenêtre d'assimilation  
Un sous-sous-dossier par particules :  
Fichiers   
param_particules_n_step_k_avant.jon  
param_particules_n_step_k_apres.json  
trajectoire_particules_n_step_k_avant.json  
trajectoire_particules_n_step_k_apres.json  

**Dossier : output_sauvegarde_plots**  
Un sous dossier par fenêtre d'assimilation  
Fichiers :  
sauvegarde_parametres_physiques_apres_k.json  
sauvegarde_parametres_physiques_avant_k.json  
sauvegarde_traj_avant_pas_k.json  
sauvegarde_traj_apres_pas_k.json  
sauvegarde_traj_ref_pas_k.json  
sauvegarde_deplacements_apres_pas_k.json  
sauvegarde_deplacements_avant_pas_k.json  
sauvegarde_deplacements_ref_pas_k.json  
sauvegarde_donnees_avant_pas_k.json  
sauvegarde_donnees_apres_pas_k.json  
sauvegarde_donnees_selection_k.json  
sauvegarde_predictions_avant_pas_k.json  
sauvegarde_predictions_apres_pas_k.json  


-----------------------------------------------------------------------------------------------
# Structure de la fonction principale fonction_principale :
-----------------------------------------------------------------------------------------------


-----------------------------------------------------------------------------------------------
## Gestion propagation des particules au pas i.
**dossier trajectoire_pas_i.py**

Sous-fonctions : 
*nouvelle_valeur  
vecteurs_temps_i  
trajectoire_pas_temps  
point_bas  
une_trajectoire_haute_basse*

Fonction principale appelée par fonction_principale : 
*une_trajectoire_pas_i(vec_Xt_eff, vec_Zt, vec_temps_eff, temps_courant, pas_temps, longueur)*

Entrées :  
vec_Xt_eff : vecteur des coordonnées X trajectoire complète  
vec_Zt : vecteur coordonnées Z trajectoire complète  
vec_temps_eff : vecteur du temps mis pas le dike pour atteindre (X,Z)  
temps_courant : temps au pas i  
pas_temps : pas temps  
longueur : longueur maximale du dike  

Sorties :
trajectoire_x_i : vecteur coordonnées X du dike au moment i  
trajectoire_z_i : vecteur coordonnées X du dike au moment i  
X_haut : coordonnées X du front du dike  
Z_haut : coordonnées Z du front du dike  
X_bas :  coordonnées X du bas du dike  
Z_bas : coordonnées Z du bas du dike  
distance_parcourue_i : distance parcourue par le dike depuis le point de départ  

-----------------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------------

## Gestion des déplacements :
**dossier deplacement_pas_i.py**

Sous-fonctions : 
*forward, chinnery, ux_ss, uy_ss, uz_ss, ux_ds, uy_ds, uz_ds, ux_tf, uy_tf, uz_tf, I1, I2, I3, I4, I5  (Okada85)  
calcul_parametres_effectifs_okada  
liaison_modele_okada*

Fonction appelée dans fonction_principale :   
*une_particule_deplacement_pas_i(X_haut, Z_haut, X_bas, Z_bas, longueur, ouverture, distance_parcourue_i, R, G, mu, E, vol, P_load, pas_okada, xmin, xmax, type_observations)*

Entrées :  
X_haut, Z_haut : coordonnées (X,Z) du front du dike  
X_bas, Z_bas : coordonnées (X,Z) du bas du dike  
longueur, ouverture : longueur et ouverture théorique du dike  
distance_parcourue_i : distance parcourue par le dike depuis le point de départ  
R, G, mu, E, vol : paramètres de la particule  
P_load : valeur charge ou décharge  
xmin, xmax : étendue axe des abscisses  
pas_okada : discrétisation de la grille  
type_observations :  'regulier','normal','random'  

Sorties :  
ux : vecteur déplacement ux associé au dike  
uz: vecteur déplacement uz associé au dike  
grille_x : liste des points observés à la surface   
vec_x_c, vec_z_c : point central entre le haut et le bas du dike  
longueur_okada, ouverture_okada : longueur et ouverture pour okada  
vec_dip : valeur du dip  
strike : 180 si X_bas >=0; 0 sinon  

-----------------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------------
## Gestion du calcul de la vraisemblance
**dossier calcul_vraisemblance_poids.py**

Sous-fonctions :  
*calcul_vraisemblances  
calcul_poids*

Fonction principale appelée par fonction_principale :  
*fonction_poids*

Entrées :   
Nb_particules : nombre de particules  
ux_ref : vecteur déplacements ux observés  
uz_ref : vecteur déplacements uz observés  
deplacements_ux : matrice des déplacements ux des n particules  
deplacements_uz : matrice des déplacements uz des n particules   
matrice_covariance : bruit associés aux déplacements  
vec_index_p_surface : dike ayant atteint la surface  

Sorties :  
vec_vraisemblances : vecteur des vraisemblances   
vecteur_poids : vecteur des poids  
coef_particules : terme ln(exp(vraisemblance))  
vraisemblances_classiques : vraisemblances   
vraisemblances_normalisee : vraisemblances normalisées par la valeur de déplacement observés  
coef_normalises : terme ln(exp(vraisemblance)) normalisés  

-----------------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------------
## Gestion du type de sélection
**dossier type_selection.py**

Selon le type de sélection choisie, une fonction est appelée dans fonction_principale :  
*fonction_selection_systematique(Nb_particules, poids_i)  
fonction_selection_stratified(Nb_particules, poids_i)  
fonction_selection_multinomial(Nb_particules, poids_i)  
fonction_selection_residual(Nb_particules, poids_i)*

Entrées :  
Nb_particules : nombre de particules  
poids_i : vecteur des poids  

Sorties :  
liste des particules sélectionnées  

-----------------------------------------------------------------------------------------------


-----------------------------------------------------------------------------------------------
## Gestion du rééchantillonnage

**dossier grilles_R_G pour le rééchantillonnage de R et G (399 combinaisons possibles)**

**dossier algorithme_reechantillonnage.py**

Sous-fonctions :  
*reech_R_G_aleatoire  
trouve_point_haut 
vec_eff_x_z_temps  
traiter_particule*


Fonctions appelées dans fonction_principale :   
*resampling_bruit(Nb_particules, vec_index_p_selec, points_centraux_x_i, points_centraux_z_i, longueurs_okada_i, ouvertures_okada_i, dips_i, strikes_i)*

Entrées :  
Nb_particules : nombre de particules  
vec_index_p_selec : liste des particules sélectionnées  
points_centraux_x_i : liste des coordonnées x des points centraux des dike pour les n particules  
points_centraux_z_i : liste des coordonnées z des points centraux des dike pour les n particules  
longueurs_okada_i : liste des longueurs okada pour les n particules  
ouvertures_okada_i : liste des ouvertures okada les n particules  
dips_i : liste des dips les n particules  
strikes_i : liste des strikes okada les n particules  

Sorties :  
points_centraux_x_i_re : liste des coordonnées x réechantillonnées des points centraux des dike pour les n particules  
points_centraux_z_i_re : liste des coordonnées z réechantillonnées des points centraux des dike pour les n particules  
longueurs_okada_i_re : liste des longueurs rééchantillonnées pour les n particules  
ouvertures_okada_i_re : liste des ouvertures rééchantillonnées pour les n particules  
dips_i_re : liste des dips rééchantillonnées pour les n particules  
strikes_i_re : liste des strikes rééchantillonnées pour les n particules  





*okada_vers_physique(Nb_particules, vec_index_particules, points_centraux_x_i_re,
                        points_centraux_z_i_re, longueurs_okada_i_re, ouvertures_okada_i_re,
                        dips_i_re, strikes_re, temps_courant, i)*


Entrées :  
Nb_particules  
vec_index_particules : vecteur des particules sélectionnées  
points_centraux_x_i_re  
points_centraux_z_i_re  
longueurs_okada_i_re  
ouvertures_okada_i_re  
dips_i_re  
strikes_re  
temps_courant  
i : fênetre d'assimilation  

Sorties :  
sauvegarde_parametres_apres : nouveau set de paramètres  
sauvegarde_trajectoires_apres : nouvelles trajectoires associées aux paramètres  







**dossier trajectoire_particules_reechanti.py**

Sous-fonctions :
*calcul_parametres_reech  
fonction_watanabe_reech  
calcul_sig1_sig3_reech  
streamplot_trajectoire_reech*

Fonction de création des trajectoires des particules rééchantillonnées : 
*trajectoire_une_particule_reech(start_point, R, G, mu, E, vol, temps_courant, pas_temps, xmin, xmax, zmin, zmax, pas_trajectoire, P_load,
                                                               rayon_load, pas_vect, norme)*
                                                               
Entrées :  
start_point : point de départ de la particule  
R, G, mu, E, vol : paramètres de la particules   
temps_courant : temps considéré  
pas_temps : pas de temps   
xmin, xmax : étendue axe abscisses (m)  
zmin, zmax : étendue axe ordonnées (m)  
pas_trajectoire : discrétisation de l'espace (m)   
P_load, rayon_load : données charge ou décharge  
pas_vect : discrétisation de représentation des vecteurs sigma1 (si besoin)  
norme : passage du mètres à une autre unité (si besoin). Pas défaut norme = 1  

Sorties :  
vec_Xt_eff, vec_Zt, vec_temps_eff, temps_courant, pas_temps, longueur, ouverture, vitesse_reech, start_x_retroprop, start_z_retroprog  
vec_Xt_eff : coordonnées X trajectoire complète  
vec_Zt : coordonnée Z trajectoire complète  
vec_temps_eff : temps depuis le début de l'expérience quand (X,Z) est atteint  
longueur : longueur maximale du dike  
ouverture : ouverture du dike  
vitesse_reech : vitesse de propagation  
start_x_retroprop : coordonnée du point de départ x de la trajectoire au temps 0  
start_z_retroprog : coordonnée du point de départ x de la trajectoire au temps 0  

-----------------------------------------------------------------------------------------------
