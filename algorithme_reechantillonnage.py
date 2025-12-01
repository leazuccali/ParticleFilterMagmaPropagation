import os
import json
from joblib import Parallel, delayed
import bisect
#from essai_fonction_point_haut_reech import *
from trajectoire_particules_reechanti import *
#import random as rd


def resampling_bruit(Nb_particules, vec_index_p_selec, points_centraux_x_i, points_centraux_z_i, longueurs_okada_i, ouvertures_okada_i, dips_i, strikes_i):
    '''

    Args:
        Nb_particules : nombre de particules
        vec_index_p_selec : liste des particules sélectionnées
        points_centraux_x_i : liste des coordonnées x des points centraux des dike pour les n particules
        points_centraux_z_i : liste des coordonnées z des points centraux des dike pour les n particules
        longueurs_okada_i : liste des longueurs okada pour les n particules
        ouvertures_okada_i : liste des ouvertures okada les n particules
        dips_i : liste des dips les n particules
        strikes_i : liste des strikes okada les n particules

    Returns:
        points_centraux_x_i_re : liste des coordonnées x réechantillonnées des points centraux des dike pour les n particules
        points_centraux_z_i_re : liste des coordonnées z réechantillonnées des points centraux des dike pour les n particules
        longueurs_okada_i_re : liste des longueurs rééchantillonnées pour les n particules
        ouvertures_okada_i_re : liste des ouvertures rééchantillonnées pour les n particules
        dips_i_re : liste des dips rééchantillonnées pour les n particules
        strikes_i_re : liste des strikes rééchantillonnées pour les n particules

    '''


    #Construction de nouveaux parametres Okada
    points_centraux_x_i_re = []
    points_centraux_z_i_re = []
    longueurs_okada_i_re = []
    ouvertures_okada_i_re = []
    dips_i_re = []
    strikes_i_re = []
    for i in range(Nb_particules):

        if vec_index_p_selec[i] != i:   #si la valeur du vecteur de reech est différente de l'indice en question
            ind_reech = vec_index_p_selec[i]
            pts_central_x_re = points_centraux_x_i[ind_reech] + np.random.randint(-500, 500)
            # print(pts_central_x_re)
            if pts_central_x_re < -30000:
                pts_central_x_re = -30000
            if pts_central_x_re > 30000:
                pts_central_x_re = 30000
            points_centraux_x_i_re.append(pts_central_x_re)
            # print(points_centraux_x_i_re)

            pts_central_z_re = points_centraux_z_i[ind_reech] + np.random.randint(-200, 200)
            if pts_central_z_re < 100:
                pts_central_z_re = 100
            if pts_central_z_re > 15000:
                pts_central_z_re = 15000
            points_centraux_z_i_re.append(pts_central_z_re)

            val_rd_lg = np.random.uniform(0.5,2)
            lg_okada_i_re = val_rd_lg * longueurs_okada_i[ind_reech]
            #lg_okada_i_re = longueurs_okada_i[ind_reech] + np.random.randint(-500, 500)
            if lg_okada_i_re < 0:
                lg_okada_i_re = longueurs_okada_i[ind_reech]
            longueurs_okada_i_re.append(lg_okada_i_re)


            ouv_okada_i_re = np.random.uniform(0.5,2) * ouvertures_okada_i[ind_reech]
            #ouv_okada_i_re = ouvertures_okada_i[ind_reech] + np.random.uniform(-1, 1)
            if ouv_okada_i_re < 0:
                ouv_okada_i_re = ouvertures_okada_i[ind_reech]
            ouvertures_okada_i_re.append(ouv_okada_i_re)

            dip_i_re = dips_i[ind_reech] + np.random.uniform(-0.25, 0.25)   ## +-0.5 avant
            if dip_i_re > np.pi / 2:
                dip_i_re = np.random.uniform(0.5,0.8) * dips_i[ind_reech]
            dips_i_re.append(dip_i_re)

            stri_i_re = strikes_i[ind_reech]
            strikes_i_re.append(stri_i_re)

        else:
            pts_central_x_re = points_centraux_x_i[i]
            points_centraux_x_i_re.append(pts_central_x_re)

            pts_central_z_re = points_centraux_z_i[i]
            points_centraux_z_i_re.append(pts_central_z_re)

            lg_okada_i_re = longueurs_okada_i[i]
            longueurs_okada_i_re.append(lg_okada_i_re)

            ouv_okada_i_re = ouvertures_okada_i[i]
            ouvertures_okada_i_re.append(ouv_okada_i_re)

            dip_i_re = dips_i[i]
            dips_i_re.append(dip_i_re)

            stri_i_re = strikes_i[i]
            strikes_i_re.append(stri_i_re)

    return points_centraux_x_i_re, points_centraux_z_i_re, longueurs_okada_i_re, ouvertures_okada_i_re, dips_i_re, strikes_i_re


def reech_R_G(nb_grilles, pt_centre_x_p, pt_centre_z_p, dip_p):
    # Paramètres de la grille
    xmin = -30000  # m
    xmax = 30000  # m
    zmin = -11000  # m
    zmax = -1  # m
    pas_trajectoire = 100

    vec_X = np.arange(xmin, xmax, pas_trajectoire)
    vec_Z = np.arange(zmin, zmax, pas_trajectoire)
    mesh_X, mesh_Z = np.meshgrid(vec_X, vec_Z)

    # Plus proche voisin de (x_c, z_c) dans le grille mesh_X mesh_Z
    # Algorithme de plus proches voisins avec une ditance au carré
    distances_squared = (mesh_X - pt_centre_x_p) ** 2 + (mesh_Z - pt_centre_z_p) ** 2
    # On cherche l'indie linéaire de la distance minimale (unravel convertit l'ind linéaire en indice 2D)
    index_min = np.unravel_index(np.argmin(distances_squared), mesh_X.shape)
    #print(index_min)
    #print(index_min[0])
    #print(index_min[1])
    # Point le plus proche de (x_c, z_c) appartenant à la grille
    pt_proche = (mesh_X[index_min], mesh_Z[index_min])
    #print('pt_proche', pt_proche)

    # Recherche parmi les dip des grilles disponibles le dip qui se rapproche le plus de la valeur d'entrée
    directory = "grilles_R_G"
    os.makedirs(directory, exist_ok=True)

    # On veut prendre les dips correspondants
    # Dictionnaire pour récolter les dips
    dico_dips = []

    for num_grille in range(1, nb_grilles + 1):
        # directory_dico = os.path.join(directory, f"grilles_{num_grille}")
        # os.makedirs(directory_dico, exist_ok=True)

        # Liste des angles correspondants au centre (x_c, z_c)
        with open(os.path.join(directory, f"grille_{num_grille}.json"), "r") as f:
            dico_i = json.load(f)

        angles = dico_i['Angles']
        angle_point = angles[index_min[0]][index_min[1]]
        #print(angle_point)
        dico_dips.append(angle_point)

    print(dico_dips)
    # Recherche de l'angle le plus proche de dip
    plus_proche_indice_m, plus_proche_valeur = min(enumerate(dico_dips), key=lambda x: abs(x[1] - dip_p))
    plus_proche_indice = plus_proche_indice_m + 1  # On rajoute 1 pour que la numérotation corresponde aux grilles générées avant


    # On retrouve les R et G correspondants
    with open(os.path.join(directory, f"grille_{plus_proche_indice}.json"), "r") as f:
        dico_i = json.load(f)

    find_R = dico_i['R']
    find_G = dico_i['G']

    return find_R, find_G



def reech_R_G_aleatoire(nb_grilles, pt_centre_x_p, pt_centre_z_p, dip_p):
    # Paramètres de la grille
    xmin = -30000  # m
    xmax = 30000  # m
    zmin = -11000  # m
    zmax = -1  # m
    pas_trajectoire = 100

    vec_X = np.arange(xmin, xmax, pas_trajectoire)
    vec_Z = np.arange(zmin, zmax, pas_trajectoire)
    mesh_X, mesh_Z = np.meshgrid(vec_X, vec_Z)

    # Plus proche voisin de (x_c, z_c) dans le grille mesh_X mesh_Z
    # Algorithme de plus proches voisins avec une ditance au carré
    distances_squared = (mesh_X - pt_centre_x_p) ** 2 + (mesh_Z - pt_centre_z_p) ** 2
    # On cherche l'indice linéaire de la distance minimale (unravel convertit l'ind linéaire en indice 2D)
    index_min = np.unravel_index(np.argmin(distances_squared), mesh_X.shape)

    # Point le plus proche de (x_c, z_c) appartenant à la grille
    pt_proche = (mesh_X[index_min], mesh_Z[index_min])

    # Recherche parmi les dip des grilles disponibles le dip qui se rapproche le plus de la valeur d'entrée
    directory = "grilles_R_G"
    os.makedirs(directory, exist_ok=True)

    # On veut prendre les dips correspondants
    # Dictionnaire pour récolter les dips
    dico_dips = []

    for num_grille in range(0, nb_grilles):
        # directory_dico = os.path.join(directory, f"grilles_{num_grille}")
        # os.makedirs(directory_dico, exist_ok=True)

        # Liste des angles correspondants au centre (x_c, z_c)
        with open(os.path.join(directory, f"grille_{num_grille}.json"), "r") as f:
            dico_i = json.load(f)

        angles = dico_i['Angles']
        angle_point = angles[index_min[0]][index_min[1]]
        #print(angle_point)
        dico_dips.append(angle_point)

    print(dico_dips)



    ##########Une fois qu'on a le dico dips, on veut choisir aleatoirement une valeur d'angle puis remonter au R G
    # Paramètre : pourcentage des meilleurs angles
    x_percent = 10  # Par exemple, 10% des meilleurs angles
    # Calcul des différences absolues
    differences_absolues = [abs(valeur - dip_p) for valeur in dico_dips]
    # Trier les indices en fonction des différences absolues
    indices_triees = np.argsort(differences_absolues)
    # Nombre d'éléments correspondant à x%
    nb_meilleurs = max(1, int(len(differences_absolues) * (x_percent / 100)))
    # Sélectionner les x% meilleurs indices
    meilleurs_indices = indices_triees[:nb_meilleurs]

    # Sélectionner un angle aléatoirement parmi les meilleurs
    indice_aleatoire = np.random.choice(meilleurs_indices)
    print('indice aléatoire', indice_aleatoire)
    #angle_aleatoire = dico_dips[indice_aleatoire]
    #indices_aleatoires_2d = np.unravel_index(indice_aleatoire, mat_angles.shape)
    with open(os.path.join(directory, f"grille_{indice_aleatoire}.json"), "r") as f:
        dico_i = json.load(f)

    find_R = dico_i['R']
    find_G = dico_i['G']
    print(f"Les paramètres correspondants sont R = {find_R} et G = {find_G}")



    return find_R, find_G








def trouve_point_haut(x_bas, z_bas, vec_Xt, vec_Zt, longueur_okada):
    ## Calcul distance cumulées
    diff = (np.diff(vec_Xt) ** 2 + np.diff(vec_Zt) ** 2) ** 0.5
    #print('diff', diff)
    diff_zero = np.insert(diff, 0, 0)
    somme_cumu = np.cumsum(diff_zero)
    #print('somme_cumu', somme_cumu)

    # print(len(somme_cumu))
    ## Intercalle la longueur_okada sur le vecteur des distances cumulées
    index_longueur = bisect.bisect(somme_cumu, longueur_okada)
    #print('index_longueur', index_longueur)

    if index_longueur >= len(vec_Xt):
        x_haut = vec_Xt[-1]  # peut etre mieux de mettre le max de la grille pour reech directement après ?
        z_haut = vec_Zt[-1]

    else:

        # Redef indice moins et plus
        indice_moins = index_longueur - 1
        indice_plus = index_longueur

        '''if indice_moins >= index_longueur:       #  or indice_plus >= len(somme_cumu):
            indice_moins = index_longueur - 1
            indice_plus = index_longueur'''

        # Def des nouvelles valeurs
        somme_c_moins = somme_cumu[indice_moins]
        somme_c_plus = somme_cumu[indice_plus]


        x_moins = vec_Xt[indice_moins]
        z_moins = vec_Zt[indice_moins]
        x_plus = vec_Xt[indice_plus]
        z_plus = vec_Zt[indice_plus]

        reste_a_parcourir = longueur_okada - somme_c_moins
        #print('rest a parcourir', reste_a_parcourir)
        distance_2_points = diff[indice_moins]
        #print('distance entre 2 points', distance_2_points)
        coef = reste_a_parcourir / distance_2_points
        #print('coef', coef)

        x_haut = (1 - coef) * x_moins + coef * x_plus
        z_haut = (1 - coef) * z_moins + coef * z_plus

        #print("xh zh ", x_haut, z_haut)

        '''plt.figure()
        plt.plot(vec_Xt[:5000], vec_Zt[:5000])
        plt.plot(x_haut, z_haut, 'bo')
        plt.show()'''



    return x_haut, z_haut, index_longueur


def vec_eff_x_z_temps(vec_Xt, vec_Zt, x_haut, z_haut, indice_insertion, temps_courant, vitesse):
    ## On veut insérer x_haut et z_haut dans vec_Xt et vec_Zt puis couper ce qu'il y a en dessous
    vec_Xt.insert(indice_insertion, x_haut)
    vec_Zt.insert(indice_insertion, z_haut)
    #print('taille vec Xt avec insertion', len(vec_Xt))
    vec_Xt_haut = vec_Xt[indice_insertion:]
    vec_Zt_haut = vec_Zt[indice_insertion:]

    ## Calcul des distances entre les points
    diff = (np.diff(vec_Xt_haut) ** 2 + np.diff(vec_Zt_haut) ** 2) ** 0.5
    diff_zero = np.insert(diff, 0, 0)

    taille_vec_Xt_h = len(vec_Xt_haut)
    #print('taille vec Xh', taille_vec_Xt_h)

    vec_temps_haut = []
    vec_temps_haut.append(temps_courant)
    for i in range(1, len(diff_zero)):
        new_temps = vec_temps_haut[i-1] + diff_zero[i] / vitesse
        vec_temps_haut.append(new_temps)

    #print("taille de vec temps", len(vec_temps_haut))

    vec_X_bas = vec_Xt[:indice_insertion]
    taille_vec_Xt_bas = len(vec_X_bas)
    vec_temps_bas = [0] * taille_vec_Xt_bas


    #print("taille vec_X_bas", len(vec_X_bas))
    #print(vec_X_bas)

    vec_temps = vec_temps_bas + vec_temps_haut
    #print("vec temps", vec_temps)
    #print("taille vec temps", len(vec_temps))


    return vec_Xt, vec_Zt, vec_temps



def traiter_particule(p, vec_index_particules, points_centraux_x_i_re, points_centraux_z_i_re,
                      longueurs_okada_i_re, ouvertures_okada_i_re,
                      dips_i_re, strikes_re, temps_courant, i):

    directory = 'output_data_and_figures'
    os.makedirs(directory, exist_ok=True)
    step_directory = os.path.join(directory, f'step_{i}')
    directory_particule = os.path.join(step_directory, f'particule_{p}')
    os.makedirs(step_directory, exist_ok=True)
    os.makedirs(directory_particule, exist_ok=True)
    data_filename_param = os.path.join(directory_particule, f'param_particules_{p}_step_{i}_apres.json')
    data_filename_traj = os.path.join(directory_particule, f'trajectoire_particules_{p}_step_{i}_apres.json')

    sauvegarde_parametres_apres = {}
    sauvegarde_trajectoires_apres = {}

    if (vec_index_particules[p] != p):

        ind_reech = vec_index_particules[p] #indice de reechantillonnage

        directory = 'output_data_and_figures'
        step_directory = os.path.join(directory, f'step_{i}')

        directory_particule = os.path.join(step_directory, f'particule_{ind_reech}')    ## au lieu de p
        data_filename = os.path.join(directory_particule, f'param_particules_{ind_reech}_step_{i}_avant.json')   ## au lieu de p
        with open(data_filename, 'r') as fichier:
            param_physiques_p_selec = json.load(fichier)

        volume_compare = param_physiques_p_selec['vol_0']
        G_compare = param_physiques_p_selec['G_0']
        mu_compare = param_physiques_p_selec['mu_0']
        constante = 5 * 10 ** (-7)
        vitesse_compare = (constante * volume_compare * G_compare) / mu_compare
        R_compare = param_physiques_p_selec['R_0']








        xmin = -30000  # m
        xmax = 30000  # m
        zmin = -15000  # m
        zmax = -1  # m
        pas_trajectoire = 100
        pas_vect = 2000
        P_load = -15000000  # MPa
        rayon_load = 10000  # m
        norme = 1  # Si norme = 1, alors l'unité est le mètre. Si norme = 1000, alors l'unité est le km.
        pas_okada = 100
        pas_temps = 60

        pt_centre_x_p = points_centraux_x_i_re[p]
        pt_centre_z_p = - points_centraux_z_i_re[p]
        dip_p = dips_i_re[p]
        longueur_okada_p = longueurs_okada_i_re[p] / 2

        l_sin = longueur_okada_p * math.sin(dip_p)
        l_cos = longueur_okada_p * math.cos(dip_p)

        z_bas = pt_centre_z_p - l_sin
        if pt_centre_x_p >= 0:
            x_bas = pt_centre_x_p - l_cos
        else:
            x_bas = pt_centre_x_p + l_cos

        if z_bas > -100:
            z_bas = -1000
        if z_bas < -15000:
            z_bas = -10000
        if x_bas < -30000:
            x_bas = -29000
        if x_bas > 30000:
            x_bas = 29000



        start_point = np.array([[x_bas, z_bas]])

        ###############################################################################################################
        nouveau_R, nouveau_G = reech_R_G_aleatoire(399, pt_centre_x_p, pt_centre_z_p, dip_p)
        ###############################################################################################################

        ouv_okada_p = ouvertures_okada_i_re[p]
        long_okada_p = longueurs_okada_i_re[p]

        nouveau_vol_1 = 2 * ouv_okada_p * (long_okada_p ** 2)
        nouveau_vol = nouveau_vol_1 / 3.14

        vitesse = vitesse_compare + np.random.uniform(-0.1, 0.1)
        #vitesse = max(vitesse, 0.01)
        if vitesse < 0:
            vitesse = vitesse_compare


        C = 5 * 10 ** (-7)
        nouveau_mu = (C * nouveau_G * nouveau_vol) / vitesse

        nouveau_E = (3.14 * long_okada_p * nouveau_G * (1 - 0.25 ** 2) * abs(P_load)) / (2 * ouv_okada_p)

        #sauvegarde_parametres_apres[p] = {}
        sauvegarde_parametres_apres['Numero'] = p
        sauvegarde_parametres_apres['x_0'] = x_bas
        sauvegarde_parametres_apres['z_0'] = z_bas
        sauvegarde_parametres_apres['R_0'] = nouveau_R
        sauvegarde_parametres_apres['G_0'] = nouveau_G
        sauvegarde_parametres_apres['mu_0'] = nouveau_mu
        sauvegarde_parametres_apres['E_0'] = nouveau_E
        sauvegarde_parametres_apres['vol_0'] = nouveau_vol

        vec_Xt_eff, vec_Zt, vec_temps_eff, temps_courant, pas_temps, longueur_p, ouverture_p, vitesse_reech, start_x_retropropag, start_z_retropropag= trajectoire_une_particule_reech(
            start_point, nouveau_R, nouveau_G, nouveau_mu, nouveau_E, nouveau_vol, temps_courant, pas_temps, xmin, xmax,
            zmin, zmax, pas_trajectoire, P_load, rayon_load, pas_vect, norme)

        x_haut_exp, z_haut_exp, indice_insertion = trouve_point_haut(x_bas, z_bas, vec_Xt_eff, vec_Zt, long_okada_p)

        vec_Xt_complet, vec_Zt_complet, vec_temps = vec_eff_x_z_temps(vec_Xt_eff, vec_Zt, x_haut_exp, z_haut_exp,
                                                                      indice_insertion, temps_courant, vitesse_reech)





        #sauvegarde_trajectoires_apres[p] = {}
        sauvegarde_trajectoires_apres['Numero'] = p
        sauvegarde_trajectoires_apres['X'] = vec_Xt_complet
        sauvegarde_trajectoires_apres['Z'] = vec_Zt_complet
        sauvegarde_trajectoires_apres['Temps effectif'] = vec_temps
        sauvegarde_trajectoires_apres['Longueur remontee'] = long_okada_p
        sauvegarde_trajectoires_apres['Ouverture'] = ouv_okada_p
        sauvegarde_trajectoires_apres['Vitesse'] = vitesse_reech

        sauvegarde_parametres_apres['Start x retropropag'] = start_x_retropropag
        sauvegarde_parametres_apres['Start z retropropag'] = start_z_retropropag

        #### Sauvegarde des fichiers des données apres selection
        with open(data_filename_param, 'w') as fichier :
            json.dump(sauvegarde_parametres_apres, fichier)

        with open(data_filename_traj, 'w') as fichier :
            json.dump(sauvegarde_trajectoires_apres, fichier)


    else:
        data_filename_param_avant = os.path.join(directory_particule, f'param_particules_{p}_step_{i}_avant.json')
        data_filename_traj_avant = os.path.join(directory_particule, f'trajectoire_particules_{p}_step_{i}_avant.json')



        with open(data_filename_param_avant, 'r') as fichier:
            param_physiques_p = json.load(fichier)
        with open(data_filename_traj_avant, 'r') as fichier:
            trajectoire_p = json.load(fichier)

        #sauvegarde_parametres_apres[p] = {}
        sauvegarde_parametres_apres['Numero'] = p
        sauvegarde_parametres_apres['x_0'] = param_physiques_p['x_0']
        sauvegarde_parametres_apres['z_0'] = param_physiques_p['z_0']
        sauvegarde_parametres_apres['R_0'] = param_physiques_p['R_0']
        sauvegarde_parametres_apres['G_0'] = param_physiques_p['G_0']
        sauvegarde_parametres_apres['mu_0'] = param_physiques_p['mu_0']
        sauvegarde_parametres_apres['E_0'] = param_physiques_p['E_0']
        sauvegarde_parametres_apres['vol_0'] = param_physiques_p['vol_0']


        cle_x = param_physiques_p.get('Start x retropropag')
        if cle_x is not None:
            sauvegarde_parametres_apres['Start x retropropag'] = param_physiques_p['Start x retropropag']
            sauvegarde_parametres_apres['Start z retropropag'] = param_physiques_p['Start z retropropag']


        #sauvegarde_trajectoires_apres[p] = {}
        sauvegarde_trajectoires_apres['Numero'] = p
        sauvegarde_trajectoires_apres['X'] = trajectoire_p['X']
        sauvegarde_trajectoires_apres['Z'] = trajectoire_p['Z']
        sauvegarde_trajectoires_apres['Temps effectif'] = trajectoire_p['Temps effectif']
        sauvegarde_trajectoires_apres['Longueur remontee'] = trajectoire_p['Longueur remontee']
        sauvegarde_trajectoires_apres['Ouverture'] = trajectoire_p['Ouverture']
        sauvegarde_trajectoires_apres['Vitesse'] = trajectoire_p['Vitesse']

        #### Sauvegarde des fichiers des données apres selection
        with open(data_filename_param, 'w') as fichier:
            json.dump(sauvegarde_parametres_apres, fichier)

        with open(data_filename_traj, 'w') as fichier:
            json.dump(sauvegarde_trajectoires_apres, fichier)


    return sauvegarde_parametres_apres, sauvegarde_trajectoires_apres


def okada_vers_physique(Nb_particules, vec_index_particules, points_centraux_x_i_re,
                        points_centraux_z_i_re, longueurs_okada_i_re, ouvertures_okada_i_re,
                        dips_i_re, strikes_re, temps_courant, i):

    '''

    Args:
        points_centraux_x_i_re : liste des coordonnées x réechantillonnées des points centraux des dike pour les n particules
        points_centraux_z_i_re : liste des coordonnées z réechantillonnées des points centraux des dike pour les n particules
        longueurs_okada_i_re : liste des longueurs rééchantillonnées pour les n particules
        ouvertures_okada_i_re : liste des ouvertures rééchantillonnées pour les n particules
        dips_i_re : liste des dips rééchantillonnées pour les n particules
        strikes_i_re : liste des strikes rééchantillonnées pour les n particules
        i: fenêtre d'assimilation

    Returns:
        sauvegarde_parametres_apres : nouveau set de paramètres
    sauvegarde_trajectoires_apres : nouvelles trajectoires associées aux paramètres

    '''
    # Utilisation de joblib pour paralléliser la boucle
    results = Parallel(n_jobs=-1)(
        delayed(traiter_particule)(
            p, vec_index_particules, points_centraux_x_i_re, points_centraux_z_i_re, longueurs_okada_i_re,
            ouvertures_okada_i_re,
            dips_i_re, strikes_re, temps_courant, i
        ) for p in range(Nb_particules)
    )

    # Collecter les résultats
    for result in results:
        sauvegarde_parametres, sauvegarde_trajectoires = result


    return sauvegarde_parametres, sauvegarde_trajectoires








