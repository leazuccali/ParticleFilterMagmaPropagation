import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
#from affichage_une_particule import affichage_trajectoire, affichage_champs
def calcul_parametres_reech(R,G,mu,E,vol, P_load):
    '''
    Calcul des paramètres extension, longueur de la fissure magmatique, et vitesse du magma à partir des rapports R, G,
    de la viscosité du magma, de l'élasticité de la croute et du volume de magma.
    :param R
    :param G
    :param viscosite
    :param rigidite
    :param volume
    :return: extension, longueur, vitesse
    '''
    #print("Début calcul paramètres")
    #Rapport entre la décharge et l'extension
    extension = abs(P_load) * R
    #G : pas de calcul pour l'instant. Pression exercée par le magma sur la roche

    #Longueur de la remontée
    longueur3 = (E * vol) / ((1 - 0.25**2)*G*abs(P_load))
    #print(longueur3)
    longueur = longueur3**(1/3)

    #Vcalitesse
    C = 5 * 10**(-7)
    vitesse = (C * vol * G)/mu

    #Ouverture
    ouverture = (math.pi * vol) / ( 2 * longueur**2)

    #print("L'extension est de ", extension)
    #print("La longueur maximale théorique de la remontée magmatique est ", longueur, "m")
    #print("La vitesse est de ", vitesse, "et la constante C est ", C)
    #print("L'ouverture moyenne de la fissure magmatique est de ", ouverture, "m")



    return extension, longueur, vitesse, ouverture


def fonction_watanabe_reech(x_min, x_max, z_min, z_max, pas, P_load, rayon_load, extension):
    '''
    Calcule les champs de contrainte selon les axes xx, zz et xz selon les équations données par watanabe.

    :param x_min: abscisse minimale
    :param x_max: abscisse maximale
    :param z_min: ordonnée minimale
    :param z_max: ordonnée maximale
    :param pas: distance entre chaque point abscisse/ordonnée
    :param P_load: charge appliquée sur la surface, symétrique par rapport à l'axe des ordonnées
    :param rayon_load: rayon de la charge appliquée
    :param extension : valeur de l'extension de la caldera
    :return: mesh_X, mesh_Z, et les champs de contrainte selon les axes xx, zz et xz
    '''

    #extension = abs(P_load) * R
    vec_X = np.arange(x_min, x_max, pas)
    vec_Z = np.arange(z_min, z_max, pas)
    mesh_X, mesh_Z = np.meshgrid(vec_X,vec_Z)

    theta_1 = np.arctan2(-mesh_Z, mesh_X - rayon_load)
    theta_2 = np.arctan2(-mesh_Z, mesh_X + rayon_load)
    diff_angle = theta_1 - theta_2

    rapport_plus = ((mesh_X + rayon_load) * (-mesh_Z)) / ((mesh_X + rayon_load) ** 2 + mesh_Z ** 2)
    rapport_moins = ((mesh_X - rayon_load) * (-mesh_Z)) / ((mesh_X - rayon_load) ** 2 + mesh_Z ** 2)

    r1_2 = (mesh_X - rayon_load) ** 2 + mesh_Z ** 2
    r2_2 = (mesh_X + rayon_load) ** 2 + mesh_Z ** 2

    sig_xx = (P_load / np.pi) * (diff_angle - rapport_plus + rapport_moins) - extension
    sig_zz = (P_load / np.pi) * (diff_angle + rapport_plus - rapport_moins)
    sig_xz = - (P_load / np.pi) * (mesh_Z ** 2 * (r2_2 - r1_2)) / (r1_2 * r2_2)


    return mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz


def calcul_sig1_sig3_reech(mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz, pas_vec, G):
    '''
    Calcule les champs de contrainte maximum et minimum Sigma1 et Sigma3 à partir de Sigma_xx, Sigma_zz,
    Sigma_xz. Affiche ces champs de vecteurs.

    :param mesh_X: meshgrid à partir du vecteur xmin-xmax divisé en un certain nombre de pas
    :param mesh_Z: meshgrid à partir du vecteur ymin-ymax divisé en un certain nombre de pas
    :param sig_xx: champs de contrainte selon l'axe xx
    :param sig_zz: champs de contrainte selon l'axe zz
    :param sig_xz: champs de contrainte selon l'axe xz
    :param pas_vec: espacement des vecteurs représentés pour sigma1
    :return: sig_1, sig_3, u_sig_1, v_sig_1, u_sig_3, v_sig_3, I (les valeurs Sigma_1, Sigma_3, et les coordonnées des vecteurs associés).
    '''


    # Calcul du champs de contrainte maximum sigma1 et du champs de contrainte minimum sigma3
    nb_l = np.shape(mesh_X)[0]
    nb_c = np.shape(mesh_X)[1]

    sig_1 = np.zeros((nb_l, nb_c))
    sig_3 = np.zeros((nb_l, nb_c))
    u_sig_1 = np.zeros((nb_l, nb_c))
    v_sig_1 = np.zeros((nb_l, nb_c))
    u_sig_3 = np.zeros((nb_l, nb_c))
    v_sig_3 = np.zeros((nb_l, nb_c))

    for i in range(nb_l):
        for j in range(nb_c):
            matrice_contraintes = [[sig_xx[i][j], sig_xz[i][j]],
                                   [sig_xz[i][j], sig_zz[i][j]]]

            valp, vecp = np.linalg.eig(matrice_contraintes)
            sorted_indexes = np.argsort(valp)
            val_propre_ordonne = valp[sorted_indexes]
            vect_propre_ordonne = vecp[:, sorted_indexes]

            sig_1[i][j] = val_propre_ordonne[1]
            sig_3[i][j] = val_propre_ordonne[0]

            u_sig_3[i][j] = vect_propre_ordonne[0][0]
            v_sig_3[i][j] = vect_propre_ordonne[0][1]

            u_sig_1[i][j] = vect_propre_ordonne[0][1]
            v_sig_1[i][j] = vect_propre_ordonne[1][1]


            if v_sig_1[i][j] < 0:
                u_sig_1[i][j] = - u_sig_1[i][j]
                v_sig_1[i][j] = - v_sig_1[i][j]

            if v_sig_3[i][j] < 0:
                u_sig_3[i][j] = - u_sig_3[i][j]
                v_sig_3[i][j] = - v_sig_3[i][j]




    # Affichage du champs de contrainte minimal sigma3 et de la direction du champs de vecteur maximal sigma1
    u_sig_1_eff = (1 - G) * u_sig_1
    v_sig_1_eff = (1 - u_sig_1_eff**2)**(1/2)

    mesh_vec = ((mesh_X%pas_vec == 0) & (mesh_Z%pas_vec == 0))

    return sig_1, sig_3, u_sig_1, v_sig_1, u_sig_1_eff, v_sig_1_eff, u_sig_3, v_sig_3, mesh_vec


def streamplot_trajectoire_reech(mesh_X, mesh_Z, u_sig_1, v_sig_1, start_point, norme, R, G, vitesse,temps_i, longueur):
    '''
    Affichage de la trajectoire du magma le long du champs de contrainte maximal sigma_1.
    Extraction des coordonnées des points de la trajectoire et calcul de la longueur maximale de la trajectoire.
    :param mesh_X: meshgrid à partir du vecteur xmin-xmax divisé en un certain nombre de pas
    :param mesh_Z: meshgrid à partir du vecteur ymin-ymax divisé en un certain nombre de pas
    :param u_sig_1: coordonnée u du champs de contrainte maximal
    :param v_sig_1: coordonnée v du champs de contrainte maximal
    :param start_point: point de départ de la trajectoire considérée :np.array([[2,3]]). Ordonnée minimale de l'affichage.
    :param R: valeurs du rapport entre l'extension et la décharge
    :param G: valeur du rapport entre la pression exercée par le magma et la décharge
    :param norme: "normalisation" des distances (passage du m au km)
    :param R : rapport entre l'extension et la décharge
    :param G : rapport entre la poussée exercée par le magma et la décharge
    :param vitesse : vitesse de remontée du magma
    :param pas_temps: Pas de temps voulu entre deux points
    :return: vec_Xt, vec_Zt, coordonnées des points de la trajectoire du magma.
    '''

    fig1, ax1 = plt.subplots(figsize=(8,7))
    start = start_point/norme
    mesh_X_new = mesh_X/norme
    mesh_Z_new = mesh_Z/norme
    strs = ax1.streamplot(mesh_X_new, mesh_Z_new, u_sig_1, v_sig_1, start_points=start, density=1000)
    plt.title("Trajectoire du magma pour R={}".format(R))
    plt.xlabel("X (km)")
    plt.ylabel("Z (km)")
    plt.close(fig1)

    segments = strs.lines.get_segments()
    #print("taille segmeents", len(segments))

    vec_coords = []
    for seg in segments:
        t_temp = np.concatenate(seg)
        nt_temps = np.reshape(t_temp, (2, 2))

        vec_coords.append(nt_temps[0])
        vec_coords.append(nt_temps[1])

    df_coords = pd.DataFrame(vec_coords)
    #print('coord', df_coords)
    df_coords_new = df_coords.drop(df_coords[(df_coords[1] < start[0][1])].index)
    df_coords_simple = df_coords_new.drop_duplicates()

    vec_coords_ok = df_coords_simple.to_numpy()

    vec_Xt = vec_coords_ok[:, 0]
    vec_Zt = vec_coords_ok[:, 1]

    #vec_Xt_eff = G * start[0][0] + (1 - G) * vec_Xt
    vec_Xt_eff = vec_Xt
    #print("nombre de points", len(vec_Xt_eff))

    #print(vec_Xt_eff)
    #print(vec_Zt)

    s = (np.diff(vec_Xt_eff) ** 2 + np.diff(vec_Zt) ** 2) ** 0.5
    # print("s=", s)
    #print('taille de s =', len(s))
    s_zero = np.insert(s, 0,0)
    s_c = np.cumsum(s_zero)
    #s_c = np.insert(s_c, 0, 0)
    #print(s_c)
    #print("taille de s_c=", len(s_c))

    long_magma = np.max(s_c)
    print("La longueur de la trajectoire streamplot du magma est de", long_magma, "m.")
    #print("nombre de différences =", len(s))


    vec_temps_eff = []
    vec_temps_eff.append(temps_i)
    for i in range(1, len(s_zero)):
        new_t = vec_temps_eff[i-1] + s_zero[i] / vitesse
        vec_temps_eff.append(new_t)

    #print("tps=", vec_temps_eff)
    #print("taille du vec temps", len(vec_t_eff))
    #print("Temps final=", vec_temps_eff[-1], "secondes")
    heure,minute,seconde = sec2hms(vec_temps_eff[-1])
    #print("Le magma met {:7.1f}".format(vec_temps_eff[-1]) + " secondes pour parvenir à la surface, soit {:3.0f} heures {:2.0f} minutes et {:2.0f} secondes. ".format(heure,minute,seconde) + ".")



    vec_Xt_eff_l = list(vec_Xt_eff)
    vec_Zt_l = list(vec_Zt)
    vec_temps_l = list(vec_temps_eff)




    ## Retroprapagation : on cherche le point de départ de la particule
    df_coords = pd.DataFrame(vec_coords)
    # print('coord', df_coords)
    df_coords_inf = df_coords.drop(df_coords[(df_coords[1] > start[0][1])].index)
    df_coords_simple_inf = df_coords_inf.drop_duplicates()

    vec_coords_inf = df_coords_simple_inf.to_numpy()

    vec_Xt_inf = vec_coords_inf[:, 0]
    vec_Zt_inf = vec_coords_inf[:, 1]

    #print(vec_Xt_inf[20:], vec_Zt_inf[20:])
    #Slicing pour inverser l'ordre
    vec_Xt_inf_inv = vec_Xt_inf[::-1]
    vec_Zt_inf_inv = vec_Zt_inf[::-1]

    # Calcul des distances
    s_inf = (np.diff(vec_Xt_inf) ** 2 + np.diff(vec_Zt_inf) ** 2) ** 0.5
    # print("s=", s)
    # print('taille de s =', len(s))
    s_zero_inf = np.insert(s_inf, 0, 0)
    #s_c_inf = np.cumsum(s_zero)
    somme_cumulee_inf = np.cumsum(s_zero_inf)

    # Construction du vecteur temps
    vec_temps_eff_inf = []
    vec_temps_eff_inf.append(temps_i)
    for i in range(1, len(s_zero_inf)):
        new_t = vec_temps_eff_inf[i - 1] - s_zero_inf[i] / vitesse
        if new_t < 0:
            break
        vec_temps_eff_inf.append(new_t)

    #print("essai vec temps inf", vec_temps_eff_inf)

    len_vec_temps_inf = len(vec_temps_eff_inf)
    vec_Xt_inf_inv_coupe = vec_Xt_inf_inv[:len_vec_temps_inf]
    vec_Zt_inf_inv_coupe = vec_Zt_inf_inv[:len_vec_temps_inf]
    somme_cumulee_inf_coupe = somme_cumulee_inf[:len_vec_temps_inf]
    #print("vec_Xt inv", vec_Xt_inf_inv_coupe)
    #print("vec_Zt inv", vec_Zt_inf_inv_coupe)

    # #point_depart_x_retropropag = vec_Xt_inf_inv_coupe[-1]
    # #point_depart_z_retropropag = vec_Zt_inf_inv_coupe[-1]
    # #print(point_depart_x_retropropag, point_depart_z_retropropag)
    # ### On détermine le point x de ec_Xt_inf_inv_coupe le plus proche de 0
    #
    # # positives = [(i, x) for i, x in enumerate(vec_Xt_inf_inv_coupe) if x > 0]
    # # indice_min = min(positives, key=lambda t: t[1])[0]
    # indice_min = min(enumerate(vec_Xt_inf_inv_coupe), key=lambda t: abs(t[1]))[0]
    # print("L'indice de la plus petite valeur positive est : ", indice_min)
    # point_depart_x_retropropag = vec_Xt_inf_inv_coupe[indice_min]
    # point_depart_z_retropropag = vec_Zt_inf_inv[indice_min]

    point_depart_x_retropropag = vec_Xt_inf_inv_coupe[-1]
    point_depart_z_retropropag = vec_Zt_inf_inv_coupe[-1]
    vitesse_reech = vitesse

    if start[0][0] >= 0:
        if point_depart_x_retropropag < 0:
            ## Recherche de l'indice positif du x le plus proche de 0
            positives = [(i, x) for i, x in enumerate(vec_Xt_inf_inv_coupe) if x > 0]
            indice_min = min(positives, key=lambda t: t[1])[0]
            ### on va chercher le nouveau point de départ correspondant
            point_depart_x_retropropag = vec_Xt_inf_inv_coupe[indice_min]
            point_depart_z_retropropag = vec_Zt_inf_inv_coupe[indice_min]
            distance_parcourue = somme_cumulee_inf_coupe[indice_min]
            vitesse_reech = distance_parcourue / temps_i

    else:
        if point_depart_x_retropropag > 0:
            negatives = [(i, x) for i, x in enumerate(vec_Xt_inf_inv_coupe) if x < 0]
            indice_min = min(negatives, key=lambda t: t[1])[0]
            ### on va chercher le nouveau point de départ correspondant
            point_depart_x_retropropag = vec_Xt_inf_inv_coupe[indice_min]
            point_depart_z_retropropag = vec_Zt_inf_inv_coupe[indice_min]
            distance_parcourue = somme_cumulee_inf_coupe[indice_min]
            vitesse_reech = distance_parcourue / temps_i



    #On retourne pour concatener les valeurs inf et supérieures
    vec_Xt_inf_coupe = vec_Xt_inf_inv_coupe[::-1]
    vec_Zt_inf_coupe = vec_Zt_inf_inv_coupe[::-1]
    vec_temps_inf = vec_temps_eff_inf[::-1]

    vec_Xt_inf_final = vec_Xt_inf_coupe[:-1]
    vec_Zt_inf_final = vec_Zt_inf_coupe[:-1]
    vec_temps_inf_final = vec_temps_inf[:-1]

    #print("vec_Xt coupe", vec_Xt_inf_final, len(vec_Xt_inf_final))
    #print("vec_Zt coupe", vec_Zt_inf_final, len(vec_Zt_inf_final))
    #print("vec_temps_eff_inf coupe", vec_temps_inf_final, len(vec_temps_inf_final))

    vec_Xt_entier = np.concatenate((vec_Xt_inf_final,vec_Xt_eff))
    vec_Zt_entier = np.concatenate((vec_Zt_inf_final, vec_Zt))
    vec_temps_entier = np.concatenate((vec_temps_inf_final,vec_temps_eff))

    #print("vec Xt entier", vec_Xt_entier, len(vec_Xt_entier))
    #print("vec_Zt entier", vec_Zt_entier, len(vec_Zt_entier))
    #print("vec_temps_eff entier", vec_temps_entier, len(vec_temps_entier))

    # Passage en liste
    vec_Xt_entier_l = list(vec_Xt_entier)
    vec_Zt_entier_l = list(vec_Zt_entier)
    vec_temps_entier_l = list(vec_temps_entier)

    #vec_Xt, vec_Xt_entier_l, vec_Zt_entier_l, long_magma, vec_temps_entier_l, point_depart_x_retropropag, point_depart_z_retropropag

    return vec_Xt_eff_l, vec_Zt_l, long_magma, vec_temps_l, point_depart_x_retropropag, point_depart_z_retropropag, vitesse_reech




def trajectoire_une_particule_reech(start_point, R, G, mu, E, vol, temps_courant, pas_temps, xmin, xmax, zmin, zmax, pas_trajectoire, P_load,
                                                               rayon_load, pas_vect, norme):
    '''

    Args:
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

    Returns:
        vec_Xt_eff, vec_Zt, vec_temps_eff, temps_courant, pas_temps, longueur, ouverture, vitesse_reech, start_x_retroprop, start_z_retroprog
        vec_Xt_eff : coordonnées X trajectoire complète
        vec_Zt : coordonnée Z trajectoire complète
        vec_temps_eff : temps depuis le début de l'expérience quand (X,Z) est atteint
        longueur : longueur maximale du dike
        ouverture : ouverture du dike
        vitesse_reech : vitesse de propagation
    start_x_retroprop : coordonnée du point de départ x de la trajectoire au temps 0
    start_z_retroprog : coordonnée du point de départ x de la trajectoire au temps 0

    '''
    x_min = -30000  # m
    x_max = 30000  # m
    z_min = -15000  # m
    z_max = -1  # m

    extension, longueur, vitesse, ouverture = calcul_parametres_reech(R, G, mu, E, vol, P_load)
    #print("La longueur de la remontée magmatique est de :", longueur, "m.")

    mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz = fonction_watanabe_reech(xmin, xmax, zmin, zmax, pas_trajectoire, P_load,
                                                               rayon_load, extension)

    #aff1 = affichage_champs(x_min, x_max, z_min, z_max, mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz, norme)

    sig_1, sig_3, u_sig_1, v_sig_1, u_sig_1_eff, v_sig_1_eff, u_sig_3, v_sig_3, mesh_vec = calcul_sig1_sig3_reech(
        mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz, pas_vect, G)

    vec_Xt_eff, vec_Zt, long_magma, vec_temps_eff, start_x_retroprop, start_z_retroprog, vitesse_reech = streamplot_trajectoire_reech(
        mesh_X, mesh_Z, u_sig_1_eff, v_sig_1_eff, start_point, norme, R, G, vitesse, temps_courant, longueur)



    return vec_Xt_eff, vec_Zt, vec_temps_eff, temps_courant, pas_temps, longueur, ouverture, vitesse_reech, start_x_retroprop, start_z_retroprog



#G, E, P_load, pas_okada, xmin, xmax

def sec2hms(ss):
	(hh, ss)=divmod(ss, 3600)
	(mm, ss)=divmod(ss, 60)
	return (hh, mm, ss)

