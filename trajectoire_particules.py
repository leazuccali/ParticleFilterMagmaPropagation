import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
#from affichage_une_particule import affichage_trajectoire, affichage_champs
def calcul_parametres(R,G,mu,E,vol, P_load):
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

    return extension, longueur, vitesse, ouverture


def fonction_watanabe(x_min, x_max, z_min, z_max, pas, P_load, rayon_load, extension):
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


def calcul_sig1_sig3(mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz, pas_vec, G):
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


def streamplot_trajectoire(mesh_X, mesh_Z, u_sig_1, v_sig_1, start_point, norme, R, G, vitesse,temps_i):
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

    vec_Xt_eff = vec_Xt

    s = (np.diff(vec_Xt_eff) ** 2 + np.diff(vec_Zt) ** 2) ** 0.5

    s_zero = np.insert(s, 0,0)
    s_c = np.cumsum(s_zero)


    long_magma = np.max(s_c)

    vec_temps_eff = []
    vec_temps_eff.append(temps_i)
    for i in range(1, len(s_zero)):
        new_t = vec_temps_eff[i-1] + s_zero[i] / vitesse
        vec_temps_eff.append(new_t)

    heure,minute,seconde = sec2hms(vec_temps_eff[-1])


    vec_Xt_eff_l = list(vec_Xt_eff)
    vec_Zt_l = list(vec_Zt)
    vec_temps_l = list(vec_temps_eff)

    return vec_Xt, vec_Xt_eff_l, vec_Zt_l, long_magma, vec_temps_l




def trajectoire_une_particule(start_point, R, G, mu, E, vol, temps_courant, pas_temps, xmin, xmax, zmin, zmax, pas_trajectoire, P_load,
                                                               rayon_load, pas_vect, norme):
    '''

    Args:
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

    Returns:
        vec_Xt_eff : coordonnées X trajectoire complète
        vec_Zt : coordonnée Z trajectoire complète
        vec_temps_eff : temps depuis le début de l'expérience quand (X,Z) est atteint
        longueur : longueur maximale du dike
        ouverture : ouverture du dike
        vitesse : vitesse de propagation

    '''
    x_min = -30000  # m
    x_max = 30000  # m
    z_min = -15000  # m
    z_max = -1  # m

    extension, longueur, vitesse, ouverture = calcul_parametres(R, G, mu, E, vol, P_load)


    mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz = fonction_watanabe(xmin, xmax, zmin, zmax, pas_trajectoire, P_load,
                                                               rayon_load, extension)


    sig_1, sig_3, u_sig_1, v_sig_1, u_sig_1_eff, v_sig_1_eff, u_sig_3, v_sig_3, mesh_vec = calcul_sig1_sig3(
        mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz, pas_vect, G)

    vec_Xt, vec_Xt_eff, vec_Zt, long_magma, vec_temps_eff = streamplot_trajectoire(
        mesh_X, mesh_Z, u_sig_1_eff, v_sig_1_eff, start_point, norme, R, G, vitesse, temps_courant)



    return vec_Xt_eff, vec_Zt, vec_temps_eff, temps_courant, pas_temps, longueur, ouverture, vitesse



def sec2hms(ss):
	(hh, ss)=divmod(ss, 3600)
	(mm, ss)=divmod(ss, 60)
	return (hh, mm, ss)

