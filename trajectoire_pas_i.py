import numpy as np
import pandas as pd

#from trajectoire_particules import *
import bisect
def nouvelle_valeur(vec_Xt, vec_Zt, vec_temps, valeur_temps):
    '''
    Computes the point located on the line between two already computed points, based on the distance at which the
    point is located, for a defined time step.

    Args:
        vec_Xt:Vector of X values
        vec_Zt:Vector of Z values
        vec_temps:Vector of the times corresponding to the (X,Z) coordinates
        valeur_temps:Time value for which the coordinates are sought.

    Returns:
        nouveau_x, nouveau_z: the (X,Z) coordinates of a point at a given time.

    '''
    index_temps = bisect.bisect(vec_temps, valeur_temps)


    if (index_temps == 0):
        nouveau_x = vec_Xt[index_temps]
        nouveau_z = vec_Zt[index_temps]

        insertion_x = vec_Xt.insert(index_temps, nouveau_x)
        insertion_z = vec_Zt.insert(index_temps, nouveau_z)
        insertion_temps = vec_temps.insert(index_temps, valeur_temps)

    elif (index_temps == len(vec_temps)):
        nouveau_x = vec_Xt[index_temps - 1]
        nouveau_z = vec_Zt[index_temps - 1]

        insertion_x = vec_Xt.insert(index_temps, nouveau_x)
        insertion_z = vec_Zt.insert(index_temps, nouveau_z)
        insertion_temps = vec_temps.insert(index_temps, valeur_temps)

    else:
        index_moins = index_temps - 1
        index_plus = index_temps + 1

        diff_moins = valeur_temps - vec_temps[index_moins]
        diff_plus = vec_temps[index_plus - 1] - valeur_temps
        diff_totale = vec_temps[index_plus - 1] - vec_temps[index_moins]

        coef1 = (diff_totale - diff_moins) / diff_totale
        coef2 = (diff_totale - diff_plus) / diff_totale

        nouveau_x = coef1 * vec_Xt[index_moins] + coef2 * vec_Xt[index_plus - 1]
        nouveau_z = coef1 * vec_Zt[index_moins] + coef2 * vec_Zt[index_plus - 1]

        insertion_x = vec_Xt.insert(index_temps, nouveau_x)
        insertion_z = vec_Zt.insert(index_temps, nouveau_z)
        insertion_temps = vec_temps.insert(index_temps, valeur_temps)



    return nouveau_x, nouveau_z, vec_Xt, vec_Zt, vec_temps


def vecteurs_temps_i(vec_Xt, vec_Zt, vec_temps, nouveau_pas):
    '''
    Function to compute the trajectory points below the front at the new_step time step.

    Args:
        vec_Xt: Vectors of X values
        vec_Zt: Vector of Z values
        vec_temps: Vector of the times corresponding to the (X,Z) coordinates
        nouveau_pas: Time considered

    Returns:
        vec_Xt_drop_haut: X coordinates of the trajectory up to the front
        vec_Zt_drop_haut: Z coordinates of the trajectory up to the front
        vec_temps_drop_haut: times of the trajectory up to the front
        vec_distance_drop_haut: vector of distances between the trajectory points up to the front
        vec_distance_parc_drop: vector of cumulative distances between the trajectory points up to the front
        distance_parcourue_i: distance traveled by the magma since the source
    '''
    #Calcul des distances entre chaque point et des distance cumulées
    vec_distance_i = (np.diff(vec_Xt) ** 2 + np.diff(vec_Zt) ** 2) ** 0.5
    vec_distance_pas_i = np.insert(vec_distance_i, 0, 0)

    vec_DC_pas_i = np.cumsum(vec_distance_pas_i)
    distance_parcourue = np.max(vec_DC_pas_i)



    # Insertion du point haut dans vec_Xt, drop ce qui est au dessus et calcul des distance cumulées jusqu'à temps_i + pas
    new_tab1 = [vec_Xt, vec_Zt, vec_temps, vec_distance_pas_i, vec_DC_pas_i]
    new_tab = np.transpose(new_tab1)
    df_new_tab = pd.DataFrame(new_tab)
    df_new_tab_drop_haut = df_new_tab.drop(df_new_tab[(df_new_tab[2] > nouveau_pas)].index)


    vec_new_tab_drop_haut = df_new_tab_drop_haut.to_numpy()
    vec_tab_drop_haut = np.transpose(vec_new_tab_drop_haut)

    vec_Xt_drop_haut = vec_tab_drop_haut[0]
    vec_Zt_drop_haut = vec_tab_drop_haut[1]
    vec_temps_drop_haut = vec_tab_drop_haut[2]
    vec_distance_drop_haut = vec_tab_drop_haut[3]
    vec_distance_parc_drop = vec_tab_drop_haut[4]

    distance_parcourue_i = vec_distance_parc_drop[-1]

    return vec_Xt_drop_haut, vec_Zt_drop_haut, vec_temps_drop_haut, vec_distance_drop_haut, vec_distance_parc_drop, distance_parcourue_i




def point_bas(x, z, vec_X, vec_Z, vec_D, vec_DC, L):
    '''
    Computes the minimum point of a trajectory based on the (x,z) coordinates of the front and the length L of the
    magmatic ascent.

    Args:
        x: X coordinate of the front.
        z: Z coordinate of the front
        vec_X: Vector of the trajectory's X coordinates
        vec_Z: Vector of the trajectory's Z coordinates
        vec_D: Vector of distances between two points of the trajectory
        vec_DC: Vector of the trajectory's cumulative distances
        L: Length of the magmatic ascent

    Returns:
        (x_b, z_b): lower coordinates of the magmatic ascent.
    '''


    if vec_DC[-1] <= L:
        x_b = vec_X[0]
        z_b = vec_Z[0]

    else:
        k = vec_DC[-1] - L    # k différence entre Lon cumulée et L. Distance cumulée " basse " DCB

        #On intercale DCB sur la trajectoire vec_X vec_Z.
        #Pour cela on la place déjà sur vec_DC le vecteur des trajectoires cumulées
        index_k = bisect.bisect(vec_DC, k)

        indice_moins = index_k - 1
        indice_plus = index_k

        reste_a_parcourir = vec_DC[indice_plus] - k

        coef = reste_a_parcourir / vec_D[indice_plus]

        x_b = coef * vec_X[indice_moins] + (1 - coef) * vec_X[indice_plus]
        z_b = coef * vec_Z[indice_moins] + (1 - coef) * vec_Z[indice_plus]


    return x_b,z_b


def une_trajectoire_haute_basse(vec_X_haut, vec_Z_haut, vec_X_bas, vec_Z_bas, x):
    '''
    Extracts the coordinates of a magmatic ascent from the coordinates of its front and its lowest point.

    Args:
        vec_X_haut: X coordinates of the trajectory - of the front.
        vec_Z_haut: Z coordinates of the trajectory - of the front.
        vec_X_bas: X coordinate of the bottom of the ascent.
        vec_Z_bas: Z coordinate of the bottom of the ascent.
        x: X coordinate of the considered front.

    Returns:
        vecX_coupe, vecZ_coupe: X and Z coordinates of the ascent, between the front and the bottom.
    '''
    # Indice de la valeur Xhaute que l'on considère.
    indice_X_haut = vec_X_haut.index(x)
    #Recupération des valeurs Xbasse et Zbasse correspondantes


    X_bas = vec_X_bas
    Z_bas = vec_Z_bas

    #Attribution d'une place pour la valeur Xbasse dans le vecteur trajectoire

    if x >=0:

        indice_X_bas = bisect.bisect(vec_X_haut, X_bas)
        #création d'un vecteur avec la valeur Xbasse considérée et les autres valeurs Xhautes
        vec_X_HB_ar = np.insert(vec_X_haut, indice_X_bas, X_bas)
        vec_Z_HB_ar = np.insert(vec_Z_haut, indice_X_bas, Z_bas)

        vec_X_HB = list(vec_X_HB_ar)
        vec_Z_HB = list(vec_Z_HB_ar)

        vecX_coupe = []
        vec_indice_X_coupe = []
        for valX in vec_X_HB:
            if valX >= X_bas and valX <= x:
                vecX_coupe.append(valX)
                indiceX = vec_X_HB.index(valX)
                vec_indice_X_coupe.append(indiceX)


        vecZ_coupe = []
        for indice in vec_indice_X_coupe:
            val = vec_Z_HB[indice]
            vecZ_coupe.append(val)


    else:
        ##### Pour x négatif
        vec_X_haut_liste = list(vec_X_haut)

        vec_X_haut_croissant = sorted(vec_X_haut_liste)

        vec_X_haut_croissant_tri = np.array(vec_X_haut_croissant)

        indice_X_bas_dans_X_haut = bisect.bisect(vec_X_haut_croissant_tri, X_bas)

        vec_X_haut_croissant_tri_ar = np.insert(vec_X_haut_croissant_tri, indice_X_bas_dans_X_haut, X_bas)

        vec_X_HB_croissant = list(vec_X_haut_croissant_tri_ar)

        vec_X_HB = sorted(vec_X_HB_croissant, reverse=True)

        indice_X_HB = len(vec_X_HB) - indice_X_bas_dans_X_haut - 1
        vec_Z_HB_ar = np.insert(vec_Z_haut, indice_X_HB, Z_bas)
        vec_Z_HB = list(vec_Z_HB_ar)



        vecX_coupe = []
        vec_indice_X_coupe = []
        for valX in vec_X_HB:
            if valX >= x and valX <= X_bas:
                vecX_coupe.append(valX)
                indiceX = vec_X_HB.index(valX)
                vec_indice_X_coupe.append(indiceX)

        vecZ_coupe = []
        for indice in vec_indice_X_coupe:
            val = vec_Z_HB[indice]
            vecZ_coupe.append(val)


    return vecX_coupe, vecZ_coupe




def une_trajectoire_pas_i(vec_Xt_eff, vec_Zt, vec_temps_eff, temps_courant, pas_temps, longueur):
    '''

    Args:
        vec_Xt_eff: vecteur des coordonnées X trajectoire complète
        vec_Zt: vecteur coordonnées Z trajectoire complète
        vec_temps_eff: vecteur du temps mis pas le dike pour atteindre (X,Z)
        temps_courant: temps au pas i
        pas_temps: pas temps
        longueur: longueur maximale du dike

    Returns:
        trajectoire_x_i : vecteur coordonnées X du dike au moment i
        trajectoire_z_i : vecteur coordonnées X du dike au moment i
        X_haut : coordonnées X du front du dike
        Z_haut : coordonnées Z du front du dike
        X_bas :  coordonnées X du bas du dike
        Z_bas : coordonnées Z du bas du dike
        distance_parcourue_i : distance parcourue par le dike depuis le point de départ


    '''
    nouveau_pas = temps_courant + pas_temps

    X_haut, Z_haut, new_vec_Xt, new_vec_Zt, new_vec_temps = nouvelle_valeur(vec_Xt_eff, vec_Zt, vec_temps_eff,
                                                                            nouveau_pas)


    vec_Xt_drop_haut, vec_Zt_drop_haut, vec_temps_drop_haut, vec_distance_drop_haut, vec_distance_parc_drop, distance_parcourue_i = vecteurs_temps_i(
        new_vec_Xt, new_vec_Zt, new_vec_temps, nouveau_pas)

    ## Maintenant calcul du point bas à au temps considéré

    X_bas, Z_bas = point_bas(X_haut, 0, vec_Xt_drop_haut, vec_Zt_drop_haut, vec_distance_drop_haut,
                             vec_distance_parc_drop, longueur)

    ## Maintenant trajectoire : tous les points de la trajectoire au moment voulu. Il faut intercaler xbas zbas

    vecXt = list(vec_Xt_drop_haut)
    vecZt = list(vec_Zt_drop_haut)

    trajectoire_x_i, trajectoire_z_i = une_trajectoire_haute_basse(vecXt, vecZt, X_bas, Z_bas, X_haut)

    return vecXt, vecZt, trajectoire_x_i, trajectoire_z_i, X_haut, Z_haut, X_bas, Z_bas, distance_parcourue_i
