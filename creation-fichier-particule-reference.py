import numpy as np
import json
from trajectoire_particules import *
import os
def creation_trajectoire_reference():

    #### Paramètres de la grille
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

    temps_courant = 0  # Initialisation au temps 0
    pas_temps = 60


    ## Direction du stockage

    directory = "data_initialisation"
    os.makedirs(directory, exist_ok=True)


    # ## Paramètres de référence
    x_0_ref = 2000
    z_0_ref = -8000
    x_ref = x_0_ref
    z_ref = z_0_ref
    R_ref = 0.5
    G_ref = 0.5
    mu_ref = 100
    E_ref = 5000000000
    vol_ref = 100000000

    # Particule hors gamme temps long
    # x_0_ref = 18000
    # z_0_ref = -6000
    # x_ref = x_0_ref
    # z_ref = z_0_ref
    # R_ref = 1
    # G_ref = 1
    # mu_ref = 75
    # E_ref = 1e8
    # vol_ref = 1e7

    #Initialisation des dictionnaires

    trajectoire_reference = {}


    # dico_X_ref = {}
    # dico_Z_ref = {}
    # dico_tempseff_ref = {}
    # dico_longueur_ref = {}
    # dico_ouverture_ref = {}
    #
    # #Constitution des dictionnaires
    # dico_X_ref['Numero'] = 'ref'
    # dico_Z_ref['Numero'] = 'ref'
    # dico_tempseff_ref['Numero'] = 'ref'
    # dico_longueur_ref['Numero'] = 'ref'
    # dico_ouverture_ref['Numero'] = 'ref'

    trajectoire_reference['Numero'] = 'ref'



    #Lancement de la fonction
    vec_Xt_ref, vec_Zt_ref, vec_temps_ref, temps_courant, pas_temps, longueur_ref, ouverture_ref, vitesse_ref = trajectoire_une_particule(
        np.array([[x_0_ref, z_0_ref]]), R_ref, G_ref, mu_ref, E_ref, vol_ref, temps_courant, pas_temps,
        xmin, xmax, zmin, zmax, pas_trajectoire, P_load, rayon_load, pas_vect, norme)


    # dico_X_ref['X'] = vec_Xt_ref
    # dico_Z_ref['Z'] = vec_Zt_ref
    # dico_tempseff_ref['Temps effectif'] = vec_temps_ref
    # dico_longueur_ref['Longueur remontee'] = longueur_ref
    # dico_ouverture_ref['Ouverture'] = ouverture_ref
    #
    trajectoire_reference['X'] = vec_Xt_ref
    trajectoire_reference['Z'] = vec_Zt_ref
    trajectoire_reference['Temps effectif'] = vec_temps_ref
    trajectoire_reference['Longueur remontee'] = longueur_ref
    trajectoire_reference['Ouverture'] = ouverture_ref
    trajectoire_reference['Vitesse'] = vitesse_ref



    # with open('sauvegarde_trajectoires_X_ref.json', 'w') as fichier :
    #     json.dump(dico_X_ref, fichier, indent = 4)
    #
    # with open('sauvegarde_trajectoires_Z_ref.json', 'w') as fichier :
    #     json.dump(dico_Z_ref, fichier, indent = 4)
    #
    # with open('sauvegarde_trajectoires_tempseff_ref.json', 'w') as fichier :
    #     json.dump(dico_tempseff_ref, fichier, indent = 4)
    #
    # with open('sauvegarde_trajectoires_longueurs_ref.json', 'w') as fichier :
    #     json.dump(dico_longueur_ref, fichier, indent = 4)
    #
    # with open('sauvegarde_trajectoires_ouvertures_ref.json', 'w') as fichier :
    #     json.dump(dico_ouverture_ref, fichier, indent = 4)

    data_filename = os.path.join(directory, f'sauvegarde_trajectoire_ref.json')
    with open(data_filename, 'w') as fichier :
        json.dump(trajectoire_reference, fichier, indent = 4)


if __name__ == '__main__':
    creation_trajectoire_reference()

