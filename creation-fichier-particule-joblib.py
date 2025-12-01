from trajectoire_particules import *
import json
import numpy as np
from joblib import Parallel, delayed

import os

def creation_particules(Nb_particules):
    # Tirage des vecteurs aléatoires
    #x_0 = np.random.randint(-15000, 15000, Nb_particules)
    #z_0 = np.random.randint(-13000, -1000, Nb_particules)
    #R_0 = np.random.uniform(0, 0.6, Nb_particules)
    # G_0 = np.random.normal(0.5, 0.1, Nb_particules)
    #G_0 = np.random.uniform(0.2, 0.8, Nb_particules)
    #mu_0 = np.random.randint(100, 1000, Nb_particules)
    #E_0 = np.random.randint(10 ** 9, 10 ** 10, Nb_particules)
    #vol_0 = np.random.randint(5 * 10 ** 6, 5 * 10 ** 8, Nb_particules)

    '''np.savetxt("fichier_parametres_X_init.txt", x_0)
    np.savetxt('fichier_parametres_Z_init.txt', z_0)
    np.savetxt('fichier_parametres_R_init.txt', R_0)
    np.savetxt('fichier_parametres_G_init.txt', G_0)
    np.savetxt('fichier_parametres_mu_init.txt', mu_0)
    np.savetxt('fichier_parametres_E_init.txt', E_0)
    np.savetxt('fichier_parametres_vol_init.txt', vol_0)'''

    ## Direction du stockage

    directory = "data_initialisation"
    os.makedirs(directory, exist_ok=True)



    sauvegarde_parametres = []

    for p in range(Nb_particules):

        directory_particule = os.path.join(directory, f'particule_{p}')
        os.makedirs(directory_particule, exist_ok=True)


        parametres_p = {}
        parametres_p['Numero'] = p
        parametres_p['x_0'] = np.random.randint(-15000, 15000)
        parametres_p['z_0'] = np.random.randint(-10000, -1000)
        parametres_p['R_0'] = np.random.uniform(0, 0.6)
        parametres_p['G_0'] = np.random.uniform(0.2, 0.8)
        parametres_p['mu_0'] = np.random.randint(100, 1000)
        parametres_p['E_0'] = np.random.randint(10 ** 9, 10 ** 10)
        parametres_p['vol_0'] = np.random.randint(5 * 10 ** 6, 5 * 10 ** 8)

        data_filename = os.path.join(directory_particule, f'param_init_particule_{p}.json')
        with open(data_filename, 'w') as f:
            json.dump(parametres_p, f)

        #print(parametres_p)

        sauvegarde_parametres.append(parametres_p)

    #x_0, z_0, R_0, G_0, mu_0, E_0, vol_0

    with open('sauvegarde_param_physiques_init.json', 'w') as fichier :
        json.dump(sauvegarde_parametres, fichier, indent = 4)


    return sauvegarde_parametres





def traiter_particule(p, sauvegarde_param_physiques, temps_courant, pas_temps, xmin, xmax, zmin, zmax, pas_trajectoire,
                      P_load, rayon_load, pas_vect, norme):
    print("Particule", p)

    #Création des repetoires et sous-repertoires
    directory = "data_initialisation"
    os.makedirs(directory, exist_ok=True)

    directory_particule = os.path.join(directory, f'particule_{p}')
    os.makedirs(directory_particule, exist_ok=True)

    # Ouverture de ce que l'on a besoin
    data_filename = os.path.join(directory_particule, f'param_init_particule_{p}.json')
    with open(data_filename, 'r') as fichier :
        param_physiques_p = json.load(fichier)

    # Ouverture des paramètres physiques de la particule p
    #param_physiques_p = sauvegarde_param_physiques[p]
    x_i = param_physiques_p['x_0']
    z_i = param_physiques_p['z_0']
    R_i = param_physiques_p['R_0']
    G_i = param_physiques_p['G_0']
    mu_i = param_physiques_p['mu_0']
    E_i = param_physiques_p['E_0']
    vol_i = param_physiques_p['vol_0']

    # Trajectoire
    vec_Xt_eff_p, vec_Zt_p, vec_temps_eff_p, temps_courant, pas_temps, longueur_p, ouverture_p, vitesse_p = trajectoire_une_particule(
        np.array([[x_i, z_i]]), R_i, G_i, mu_i, E_i, vol_i, temps_courant, pas_temps,
        xmin, xmax, zmin, zmax, pas_trajectoire, P_load, rayon_load, pas_vect, norme)

    dico_trajectoire_p = {
        'Numero': p,
        'X': vec_Xt_eff_p,
        'Z': vec_Zt_p,
        'Temps effectif': vec_temps_eff_p,
        'Longueur remontee': longueur_p,
        'Ouverture': ouverture_p,
        'Vitesse': vitesse_p
    }

    data_filename = os.path.join(directory_particule, f'trajectoire_init_particule_{p}.json')
    with open(data_filename, 'w') as fichier :
        json.dump(dico_trajectoire_p, fichier, indent = 4)

    return dico_trajectoire_p


def creation_trajectoires(Nb_particules):
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

    sauvegarde_param_trajectoires = []

    with open('sauvegarde_param_physiques_init.json', 'r') as fichier:
        sauvegarde_param_physiques = json.load(fichier)

    # Utilisation de joblib pour paralléliser la boucle
    resultats = Parallel(n_jobs=-1)(
        delayed(traiter_particule)(
            p, sauvegarde_param_physiques, temps_courant, pas_temps, xmin, xmax, zmin, zmax, pas_trajectoire, P_load,
            rayon_load, pas_vect, norme
        ) for p in range(Nb_particules)
    )

    # Collecter les résultats
    sauvegarde_param_trajectoires.extend(resultats)

    with open('sauvegarde_param_trajectoires_init.json', 'w') as fichier:
        json.dump(sauvegarde_param_trajectoires, fichier, indent=4)

# Vous devez définir la fonction trajectoire_une_particule ici ou l'importer si elle est définie ailleurs.




if __name__ == '__main__':
    Nb_particules = 100
    #x_0, z_0, R_0, G_0, mu_0, E_0, vol_0 = creation_particules(Nb_particules)
    sauve_parametres = creation_particules(Nb_particules)
    #print(sauvegarde_particules)
    creation_fichier_trajectoire_X = creation_trajectoires(Nb_particules)  #, x_0, z_0, R_0, G_0, mu_0, E_0, vol_0)
    #print(sauvegarde_X)


