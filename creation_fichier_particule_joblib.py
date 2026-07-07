from trajectoire_particules import *
import json
import numpy as np
from joblib import Parallel, delayed

import os

def creation_particules(test_data):

    Nb_particules = test_data['Nb particles']
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


    with open('sauvegarde_param_physiques_init.json', 'w') as fichier :
        json.dump(sauvegarde_parametres, fichier, indent = 4)


    return sauvegarde_parametres



#temps_courant, pas_temps, xmin, xmax, zmin, zmax, pas_trajectoire,P_load, rayon_load, pas_vect

def traiter_particule(p, sauvegarde_param_physiques, global_data, grid_data, current_time ):
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
        global_data, grid_data, current_time, x_i, z_i, R_i, G_i, mu_i, E_i, vol_i)

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


def creation_trajectoires(global_data, grid_data,test_data, current_time):

    # #### Global data
    # pas_temps = global_data['step time']
    # P_load = global_data['p_load']
    # rayon_load = global_data['radius load']
    # pas_vect = global_data['step vectors']
    # pas_trajectoire = grid_data['discretisation step']
    # ### Grid Data
    # xmin = grid_data['xmin']
    # xmax = grid_data['xmax']
    # zmin = grid_data['zmin']
    # zmax = grid_data['zmax']
    ### Test data
    Nb_particules = test_data['Nb particles']


    sauvegarde_param_trajectoires = []

    with open('sauvegarde_param_physiques_init.json', 'r') as fichier:
        sauvegarde_param_physiques = json.load(fichier)

    # Utilisation de joblib pour paralléliser la boucle
    resultats = Parallel(n_jobs=-1)(
        delayed(traiter_particule)(
            p, sauvegarde_param_physiques, global_data, grid_data, current_time
        ) for p in range(Nb_particules)
    )

    # Collecter les résultats
    sauvegarde_param_trajectoires.extend(resultats)

    with open('sauvegarde_param_trajectoires_init.json', 'w') as fichier:
        json.dump(sauvegarde_param_trajectoires, fichier, indent=4)




if __name__ == '__main__':

    ### Global data/parameters - Units : step time (s), Pload (Pa), radius (m), step vectors (m)
    global_data = {'step time': 60, 'p_load': -15000000, 'radius load': 10000, 'step vectors': 400,
                   'nb free propag': 60}

    ### Grid data/parameters - Units : m
    grid_data = {'xmin': -30000, 'xmax': 30000, 'zmin': -15000, 'zmax': -1, 'discretisation step': 100}

    test_data = {'Nb particles': 100, 'Assim. window': 5, 'Selection type': 'systematic', 'Observation type': 'regulier',
                 'Observation step (regular case)': 100, 'Nb observations (other cases)': 0,
                 'dossier output': "/home/zuccalil/WS1-NAS-colddata/zuccalil/output_30.06.25_systematic_100p_fenetre5_600obsregulieres_particulereference_30pourcent"}

    #current_time=0

    sauve_parametres = creation_particules(test_data)

    creation_fichier_trajectoire_X = creation_trajectoires(global_data, grid_data, test_data, 0)

