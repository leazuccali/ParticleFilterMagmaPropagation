import numpy as np
import json
from trajectoire_particules import *
import os

def creation_trajectoire_reference(global_data, grid_data,p_reference_data, current_time):
    ##### Récupération des paramètres
    #### Global data
    pas_temps = global_data['step time']
    P_load = global_data['p_load']
    rayon_load = global_data['radius load']
    pas_vect = global_data['step vectors']
    ### Grid Data
    xmin = grid_data['xmin']
    xmax = grid_data['xmax']
    zmin = grid_data['zmin']
    zmax = grid_data['zmax']
    pas_trajectoire = grid_data['discretisation step']
    ### Particle reference data
    x_0_ref = p_reference_data['x0_ref']
    z_0_ref = p_reference_data['z0_ref']
    R_ref = p_reference_data['R_ref']
    G_ref = p_reference_data['G_ref']
    mu_ref = p_reference_data['mu_ref']
    E_ref = p_reference_data['E_ref']
    vol_ref = p_reference_data['V_ref']


    ## Direction du stockage

    directory = "data_initialisation"
    os.makedirs(directory, exist_ok=True)


    #Initialisation des dictionnaires

    trajectoire_reference = {}


    trajectoire_reference['Numero'] = 'ref'

    #np.array([[x_0_ref, z_0_ref]]),
    #Lancement de la fonction
    #vec_Xt_ref, vec_Zt_ref, vec_temps_ref, temps_courant, pas_temps, longueur_ref, ouverture_ref, vitesse_ref = trajectoire_une_particule(
    #    np.array([[x_0_ref, z_0_ref]]), R_ref, G_ref, mu_ref, E_ref, vol_ref, current_time, pas_temps,
    #    xmin, xmax, zmin, zmax, pas_trajectoire, P_load, rayon_load, pas_vect)

    vec_Xt_ref, vec_Zt_ref, vec_temps_ref, temps_courant, pas_temps, longueur_ref, ouverture_ref, vitesse_ref = trajectoire_une_particule(
        global_data, grid_data, current_time, x_0_ref, z_0_ref, R_ref, G_ref, mu_ref, E_ref, vol_ref)


    trajectoire_reference['X'] = vec_Xt_ref
    trajectoire_reference['Z'] = vec_Zt_ref
    trajectoire_reference['Temps effectif'] = vec_temps_ref
    trajectoire_reference['Longueur remontee'] = longueur_ref
    trajectoire_reference['Ouverture'] = ouverture_ref
    trajectoire_reference['Vitesse'] = vitesse_ref



    data_filename = os.path.join(directory, f'sauvegarde_trajectoire_ref.json')
    with open(data_filename, 'w') as fichier :
        json.dump(trajectoire_reference, fichier, indent = 4)

    print("Reference particle")


if __name__ == '__main__':
    ### Global data/parameters
    global_data = {'step time': 60, 'p_load': 15000, 'radius load': 10000, 'step vectors': 400,
                   'nb free propag': 60}
    ### Grid data/parameters
    grid_data = {'xmin': -30000, 'xmax': 30000, 'zmin': -15000, 'zmax': -1, 'discretisation step': 100}
    ### Particle reference data/parameters
    p_reference_data = {'x0_ref': 2000, 'z0_ref': -8000, 'R_ref': 0.5, 'G_ref': 0.5, 'mu_ref': 100,
                        'E_ref': 5 * 10 ** 9,
                        'V_ref': 10 ** 8, 'Step max': 65}

    current_time = 0
    creation_trajectoire_reference(global_data, grid_data,p_reference_data, current_time)

