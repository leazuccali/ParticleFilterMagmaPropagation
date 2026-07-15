import numpy as np
import json
from trajectoire_particules import *
import os

def creation_trajectoire_reference(global_data, grid_data,p_reference_data, current_time):
    '''
    Creation of truth particle trajectory

    Args :
        global_data,
        grid_data,
        p_reference_data
        current_time
    Return :
    sauvegarde_trajectoire_ref.json in data_initialisation folder
    '''


    #### Paramètres de la grille
    xmin = grid_data['xmin']  # m
    xmax = grid_data['xmax']  # m
    zmin = grid_data['zmin']  # m
    zmax = grid_data['zmax']  # m
    pas_trajectoire = grid_data['discretisation step']   ##
    pas_vect = global_data['step vectors']
    P_load = global_data['p_load']  # MPa
    rayon_load = global_data['radius load']  # m
    norme = 1  # Si norme = 1, alors l'unité est le mètre. Si norme = 1000, alors l'unité est le km.
    #pas_okada = 100

    current_time = 0  # Initialisation au temps 0
    pas_temps = global_data['step time']


    ## Direction du stockage

    directory = "data_initialisation"
    os.makedirs(directory, exist_ok=True)


    # # ## Paramètres de référence
    x_0_ref = p_reference_data['x0_ref']
    z_0_ref = p_reference_data['z0_ref']
    x_ref = x_0_ref
    z_ref = z_0_ref
    R_ref = p_reference_data['R_ref']
    G_ref = p_reference_data['G_ref']
    mu_ref = p_reference_data['mu_ref']
    E_ref = p_reference_data['E_ref']
    vol_ref = p_reference_data['V_ref']


    trajectoire_reference = {}

    trajectoire_reference['Numero'] = 'ref'

    #Lancement de la fonction
    vec_Xt_ref, vec_Zt_ref, vec_temps_ref, current_time, pas_temps, longueur_ref, ouverture_ref, vitesse_ref = trajectoire_une_particule(
        np.array([[x_0_ref, z_0_ref]]), R_ref, G_ref, mu_ref, E_ref, vol_ref, current_time, pas_temps,
        xmin, xmax, zmin, zmax, pas_trajectoire, P_load, rayon_load, pas_vect, norme)


    trajectoire_reference['X'] = vec_Xt_ref
    trajectoire_reference['Z'] = vec_Zt_ref
    trajectoire_reference['Temps effectif'] = vec_temps_ref
    trajectoire_reference['Longueur remontee'] = longueur_ref
    trajectoire_reference['Ouverture'] = ouverture_ref
    trajectoire_reference['Vitesse'] = vitesse_ref




    data_filename = os.path.join(directory, f'sauvegarde_trajectoire_ref.json')
    with open(data_filename, 'w') as fichier :
        json.dump(trajectoire_reference, fichier, indent = 4)


if __name__ == '__main__':
    ### Global data/parameters - Units : step time (s), Pload (Pa), radius (m), step vectors (m)
    global_data = {'step time': 60, 'p_load': -15000000, 'radius load': 10000, 'step vectors': 400,
                   'nb free propag': 60, 'grid RG nb': 399}

    ### Grid data/parameters - Units : m
    grid_data = {'xmin': -30000, 'xmax': 30000, 'zmin': -15000, 'zmax': -1, 'discretisation step': 100}

    ### Particle reference data/parameters - Units : x0_ref, z0_ref (m) ; R,G (-) ; mu (Pa.s) ; E (Pa) ; Vol (m^3)
    p_reference_data = {'x0_ref': 2000, 'z0_ref': -8000, 'R_ref': 0.5, 'G_ref': 0.5, 'mu_ref': 100,
                        'E_ref': 5 * 10 ** 9,
                        'V_ref': 10 ** 8, 'Step max': 550}

    current_time = 0
    creation_trajectoire_reference(global_data, grid_data,p_reference_data, current_time)

