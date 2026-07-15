import os
import json
import itertools
import numpy as np
from joblib import Parallel, delayed


def fonction_watanabe(x_min, x_max, z_min, z_max, pas, P_load, rayon_load, extension):
    '''
    Calcule les champs de contrainte selon les axes xx, zz et xz selon les équations données par watanabe.
    ...
    '''
    vec_X = np.arange(x_min, x_max, pas)
    vec_Z = np.arange(z_min, z_max, pas)
    mesh_X, mesh_Z = np.meshgrid(vec_X, vec_Z)

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
    Calcule les champs de contrainte maximum et minimum Sigma1 et Sigma3.
    ...
    '''
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

    u_sig_1_eff = (1 - G) * u_sig_1
    v_sig_1_eff = (1 - u_sig_1_eff ** 2) ** (1 / 2)

    mesh_vec = ((mesh_X % pas_vec == 0) & (mesh_Z % pas_vec == 0))

    return sig_1, sig_3, u_sig_1, v_sig_1, u_sig_1_eff, v_sig_1_eff, u_sig_3, v_sig_3, mesh_vec


def traiter_grille(num_dico, R, G, xmin, xmax, zmin, zmax, pas_trajectoire, pas_vect,
                    P_load, rayon_load, directory):
    '''
    Traite une combinaison (R, G) : calcule le champ de contrainte, les angles,
    et sauvegarde le résultat dans un fichier JSON. Fonction isolée pour permettre
    la parallélisation avec joblib (chaque appel est indépendant).
    '''
    dico_i = {}
    dico_i['Numero'] = num_dico
    dico_i['R'] = R
    dico_i['G'] = G

    extension = abs(P_load) * R
    mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz = fonction_watanabe(
        xmin, xmax, zmin, zmax, pas_trajectoire, P_load, rayon_load, extension
    )

    sig_1, sig_3, u_sig_1, v_sig_1, u_sig_1_eff, v_sig_1_eff, u_sig_3, v_sig_3, mesh_vec = calcul_sig1_sig3(
        mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz, pas_vect, G
    )

    #val_tan = v_sig_1_eff / u_sig_1_eff
    #vec_angles = np.arctan(abs(val_tan))
    vec_angles = np.arctan2(np.abs(v_sig_1_eff), np.abs(u_sig_1_eff))

    dico_i['Angles'] = vec_angles.tolist()

    data_filename = os.path.join(directory, f'grille_{num_dico}.json')
    with open(data_filename, 'w') as f:
        json.dump(dico_i, f)

    print('grille', num_dico, 'ok')
    return num_dico


if __name__ == '__main__':

    directory = "grilles_R_G"
    os.makedirs(directory, exist_ok=True)


    global_data = {'step time': 60, 'p_load': -15000000, 'radius load': 10000, 'step vectors': 400,
                   'nb free propag': 60, 'grid RG nb': 399}

    ### Grid data/parameters - Units : m
    grid_data = {'xmin': -30000, 'xmax': 30000, 'zmin': -15000, 'zmax': -1, 'discretisation step': 100}

    xmin = grid_data['xmin']
    xmax = grid_data['xmax']
    zmin = grid_data['zmin']
    zmax = grid_data['zmax']
    pas_trajectoire = grid_data['discretisation step']
    pas_vect = global_data['step vectors']
    P_load = global_data['p_load']
    rayon_load = global_data['radius load']
    norme = 1
    #pas_okada = 100
    # temps_i = 60






    vec_R = np.arange(0, 1 + 0.01, 0.05)
    vec_G = np.arange(0.05, 1, 0.05)

    # Génère la liste (num_dico, R, G) dans le même ordre que la double boucle originale
    combinaisons = [(num_dico, R, G) for num_dico, (R, G) in
                     enumerate(itertools.product(vec_R, vec_G))]

    Parallel(n_jobs=-1)(
        delayed(traiter_grille)(num_dico, R, G, xmin, xmax, zmin, zmax,
                                 pas_trajectoire, pas_vect, P_load, rayon_load, directory)
        for num_dico, R, G in combinaisons
    )