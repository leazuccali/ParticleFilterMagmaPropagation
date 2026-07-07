import os
import numpy as np
import json

def fonction_watanabe(x_min, x_max, z_min, z_max, pas, P_load, rayon_load, extension):
    '''
    Computes the stress fields along the xx, zz, and xz axes according to the equations given by Watanabe.
    :param x_min: minimum abscissa
    :param x_max: maximum abscissa
    :param z_min: minimum ordinate
    :param z_max: maximum ordinate
    :param pas: distance between each abscissa/ordinate point
    :param P_load: load applied on the surface, symmetric with respect to the ordinate axis
    :param rayon_load: radius of the applied load
    :param extension: value of the caldera extension
    :return: mesh_X, mesh_Z, and the stress fields along the xx, zz, and xz axes
    '''

    vec_X = np.arange(x_min, x_max, pas)
    vec_Z = np.arange(z_min, z_max, pas)
    mesh_X, mesh_Z = np.meshgrid(vec_X,vec_Z)


    # plt.figure(figsize=(6, 6))
    # plt.scatter(mesh_X, mesh_Z, color='red')  # Points de la grille
    # for i in range(len(vec_X)):
    #     for j in range(len(vec_Z)):
    #         plt.text(mesh_X[j, i], mesh_Z[j, i], f"({mesh_X[j, i]}, {mesh_Z[j, i]})", fontsize=10, ha='center')
    # plt.title("Grille générée par np.meshgrid")
    # plt.xlabel("X-axis")
    # plt.ylabel("Y-axis")
    # plt.grid()
    # plt.show()

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
    Computes the maximum and minimum stress fields Sigma1 and Sigma3 from Sigma_xx, Sigma_zz,
    Sigma_xz. Displays these vector fields.
    :param mesh_X: meshgrid from the xmin-xmax vector divided into a certain number of steps
    :param mesh_Z: meshgrid from the ymin-ymax vector divided into a certain number of steps
    :param sig_xx: stress field along the xx axis
    :param sig_zz: stress field along the zz axis
    :param sig_xz: stress field along the xz axis
    :param pas_vec: spacing of the vectors displayed for sigma1
    :return: sig_1, sig_3, u_sig_1, v_sig_1, u_sig_3, v_sig_3, I (the Sigma_1, Sigma_3 values, and the coordinates of the associated vectors).
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
    v_sig_1_eff = (1 - u_sig_1_eff**2)**(1/2)

    mesh_vec = ((mesh_X%pas_vec == 0) & (mesh_Z%pas_vec == 0))

    return sig_1, sig_3, u_sig_1, v_sig_1, u_sig_1_eff, v_sig_1_eff, u_sig_3, v_sig_3, mesh_vec






if __name__ == '__main__':

    directory = "grilles_R_G"
    os.makedirs(directory, exist_ok=True)

    ### Global data/parameters - Units : step time (s), Pload (Pa), radius (m), step vectors (m)
    global_data = {'step time': 60, 'p_load': 15000000, 'radius load': 10000, 'step vectors': 400,
                   'nb free propag': 60}

    ### Grid data/parameters - Units : m
    grid_data = {'xmin': -30000, 'xmax': 30000, 'zmin': -15000, 'zmax': -1, 'discretisation step': 100}

    # Parameters
    xmin = grid_data['xmin']
    xmax = grid_data['xmax']
    zmin = grid_data['zmin']   #  -11000  # m
    zmax = grid_data['zmax']
    pas_trajectoire = grid_data['discretisation step']

    pas_vect = global_data['step vectors']  #2000
    P_load = global_data['p_load']
    rayon_load = global_data['radius load']

    ## R/G values discretisation
    vec_R = np.arange(0, 1+0.01, 0.05)
    vec_G = np.arange(0.05, 1, 0.05)

    ## Grids creation
    num_dico = 0

    for R in vec_R:
        for G in vec_G:
            #num_dico += 1
            print("Grid number ", num_dico, "(R,G) = ",R,G)
            # Creation d'un nouveau dictionnaire
            dico_i = {}
            dico_i['Numero'] = num_dico
            dico_i['R'] = R
            dico_i['G'] = G
            #Creation de la grille mesh_X mesh_Z
            extension = abs(P_load) * R
            mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz = fonction_watanabe(xmin, xmax, zmin, zmax, pas_trajectoire,
                                                                       P_load, rayon_load, extension)

            # Creation des vecteurs de direction principaux en chaque point de la grille
            sig_1, sig_3, u_sig_1, v_sig_1, u_sig_1_eff, v_sig_1_eff, u_sig_3, v_sig_3, mesh_vec = calcul_sig1_sig3(mesh_X,
                                                                    mesh_Z, sig_xx, sig_zz, sig_xz, pas_vect, G)


            # Creation des dip correspondants aux orientations des vecteurs principaux en chaque point de la grille

            val_tan = v_sig_1_eff / u_sig_1_eff
            vec_angles = np.arctan(abs(val_tan))

            liste_angles = vec_angles.tolist()
            dico_i['Angles'] = liste_angles

            data_filename = os.path.join(directory, f'grille_{num_dico}.json')
            with open(data_filename, 'w') as f:
                 json.dump(dico_i, f)

            print("-> grid",num_dico, "created")

            num_dico += 1