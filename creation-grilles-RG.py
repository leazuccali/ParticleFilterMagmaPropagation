import os
import numpy as np
import matplotlib.pyplot as plt
import json
#from affichage_champs import *
#from trajectoire_une_particule import fonction_watanabe

#Dans cette fonction notre but est de créer n grilles qui ont chacunes des caractéristiques différentes puis
# de les stocker dans un dictionnaire
# Un dictionnaire : {'Numero', 'R', 'G', 'u_sig_1', 'v_sig_1', 'grille_dip'}


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






if __name__ == '__main__':

    directory = "grilles_R_G"
    os.makedirs(directory, exist_ok=True)

    #Caractéristiques des grilles et de la décharge (pression/rayon)
    xmin = -30000  # m
    xmax = 30000  # m
    zmin = -11000  # m
    zmax = -1  # m
    pas_trajectoire = 100
    pas_vect = 2000
    P_load = -15000000  # MPa
    rayon_load = 10000  # m
    norme = 1  # Si norme = 1, alors l'unité est le mètre. Si norme = 1000, alors l'unité est le km.
    pas_okada = 100
    temps_i = 60




    vec_R = np.arange(0, 1+0.01, 0.05)
    vec_G = np.arange(0.05, 1, 0.05)
    print('vec_R', vec_R)
    print('vec_G', vec_G)
    num_dico = 0
    # A mettre dans une boucle
    for R in vec_R:

        for G in vec_G:
            print(G)

            num_dico += 1
            print("num_dico", num_dico)
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
            #aff2 = affichage_trajectoire(mesh_X, mesh_Z, sig_3, u_sig_1_eff, v_sig_1_eff, mesh_vec, P_load, norme, R, G)

            # Creation des dip correspondants aux orientations des vecteurs principaux en chaque point de la grille

            #angles = np.arctan2(v_sig_1_eff, u_sig_1_eff)
            val_tan = v_sig_1_eff / u_sig_1_eff
            vec_angles = np.arctan(abs(val_tan))

            liste_angles = vec_angles.tolist()
            dico_i['Angles'] = liste_angles
            #print("dico_i", dico_i)

            data_filename = os.path.join(directory, f'grille_{num_dico}.json')
            with open(data_filename, 'w') as f:
                 json.dump(dico_i, f)


            # Calcul de la carte des angles
            # ngles_deg = np.degrees(angles)

            #print('angles', vec_angles)
            #print('vec angles', vec_angles)

            # plt.figure(figsize=(8, 6))
            # plt.contourf(mesh_X, mesh_Z, vec_angles, levels=20, cmap='twilight')
            # plt.colorbar(label="Angle (radians)")
            # plt.xlabel("X")
            # plt.ylabel("Z")
            # plt.title(f"Carte des angles : arctan(Z / X) pour R={R} et G={G}")
            # plt.grid()
            # plt.show()
