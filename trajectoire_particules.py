import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


def calcul_parametres(R,G,mu,E,vol, P_load):
    '''
    Computation of the extension, magmatic crack length, and magma velocity parameters from the R and G ratios,
    the magma viscosity, the crust elasticity, and the magma volume.

    Args :
        R
        G
        viscosite
        rigidite
        volume

    Return:
        extension
        longueur   (length)
        vitesse   (velocity)
        ouverture  (opening)
    '''

    #Ratio unloading décharge/extension
    extension = abs(P_load) * R

    #length
    longueur3 = (E * vol) / ((1 - 0.25**2)*G*abs(P_load))
    longueur = longueur3**(1/3)

    #Constant C, velocity
    C = 5 * 10**(-7)
    vitesse = (C * vol * G)/mu

    #Opening
    ouverture = (math.pi * vol) / ( 2 * longueur**2)

    return extension, longueur, vitesse, ouverture


def fonction_watanabe(x_min, x_max, z_min, z_max, pas, P_load, rayon_load, extension):
    '''
    Computes the stress fields along the xx, zz, and xz axes according to the equations given by Watanabe.

    Args :
        x_min: minimum abscissa
        x_max: maximum abscissa
        z_min: minimum ordinate
        z_max: maximum ordinate
        pas: distance between each abscissa/ordinate point
        P_load: load applied on the surface, symmetric with respect to the ordinate axis
        rayon_load: radius of the applied load
        extension: value of the caldera extension

    Return:
        mesh_X, mesh_Z, and the stress fields along the xx, zz, and xz axes
    '''

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
    Computes the maximum and minimum stress fields Sigma1 and Sigma3 from Sigma_xx, Sigma_zz,
    Sigma_xz. Displays these vector fields.

    Args :
        mesh_X: meshgrid from the xmin-xmax vector divided into a certain number of steps
        mesh_Z: meshgrid from the ymin-ymax vector divided into a certain number of steps
        sig_xx: stress field along the xx axis
        sig_zz: stress field along the zz axis
        sig_xz: stress field along the xz axis
        pas_vec: spacing of the vectors displayed for sigma1

    Return:
        sig_1, sig_3, u_sig_1, v_sig_1, u_sig_3, v_sig_3, I (the Sigma_1, Sigma_3 values, and the coordinates of the associated vectors).
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


def streamplot_trajectoire(mesh_X, mesh_Z, u_sig_1, v_sig_1, start_point, R, G, vitesse,temps_i):
    '''
    Displays the magma trajectory along the maximum stress field sigma_1.
    Extracts the coordinates of the trajectory points and computes the maximum length of the trajectory.

    Args :
        mesh_X: meshgrid from the xmin-xmax vector divided into a certain number of steps
        mesh_Z: meshgrid from the ymin-ymax vector divided into a certain number of steps
        u_sig_1: u coordinate of the maximum stress field
        v_sig_1: v coordinate of the maximum stress field
        start_point: starting point of the considered trajectory: np.array([[2,3]]). Minimum ordinate of the display.
        R: values of the ratio between the extension and the discharge
        G: value of the ratio between the pressure exerted by the magma and the discharge
        vitesse: magma ascent velocity
        pas_temps: desired time step between two points

    Return:
        vec_Xt, vec_Zt, coordinates of the magma trajectory points.
    '''

    #mesh_X, mesh_Z, u_sig_1_eff, v_sig_1_eff, start_point, R, G, vitesse, temps_courant)
    fig1, ax1 = plt.subplots(figsize=(8,7))
    start = start_point
    mesh_X_new = mesh_X
    mesh_Z_new = mesh_Z
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




def trajectoire_une_particule(global_data, grid_data, current_time, x0, z0, R, G, mu, E, vol):
    '''
    Global computation of a particle's trajectory over a discretized grid, using the given parameters, the Watanabe
    function, and the streamplot fields.

    Args:
        global_data, grid_data, current_time : description in main
        x0,z0,R,G,mu,E,vol : particle parameters
    Returns:
        vec_Xt_eff: X coordinates of the complete trajectory
        vec_Zt: Z coordinate of the complete trajectory
        vec_temps_eff: time since the start of the experiment when (X,Z) is reached
        longueur: maximum dike length
        ouverture: dike opening
        vitesse: propagation velocity
    '''

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


    start_point = np.array([[x0, z0]])
    extension, longueur, vitesse, ouverture = calcul_parametres(R, G, mu, E, vol, P_load)


    mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz = fonction_watanabe(xmin, xmax, zmin, zmax, pas_trajectoire, P_load,
                                                               rayon_load, extension)


    sig_1, sig_3, u_sig_1, v_sig_1, u_sig_1_eff, v_sig_1_eff, u_sig_3, v_sig_3, mesh_vec = calcul_sig1_sig3(
        mesh_X, mesh_Z, sig_xx, sig_zz, sig_xz, pas_vect, G)

    vec_Xt, vec_Xt_eff, vec_Zt, long_magma, vec_temps_eff = streamplot_trajectoire(
        mesh_X, mesh_Z, u_sig_1_eff, v_sig_1_eff, start_point, R, G, vitesse, current_time)



    return vec_Xt_eff, vec_Zt, vec_temps_eff, current_time, pas_temps, longueur, ouverture, vitesse



def sec2hms(ss):
	(hh, ss)=divmod(ss, 3600)
	(mm, ss)=divmod(ss, 60)
	return (hh, mm, ss)

