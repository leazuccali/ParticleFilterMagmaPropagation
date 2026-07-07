import numpy as np
#from pymc3 import logsumexp
from scipy.special import logsumexp

def calcul_vraisemblances(Nb_particules, ux_ref, uz_ref, deplacements_ux, deplacements_uz, deplacements_mY, incertitudes, vec_index_p_surface):
    '''

    Args:
        Nb_particules:
        ux_ref: UX displacements of the reference particle
        uz_ref: UZ displacements of the reference particle
        deplacements_ux: list of UX displacements of the particles
        deplacements_uz: list of UZ displacements of the particles
        deplacements_mY:
        incertitudes: uncertainties
        vec_index_p_surface: index of the particles that reached the surface

    Returns:

    '''
    vraisemblances_particules = []

    coefficients_particules = []

    coefficients_normalise = []

    vraisemblances_normalisee = []

    #erreurs_particules = []
    Nb_points_x = len(ux_ref)
    Nb_points_z = len(uz_ref)

    #### Construction de la matrice de covariance
    erreur_ux = incertitudes[0]
    erreur_uz = incertitudes[1]

    diago = [erreur_ux]*Nb_points_x + [erreur_uz]*Nb_points_z
    mat_coefs = np.diag(diago)
    mat_cov = np.dot(mat_coefs, mat_coefs)

    cov_inv = np.linalg.inv(mat_cov)  # Inverse de la matrice de covariance
    cov_det = np.linalg.det(mat_cov)  # Déterminant pour la normalisation

    ### Calcul de la vraisemblance de l'hypothèse nulle
    ux_ref_ar_nulle = np.array(ux_ref)
    uz_ref_ar_nulle = np.array(uz_ref)
    vec_dep_ref_nulle = np.concatenate((ux_ref_ar_nulle, uz_ref_ar_nulle))

    erreur_nulle1 = np.dot(vec_dep_ref_nulle, cov_inv)
    erreur_nulle2 = np.dot(erreur_nulle1, vec_dep_ref_nulle.T)

    ## Vraisemblance de l'erreur classique
    vraisemb_nulle = np.exp(-0.5 * erreur_nulle2)

    ## Vecteur du bruit de l'observation. Ici un point = un bruit.
    vec_bruit_ux = np.random.normal(0, erreur_ux, Nb_points_x)
    vec_bruit_uz = np.random.normal(0, erreur_uz, Nb_points_z)


    #coefficient = round(1 / ((2 * np.pi) ** (Nb_points) * cov_det ** (Nb_points / 2)), 5)
    #coefficient = 1 / (np.sqrt((2* np.pi)**2  * cov_det))
    for p in range(Nb_particules):
        if p in vec_index_p_surface:
            vraisemb = 0
            #vraisemb = np.array([0])

            vraisemblances_particules.append(vraisemb)

            ## calcul de la min-max vraisemblance
            coefficients_particules.append('S')
            ## coef

            coefficients_normalise.append(10000)
            vraisemblances_normalisee.append(0)

        else:
            ind_0 = np.nonzero(deplacements_mY[p] == 0)
            dep_x = deplacements_ux[p]
            dep_x_p = dep_x[ind_0]
            dep_z = deplacements_uz[p]
            dep_z_p = dep_z[ind_0]


            #Gaussienne normalisee carree
            ux_ref_ar = np.array(ux_ref) + vec_bruit_ux   #obs bruitée
            uz_ref_ar = np.array(uz_ref) + vec_bruit_uz   #obs bruitéé
            vec_dep_ref = np.concatenate((ux_ref_ar, uz_ref_ar))

            dep_x_ar = np.array(dep_x_p)
            dep_z_ar = np.array(dep_z_p)
            vec_dep_pred = np.concatenate((dep_x_ar, dep_z_ar))

            #erreur = (vec_dep_ref - vec_dep_pred)/vec_dep_ref   ##########erreur normalisee
            erreur = vec_dep_ref - vec_dep_pred
            erreur_ponderee = np.dot(erreur, cov_inv)
            erreur_carree = np.dot(erreur_ponderee, erreur.T)

            ## Coefficient de l'erreur classique
            coefficients_particules.append(0.5 * erreur_carree)
            ## Vraisemblance de l'erreur classique
            vraisemb = np.exp(-0.5 * erreur_carree)
            vraisemblances_particules.append(vraisemb)

            ###### coefficient normalise
            erreur_normalise1 = erreur / vec_dep_ref
            erreur_normalise2 = np.dot(erreur_normalise1, cov_inv)
            erreur_normalise3 = np.dot(erreur_normalise2, erreur_normalise1.T)
            coefficients_normalise.append(0.5* erreur_normalise3)
            vraisemb_normal = np.exp(-0.5 * erreur_normalise3)
            vraisemblances_normalisee.append(vraisemb_normal)


    ###############Normalisation entre 0 et 3 du vecteur coefficients particules

    # Extraction des valeurs numériques
    numerical_values = np.array([x for x in coefficients_particules if x != 'S'], dtype=float)
    # Min et Max des valeurs numériques
    X_min, X_max = numerical_values.min(), numerical_values.max()
    # Min-Max Normalization
    normalized_values = 3*(numerical_values - X_min) / (X_max - X_min)
    # Création d'un nouveau vecteur avec les valeurs normalisées
    #normalized_vec = np.array(coefficients_particules, dtype=object)
    normalized_vec = np.array(coefficients_particules, dtype=object)
    num_idx = 0  # Index pour parcourir les valeurs numériques normalisées

    for i in range(len(coefficients_particules)):
        if coefficients_particules[i] != 'S':
            #normalized_vec.append(normalized_values[i])
            normalized_vec[i] = normalized_values[num_idx]  # Remplace par la valeur normalisée
            num_idx += 1
        else:
            #normalized_vec.append(3)
            normalized_vec[i] = 3  # Remplace 'S' par 2

    # Conversion explicite en float pour éviter l'erreur
    normalized_vec = normalized_vec.astype(float)

    #vec_norm = np.array(normalized_vec)
    norm_vraisemblance = np.exp(-normalized_vec)

    return norm_vraisemblance, coefficients_particules, vraisemblances_particules, vraisemblances_normalisee, coefficients_normalise, vraisemb_nulle





def fonction_poids(Nb_particules, ux_ref, uz_ref, deplacements_ux, deplacements_uz, deplacements_mY, matrice_covariance, vec_index_p_surface):
    '''

    Args:
        Nb_particules : nombre de particules
        ux_ref : vecteur déplacements ux observés
        uz_ref : vecteur déplacements uz observés
        deplacements_ux : matrice des déplacements ux des n particules
        deplacements_uz : matrice des déplacements uz des n particules
        matrice_covariance : bruit associés aux déplacements
        vec_index_p_surface : dike ayant atteint la surface

    Returns:
        vec_vraisemblances : vecteur des vraisemblances
        vecteur_poids : vecteur des poids
        coef_particules : terme ln(exp(vraisemblance))
        vraisemblances_classiques : vraisemblances
        vraisemblances_normalisee : vraisemblances normalisées par la valeur de déplacement observés
        coef_normalises : terme ln(exp(vraisemblance)) normalisés
    '''

    vec_vraisemblances, coef_particules, vraisemblances_classiques, vraisemblances_normalisee, coef_normalises, vraisemb_hpnulle = calcul_vraisemblances(Nb_particules, ux_ref, uz_ref, deplacements_ux, deplacements_uz, deplacements_mY, matrice_covariance, vec_index_p_surface)

    vecteur_poids = vec_vraisemblances/np.sum(vec_vraisemblances)



    return vec_vraisemblances, vecteur_poids, coef_particules, vraisemblances_classiques, vraisemblances_normalisee, coef_normalises, vraisemb_hpnulle

