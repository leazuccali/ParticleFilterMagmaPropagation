import numpy as np

def fonction_selection_systematique(Nb_particules, poids_i):
    '''
    Systematic resampling algorithm

    Args:
        Nb_particules: number of particles
        poids_i: vector of weights
    Returns:
        list of selected particles
    '''


    L1 = [x for x in range(Nb_particules)]
    #print('L1', L1)


    # Définition distance entre les pointeurs
    distance_pointeur = 1/Nb_particules
    # Point de départ du pointeur
    pt_depart = np.random.uniform(0, distance_pointeur)
    #print("Le nombre de depart est : ", pt_depart)
    # Positions des pointeurs
    positions_pointeurs = [pt_depart + i * distance_pointeur for i in range(Nb_particules)]

    # Poids cumulés
    poids_cumules = np.cumsum(poids_i)

    ################################################################################################
    # Selection
    ################################################################################################
    vec_index_p_selec = []
    pointeur_index = 0

    for index_particule, pds_cum in enumerate(poids_cumules):
        #print('ind particule', index_particule)
        while pointeur_index < Nb_particules and positions_pointeurs[pointeur_index] < pds_cum:
            #print(vec_index_p_selec)
            vec_index_p_selec.append(index_particule)
            pointeur_index += 1
            if pointeur_index >= Nb_particules:
                break
    #print('L2', vec_index_p_selec)


    # Étape 1 : Créer L3
    L3 = ['_' for _ in L1]  # Liste initiale avec des '_'
    # Dictionnaire pour compter les occurrences de chaque élément
    count_L2 = {x: vec_index_p_selec.count(x) for x in vec_index_p_selec}
    # Remplir L3 avec les valeurs uniques de L2
    for num in L1:
        if count_L2.get(num, 0) > 0:
            L3[L1.index(num)] = num
            count_L2[num] -= 1
    #print("L3:", L3)  # Affiche L3


    # Étape 2 : Compléter L3 pour obtenir L4
    L4 = L3[:]  # Copier L3 pour créer L4
    # Liste des éléments restants de L2 (ceux en trop)
    remaining_elements = [num for num, count in count_L2.items() for _ in range(count)]
    # Remplacer les '_' par les éléments restants
    for i in range(len(L4)):
        if L4[i] == '_':
            L4[i] = remaining_elements.pop(0)
    #print("L4:", L4)  # Affiche L4

    return L4



def fonction_selection_stratified(Nb_particules, poids_i):
    '''
    Stratified resampling algorithm

    Args:
        Nb_particules: number of particles
        poids_i: vector of weights
    Returns:
        list of selected particles
    '''


    L1 = [x for x in range(Nb_particules)]
    #print('L1', L1)


    # on tire un U ~ Unif(i/N, (i+1)/N)
    positions_pointeurs = (np.arange(Nb_particules) + np.random.uniform(0, 1, Nb_particules)) / Nb_particules
    #print("positions pointeurs", positions_pointeurs)

    poids_cumules = np.cumsum(poids_i)


    vec_index_p_selec = []
    pointeur_index = 0


    for index_particule, pds_cum in enumerate(poids_cumules):
        while pointeur_index < Nb_particules and positions_pointeurs[pointeur_index] < pds_cum:
            vec_index_p_selec.append(index_particule)
            pointeur_index += 1
            if pointeur_index >= Nb_particules:
                break

    #print("vec_index_p_selec", vec_index_p_selec)

    if len(vec_index_p_selec) < Nb_particules:
        vec_index_p_selec += [len(poids_i) - 1] * (Nb_particules - len(vec_index_p_selec))
    elif len(vec_index_p_selec) > Nb_particules:
        vec_index_p_selec = vec_index_p_selec[:Nb_particules]




    # Étape 1 : Créer L3
    L3 = ['_' for _ in L1]

    count_L2 = {x: vec_index_p_selec.count(x) for x in vec_index_p_selec}
    # Remplir L3 avec les valeurs uniques de L2
    for num in L1:
        if count_L2.get(num, 0) > 0:
            L3[L1.index(num)] = num
            count_L2[num] -= 1
    #print("L3:", L3)  # Affiche L3


    # Étape 2 : Compléter L3 pour obtenir L4
    L4 = L3[:]
    remaining_elements = [num for num, count in count_L2.items() for _ in range(count)]
    # Remplacer les '_' par les éléments restants
    for i in range(len(L4)):
        if L4[i] == '_':
            L4[i] = remaining_elements.pop(0)


    return L4



def fonction_selection_multinomial(Nb_particules, poids_i):
    '''
        Multinomial resampling algorithm

        Args:
            Nb_particules: number of particles
            poids_i: vector of weights
        Returns:
            list of selected particles
    '''


    poids_i = np.asarray(poids_i, dtype=float)
    somme = poids_i.sum()
    if somme <= 0:
        raise ValueError("La somme des poids doit être strictement positive.")
    poids_i = poids_i / somme


    poids_cumules = np.cumsum(poids_i)

    poids_cumules[-1] = 1.0


    # N uniforms dans [0,1)
    u = np.random.uniform(0.0, 1.0, Nb_particules)
    # Pour chaque u, on trouve l'indice de la march
    vec_index_p_selec = np.searchsorted(poids_cumules, u, side="right").tolist()
    #print('vec index_p_selec:', vec_index_p_selec)


    if len(vec_index_p_selec) < Nb_particules:
        vec_index_p_selec += [len(poids_i) - 1] * (Nb_particules - len(vec_index_p_selec))
    elif len(vec_index_p_selec) > Nb_particules:
        vec_index_p_selec = vec_index_p_selec[:Nb_particules]

    L1 = [x for x in range(Nb_particules)]


    L3 = ['_' for _ in L1]

    count_L2 = {x: vec_index_p_selec.count(x) for x in vec_index_p_selec}

    for num in L1:
        if count_L2.get(num, 0) > 0:
            L3[L1.index(num)] = num
            count_L2[num] -= 1
    #print("L3:", L3)  # Affiche L3


    L4 = L3[:]

    remaining_elements = [num for num, count in count_L2.items() for _ in range(count)]

    for i in range(len(L4)):
        if L4[i] == '_':
            L4[i] = remaining_elements.pop(0)
    #print("L4:", L4)

    return L4




def fonction_selection_residual(Nb_particules, poids_i):
    '''
        Residual resampling algorithm

        Args:
            Nb_particules: number of particles
            poids_i: vector of weights
        Returns:
            list of selected particles
    '''


    L1 = [x for x in range(Nb_particules)]
    poids_i = np.asarray(poids_i, dtype=float)
    somme = poids_i.sum()
    if somme <= 0:
        raise ValueError("La somme des poids doit être strictement positive.")
    poids_i = poids_i / somme  # normalisation

    # ---------- Partie déterministe ----------
    copies_det = np.floor(Nb_particules * poids_i).astype(int)
    somme_det = copies_det.sum()
    R = Nb_particules - somme_det  # copies restantes à attribuer
    #print("copies_det:", copies_det, " | somme_det:", somme_det, " | R:", R)

    # ---------- Distribution résiduelle (multinomial) ----------
    vec_index_p_selec = []

    # 1) Ajouter la partie déterministe
    for idx, c in enumerate(copies_det):
        if c > 0:
            vec_index_p_selec.extend([idx] * c)

    # 2) S'il reste des copies, on les distribue selon les poids résiduels
    if R > 0:
        residus = Nb_particules * poids_i - copies_det
        somme_residus = residus.sum()
        if somme_residus > 0:
            p_res = residus / somme_residus
        else:
            p_res = np.ones_like(residus) / len(residus)

        # On tire R indices selon p_res
        cdf_res = np.cumsum(p_res)
        cdf_res[-1] = 1.0
        u = np.random.uniform(0.0, 1.0, R)
        extra = np.searchsorted(cdf_res, u, side="right").tolist()

        #print("poids_residuels_norm:", p_res)
        #print("u (residual):", u)
        #print("extra (indices résiduels):", extra)

        vec_index_p_selec.extend(extra)

    # ---------- Sécurité : taille exacte N ----------
    if len(vec_index_p_selec) < Nb_particules:
        vec_index_p_selec += [len(poids_i) - 1] * (Nb_particules - len(vec_index_p_selec))
    elif len(vec_index_p_selec) > Nb_particules:
        vec_index_p_selec = vec_index_p_selec[:Nb_particules]

    #print("vec_index_p_selec:", vec_index_p_selec)


    L3 = ['_' for _ in L1]
    count_L2 = {x: vec_index_p_selec.count(x) for x in vec_index_p_selec}
    # Remplir L3 avec les valeurs uniques de L2
    for num in L1:
        if count_L2.get(num, 0) > 0:
            L3[L1.index(num)] = num
            count_L2[num] -= 1
    #print("L3:", L3)  # Affiche L3

    L4 = L3[:]  # Copier L3 pour créer L4

    remaining_elements = [num for num, count in count_L2.items() for _ in range(count)]

    for i in range(len(L4)):
        if L4[i] == '_':
            L4[i] = remaining_elements.pop(0)
    #print("L4:", L4)

    return L4