import numpy as np
import matplotlib.pyplot as plt
import json
import matplotlib as mpl
import os
import math



plt.style.use("seaborn-v0_8-white")   # ou: seaborn-v0_8, -darkgrid, -dark, -white, -ticks
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm",    # math en Computer Modern
    "axes.unicode_minus": False, # joli signe moins
})


########################################################################################################################
## Display functions
########################################################################################################################
# ============================================================
# GEOMETRIE
# ============================================================

def _cumulative_arclength(points):
    pts = np.asarray(points, float)
    seg = np.diff(pts, axis=0)
    ds = np.hypot(seg[:, 0], seg[:, 1])
    return np.concatenate([[0.0], np.cumsum(ds)])

def _resample_polyline(points, n_samples=400):
    """Rééchantillonne la polyligne à arclength uniforme (n_samples points)."""
    points = np.asarray(points, float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points doit être de forme (N,2).")
    s_old = _cumulative_arclength(points)
    L = s_old[-1]
    if L <= 0:
        raise ValueError("Polyligne de longueur nulle.")
    s_new = np.linspace(0.0, L, int(n_samples))
    x = np.interp(s_new, s_old, points[:, 0])
    y = np.interp(s_new, s_old, points[:, 1])
    return np.column_stack([x, y]), s_new, L

def _unit_tangent_and_normal(pts, s):
    """Tangente et normale unitaires le long de la ligne centrale."""
    dx = np.gradient(pts[:, 0], s)
    dy = np.gradient(pts[:, 1], s)
    t = np.column_stack([dx, dy])
    nrm = np.linalg.norm(t, axis=1)
    nrm[nrm == 0] = 1.0
    t = t / nrm[:, None]
    n = np.column_stack([-t[:, 1], t[:, 0]])
    return t, n

# ============================================================
# POLYGONE
# ============================================================
def polygon_from_centerline_and_opening(points, w, exaggeration=1.0, n_samples=400):
    """
    points : (N,2) centre du dyke (X,Z)
    w      : ouverture TOTALE le long de la polyligne
             (scalaire | vecteur long N | fonction f(s_frac))
    """
    center, s_res, L = _resample_polyline(points, n_samples=n_samples)
    _, n = _unit_tangent_and_normal(center, s_res)
    s_frac = s_res / L

    # fabrique w sur s_res (forme fixe)
    if np.isscalar(w):
        w_res = np.full_like(s_res, float(w))
    elif callable(w):
        w_res = np.asarray(w(s_frac), float)
        if w_res.shape != s_res.shape:
            raise ValueError("w(s_frac) doit renvoyer un vecteur de même taille que s.")
    else:
        w = np.asarray(w, float)
        if w.ndim != 1:
            raise ValueError("w doit être 1D.")
        s_old = _cumulative_arclength(points)
        w_res = np.interp(s_res, s_old, w)

    w_res = np.clip(w_res, 0.0, None) * float(exaggeration)

    upper = center + 0.5 * w_res[:, None] * n
    lower = center - 0.5 * w_res[:, None] * n
    poly = np.vstack([upper, lower[::-1], upper[0:1]])
    xi = 2.0 * s_frac - 1.0
    return poly, w_res, s_res, xi

# ============================================================
# OPENING PROFILES
# ============================================================
def beta_drop(a, b, wmax):
    """
    Profil type Beta : f(s) ∝ s^a (1-s)^b, normalisé à max=1 puis * wmax.
    s est la fraction curviligne s_frac ∈ [0,1].
    """
    def f(s_frac):
        s = np.clip(np.asarray(s_frac, float), 0.0, 1.0)
        val = np.power(s, a) * np.power(1.0 - s, b)
        m = np.max(val) if val.size > 0 else 0.0
        if m > 0:
            val = val / m
        return wmax * val
    return f

# (Optionnel) profil type Weertman
def weertman_opening_profile(C, L, s, skew=0.0, power=0.5):
    xi = 2.0 * (s / L) - 1.0
    xi = np.clip(xi - float(skew), -1.0, 1.0)
    w = float(C) * np.power(np.maximum(0.0, 1.0 - xi**2), power)
    return w

# ============================================================
# ============================================================
def _clip_polyline_to_top_by_zcut(points, z_cut):
    """
    Coupe la polyligne par un plan horizontal z=z_cut et garde la partie AU-DESSUS (z >= z_cut).
    Insère un point d'intersection si nécessaire.
    """
    pts = np.asarray(points, float)
    if pts.ndim != 2 or pts.shape[1] != 2 or pts.shape[0] == 0:
        return pts

    kept = []
    for (x0, z0), (x1, z1) in zip(pts[:-1], pts[1:]):
        above0 = (z0 >= z_cut)
        above1 = (z1 >= z_cut)

        if above0:
            if not kept or (kept[-1][0] != x0 or kept[-1][1] != z0):
                kept.append([x0, z0])

        if above0 != above1:
            # intersection linéaire en z
            if z1 != z0:
                t = (z_cut - z0) / (z1 - z0)
            else:
                t = 0.0
            xi = x0 + t * (x1 - x0)
            kept.append([xi, z_cut])

    if pts[-1, 1] >= z_cut:
        if not kept or (kept[-1][0] != pts[-1, 0] or kept[-1][1] != pts[-1, 1]):
            kept.append([pts[-1, 0], pts[-1, 1]])

    return np.asarray(kept, float)

def _reorder_from_top(center, w_res):
    """Réordonne pour que l'index 0 soit le point le plus HAUT (Z max), et que la séquence descende globalement en Z."""
    i_top = int(np.argmax(center[:, 1]))
    def descend_score(arr):
        z = arr[:, 1]
        return np.sum(np.diff(z[:min(10, len(z))]) < 0)
    c1 = np.roll(center, -i_top, axis=0); w1 = np.roll(w_res, -i_top, axis=0)
    c2 = c1[::-1].copy();               w2 = w1[::-1].copy()
    center_ord, w_ord = (c1, w1) if descend_score(c1) >= descend_score(c2) else (c2, w2)
    return center_ord, w_ord

def find_head_tail_split(center, w_res, frac_threshold=0.30, min_run=3):
    """
    Trouve le point de séparation tête/queue : premier run de 'min_run' points
    où w <= frac_threshold*wmax en partant du sommet.
    Retourne (split_idx_local, xy_split) dans l'espace réordonné (sommet -> bas).
    """
    if len(center) == 0:
        return None, None
    center_ord, w_ord = _reorder_from_top(center, w_res)
    wmax = float(np.max(w_ord)) if w_ord.size else 0.0
    if wmax <= 0:
        return 0, tuple(center_ord[0])

    thresh = frac_threshold * wmax
    run = 0
    split_idx_local = None
    for j in range(1, len(w_ord)):
        if w_ord[j] <= thresh:
            run += 1
            if run >= min_run:
                split_idx_local = j - min_run + 1
                break
        else:
            run = 0
    if split_idx_local is None:
        split_idx_local = len(w_ord) - 1
    xy_split = tuple(center_ord[split_idx_local])
    return split_idx_local, xy_split

# ============================================================
# ============================================================

def _auto_zcut_from_geometry(pts, w, L_theoretical, p_threshold):
    """
    Calcule automatiquement le plan de coupe z_cut (dévoilement) à partir
    de la hauteur réellement atteinte par la polyligne et de L_theoretical.
    Retourne : z_cut, dict_infos
    """
    pts = np.asarray(pts, float)
    if pts.ndim != 2 or pts.shape[1] != 2 or pts.shape[0] == 0:
        return None, {}

    z_top = float(np.max(pts[:, 1]))
    z_bottom_theory = z_top - float(L_theoretical)
    z_bottom_actual = float(np.min(pts[:, 1]))
    z_bottom_reachable = min(z_bottom_theory, z_bottom_actual)

    # Partie atteignable (pour mesurer tête/queue)
    pts_full = _clip_polyline_to_top_by_zcut(pts, z_bottom_reachable)
    if pts_full.shape[0] < 2:
        pts_full = np.vstack([pts_full, pts_full])

    center_full, s_full, L_full = _resample_polyline(pts_full, n_samples=400)
    # w(s) – forme fixe
    if np.isscalar(w):
        w_full = np.full_like(s_full, float(w))
    elif callable(w):
        w_full = np.asarray(w(s_full / L_full), float)
    else:
        s_old = _cumulative_arclength(pts_full)
        w_full = np.interp(s_full, s_old, np.asarray(w, float))

    # Split tête/queue
    _, xy_split = find_head_tail_split(center_full, w_full, frac_threshold=p_threshold, min_run=3)
    z_split = float(xy_split[1])

    H_now = z_top - z_bottom_actual                      # hauteur réellement atteinte
    H_target = min(H_now, float(L_theoretical))          # ce qu'on souhaite montrer
    H_head = max(0.0, z_top - z_split)
    H_tail = max(0.0, z_split - z_bottom_reachable)

    # Deux phases de dévoilement (auto)
    if H_head <= 1e-12:   # pas de tête détectable -> tout est queue
        frac_tail = 0.0 if H_tail <= 0 else H_target / H_tail
        z_cut = z_split - frac_tail * H_tail
    elif H_target <= H_head:
        # phase "tête" uniquement
        frac_head = H_target / H_head
        z_cut = z_top - frac_head * H_head
    else:
        # tête complète + fraction de queue
        remain = H_target - H_head
        frac_tail = 0.0 if H_tail <= 0 else remain / H_tail
        z_cut = z_split - frac_tail * H_tail

    # bornage
    z_cut = max(z_bottom_reachable, min(z_top, z_cut))
    info = dict(z_top=z_top, z_split=z_split, z_bottom_reachable=z_bottom_reachable,
                H_now=H_now, H_target=H_target, H_head=H_head, H_tail=H_tail, z_cut=z_cut)
    return z_cut, info

# ============================================================
# ============================================================
def _as_list(value, n, name="param"):
    """Si value est scalaire -> [value]*n, si liste de bonne longueur -> telle quelle."""
    if value is None:
        return [None] * n
    if np.isscalar(value):
        return [value] * n
    value = list(value)
    if len(value) != n:
        raise ValueError(f"{name} doit être un scalaire ou une liste de longueur {n} (reçu {len(value)}).")
    return value

# ============================================================
# Tracé
# ============================================================
def plot_openings_progressive_threshold(points_list, w,
                                        L_theoretical,            # scalaire OU liste
                                        p_threshold=0.30,
                                        exaggeration=1.0,
                                        alpha_len=2.0,            # <-- ton α
                                        opening_ref_mode="max",   # 'max' ou 'mean'
                                        scatter_kwargs=None,      # style des points en mode "points"
                                        ax=None,
                                        **fill_kwargs):
    """
    Si la longueur visible L_temps est < alpha_len * W_ref (ouverture de référence),
    on affiche uniquement des points le long de la trajectoire visible.
    Sinon, on trace le polygone (comportement normal).

    Retourne (results, ax) avec pour chaque dyke : mode ('points'/'polygon'), L_temps, W_ref, z_cut, L_theoretical.
    """
    n = len(points_list)
    L_th_list = _as_list(L_theoretical, n, name="L_theoretical")

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 10))

    if scatter_kwargs is None:
        scatter_kwargs = dict(s=16)  # taille des marqueurs par défaut

    results = []

    for i, pts in enumerate(points_list):
        pts = np.asarray(pts, float)
        if pts.ndim != 2 or pts.shape[1] != 2:
            raise ValueError("Chaque entrée de points_list doit être un tableau (N,2).")

        # 1) Plan de coupe auto pour la longueur théorique de CE dyke
        z_cut, auto_info = _auto_zcut_from_geometry(pts, w, L_th_list[i], p_threshold)

        # 2) Partie visible (au-dessus de z_cut)
        pts_vis = _clip_polyline_to_top_by_zcut(pts, z_cut)
        if pts_vis.shape[0] < 2:
            pts_vis = np.vstack([pts_vis, pts_vis])

        # 3) Rééchantillonnage de la partie visible
        center, s_res, L_vis = _resample_polyline(pts_vis, n_samples=200)  # L_vis = L_temps

        # 4) Profil d'ouverture sur la partie visible (forme fixe)
        if np.isscalar(w):
            w_res = np.full_like(s_res, float(w))
        elif callable(w):
            # on normalise par la longueur visible pour évaluer w(s_frac)
            w_res = np.asarray(w(s_res / s_res[-1]), float)
        else:
            s_old = _cumulative_arclength(pts_vis)
            w_res = np.interp(s_res, s_old, np.asarray(w, float))

        # 5) Choix de W_ref (max ou moyenne)
        W_ref = float(np.max(w_res)) if opening_ref_mode == "max" else float(np.mean(w_res))

        # 6) Test de seuil
        L_temps = float(L_vis)
        show_points_only = (L_temps < float(alpha_len) * W_ref)

        # 7) Tracé selon le mode
        if show_points_only:
            n_skip = 5  # ou 5, selon ton besoin
            center_sparse = center[::n_skip]
            ax.scatter(center_sparse[:, 0] / 1000.0, center_sparse[:, 1] / 1000.0, **scatter_kwargs)
            #ax.scatter(center[:, 0]/1000.0, center[:, 1]/1000.0, **scatter_kwargs)
            mode = "points"
        else:
            poly, w_plot, s_plot, xi = polygon_from_centerline_and_opening(
                center, w_res, exaggeration=exaggeration, n_samples=len(s_res)
            )
            ax.fill(poly[:, 0]/1000.0, poly[:, 1]/1000.0, **fill_kwargs)
            mode = "polygon"

        # 8) Sauvegarde des infos
        results.append(dict(
            mode=mode,
            L_temps=L_temps,
            W_ref=W_ref,
            z_cut=z_cut,
            L_theoretical=L_th_list[i]
        ))

    # Mise en forme globale
    plt.xlim(-15, 15)
    ticks = [-14.6,-10,-5,0,5,10,15]
    labels = ["-15", "-10", "-5", "0", "5","10","15"]
    plt.xticks(ticks, labels)




    plt.ylim(-10, 0)

    ax.set_aspect("equal")
    ax.set_xlabel("X (km)", fontsize=30, labelpad=10)
    ax.set_ylabel("Z (km)", fontsize=30, labelpad=-10)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')  # si tu as un xlabel
    ax.tick_params(axis='x', which='major', length=8, width=1.2)  # tirets plus longs/épais
    ax.tick_params(axis='y', which='major', length=8, width=1.2)  # tirets plus longs/épais

    return results, ax


if __name__ == "__main__":

    #dossier_source = 'output_6.07.25_test_github_2'
    directory_images = 'dossier_images_github_test'


    global_data = {'step time': 60, 'p_load': -15000000, 'radius load': 10000, 'step vectors': 400,
                   'nb free propag': 60, 'grid RG nb': 399}

    p_reference_data = {'x0_ref': 2000, 'z0_ref': -8000, 'R_ref': 0.5, 'G_ref': 0.5, 'mu_ref': 100,
                        'E_ref': 5 * 10 ** 9,
                        'V_ref': 10 ** 8, 'Step max': 550}

    test_data = {'Nb particles': 100, 'Assim. window': 5, 'Selection type': 'systematic',
                 'Observation type': 'regulier',
                 'Observation step (regular case)': 100, 'Nb observations (other cases)': 0,
                 'dossier output': "output_github_test"}


    Nb_particules = test_data['Nb particles']
    P_load = global_data['p_load']
    nb_step_min = global_data['nb free propag']
    nb_step_max = p_reference_data['Step max']
    pas_image = 20       #test_data['Assim. window']
    dossier_source = test_data['dossier output']
    ######### STEPS
    vec_pas = np.arange(nb_step_min, nb_step_min + 1, pas_image)
    ########



    # Chemin vers tes fichiers (utilise {pas} pour injecter le step)
    with open('sauvegarde_trajectoire_ref.json', 'r') as f:
        donnees_ref = json.load(f)
    longueur = donnees_ref[('Longueur remontee')]




    #### Parameters
    # Profil d'ouverture fixe (identique pour tous)
    wmax_m = 9.1515551378438
    w_fixed = beta_drop(a=3.0, b=0.6, wmax=wmax_m)
    # Longueur verticale théorique du dyke (en mètres)
    L_theoretical = longueur

    # Paramètres
    p_seuil = 0.30  # seuil tête/queue en fraction de w_max
    exaggeration = 50.0
    facecolor = "steelblue"
    facealpha = 0.5


    for pas in vec_pas:

        print(f"--- Step {pas} ---")
        with open(f'{dossier_source}/step_{pas}/sauvegarde_traj_ref_pas_{pas}.json', 'r') as f:
            data = json.load(f)
        coordsX = np.asarray(data['Trajectoire X'], float)
        coordsZ = np.asarray(data['Trajectoire Z'], float)
        M = np.column_stack((coordsX, coordsZ))

        with open(f'{dossier_source}/step_{pas}/sauvegarde_traj_ref_pas_{pas}.json', 'r') as fichier:
            trajectoires_ref_pas = json.load(fichier)
        memoire_x_ref = trajectoires_ref_pas['Memoire X']
        memoire_z_ref = trajectoires_ref_pas['Memoire Z']
        pt_bas_x_ref = trajectoires_ref_pas['Point bas X']
        pt_bas_z_ref = trajectoires_ref_pas['Point bas Z']
        memoire_x_ref_bas = [x for x in memoire_x_ref if x <= pt_bas_x_ref]
        memoire_z_ref_bas = [x for x in memoire_z_ref if x <= pt_bas_z_ref]


        with open(f'{dossier_source}/step_{pas}/sauvegarde_traj_apres_pas_{pas}.json', 'r') as fichier:
            trajectoires_pas = json.load(fichier)

        matrice_trajectoires = []
        matrice_memoire_x = []
        matrice_memoire_z = []
        matrice_bas_x = []
        matrice_bas_z = []
        for p in range(Nb_particules):
            trajectoires_p_x = trajectoires_pas[p]["Trajectoire X"]
            trajectoires_p_z = trajectoires_pas[p]["Trajectoire Z"]
            M_particule = np.column_stack((trajectoires_p_x, trajectoires_p_z))
            matrice_trajectoires.append(M_particule)

            memoire_p_x = trajectoires_pas[p]["Memoire X"]
            memoire_p_z = trajectoires_pas[p]["Memoire Z"]
            pt_bas_p_x = trajectoires_pas[p]["Point bas X"]
            pt_bas_p_z = trajectoires_pas[p]["Point bas Z"]

            matrice_memoire_x.append(memoire_p_x)
            matrice_memoire_z.append(memoire_p_z)
            matrice_bas_x.append(np.array(pt_bas_p_x))
            matrice_bas_z.append(np.array(pt_bas_p_z))



        with open(f'{dossier_source}/step_{pas}/sauvegarde_parametres_physiques_apres_{pas}.json',
                  'r') as fichier:
            parametres_physiques_apres = json.load(fichier)
        volumes = parametres_physiques_apres['vol']
        liste_G = parametres_physiques_apres['G']
        liste_E = parametres_physiques_apres['E']
        liste_R = parametres_physiques_apres['R']
        liste_mu = parametres_physiques_apres['mu']



        liste_X0_non = parametres_physiques_apres['X_0']
        liste_Z0_non = parametres_physiques_apres['Z_0']

        with open(f'{dossier_source}/step_{pas}/sauvegarde_donnees_apres_pas_{pas}.json', 'r') as fichier:
            d = json.load(fichier)

        liste_X0 = []
        liste_Z0 = []
        for p in range(Nb_particules):
            dico = d[p]
            X0 = dico["Retro x_0"] if "Retro x_0" in dico else dico["Position x"]
            Z0 = dico["Retro z_0"] if "Retro z_0" in dico else dico["Position z"]
            liste_X0.append(X0)
            liste_Z0.append(Z0)


        longueurs_theoriques = []
        ouvertures_theoriques = []
        for p in range(Nb_particules):
            r3_long = (liste_E[p] * volumes[p]) / ((1 - 0.25 ** 2) * liste_G[p] * abs(P_load))
            longueur = r3_long ** (1 / 3)
            longueurs_theoriques.append(longueur)
            ouverture = (math.pi * volumes[p]) / (2 * longueur ** 2)
            ouvertures_theoriques.append(ouverture)

        print(longueurs_theoriques)

        # Display

        fig, ax = plt.subplots(figsize=(12, 10))

        n = len(matrice_trajectoires)
        cmap = plt.get_cmap("RdPu")
        # on évite la zone [cut, 1]
        cut1 = 0.15
        cut = 0.6
        colors = cmap(np.linspace(cut1, cut, n))

        for i, pts in enumerate(matrice_trajectoires):
            pts = np.asarray(pts, float)
            print(i, pts.shape, pts[:5])

        for i, (pts, color) in enumerate(zip(matrice_trajectoires, colors)):
            res, _ = plot_openings_progressive_threshold(
                points_list=[pts],
                w=w_fixed,
                L_theoretical=float(longueurs_theoriques[i]),
                p_threshold=p_seuil,
                exaggeration=exaggeration,
                color=color, alpha=0.25,     #alpha=0.3
                ax=ax,
                alpha_len=25.0,  # <-- ton α (à ajuster)
                opening_ref_mode="max",  # 'max' ou 'mean'
                scatter_kwargs=dict(s=65, alpha=0.01, marker=".", color=color),
            )
            #print(res[0]["mode"], res[0]["L_temps"], res[0]["W_ref"])

        # 2) Dyke théorique sur le même ax (couleur différente)
        results_theo, ax = plot_openings_progressive_threshold(
            points_list=[M],  # polyligne théorique
            w=w_fixed,
            L_theoretical=L_theoretical,  # scalaire ou liste de longueur 1
            p_threshold=p_seuil,
            exaggeration=exaggeration,
            label="Vérité synthétique",   #Synthetic truth
            color='red', alpha=1,
            ax=ax,  # <-- important
            alpha_len=2.0,  # <-- ton α (à ajuster)
            opening_ref_mode="max",  # 'max' ou 'mean'
            scatter_kwargs=dict(s=20, alpha=0.9),
        )

        leg = ax.legend(title=f"t = {round((pas * 100) / 550)}%", title_fontsize=16, fontsize=14,
                        loc="upper left", frameon=True)
        leg.get_frame().set_alpha(1.0)


        #directory_images = 'dossier_images_7.7.26'
        os.makedirs(directory_images, exist_ok=True)
        plot_filename = os.path.join(directory_images, f'position_pas_{pas}.png')
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')

