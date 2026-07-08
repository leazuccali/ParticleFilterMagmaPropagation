import numpy as np
from scipy.interpolate import make_interp_spline  # Pour lisser les courbes
import json
import os
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
import matplotlib as mpl


########################################################################################################################
### Directories ########################################################################################################

directory = 'display_boxplots'
os.makedirs(directory, exist_ok=True)
directory_RMSE = f'RMSE_files'
os.makedirs(directory_RMSE, exist_ok=True)
date = '8.07.26'

########################################################################################################################
########################################################################################################################

### Global data/parameters - Units : step time (s), Pload (Pa), radius (m), step vectors (m)
global_data = {'step time' : 60, 'p_load' : -15000000, 'radius load' : 10000, 'step vectors' : 400,
                   'nb free propag' : 60, 'grid RG nb':399}

### Grid data/parameters - Units : m
grid_data = {'xmin': -30000,'xmax': 30000,'zmin': -15000,'zmax': -1,'discretisation step': 100}

### Particle reference data/parameters - Units : x0_ref, z0_ref (m) ; R,G (-) ; mu (Pa.s) ; E (Pa) ; Vol (m^3)
p_reference_data = {'x0_ref': 2000,'z0_ref': -8000,'R_ref': 0.5,'G_ref': 0.5,'mu_ref': 100, 'E_ref': 5*10**9,
                        'V_ref': 10**8, 'Step max': 550}

### Test data/parameters
# Selection type : 'systematic', 'stratified','multinomial', 'residual'
# Observation type : 'regulier', 'normal', 'random'
# If 'regulier' => Observation step.  If 'normal' or 'random' => Nb observations.
test_data = {'Nb particles':100, 'Assim. window': 10, 'Selection type': 'systematic', 'Observation type': 'regulier',
                 'Observation step (regular case)': 100, 'Nb observations (other cases)' : 0,
                   'dossier output': "output_github_test"}

### Resampling data
resampling_data = {'rd_xc_min': -200, 'rd_xc_max': 200,
                       'rd_zc_min': -200, 'rd_zc_max': 200,
                       'rd_lg_min': 0.5, 'rd_lg_max': 0.5,
                       'rd_open_min': 0.5, 'rd_open_max': 2,
                       'rd_dip_min': -0.25, 'rd_dip_max': 0.25,
                       'rd_veloc_min': -0.1, 'rd_veloc_max': 0.1,
                       'percentage':10}




### Data recuperation

pas_temps = global_data['step time']
P_load = global_data['p_load']
nb_it_depart = global_data['nb free propag']

xmin = grid_data['xmin']
xmax = grid_data['xmax']



x_0_ref = p_reference_data['x0_ref']
z_0_ref = p_reference_data['z0_ref']
R_ref = p_reference_data['R_ref']
G_ref = p_reference_data['G_ref']
mu_ref = p_reference_data['mu_ref']
E_ref = p_reference_data['E_ref']
vol_ref = p_reference_data['V_ref']
pas_max = p_reference_data['Step max']



Nb_particules = test_data['Nb particles']
pas_assimilation = test_data['Assim. window']
sorte_test = 'regulier'
if sorte_test == 'regulier':
    pts_observed = np.arange(xmin,xmax,100)
    nb_obs_points = len(pts_observed)
else:
    nb_obs_points=test_data['Nb observations (other cases)']

data_file = test_data['dossier output']


vec_pas = np.arange(nb_it_depart, pas_max, pas_assimilation)













plt.style.use("seaborn-v0_8-white")   # ou: seaborn-v0_8, -darkgrid, -dark, -white, -ticks
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm",    # math en Computer Modern
    "axes.unicode_minus": False, # joli signe moins
})


######################################################################################################################
### TEMPS & OBSERVATIONS
######################################################################################################################

print('vec_pas', vec_pas)
n_steps = len(vec_pas)
print('n_steps', n_steps)

temps_depart = pas_temps*nb_it_depart
temps_fin = pas_temps*pas_max
temps_assimilation = pas_temps*pas_assimilation
#time = np.arange(temps_depart, temps_fin+temps_assimilation, temps_assimilation)
time = np.arange(temps_depart, temps_fin, temps_assimilation)
#Temps lissé
nb_pts_lisses = 600
smooth_time = np.linspace(time.min(), time.max(), nb_pts_lisses)
time_pourcentage = (smooth_time/ (pas_max*pas_temps)) * 100
time_pourcentage_RMSE = (time/ (pas_max*pas_temps)) * 100
#fonction de lissage des courbes
def smooth(y): return make_interp_spline(time, y)(smooth_time)

###########################################################################################################################################################################################
###########################################################################################################################################################################################
def safe_smooth(y):
    spline = make_interp_spline(time, y)
    y_smooth = spline(smooth_time)
    return np.clip(y_smooth, 0, None)  # Force les valeurs positives




vec_pas_liste = vec_pas.tolist()
time_liste = time.tolist()
time_pourcentage_RMSE_liste = time_pourcentage_RMSE.tolist()

######  Création du dictionnaire des RMSE
dico_RMSE = {}
dico_RMSE['Vec pas'] = vec_pas_liste
dico_RMSE['Time'] = time_liste
dico_RMSE['Time pourcentage RMSE'] = time_pourcentage_RMSE_liste



#### Couleurs
forecast_color = 'pink'
analysis_color = 'darkorange'
true_color = '#8B0000'






######################################################################################################################
##import des données pour les prévisions des location/temps d'arrivée
######################################################################################################################

### paramètres de la particule de référence
data_filename = f'sauvegarde_trajectoire_ref.json'
with open(data_filename, 'r') as fichier:
    sauvegarde_trajectoire_ref = json.load(fichier)
vec_Xt_ref = sauvegarde_trajectoire_ref['X']
vec_temps_ref = sauvegarde_trajectoire_ref['Temps effectif']
## Temps final
position_final_ref = vec_Xt_ref[-1]
temps_final_ref = vec_temps_ref[-1]





#Constrution des vecteurs pour les paramètres
sauv_vecteurs_repartition_x_apres = []
sauv_vecteurs_repartition_temps_apres = []
sauv_vecteurs_repartition_x_avant = []
sauv_vecteurs_repartition_temps_avant = []

for pas in vec_pas:

    with open(f'{data_file}/step_{pas}/sauvegarde_predictions_avant_pas_{pas}.json', 'r') as fichier:
        param_trajectoires_i_avant = json.load(fichier)

    predictions_location_pas_avant = param_trajectoires_i_avant['Location']
    predictions_temps_final_avant = param_trajectoires_i_avant['Temps']
    sauv_vecteurs_repartition_x_avant.append(predictions_location_pas_avant)
    sauv_vecteurs_repartition_temps_avant.append(predictions_temps_final_avant)

    with open(f'{data_file}/step_{pas}/sauvegarde_predictions_apres_pas_{pas}.json', 'r') as fichier:
        param_trajectoires_i_apres = json.load(fichier)

    predictions_location_pas_apres = param_trajectoires_i_apres['Location']
    predictions_temps_final_apres = param_trajectoires_i_apres['Temps']
    sauv_vecteurs_repartition_x_apres.append(predictions_location_pas_apres)
    sauv_vecteurs_repartition_temps_apres.append(predictions_temps_final_apres)


print(sauv_vecteurs_repartition_temps_apres)



vec_X0_avant = []
vec_Z0_avant = []
vec_R_avant = []
vec_G_avant = []
vec_mu_avant = []
vec_E_avant= []
vec_vol_avant = []
vec_X0_retro_avant = []
vec_Z0_retro_avant = []

for pas in vec_pas:
    vec_X0_pas = []
    vec_Z0_pas = []
    vec_R_pas = []
    vec_G_pas = []
    vec_mu_pas = []
    vec_E_pas = []
    vec_vol_pas = []
    vec_X0_retro_pas = []
    vec_Z0_retro_pas = []
    with open(f'{data_file}/step_{pas}/sauvegarde_donnees_avant_pas_{pas}.json', 'r') as fichier:
        parametres_pas = json.load(fichier)
    for p in range(Nb_particules):
        dico_pas_p = parametres_pas[p]
        vec_X0_pas.append(dico_pas_p['Position x'])
        vec_Z0_pas.append(dico_pas_p['Position z'])
        vec_R_pas.append(dico_pas_p['R_0'])
        vec_G_pas.append(dico_pas_p['G_0'])
        vec_mu_pas.append(dico_pas_p['mu_0'])
        vec_E_pas.append(dico_pas_p['E_0'])
        vec_vol_pas.append(dico_pas_p['vol_0'])

        if 'Retro x_0' in dico_pas_p:
            vec_X0_retro_pas.append(dico_pas_p['Retro x_0'])
            vec_Z0_retro_pas.append(dico_pas_p['Retro z_0'])
        else:
            vec_X0_retro_pas.append(dico_pas_p['Position x'])
            vec_Z0_retro_pas.append(dico_pas_p['Position z'])

    vec_X0_avant.append(vec_X0_pas)
    vec_Z0_avant.append(vec_Z0_pas)
    vec_R_avant.append(vec_R_pas)
    vec_G_avant.append(vec_G_pas)
    vec_mu_avant.append(vec_mu_pas)
    vec_E_avant.append(vec_E_pas)
    vec_vol_avant.append(vec_vol_pas)
    vec_X0_retro_avant.append(vec_X0_retro_pas)
    vec_Z0_retro_avant.append(vec_Z0_retro_pas)



vec_X0_apres = []
vec_Z0_apres = []
vec_R_apres = []
vec_G_apres = []
vec_mu_apres = []
vec_E_apres= []
vec_vol_apres = []
vec_X0_retro_apres = []
vec_Z0_retro_apres = []
for pas in vec_pas:
    vec_X0_pas = []
    vec_Z0_pas = []
    vec_R_pas = []
    vec_G_pas = []
    vec_mu_pas = []
    vec_E_pas = []
    vec_vol_pas = []
    vec_X0_retro_pas = []
    vec_Z0_retro_pas = []
    with open(f'{data_file}/step_{pas}/sauvegarde_donnees_apres_pas_{pas}.json', 'r') as fichier:
        parametres_pas = json.load(fichier)

    for p in range(Nb_particules):
        dico_pas_p = parametres_pas[p]

        vec_X0_pas.append(dico_pas_p['Position x'])
        vec_Z0_pas.append(dico_pas_p['Position z'])
        vec_R_pas.append(dico_pas_p['R_0'])
        vec_G_pas.append(dico_pas_p['G_0'])
        vec_mu_pas.append(dico_pas_p['mu_0'])
        vec_E_pas.append(dico_pas_p['E_0'])
        vec_vol_pas.append(dico_pas_p['vol_0'])

        if 'Retro x_0' in dico_pas_p:
            vec_X0_retro_pas.append(dico_pas_p['Retro x_0'])
            vec_Z0_retro_pas.append(dico_pas_p['Retro z_0'])
        else:
            vec_X0_retro_pas.append(dico_pas_p['Position x'])
            vec_Z0_retro_pas.append(dico_pas_p['Position z'])

    vec_X0_apres.append(vec_X0_pas)
    vec_Z0_apres.append(vec_Z0_pas)
    vec_R_apres.append(vec_R_pas)
    vec_G_apres.append(vec_G_pas)
    vec_mu_apres.append(vec_mu_pas)
    vec_E_apres.append(vec_E_pas)
    vec_vol_apres.append(vec_vol_pas)
    vec_X0_retro_apres.append(vec_X0_retro_pas)
    vec_Z0_retro_apres.append(vec_Z0_retro_pas)




######################################################################################################################
## PLOTS
######################################################################################################################

##########################################################################################################################################################################################
### LOCATION ARRIVAL
###########################################################################################################################################################################################
# Forecats et analyse : récupération des données et calcul des statistiques
ensemble_avant = sauv_vecteurs_repartition_x_avant
ensemble_apres = sauv_vecteurs_repartition_x_apres
print('en apres', ensemble_apres, len(ensemble_apres))
# === True value ===
true_value = position_final_ref
print(true_value)
# === Forecast ===
ensemble_forecast = np.array(ensemble_avant).T
median_forecast = np.mean(ensemble_forecast, axis=0)
std_forecast = np.std(ensemble_forecast, axis=0)
# === Analyse ===
ensemble_analysis = np.array(ensemble_apres).T
median_analysis = np.mean(ensemble_analysis, axis=0)
std_analysis = np.std(ensemble_analysis, axis=0)



print('smooth time', smooth_time.shape)
mf = smooth(median_forecast)
print('smooth median', mf.shape)

############# smooth_time/1000
###Plot######
#############
#forecast_color = '#E74C3C'
#analysis_color = '#1F77B4'
# Forecast
fig, ax = plt.subplots(figsize=(15, 5))
plt.plot(time_pourcentage, smooth(median_forecast)/1000, color= 'grey', label="Forecast mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_forecast - std_forecast)/1000,
                 smooth(median_forecast + std_forecast)/1000,
                 color=forecast_color, alpha=0.4, label="Forecast spread (1σ)")
# linestyle='--'
# Analyse
plt.plot(time_pourcentage, smooth(median_analysis)/1000, color=analysis_color, label="Analysis mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_analysis-std_analysis)/1000,
                 smooth(median_analysis+std_analysis)/1000,
                 color=analysis_color, alpha=0.3, label="Analysis spread (1σ)")

# Theoretical value
plt.axhline(true_value/1000, color=true_color, linestyle='--', linewidth=2, label="Synthetic truth")

# Layout
plt.xlabel(r"Time (%)", fontsize=25)
plt.ylabel("Vents locations (km)", fontsize=25)

plt.yticks(fontsize=25)
plt.ylim(-15, 15)
plt.legend(ncol=3, fontsize=15, loc='lower right', frameon=True,  # active le cadre
                     fancybox=True,  # coins arrondis
                     framealpha=1.0,  # opacité du fond
                    edgecolor='black',  # couleur du bord
                     facecolor='white', bbox_to_anchor=(1.005, -0.02))

plt.xlim(11, 100)
plt.xticks(fontsize=25)

plt.tight_layout()

plt.annotate('10',
              xy=(0.01, -0.1),  # 0.5 = centre en x, 1.02 = légèrement au-dessus du haut du plot
              xycoords='axes fraction',  # coordonnées relatives au plot
              ha='center', va='bottom',  # alignement horizontal et vertical
              fontsize=25)  # optionnel : mettre en gras

plot_filename = os.path.join(directory, f'boxplots_evolution_location_new_{date}.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()





# === Calcul RMSE brute ===
erreur = true_value - np.array(ensemble_apres)  # shape = (n_steps, c)
rmse_instant = np.sqrt(np.mean((erreur)**2, axis=1))  # shape = (n_steps,)
rmse_instant_liste = rmse_instant.tolist()
dico_RMSE['RMSE Location'] = rmse_instant_liste

rmse_instant_normalisee = np.sqrt(np.mean( erreur**2 / true_value**2 , axis=1))
rmse_instant_normalisee_liste = rmse_instant_normalisee.tolist()
dico_RMSE['RMSE Location normalisee'] = rmse_instant_normalisee_liste













##########################################################################################################################################################################################
### TIMING ARRIVAL
###########################################################################################################################################################################################
# Forecats et analyse : récupération des données et calcul des statistiques
ensemble_avant = sauv_vecteurs_repartition_temps_avant
ensemble_apres = sauv_vecteurs_repartition_temps_apres
# === True value ===
true_value = temps_final_ref
# === Forecast ===
ensemble_forecast = np.array(ensemble_avant).T
median_forecast = np.mean(ensemble_forecast, axis=0)
std_forecast = np.std(ensemble_forecast, axis=0)
# === Analyse ===
ensemble_analysis = np.array(ensemble_apres).T
median_analysis = np.mean(ensemble_analysis, axis=0)
std_analysis = np.std(ensemble_analysis, axis=0)

print("Valeurs brutes forecast min/max:", np.min(ensemble_forecast), np.max(ensemble_forecast))
print("Valeurs brutes analysis min/max:", np.min(ensemble_analysis), np.max(ensemble_analysis))
print("Valeur de référence:", true_value)

#############
###Plot######
#############
#forecast_color = '#E74C3C'
#analysis_color = '#1F77B4'
# Forecast
fig, ax = plt.subplots(figsize=(15, 5))
plt.plot(time_pourcentage, safe_smooth(median_forecast)/1000, color='grey', label="Forecast mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 safe_smooth(median_forecast - std_forecast)/1000,
                 safe_smooth(median_forecast + std_forecast)/1000,
                 color=forecast_color, alpha=0.4, label="Forecast spread (1σ)")

# Analyse
plt.plot(time_pourcentage, safe_smooth(median_analysis)/1000, color=analysis_color, label="Analysis mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 safe_smooth(median_analysis - std_analysis)/1000,
                 safe_smooth(median_analysis + std_analysis)/1000,
                 color=analysis_color, alpha=0.3, label="Analysis spread (1σ)")

# Theoretical value
plt.axhline(true_value/1000, color=true_color, linestyle='--', linewidth=2, label="Synthetic truth")

# Layout
plt.xlabel(r"Time (%)", fontsize=25)
plt.ylabel(r"Arrival time ($10^3$s)", fontsize=25)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlim(time_pourcentage[0], time_pourcentage[-1])


plt.ylim(0, 3000)


plt.tight_layout()

plt.annotate('10',
              xy=(0.01, -0.1),  # 0.5 = centre en x, 1.02 = légèrement au-dessus du haut du plot
              xycoords='axes fraction',  # coordonnées relatives au plot
              ha='center', va='bottom',  # alignement horizontal et vertical
              fontsize=25)  # optionnel : mettre en gras


plot_filename = os.path.join(directory, f'boxplots_evolution_timing_new_{date}.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()





# === Calcul RMSE brute ===
erreur = true_value - np.array(ensemble_apres)  # shape = (n_steps, c)
rmse_instant = np.sqrt(np.mean((erreur)**2, axis=1))  # shape = (n_steps,)


rmse_instant_liste = rmse_instant.tolist()
dico_RMSE['RMSE Timing'] = rmse_instant_liste


rmse_instant_normalisee = np.sqrt(np.mean( erreur**2 / true_value**2 , axis=1))
rmse_instant_normalisee_liste = rmse_instant_normalisee.tolist()
dico_RMSE['RMSE Timing normalisee'] = rmse_instant_normalisee_liste




###########################################################################################################################################################################################
###########################################################################################################################################################################################

#####ZOOM Z0

#import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
#import os

fig, ax = plt.subplots(figsize=(15, 5))

# --- Forecast ---
ax.plot(time_pourcentage, safe_smooth(median_forecast)/1000,
        color='grey', label="Forecast mean", linewidth=3)

ax.fill_between(time_pourcentage,
                safe_smooth(median_forecast - std_forecast)/1000,
                safe_smooth(median_forecast + std_forecast)/1000,
                color=forecast_color, alpha=0.4, label="Forecast spread (1σ)")

# --- Analysis ---
ax.plot(time_pourcentage, safe_smooth(median_analysis)/1000,
        color=analysis_color,  label="Analysis mean", linewidth=3)

ax.fill_between(time_pourcentage,
                safe_smooth(median_analysis - std_analysis)/1000,
                safe_smooth(median_analysis + std_analysis)/1000,
                color=analysis_color, alpha=0.3, label="Analysis spread (1σ)")

# --- Valeur théorique ---
ax.axhline(true_value/1000, color=true_color, linestyle='--',
           linewidth=2, label="Synthetic truth")

# --- Layout ---
ax.set_xlabel(r"Time (%)", fontsize=25)
ax.set_ylabel(r"Arrival time ($10^3$s)", fontsize=25)
ax.tick_params(axis='x', labelsize=25)
ax.tick_params(axis='y', labelsize=25)
ax.set_xlim(time_pourcentage[0], time_pourcentage[-1])
ax.set_ylim(0, 3000)

plt.tight_layout()

# --- EXEMPLE ANNOTATION ---
plt.annotate('10',
             xy=(0.01, -0.1),
             xycoords='axes fraction',
             ha='center', va='bottom',
             fontsize=25)

# ==============================================
#                 ZOOM (INSET)
# ==============================================
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# 1) Créer un inset dans un coin
#axins = zoomed_inset_axes(ax, zoom=5.0, loc='upper right')
# loc = 'upper left', 'upper right', 'lower left', 'lower right'

axins = inset_axes(ax,
                   width="40%",   # largeur relative
                   height="50%",  # hauteur relative
                   loc="upper right")

# 2) Replotter ce qu’il y a dans le zoom
axins.plot(time_pourcentage, safe_smooth(median_forecast)/1000,
           color='grey', linewidth=2)

axins.fill_between(time_pourcentage,
                   safe_smooth(median_forecast - std_forecast)/1000,
                   safe_smooth(median_forecast + std_forecast)/1000,
                   color=forecast_color, alpha=0.4)

axins.plot(time_pourcentage, safe_smooth(median_analysis)/1000,
           color=analysis_color, linewidth=2)

axins.fill_between(time_pourcentage,
                   safe_smooth(median_analysis - std_analysis)/1000,
                   safe_smooth(median_analysis + std_analysis)/1000,
                   color=analysis_color, alpha=0.3)

axins.axhline(true_value/1000, color=true_color, linestyle='--', linewidth=1)

# 3) Limites du zoom — À AJUSTER SELON CE QUE TU VEUX MONTRER
x1, x2 = 80, 100      # Time (%) zoomé
y1, y2 = 0, 100    # Arrival time zoomée
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# 4) Ticks plus petits ou suppressions
axins.tick_params(labelsize=8)
l = [80, 85, 90, 95, 100]
axins.set_xticks(l)

# 5) Tracer rectangle + lignes de connexion
mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="black", linewidth=1)
plot_filename = os.path.join(directory, f'boxplots_evolution_timing_new_{date}_zoom.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()














##########################################################################################################################################################################################
### X0 retroprop
###########################################################################################################################################################################################
# Forecats et analyse : récupération des données et calcul des statistiques
ensemble_avant = vec_X0_retro_avant
ensemble_apres = vec_X0_retro_apres
# === True value ===
true_value = x_0_ref
# === Forecast ===
ensemble_forecast = np.array(ensemble_avant).T
median_forecast = np.mean(ensemble_forecast, axis=0)
std_forecast = np.std(ensemble_forecast, axis=0)
# === Analyse ===
ensemble_analysis = np.array(ensemble_apres).T
median_analysis = np.mean(ensemble_analysis, axis=0)
std_analysis = np.std(ensemble_analysis, axis=0)



#############
###Plot######
#############
#forecast_color = '#E74C3C'
#analysis_color = '#1F77B4'
# Forecast
fig, ax = plt.subplots(figsize=(15, 5))
plt.plot(time_pourcentage, smooth(median_forecast)/1000, color='grey', label="Forecast mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_forecast - std_forecast)/1000,
                 smooth(median_forecast + std_forecast)/1000,
                 color=forecast_color, alpha=0.4, label="Forecast spread (1σ)")

# Analyse
plt.plot(time_pourcentage, smooth(median_analysis)/1000, color=analysis_color, label="Analysis mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_analysis - std_analysis)/1000,
                 smooth(median_analysis + std_analysis)/1000,
                 color=analysis_color, alpha=0.3, label="Analysis spread (1σ)")

# Theoretical value
plt.axhline(true_value/1000, color=true_color, linestyle='--', linewidth=2, label="Synthetic value")

# Layout
plt.xlabel(r"Time (%)", fontsize=25)
plt.ylabel(r"Starting point $x_0$ (km)", fontsize=25)
plt.yticks(fontsize=25)
plt.ylim(-15, 15)


####
plt.xlim(10, 100)
plt.xticks([20, 30, 40, 50, 60, 70, 80, 90, 100], fontsize=25)
plt.annotate('10',
              xy=(0.011, -0.1),  # 0.5 = centre en x, 1.02 = légèrement au-dessus du haut du plot
              xycoords='axes fraction',  # coordonnées relatives au plot
              ha='center', va='bottom',  # alignement horizontal et vertical
              fontsize=25)  # optionnel : mettre en gras

param_gamme_m = -15
param_gamme_p = 15
# On trace une ligne verticale sur l’axe Y (x=0 dans les axes data)
ax.plot(
    [0, 0],               # x début et fin (0 = axe Y)
    [param_gamme_m, param_gamme_p],       # intervalle à colorer
    transform=ax.get_yaxis_transform(),  # Attache au système de l’axe Y
    color='indianred',
    linewidth=19, alpha = 0.5,
    zorder=10
)

#Ajout à la légende
handles, labels = ax.get_legend_handles_labels()
new_handle = mlines.Line2D([], [], color='indianred', alpha = 0.5, linewidth=10)
handles.append(new_handle)
labels.append('Initial distribution')

ax.legend(handles, labels, ncol=3, fontsize=15, loc='lower right', frameon=True,  # active le cadre
                    fancybox=True,  # coins arrondis
                    framealpha=1.0,  # opacité du fond
                    edgecolor='black',  # couleur du bord
                    facecolor='white', bbox_to_anchor=(1.005, -0.02))

#############
plt.tight_layout()


plot_filename = os.path.join(directory, f'boxplots_evolution_stpoint_x_new_{date}.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()





# === Calcul RMSE brute ===
erreur = true_value - np.array(ensemble_apres)  # shape = (n_steps, c)
rmse_instant = np.sqrt(np.mean((erreur)**2, axis=1))  # shape = (n_steps,)

rmse_instant_liste = rmse_instant.tolist()
dico_RMSE['RMSE X0'] = rmse_instant_liste


rmse_instant_normalisee = np.sqrt(np.mean( erreur**2 / true_value**2 , axis=1))
rmse_instant_normalisee_liste = rmse_instant_normalisee.tolist()
dico_RMSE['RMSE X0 normalisee'] = rmse_instant_normalisee_liste


###########################################################################################################################################################################################
###########################################################################################################################################################################################

## essai pour l'évolution. on choisi 5 temps
t0 = (nb_it_depart // 60 ) -1
t1 = (len(vec_pas)//2) #   * pas_assimilation
t2 = ((pas_max - nb_it_depart) // pas_assimilation) -1
temps_choisis = [t0, t1, t2] #, t3, t4, t5, t6]     #### à changer







# Paramètres subplot
n = len(temps_choisis)
fig, axes = plt.subplots(1, n, figsize=(6*n, 6), sharex=True, sharey=True)

# Si un seul subplot, axes n’est pas une liste → on le force à l’être
if n == 1:
    axes = [axes]

# Boucle sur les temps choisis
for i, t in enumerate(temps_choisis):
    ax = axes[i]
    mat_temps_x0 = np.array(vec_X0_retro_apres[t])
    mat_temps_z0 = np.array(vec_Z0_retro_apres[t])

    ax.scatter(mat_temps_x0 / 1000, mat_temps_z0 / 1000, color = 'k')
    ax.scatter(x_0_ref / 1000, z_0_ref / 1000, color='red', marker='*', s = 70, label='Theoretical value')

    # Annotations et style
    ax.set_xlim(-15, 15)
    ax.set_ylim(-16, -2)
    ax.grid(True)
    ax.set_title(f't = { (nb_it_depart*pas_temps + t*pas_temps*pas_assimilation) / (pas_max*pas_temps) * 100 :.0f} %', fontsize=16)
    ax.tick_params(labelsize=14)
    ax.legend(fontsize=14, loc='lower right')
    if i == 0:
        ax.set_ylabel("Z (km)", fontsize=20)
    ax.set_xlabel("X (km)", fontsize=20)


# Ajustement de l'espacement
plt.tight_layout(rect=[0, 0, 1, 0.95])  # Laisse de la place pour la légende globale
plt.suptitle(r"Spatial distribution of starting points ($x_0, z_0$)", fontsize=20)

plot_filename = os.path.join(directory, f'start_x_z_new_{date}.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()
plt.show()









##########################################################################################################################################################################################
### Z0 retropop
###########################################################################################################################################################################################
# Forecats et analyse : récupération des données et calcul des statistiques
ensemble_avant = vec_Z0_retro_avant
ensemble_apres = vec_Z0_retro_apres
# === True value ===
true_value = z_0_ref
# === Forecast ===
ensemble_forecast = np.array(ensemble_avant).T
median_forecast = np.mean(ensemble_forecast, axis=0)
std_forecast = np.std(ensemble_forecast, axis=0)
# === Analyse ===
ensemble_analysis = np.array(ensemble_apres).T
median_analysis = np.mean(ensemble_analysis, axis=0)
std_analysis = np.std(ensemble_analysis, axis=0)



#############
###Plot######
#############
#forecast_color = '#E74C3C'
#analysis_color = '#1F77B4'
# Forecast
fig, ax = plt.subplots(figsize=(15, 5))
plt.plot(time_pourcentage, smooth(median_forecast)/1000, color='grey', label="Forecast mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_forecast - std_forecast)/1000,
                 smooth(median_forecast + std_forecast)/1000,
                 color=forecast_color, alpha=0.4, label="Forecast spread (1σ)")

# Analyse
plt.plot(time_pourcentage, smooth(median_analysis)/1000, color=analysis_color, label="Analysis mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_analysis - std_analysis)/1000,
                 smooth(median_analysis + std_analysis)/1000,
                 color=analysis_color, alpha=0.3, label="Analysis spread (1σ)")

# Theoretical value
plt.axhline(true_value/1000, color=true_color, linestyle='--', linewidth=2, label="Theoretical value")

# Layout
plt.xlabel(r"Time (%)", fontsize=25)
plt.ylabel(r"Z-axis starting point position (km)", fontsize=25)
plt.yticks(fontsize=20)
plt.ylim(-16, -2)
plt.tight_layout()


plt.xlim(10, 100)
plt.xticks([20, 30, 40, 50, 60, 70, 80, 90, 100], fontsize=25)
plt.annotate('10',
              xy=(0.011, -0.1),  # 0.5 = centre en x, 1.02 = légèrement au-dessus du haut du plot
              xycoords='axes fraction',  # coordonnées relatives au plot
              ha='center', va='bottom',  # alignement horizontal et vertical
              fontsize=25)  # optionnel : mettre en gras

param_gamme_m = -10
param_gamme_p = -1
# On trace une ligne verticale sur l’axe Y (x=0 dans les axes data)
ax.plot(
    [0, 0],               # x début et fin (0 = axe Y)
    [param_gamme_m, param_gamme_p],       # intervalle à colorer
    transform=ax.get_yaxis_transform(),  # Attache au système de l’axe Y
    color='indianred',
    linewidth=19, alpha = 0.5,
    zorder=10
)

plot_filename = os.path.join(directory, f'boxplots_evolution_stpoint_z_new_{date}.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()





# === Calcul RMSE brute ===
erreur = true_value - np.array(ensemble_apres)  # shape = (n_steps, c)
rmse_instant = np.sqrt(np.mean((erreur)**2, axis=1))  # shape = (n_steps,)

rmse_instant_liste = rmse_instant.tolist()
dico_RMSE['RMSE Z0'] = rmse_instant_liste

rmse_instant_normalisee = np.sqrt(np.mean( erreur**2 / true_value**2 , axis=1))
rmse_instant_normalisee_liste = rmse_instant_normalisee.tolist()
dico_RMSE['RMSE Z0 normalisee'] = rmse_instant_normalisee_liste



###########################################################################################################################################################################################
###########################################################################################################################################################################################


































































##########################################################################################################################################################################################
### R
###########################################################################################################################################################################################
# Forecats et analyse : récupération des données et calcul des statistiques
ensemble_avant = vec_R_avant
ensemble_apres = vec_R_apres
# === True value ===
true_value = R_ref
# === Forecast ===
ensemble_forecast = np.array(ensemble_avant).T
median_forecast = np.mean(ensemble_forecast, axis=0)
std_forecast = np.std(ensemble_forecast, axis=0)
# === Analyse ===
ensemble_analysis = np.array(ensemble_apres).T
median_analysis = np.mean(ensemble_analysis, axis=0)
std_analysis = np.std(ensemble_analysis, axis=0)



#############
###Plot######
#############
#forecast_color = '#E74C3C'
#analysis_color = '#1F77B4'
# Forecast
fig, ax = plt.subplots(figsize=(15, 5))
plt.plot(time_pourcentage, smooth(median_forecast)/1, color='grey', label="Forecast mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_forecast - std_forecast)/1,
                 smooth(median_forecast + std_forecast)/1,
                 color=forecast_color, alpha=0.4, label="Forecast spread (1σ)")

# Analyse
plt.plot(time_pourcentage, smooth(median_analysis)/1, color=analysis_color, label="Analysis mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_analysis - std_analysis)/1,
                 smooth(median_analysis + std_analysis)/1,
                 color=analysis_color, alpha=0.3, label="Analysis spread (1σ)")

# Theoretical value
plt.axhline(true_value/1, color=true_color, linestyle='--', linewidth=2, label="Synthetic truth")

# Layout
plt.xlabel(r"Time (%)", fontsize=25)
plt.ylabel("R", fontsize=25)
plt.yticks(fontsize=25)
plt.ylim(0, 1)

plt.xlim(10, 100)
plt.xticks([20, 30, 40, 50, 60, 70, 80, 90, 100], fontsize=25)
plt.annotate('10',
              xy=(0.011, -0.1),  # 0.5 = centre en x, 1.02 = légèrement au-dessus du haut du plot
              xycoords='axes fraction',  # coordonnées relatives au plot
              ha='center', va='bottom',  # alignement horizontal et vertical
              fontsize=25)  # optionnel : mettre en gras

param_gamme_m = 0
param_gamme_p = 0.6
# On trace une ligne verticale sur l’axe Y (x=0 dans les axes data)
ax.plot(
    [0, 0],               # x début et fin (0 = axe Y)
    [param_gamme_m, param_gamme_p],       # intervalle à colorer
    transform=ax.get_yaxis_transform(),  # Attache au système de l’axe Y
    color='indianred',
    linewidth=19, alpha = 0.5,
    zorder=10
)

plt.tight_layout()

plot_filename = os.path.join(directory, f'boxplots_evolution_R_new_{date}.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()





# === Calcul RMSE brute ===
erreur = true_value - np.array(ensemble_apres)  # shape = (n_steps, c)
rmse_instant = np.sqrt(np.mean((erreur)**2, axis=1))  # shape = (n_steps,)
rmse_instant_liste = rmse_instant.tolist()
dico_RMSE['RMSE R'] = rmse_instant_liste

rmse_instant_normalisee = np.sqrt(np.mean( erreur**2 / true_value**2 , axis=1))
rmse_instant_normalisee_liste = rmse_instant_normalisee.tolist()
dico_RMSE['RMSE R normalisee'] = rmse_instant_normalisee_liste


###########################################################################################################################################################################################
###########################################################################################################################################################################################




##########################################################################################################################################################################################
### G
###########################################################################################################################################################################################
# Forecats et analyse : récupération des données et calcul des statistiques
ensemble_avant = vec_G_avant
ensemble_apres = vec_G_apres
# === True value ===
true_value = G_ref
# === Forecast ===
ensemble_forecast = np.array(ensemble_avant).T
median_forecast = np.mean(ensemble_forecast, axis=0)
std_forecast = np.std(ensemble_forecast, axis=0)
# === Analyse ===
ensemble_analysis = np.array(ensemble_apres).T
median_analysis = np.mean(ensemble_analysis, axis=0)
std_analysis = np.std(ensemble_analysis, axis=0)



#############
###Plot######
#############
#forecast_color = '#E74C3C'
#analysis_color = '#1F77B4'
# Forecast
fig, ax = plt.subplots(figsize=(15, 5))
plt.plot(time_pourcentage, smooth(median_forecast)/1, color='grey', label="Forecast mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_forecast - std_forecast)/1,
                 smooth(median_forecast + std_forecast)/1,
                 color=forecast_color, alpha=0.4, label="Forecast spread (1σ)")

# Analyse
plt.plot(time_pourcentage, smooth(median_analysis)/1, color=analysis_color, label="Analysis mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_analysis - std_analysis)/1,
                 smooth(median_analysis + std_analysis)/1,
                 color=analysis_color, alpha=0.3, label="Analysis spread (1σ)")

# Theoretical value
plt.axhline(true_value/1, color=true_color, linestyle='--', linewidth=2, label="Synthetic truth")

# Layout
plt.xlabel(r"Time (%)", fontsize=20)
plt.ylabel(r"$\gamma$", fontsize=20)
plt.yticks(fontsize=25)
plt.ylim(0, 1)

plt.xlim(10, 100)
plt.xticks([20, 30, 40, 50, 60, 70, 80, 90, 100], fontsize=25)
plt.annotate('10',
              xy=(0.011, -0.1),  # 0.5 = centre en x, 1.02 = légèrement au-dessus du haut du plot
              xycoords='axes fraction',  # coordonnées relatives au plot
              ha='center', va='bottom',  # alignement horizontal et vertical
              fontsize=25)  # optionnel : mettre en gras

param_gamme_m = 0.2
param_gamme_p = 0.8
# On trace une ligne verticale sur l’axe Y (x=0 dans les axes data)
ax.plot(
    [0, 0],               # x début et fin (0 = axe Y)
    [param_gamme_m, param_gamme_p],       # intervalle à colorer
    transform=ax.get_yaxis_transform(),  # Attache au système de l’axe Y
    color='indianred',
    linewidth=19, alpha = 0.5,
    zorder=10
)

plt.tight_layout()

plot_filename = os.path.join(directory, f'boxplots_evolution_G_new_{date}.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()





# === Calcul RMSE brute ===
erreur = true_value - np.array(ensemble_apres)  # shape = (n_steps, c)
rmse_instant = np.sqrt(np.mean((erreur)**2, axis=1))  # shape = (n_steps,)

rmse_instant_liste = rmse_instant.tolist()
dico_RMSE['RMSE G'] = rmse_instant_liste

rmse_instant_normalisee = np.sqrt(np.mean( erreur**2 / true_value**2 , axis=1))
rmse_instant_normalisee_liste = rmse_instant_normalisee.tolist()
dico_RMSE['RMSE G normalisee'] = rmse_instant_normalisee_liste




###########################################################################################################################################################################################
###########################################################################################################################################################################################




##########################################################################################################################################################################################
### mu
###########################################################################################################################################################################################
# Forecats et analyse : récupération des données et calcul des statistiques
ensemble_avant = vec_mu_avant
ensemble_apres = vec_mu_apres
# === True value ===
true_value = mu_ref
# === Forecast ===
ensemble_forecast = np.array(ensemble_avant).T
median_forecast = np.mean(ensemble_forecast, axis=0)
std_forecast = np.std(ensemble_forecast, axis=0)
# === Analyse ===
ensemble_analysis = np.array(ensemble_apres).T
median_analysis = np.mean(ensemble_analysis, axis=0)
std_analysis = np.std(ensemble_analysis, axis=0)



#############
###Plot######
#############
#forecast_color = '#E74C3C'
#analysis_color = '#1F77B4'
# Forecast
fig, ax = plt.subplots(figsize=(15, 5))
plt.plot(time_pourcentage, smooth(median_forecast)/1, color='grey', label="Forecast mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_forecast - std_forecast)/1,
                 smooth(median_forecast + std_forecast)/1,
                 color=forecast_color, alpha=0.4, label="Forecast spread (1σ)")

# Analyse
plt.plot(time_pourcentage, smooth(median_analysis)/1, color=analysis_color, label="Analysis mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_analysis - std_analysis)/1,
                 smooth(median_analysis + std_analysis)/1,
                 color=analysis_color, alpha=0.3, label="Analysis spread (1σ)")

# Theoretical value
plt.axhline(true_value/1, color=true_color, linestyle='--', linewidth=2, label="Synthetic truth")

# Layout
plt.xlabel(r"Time (%)", fontsize=25)
plt.ylabel(r"$\mu$ (Pas)", fontsize=25)
plt.yticks(fontsize=25)
plt.ylim(0, 2000)



plt.xlim(10, 100)
plt.xticks([20, 30, 40, 50, 60, 70, 80, 90, 100], fontsize=25)
plt.annotate('10',
              xy=(0.011, -0.1),  # 0.5 = centre en x, 1.02 = légèrement au-dessus du haut du plot
              xycoords='axes fraction',  # coordonnées relatives au plot
              ha='center', va='bottom',  # alignement horizontal et vertical
              fontsize=25)  # optionnel : mettre en gras

param_gamme_m = 100
param_gamme_p = 1000
# On trace une ligne verticale sur l’axe Y (x=0 dans les axes data)
ax.plot(
    [0, 0],               # x début et fin (0 = axe Y)
    [param_gamme_m, param_gamme_p],       # intervalle à colorer
    transform=ax.get_yaxis_transform(),  # Attache au système de l’axe Y
    color='indianred',
    linewidth=19, alpha = 0.5,
    zorder=10
)

plt.tight_layout()

plot_filename = os.path.join(directory, f'boxplots_evolution_mu_new_{date}.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()





# === Calcul RMSE brute ===
erreur = true_value - np.array(ensemble_apres)  # shape = (n_steps, c)
rmse_instant = np.sqrt(np.mean((erreur)**2, axis=1))  # shape = (n_steps,)


rmse_instant_liste = rmse_instant.tolist()
dico_RMSE['RMSE mu'] = rmse_instant_liste

rmse_instant_normalisee = np.sqrt(np.mean( erreur**2 / true_value**2 , axis=1))
rmse_instant_normalisee_liste = rmse_instant_normalisee.tolist()
dico_RMSE['RMSE mu normalisee'] = rmse_instant_normalisee_liste


###########################################################################################################################################################################################
###########################################################################################################################################################################################




##########################################################################################################################################################################################
### E
###########################################################################################################################################################################################
# Forecats et analyse : récupération des données et calcul des statistiques
ensemble_avant = vec_E_avant
ensemble_apres = vec_E_apres
# === True value ===
true_value = E_ref
# === Forecast ===
ensemble_forecast = np.array(ensemble_avant).T
median_forecast = np.mean(ensemble_forecast, axis=0)
std_forecast = np.std(ensemble_forecast, axis=0)
# === Analyse ===
ensemble_analysis = np.array(ensemble_apres).T
median_analysis = np.mean(ensemble_analysis, axis=0)
std_analysis = np.std(ensemble_analysis, axis=0)



#############
###Plot######
#############
#forecast_color = '#E74C3C'
#analysis_color = '#1F77B4'
# Forecast
fig, ax = plt.subplots(figsize=(15, 5))
plt.plot(time_pourcentage, safe_smooth(median_forecast)/1e9, color='grey', label="Forecast mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 safe_smooth(median_forecast - std_forecast)/1e9,
                 safe_smooth(median_forecast + std_forecast)/1e9,
                 color=forecast_color, alpha=0.4, label="Forecast spread (1σ)")

# Analyse
plt.plot(time_pourcentage, safe_smooth(median_analysis)/1e9, color=analysis_color, label="Analysis mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 safe_smooth(median_analysis - std_analysis)/1e9,
                 safe_smooth(median_analysis + std_analysis)/1e9,
                 color=analysis_color, alpha=0.3, label="Analysis spread (1σ)")

# Theoretical value
plt.axhline(true_value/1e9, color=true_color, linestyle='--', linewidth=2, label="Synthetic value")

# Layout
plt.xlabel(r"Time (%)", fontsize=20)
plt.ylabel(r"E (GPa)", fontsize=20)
plt.yticks(fontsize=25)
plt.ylim(0, 30)
plt.tight_layout()

plt.xlim(10, 100)
plt.xticks([20, 30, 40, 50, 60, 70, 80, 90, 100], fontsize=25)
plt.annotate('10',
              xy=(0.011, -0.1),  # 0.5 = centre en x, 1.02 = légèrement au-dessus du haut du plot
              xycoords='axes fraction',  # coordonnées relatives au plot
              ha='center', va='bottom',  # alignement horizontal et vertical
              fontsize=25)  # optionnel : mettre en gras

param_gamme_m = 10**9 /1e9
param_gamme_p = 10**10 /1e9
# On trace une ligne verticale sur l’axe Y (x=0 dans les axes data)
ax.plot(
    [0, 0],               # x début et fin (0 = axe Y)
    [param_gamme_m, param_gamme_p],       # intervalle à colorer
    transform=ax.get_yaxis_transform(),  # Attache au système de l’axe Y
    color='indianred',
    linewidth=19, alpha = 0.5,
    zorder=10
)

plot_filename = os.path.join(directory, f'boxplots_evolution_E_new_{date}.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()





# === Calcul RMSE brute ===
erreur = true_value - np.array(ensemble_apres)  # shape = (n_steps, c)
rmse_instant = np.sqrt(np.mean((erreur)**2, axis=1))  # shape = (n_steps,)
rmse_instant_liste = rmse_instant.tolist()
dico_RMSE['RMSE E'] = rmse_instant_liste
true_value_i = np.asarray(true_value)
rmse_instant_normalisee = np.sqrt(np.mean((erreur**2) / (true_value_i**2), axis=1))  # retourne (98,)


dico_RMSE['RMSE E normalisee'] = rmse_instant_normalisee_liste

###########################################################################################################################################################################################
###########################################################################################################################################################################################





##########################################################################################################################################################################################
### Volume
###########################################################################################################################################################################################
# Forecats et analyse : récupération des données et calcul des statistiques
ensemble_avant = vec_vol_avant
ensemble_apres = vec_vol_apres
# === True value ===
true_value = vol_ref
# === Forecast ===
ensemble_forecast = np.array(ensemble_avant).T
median_forecast = np.mean(ensemble_forecast, axis=0)
std_forecast = np.std(ensemble_forecast, axis=0)
# === Analyse ===
ensemble_analysis = np.array(ensemble_apres).T
median_analysis = np.mean(ensemble_analysis, axis=0)
std_analysis = np.std(ensemble_analysis, axis=0)



#############
###Plot######
#############

# Forecast
fig, ax = plt.subplots(figsize=(15, 5))
plt.plot(time_pourcentage, smooth(median_forecast)/1e9, color='grey', label="Forecast mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_forecast - std_forecast)/1e9,
                 smooth(median_forecast + std_forecast)/1e9,
                 color=forecast_color, alpha=0.4, label="Forecast spread (1σ)")

# Analyse
plt.plot(time_pourcentage, smooth(median_analysis)/1e9, color=analysis_color,  label="Analysis mean", linewidth=3)
plt.fill_between(time_pourcentage,
                 smooth(median_analysis - std_analysis)/1e9,
                 smooth(median_analysis + std_analysis)/1e9,
                 color=analysis_color, alpha=0.3, label="Analysis spread (1σ)")

# Theoretical value
plt.axhline(true_value/1e9, color=true_color, linestyle='--', linewidth=2, label="Synthetic truth")

# Layout
plt.xlabel(r"Time (%)", fontsize=20)
plt.ylabel(r"V ($km^3$)", fontsize=20)
plt.yticks(fontsize=25)
plt.ylim(0, 0.6)



plt.xlim(10, 100)
plt.xticks([20, 30, 40, 50, 60, 70, 80, 90, 100], fontsize=25)
plt.annotate('10',
              xy=(0.011, -0.1),  # 0.5 = centre en x, 1.02 = légèrement au-dessus du haut du plot
              xycoords='axes fraction',  # coordonnées relatives au plot
              ha='center', va='bottom',  # alignement horizontal et vertical
              fontsize=25)  # optionnel : mettre en gras

param_gamme_m = 5*10**6 / 1e9
param_gamme_p = 5*10**8 /1e9

# On trace une ligne verticale sur l’axe Y (x=0 dans les axes data)
ax.plot(
    [0, 0],               # x début et fin (0 = axe Y)
    [param_gamme_m, param_gamme_p],       # intervalle à colorer
    transform=ax.get_yaxis_transform(),  # Attache au système de l’axe Y
    color='indianred',
    linewidth=19, alpha = 0.5,
    zorder=10
)

plt.tight_layout()


plot_filename = os.path.join(directory, f'boxplots_evolution_vol_new_{date}.png')
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
plt.show()

# === Calcul RMSE brute ===
erreur = true_value - np.array(ensemble_apres)  # shape = (n_steps, c)
rmse_instant = np.sqrt(np.mean((erreur)**2, axis=1))  # shape = (n_steps,)
print("rmse_instant shape :", rmse_instant.shape)
print("time shape :", time.shape)

rmse_instant_liste = rmse_instant.tolist()
dico_RMSE['RMSE Vol'] = rmse_instant_liste

rmse_instant_normalisee = np.sqrt(np.mean( erreur**2 / true_value**2 , axis=1))
rmse_instant_normalisee_liste = rmse_instant_normalisee.tolist()
dico_RMSE['RMSE Vol normalisee'] = rmse_instant_normalisee_liste







###########################################################################################################################################################################################
###########################################################################################################################################################################################
data_filename = os.path.join(directory_RMSE, f'RMSE_{pas_assimilation}_min_{nb_obs_points}_obs_{sorte_test}.json')
with open(data_filename, 'w') as f:
     json.dump(dico_RMSE, f, indent=4)