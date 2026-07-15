from propagation_steps import *
from creation_fichier_particule_reference import *
from creation_fichier_particule_joblib import *

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

    ### Test data/parameters
    # Selection type : 'systematic', 'stratified','multinomial', 'residual'
    # Observation type : 'regulier', 'normal', 'random'
    # If 'regulier' => Observation step.  If 'normal' or 'random' => Nb observations.
    test_data = {'Nb particles': 100, 'Assim. window': 10, 'Selection type': 'systematic',
                 'Observation type': 'regulier',
                 'Observation step (regular case)': 100, 'Nb observations (other cases)': 0,
                 'dossier output': "output_test"}


    ### Resampling data
    resampling_data = {'rd_xc_min': -200, 'rd_xc_max': 200,
                       'rd_zc_min': -200, 'rd_zc_max': 200,
                       'rd_lg_min': 0.5, 'rd_lg_max': 2,
                       'rd_open_min': 0.5, 'rd_open_max': 2,
                       'rd_dip_min': -0.25, 'rd_dip_max': 0.25,
                       'rd_veloc_min': -0.1, 'rd_veloc_max': 0.1,
                       'percentage': 10}



    ####################################################################################################################
    ####################################################################################################################
    creation_trajectoire_reference(global_data, grid_data, p_reference_data, 0)

    ### Creation of random particles and associated trajectories
    sauvegarde_parametres = creation_particules(test_data)
    creation_fichier_trajectoire_X = creation_trajectoires(global_data, grid_data, test_data, 0)

    ### Lancement de la simulation de la propagation
    propagation_step_function(global_data, grid_data, p_reference_data, test_data, resampling_data)