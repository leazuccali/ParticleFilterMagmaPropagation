README at [Voir le PDF](https://github.com/leazuccali/ParticleFilterMagmaPropagation/blob/PF_V4/README.pdf)
or download at https://raw.githubusercontent.com/leazuccali/ParticleFilterMagmaPropagation/PF_V4/README.pdf

All our experiments were performed on 2 machines. Each machine is equipped with 2 Intel Xeon Gold 5220R processors (2.20 GHz), providing 96 CPU cores, and 125 GB of RAM. The operating system is Ubuntu 22.04.5 LTS.

----------------------------------------------------------------------------------------------------------------

Structure of one simulation

BEFORE SIMULATION

**0. Definition of simulation parameters.** global_data, grid_data, p_reference_data, resampling_data are 5 dictionnaries that defines the simulation. Each dictionnary is describes in the pdf.

**1. Generation of various elements**

* (R, $\gamma$) grids

-> To generate independently with creation_grilles_RG.py, $\approx 3$ minutes.

* truth particle : creation_fichier_particule_reference.py
* a set of n particles : creation_fichier_particule_joblib.py $\approx 3$ minutes.

-> Can be generate with the independently of directly by running main.py

SIMULATION

**2. Simulation of propagation**

By running main.py 

Outputs are stored step by step (step_X) in a general folder output_test. 

AFTER SIMULATION

**3. Results of simulation**

Particle evolution (in output_test_github) is illustrated with plot_particles_positions.py, genere_animation_gif.py and plot_distribution_values.py.


-----------------------------------------------------------------------------------------------
This code is ready to run the simulation case:

- Truth particle parameters ($x_0=2000, z_0=-8000, R=0.5, \gamma=0.5, \mu=100, E=5.10^9, V=10^8$).

- Test with 100 particles, assimilation window each 5 time_step (5x60=300 seconds), 600 observations regularly spaced, systematic selection and 399 $(R,\gamma )$ grid possibilities.

- Time to simulate : $\approx 1$ hour. Memory : $\approx$ 4-5 Go.


Simulation parameters definition : 

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
    
    test_data = {'Nb particles':100, 'Assim. window': 5, 'Selection type': 'systematic', 'Observation type': 'regulier',
                 'Observation step (regular case)': 100, 'Nb observations (other cases)' : 0,
                   'dossier output': "output_test }

### Resampling data
    resampling_data = {'rd_xc_min': -200, 'rd_xc_max': 200,
                       'rd_zc_min': -200, 'rd_zc_max': 200,
                       'rd_lg_min': 0.5, 'rd_lg_max': 0.5,
                       'rd_open_min': 0.5, 'rd_open_max': 2,
                       'rd_dip_min': -0.25, 'rd_dip_max': 0.25,
                       'rd_veloc_min': -0.1, 'rd_veloc_max': 0.1,
                       'percentage':10}

