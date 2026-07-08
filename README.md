README at [Voir le PDF](https://github.com/leazuccali/ParticleFilterMagmaPropagation/blob/PF_V4/README.pdf)
or download with https://raw.githubusercontent.com/leazuccali/ParticleFilterMagmaPropagation/PF_V4/README.pdf

All our experiments were performed on 2 machines. Each machine is equipped with 2 Intel Xeon Gold 5220R processors (2.20 GHz), providing 96 CPU cores, and 125 GB of RAM. The operating system is Ubuntu 22.04.5 LTS.

----------------------------------------------------------------------------------------------------------------

Structure of one simulation

BEFORE SIMULATION

**0. Definition of simulation parameters.** global_data, grid_data, p_reference_data, resampling_data are 5 dictionnaries that defines the simulation. Each dictionnary is describes in the pdf.

**1. Generation of various elements**

* (R, $\gamma$) grids

-> To generate independently with creation_grilles_RG.py

* truth particle : creation_fichier_particule_reference.py
* a set of n particles : creation_fichier_particule_joblib.py

-> Can be generate with the independently of directly by running main.py

SIMULATION

**2. Simulation of propagation**

By running main.py 

Outputs are stored step by step (step_X) in a general folder output_test. 

AFTER SIMULATION

**3. Results of simulation**

Particle evolution (in output_test) is illustrated with plot_particles_positions.py, genere_animation_gif.py and plot_distribution_values.py.


-----------------------------------------------------------------------------------------------
This code is ready to run the simulation :

- Truth particle parameters ($x_0=2000, z_0=-8000, R=0.5, \gamma=0.5, \mu=100, E=5.10^9, V=10^8$).

- Test with 100 particles, assimilation window each 5 time_step (5x60=300 seconds), 600 observations regularly spaced, systematic selection.

- Time to simulate : $\approx 1$ hour. Memory : $\approx$ 4-5 Go.
