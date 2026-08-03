# A workflow for modelling the design of (cluster-) randomised clinical trials for malaria interventions

The code in this repository presents a modelling framework for (randomised and) cluster-randomised trials for malaria interventions which tracks individual malaria outcomes over time.

The modelling is based on malaria transmission simulations run using Imperial College London's [malariasimulation](https://github.com/mrc-ide/malariasimulation) model, for which there is much more and far better documentation available [here](https://mrc-ide.github.io/malariasimulation/index.html).
The workflow uses the model above, and using helper functions and wrappers, run simulations and analyse their outcomes in the context of simulating clinical trials.

Briefly, the two main scripts are to be used in this order:

1. `rum_sim_script.R`: runs the trial simulation for a control and an intervention arm, allowing to modulate baseline malaria transmission parameters and other parameters to adjust trial design.

2. `trial_analysis_script.R`: reads the data from the previous step; and analyses prevalence, incidence and time-to-infection in control and intervention arms.
