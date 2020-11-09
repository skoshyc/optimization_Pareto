# optimization_Pareto_list of files and folders
paretofrontierhelper.R – contains helper code used in the pareto optimization as well as to create the ordinal lists

thickness_ordinal_list.R- creates the list of controls for the original model using latin hypercube sampling and runs iDynoMiCS for the original model to obtain the mean biofilm thickness for each combination of the controls.

thickness_ordinal_list_reduced.R- uses the list of controls from the original model and then the reduced models are simulated to obtained their mean biofilm thickness. The ordinal lists of the original and reduced models are then compared using Cohen’s kappa. 

pareto_repetitions.R- Code to find the Pareto frontier. It uses the package \texttt{nsga2R} with modifications to incorporate the iDynoMiCS software. It can be used for the original and reduced models.


iDynomics_protocol folder contains all the input xml files for the original and reduced multi_heterotroph model and the gingivalis-gordonii model.
