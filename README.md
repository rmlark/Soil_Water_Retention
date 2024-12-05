# Soil_Water_Retention
Demonstration R script which uses CEPHaS project functions for inference about van Genuchten models of soil water retention.  The input data are measurements of the volumetric water content of the soil at a standard set of tensions (expressed in kPa, positive values).

The functions use the saem library for estimation of non-linear mixed effects models.  They are set up for inference about differences between the water retention curve (SWRC) of soils in different categories (e.g. experimental treatments, in this case two contrasting classes of the World Reference Base classification).  The soil data are from Kenya, extracted from the WOSIS data base, available under the CC-BY-NC licence. Details are in the script file Kenya.R.

The script also demonstrates parametric bootstrap estimation of the estimation distribution of the SWRC parameters.  It contains a function to return soil physical quality parameters (details in the script and the paper in Geoderma).  

Finally, the script contains a function to simulate infiltration of water into the soil, using a 15-minute timestep weather data file which is provided.  This is a version of the Green Ampt model of water movement into soil, as proposed for layered soil by Liu et al.


