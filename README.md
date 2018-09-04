## Simulating data acquisition timing errors in Eddy Covariance

This repo contains the code used to perform the Random Timing Error (RTE) and 
Systematic Timing Error (STE) simulation in Eddy Covariance (EC) data acquisition, 
starting from existing raw EC anemometric data.

The simulation code was built as an add-on to the source code of EddyPro 5.2.0's engine. 
The key routines for the simulation are `sync_simulate_misalignment` and those called therein. 

The simulation can be performed using any raw EC dataset containing data of 
the three wind components and of sonic temperature or speed of sound, at frequencies > 1 Hz.

For an explanation of the simulation rationale and design, and for evaluation of
corresponding results, see:

[Fratini, G., Sabbatini, S., Ediger, K., Nicolini, G., Vitale, D., Papale, D., Eddy Covariance flux errors due to random and systematic timing errors during data acquisition, Biogeosciences Discussions, pp 1-21, doi:10.5194/bg-2018-177, 2018.](https://www.biogeosciences-discuss.net/bg-2018-177/)


### EddyPro
[EddyPro](www.licor.com/eddypro) is developed, maintained and supported by [LI‑COR Biosciences](www.licor.com). 
It originates from [ECO<sub>2</sub>S](http://gaia.agraria.unitus.it/eco2s),
the Eddy COvariance COmmunity Software project, which was developed as part
of the Infrastructure for Measurement of the European Carbon Cycle (IMECC-EU)
research project. We gratefully acknowledge the
[IMECC](http://imecc.ipsl.jussieu.fr/index.html) consortium,
the ECO<sub>2</sub>S development team, the [University of Tuscia](www.unitus.it) (Italy)
and scientists around the world who assisted with development and testing of
the original version of this software.
