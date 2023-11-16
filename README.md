# gridEZ

gridEZ algorithm for generating enumeration zones (EZs) with user-defined target population and geographic size. 

EZs are generated from three input rasters that cover the study region. These are: (1) population counts or densities; (2) stratum IDs; (3) settlement type IDs.

### Instructions

1. Download 'gridEZ_running_script_public_release_v2.R' and 'gridEZ_fn_public_release_v2.R'
2. Open 'gridEZ_running_script_public_release_v2.R' in R 
3. Edit ncores and par_type to suit your computing system's parallel processing 
4. Edit the line of code that sources the gridEZ function to point to the file location of gridEZ_fn_public_release_v2.R on your computer 
5. Load in the population, strata, and settlement rasters for your project by editing the file paths in the 'set input files' section
6. Edit the gridEZ() specifications for generating EZs for your project under the 'run grid EZ code' section
7. Run the entire gridEZ_running_script_public_release_v2.R script

### Saint Kitts and Nevis example

Please try out this example to see what an output sampling frame looks like. The R script includes code for producing summaries of the output EZs; this code may be helpful for summarising other gridded sampling frames.

### Recommendations (as of 16/07/2019)

As the predefined_EZ_size functionality has been thoroughly tested, it is recommended that users take advantage of the predefined_EZ_size specification by setting this to = "small", "medium" or "large". 

Testing of the alternative process of setting EZ size—using the EZ_by_hh, target_hh_per_EZ, pop_per_hh, EZ_by_pop, target_pop_per_EZ, and max_cells_per_EZ arguments—is ongoing. 
