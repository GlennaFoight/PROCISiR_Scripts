README for PROCISiR scripts

For processing high throughput microscopy data for cell size parameters:

Run ADAPT (Barry, DJ, Durkin, CH, Abella, JV and Way, M. 2015. Open source software for quantification of cell migration, protrusions, and fluorescence intensities. J Cell Biol. 209: 163-180) on FIJI. Set Morphology statistics in FIJI to output desired cell parameters.

Use these scripts to process the output files:

1) ADAPTprocossing1.ipynb

Inputs: Documents//ADAPT/morphology.txt files
	Documents//ADAPT/*_Output/MeanFluorescenceIntensity.csv

Outputs: >*_Allmeans.csv: expt_well_site_Cell_# rows, with columns of means: Area, Perim., Circ., AR, Round, Solidity, GFPfluor


2) hist_from_dataframe.ipynb

Generates histograms of these statistics and percentiles to use for setting filtering cutoffs.

3) ADAPTprocessing2.ipynb

Runs on all input experiments/wells/sites for a given experiment batch
Input: Cutoffs for stats from previous script
	Documents//ADAPT/morphology.txt files
	Documents//ADAPT/*_Output/MeanFluorescenceIntensity.csv

Output: csvs by condition_stat with rows as frames, columns as cells

4) ADAPTplotting.ipynb

Run on full experiment. Plots means and s.e.m over time.

Variations on this script:

ADAPTplotting_DMSOsubtract.ipynb
	Subtracts signal from "condition 7" (DMSO) from others when plotting.

ADAPTplotting_KDE.ipynb
	Plots kernal denisty estimates for the change in statistics.

ADAPTplotting_KDE_AllExpts.ipynb
	Plots kernal density estimates of raw statistics.



For processing site saturation mutagenesis Illumina data: simple_countenrich.py
Script derived from one written by Jorgen Nelson

Usage:
python simple_countenrich.py [SSM name] [WT sequence] [unselected mapcounts file] [selected mapcounts file] 


For plotting time series data for membrane association by generating a 1-D interpolation of time series:
plotAverageFunction.ipynb








