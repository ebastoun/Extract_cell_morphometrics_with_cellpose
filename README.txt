This script is meant to analyse cell variables (area, perimeter, cell shape factor, circularity,...)
of cell monolayers from data generated in Fiji (plugin: Trackmate + Cellpose)
Created: 21/11/2023 by Raul Aparicio Yuste and Effie Bastounis


INPUT DATA: once you run Fiji->Plugins->Tracking (+cellpose), in the display options step, export spots to CSV.
To avoid format issues, save correctly the .CSV file as .XLSX (excel format).
Place your excel file with the spots in the same folder as the function "import_data_Trackmate_cellpose.m"
Example: we provide an example, see "my_sample.xlsx" file

"import_data_Trackmate_cellpose.m": we strongly recommend to run the code block by block to adapt features

Structure of the code: 

(1) Settings: adapt here your parameters

(2) Extract specific variables (columns) from excel file

(3) Process tracks by applying filters considering:
	(3.1) Track lenght
	(3.2) Set a ROI, do not take into account cells at the border

(4) Analyse frames [0,1,2....,N] over time, we apply a filter:
	(4.1) Set a ROI, do not take into account cells at the border

(5) Plots
	(5.1) Plot information from individual tracked cells. Select the id of cells to plot them
		IMPORTANT: choose cells that have the same track length, so we do not have problems
		of vector size when plotting
	(5.2) Plot mean values of variables over time

(6) Normalize data with respect to the first frame
	IMPORTANT: run this part of the code ONLY if all the tracked_cells have the same track length
	Example: if there are 50 frames in your data, choose min_track_length = 49 so you will get tracks
	with the same length and all the information over the whole time

(7) Compute displacements
	Given that we follow cells and we know their position, we can compute their displacement and velocities
	IMPORTANT: run this part of the code ONLY if all the tracked_cells have the same track length