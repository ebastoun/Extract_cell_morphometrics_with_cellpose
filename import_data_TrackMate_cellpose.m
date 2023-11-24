%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is meant to analyse cell variables (area, perimeter, cell s
% hape factor, circularity,...) of cell monolayers from data generated in 
% Fiji (plugin: Trackmate + Cellpose)
% Created: 21/11/2023 by Raul Aparicio Yuste and Effie Bastounis

% INPUT DATA: once you run Fiji->Plugins->Tracking (+cellpose), in the 
% display options step, export spots to CSV. To avoid format issues, save 
% correctly the .CSV file as .XLSX (excel format).
% Place your excel file with the spots in the same folder as the function 
% "import_data_Trackmate_cellpose.m". Example: we provide an example, 
% see my_sample.xlsx file

%"import_data_Trackmate_cellpose.m": we strongly recommend to run the code 
% block by block to adapt features

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure of the code: 
 
% (1) Settings: adapt here your parameters
 
% (2) Extract specific variables (columns) from excel file
 
% (3) Process tracks by applying filters considering:
% 	(3.1) Track lenght
% 	(3.2) Set a ROI, do not take into account cells at the border
 
% (4) Analyse frames [0,1,2....,N] over time, we apply a filter:
% 	(4.1) Set a ROI, do not take into account cells at the border
 
% (5) Plots
% 	(5.1) Plot information from individual tracked cells. Select the id of cells to plot them
% 		IMPORTANT: choose cells that have the same track length, so we do not have problems
% 		of vector size when plotting
% 	(5.2) Plot mean values of variables over time
 
% (6) Normalize data with respect to the first frame
% 	IMPORTANT: run this part of the code ONLY if all the tracked_cells have the same track length
% 	Example: if there are 50 frames in your data, choose min_track_length = 49 so you will get tracks
% 	with the same length and all the information over the whole time
 
% (7) Compute displacements
% 	Given that we follow cells and we know their position, we can compute their displacement and velocities
% 	IMPORTANT: run this part of the code ONLY if all the tracked_cells have the same track length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

excel_name         = 'my_sample.xlsx'                    ; % Name of the .xlsx file
numHeaderRows      = 3                                   ; % Specify the number of header rows in your excel file - 1
labels             = readcell(excel_name, 'Range', '2:2'); % Get second row in the excel file (labels)
variableNames      = string(labels)                      ;

min_track_length   = 10                                  ; % Analyse cells that are tracked at least XX steps
monolayer_size     = 450.24                              ; % Monolayer size of your image (MICRONS) (assuming is square)
crop_border_length = 20                                  ; % Set a ROI to crop the borders to avoid the analysis of cells at the border (MICRONS)
delta_time         = 20                                  ; % delta time between frames (MIN)
total_numb_frames  = 11                                  ; % total number of frames


time_min           = 0:delta_time:total_numb_frames*delta_time-delta_time; % vector of minutes
time_hou           = time_min/60                                         ; % vector of hours


meta_data                          = readtable(excel_name, 'ReadVariableNames', true, 'HeaderLines', numHeaderRows);
meta_data.Properties.VariableNames = variableNames                                                                 ; % Rename table labels






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Select variables to extract (order of columns might differ from file 
% to file)
% Make sure you write the right name of the variable (names from second row
% in the excel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extract_names{1} = 'X'          ;  % Do not remove/overwrite this one
extract_names{2} = 'Y'          ;  % Do not remove/overwrite this one
extract_names{3} = 'Frame'      ;  % Do not remove/overwrite this one
extract_names{4} = 'Track ID'   ;  % Do not remove/overwrite this one
extract_names{5} = 'Area'       ;
extract_names{6} = 'Shape index';
extract_names{7} = 'Perimeter'  ;
extract_names{8} = 'Circularity';

string_pos = zeros(1,length(extract_names));

for i = 1:length(extract_names)

    look_for_string = extract_names{i}                            ;
    string_pos(i)   = find(strcmp(variableNames, look_for_string));
    
    % Check if the label was found
    if isempty(string_pos)
        error(['String "', look_for_string, '" not found.']);
    end

end

% Assign variables, keep the same order as before. UPDATE it if you add new
% variables
X            = meta_data(:,string_pos(1)      )   ; 
Y            = meta_data(:,string_pos(2)      )   ;
frames       = meta_data(:,string_pos(3)      )   ; 
track_id     = meta_data(:,string_pos(4)      )   ; 
area         = meta_data(:,string_pos(5)      )   ;
shape_fact   = meta_data(:,string_pos(6)      )   ;
perimeter    = meta_data(:,string_pos(7)      )   ;
circularity  = meta_data(:,string_pos(8)      )   ;
clear meta_data

% Convert data table to vectors of doubles and clean table with tons of
% data for the analysis. UPDATE it if you add new variables
X            = double(table2array(X)          )   ;
Y            = double(table2array(Y)          )   ; 
frames       = double(table2array(frames)     )   ; 
track_id     = double(table2array(track_id)   )   ; 
area         = double(table2array(area)       )   ;
shape_fact   = double(table2array(shape_fact) )   ;
perimeter    = double(table2array(perimeter)  )   ;
circularity  = double(table2array(circularity))   ;

% Compute min and max track ID number
min_ID       = min(track_id)                      ;
max_ID       = max(track_id)                      ;
% Initialize the number of cells that will be tracked in the analysis
cell_counter = 0                                  ;

% Define region of interest to analyse cells
x_min        = crop_border_length                 ;
y_min        = crop_border_length                 ;
x_max        = monolayer_size - crop_border_length;
y_max        = monolayer_size - crop_border_length;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Extract data from tracks filtering:
%       (3.1) track length
%       (3.2) cells located at the border
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run over all the tracks
for k = min_ID:max_ID

    % Find positions assigned to that track
    pos     = find(track_id==k);
    
    %......................................................................
    % Filter track lenght
    %......................................................................
    if length(pos) > min_track_length

        %..................................................................
        % Filter cells at the border
        %..................................................................
        x_local   = X(pos);
        y_local   = Y(pos);
        insideROI = (x_local > x_min) & (x_local < x_max) & ...
                    (y_local > y_min) & (y_local < y_max);
        
        if any(insideROI == 0)
            % Don't do anything with cells at the border
        else
            % Analyse cells
            cell_counter                            = cell_counter + 1 ;

            tracked_cells(cell_counter).area        =        area(pos); % area um^2
            tracked_cells(cell_counter).shape_fact  =  shape_fact(pos); % perimeter to sqrt of area
            tracked_cells(cell_counter).perimeter   =   perimeter(pos); % perimeter um
            tracked_cells(cell_counter).circularity = circularity(pos); % 4*pi*area/perimeter^2 
            tracked_cells(cell_counter).track_id    =                k;

            % Save coordinates to compute displacements and velocity
            % vectors
            tracked_cells(cell_counter).x_pos       =          x_local;
            tracked_cells(cell_counter).y_pos       =          y_local;

        end             
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4) Analyze each frame over time filtering:
%       (4.1) Cells located at the border
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of cells tracked
tracks_analyzed        = cell_counter      ;
max_frame              = max(frames)       ;

% Initialize structures
mean_cell_area_time    = zeros(1,max_frame);
mean_shape_factor_time = zeros(1,max_frame);
mean_perimeter_time    = zeros(1,max_frame);
mean_circularity_time  = zeros(1,max_frame);

% Compute mean over time
for t_i = 0:max_frame

    pos = find(frames == t_i);

    %......................................................................
    % Filter cells at the border
    %......................................................................
    x_local      = X(pos);
    y_local      = Y(pos);
    insideROI    = (x_local > x_min) & (x_local < x_max) & ...
                   (y_local > y_min) & (y_local < y_max);

    % Retrieve cells within the ROI
    cells_in     = find(insideROI==1);
    cells_in_roi = pos(cells_in)     ;
   
    % (t_i+1): arrays cannot start from 0
    mean_cell_area_time   (t_i+1)  =   mean(        area(cells_in_roi)); % area um^2
    mean_shape_factor_time(t_i+1)  =   mean(  shape_fact(cells_in_roi)); % perimeter to sqrt of area
    mean_perimeter_time   (t_i+1)  =   mean(   perimeter(cells_in_roi)); % perimeter um
    mean_circularity_time (t_i+1)  =   mean( circularity(cells_in_roi)); % 4*pi*area/perimeter^2 

end





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5) Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%..........................................................................
% (5.1) Plot individual cells
%..........................................................................
% Select id tracked cells to plot
% IMPORTANT: 
%   (1) select cells that are tracked the same number of frames
%   (2) if you choose cells that are not followed over all the time frames,
%       adapt accordingly the time vector when plotting (time_hou)
cell_a = 68;
cell_b = 63;
cell_c = 9;
cell_d = 11;
%.......................
% Area
%.......................
fig = figure(1);
hold on;
plot(time_hou',tracked_cells(cell_a).area,'r');  
plot(time_hou',tracked_cells(cell_b).area,'g'); 
plot(time_hou',tracked_cells(cell_c).area,'b');
plot(time_hou',tracked_cells(cell_d).area,'m'); 
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Cell area (um)^2','Fontsize',20,'Color','k');
legend(['Cell ' num2str(cell_a)],['Cell ' num2str(cell_b)], ...
       ['Cell ' num2str(cell_c)],['Cell ' num2str(cell_d)]);
%axis([0 3.5 200 550])
hold off
Save_figure = getframe(fig); 
imwrite(Save_figure.cdata, 'area_indiv_cells.tif'); 
close all

%.......................
% Shape factor
%.......................
fig = figure(1);
hold on;
plot(time_hou',tracked_cells(cell_a).shape_fact,'r');  
plot(time_hou',tracked_cells(cell_b).shape_fact,'g'); 
plot(time_hou',tracked_cells(cell_c).shape_fact,'b');
plot(time_hou',tracked_cells(cell_d).shape_fact,'m'); 
xlabel('Time (h)','Fontsize',18);
ylabel('Shape factor (P/sqrt(A))','Fontsize',18);
legend(['Cell ' num2str(cell_a)],['Cell ' num2str(cell_b)], ...
       ['Cell ' num2str(cell_c)],['Cell ' num2str(cell_d)]);
set(gca,'Fontsize',18); grid on;
%axis([0 3.5 200 550])
hold off
Save_figure = getframe(fig); 
imwrite(Save_figure.cdata, 'shape_f_indiv_cells.tif'); 
close all

%.......................
% Perimeter
%.......................
fig = figure(1);
hold on;
plot(time_hou',tracked_cells(cell_a).perimeter,'r');  
plot(time_hou',tracked_cells(cell_b).perimeter,'g'); 
plot(time_hou',tracked_cells(cell_c).perimeter,'b');
plot(time_hou',tracked_cells(cell_d).perimeter,'m');
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Perimeter (um)','Fontsize',20,'Color','k');
legend(['Cell ' num2str(cell_a)],['Cell ' num2str(cell_b)], ...
       ['Cell ' num2str(cell_c)],['Cell ' num2str(cell_d)]);
hold off
Save_figure = getframe(fig); 
imwrite(Save_figure.cdata, 'perim_indiv_cells.tif'); 
close all

%.......................
% Circularity
%.......................
fig = figure(1);
hold on;
plot(time_hou',tracked_cells(cell_a).circularity,'r');  
plot(time_hou',tracked_cells(cell_b).circularity,'g'); 
plot(time_hou',tracked_cells(cell_c).circularity,'b');
plot(time_hou',tracked_cells(cell_d).circularity,'m'); 
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Circularity (um)','Fontsize',20,'Color','k');
legend(['Cell ' num2str(cell_a)],['Cell ' num2str(cell_b)], ...
       ['Cell ' num2str(cell_c)],['Cell ' num2str(cell_d)]);
hold off
Save_figure = getframe(fig); 
imwrite(Save_figure.cdata, 'circul_indiv_cells.tif'); 
close all





%..........................................................................
% (5.2) Plot mean values over time
%..........................................................................
% Area
fig = figure(1);
plot(time_hou',mean_cell_area_time,'r'); 
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Cell area (um)^2','Fontsize',20,'Color','k');
Save_figure = getframe(fig); 
imwrite(Save_figure.cdata, 'mean_area.tif'); 

% Shape factor
fig = figure(1);
plot(time_hou',mean_shape_factor_time,'r'); 
xlabel('Time (h)','Fontsize',18);
ylabel('Shape factor (P/sqrt(A))','Fontsize',18);
set(gca,'Fontsize',18); grid on;
Save_figure = getframe(fig); 
imwrite(Save_figure.cdata, 'mean_shape_fact.tif');

% Perimeter
fig = figure(1);
plot(time_hou',mean_perimeter_time,'r'); 
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Perimeter (um)','Fontsize',20,'Color','k');
Save_figure = getframe(fig); 
imwrite(Save_figure.cdata, 'mean_perimeter.tif'); 

% Circularity
fig = figure(1);
plot(time_hou',mean_circularity_time,'r'); 
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Circularity (um)','Fontsize',20,'Color','k');
Save_figure = getframe(fig); 
imwrite(Save_figure.cdata, 'mean_circularity.tif'); 




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6) Normalize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT: run this part of the code ONLY if all the tracked cells have
% the same length. You need to filter them, for example, by choosing on the
% settings before:
% min_track_length = your max number of frames -1


% .........................................................................
% Area 
% .........................................................................
% Normalized area with respect to the first frame
for k=1:tracks_analyzed
    tracked_cells(k).arean = tracked_cells(k).area  ./   tracked_cells(k).area(1);
end

% Extract values from tracked_cells and make a matrix so you can calculate 
% mean, SD
cell_area_matrix = zeros(tracks_analyzed,total_numb_frames);
for k=1:tracks_analyzed
    cell_area_matrix(k,:) = tracked_cells(k).arean;
end
plot(mean(cell_area_matrix(:,:),1))
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Normalized cell area','Fontsize',20,'Color','k');
saveas(gcf, 'norm_area.tif');


% .........................................................................
% Shape factor 
% .........................................................................
% Normalized shape factor with respect to the first frame
for k=1:tracks_analyzed
    tracked_cells(k).shapefn = tracked_cells(k).shape_fact  ./   tracked_cells(k).shape_fact(1);
end

% Extract values from tracked_cells and make a matrix so you can calculate 
% mean, SD
cell_shapef_matrix = zeros(tracks_analyzed,total_numb_frames);
for k=1:tracks_analyzed
    cell_shapef_matrix(k,:) = tracked_cells(k).shapefn;
end
plot(mean(cell_shapef_matrix(:,:),1))
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Normalized cell shape factor','Fontsize',20,'Color','k');
saveas(gcf, 'norm_shapef.tif');


% .........................................................................
% Perimeter
% .........................................................................
% Normalized perimeter with respect to the first frame
for k=1:tracks_analyzed
    tracked_cells(k).perimn = tracked_cells(k).perimeter  ./   tracked_cells(k).perimeter(1);
end

% Extract values from tracked_cells and make a matrix so you can calculate 
% mean, SD
cell_perim_matrix = zeros(tracks_analyzed,total_numb_frames);
for k=1:tracks_analyzed
    cell_perim_matrix(k,:) = tracked_cells(k).perimn;
end
plot(mean(cell_perim_matrix(:,:),1))
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Normalized perimeter','Fontsize',20,'Color','k');
saveas(gcf, 'norm_perimeter.tif');


% .........................................................................
% Circularity
% .........................................................................
% Normalized circularity with respect to the first frame
for k=1:tracks_analyzed
    tracked_cells(k).circn = tracked_cells(k).circularity  ./   tracked_cells(k).circularity(1);
end

% Extract values from tracked_cells and make a matrix so you can calculate 
% mean, SD
cell_circ_matrix = zeros(tracks_analyzed,total_numb_frames);
for k=1:tracks_analyzed
    cell_circ_matrix(k,:) = tracked_cells(k).circn;
end
plot(mean(cell_circ_matrix(:,:),1))
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Normalized circularity','Fontsize',20,'Color','k');
saveas(gcf, 'norm_circularity.tif');
close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (7) Compute displacement and velocity vector for each track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT: run this part of the code ONLY if all the tracked cells have
% the same length. You need to filter them, for example, by choosing on the
% settings before:
% min_track_length = your max number of frames -1

% Another alternative, you could adapt the code to get the displacement or
% the velocity of some specific cells


% Initialize structures
displacement_magn = zeros(total_numb_frames-1,tracks_analyzed);
velocity_magn     = zeros(total_numb_frames-1,tracks_analyzed);

for k = 1:tracks_analyzed

    x_track                = tracked_cells(k).x_pos                     ;
    y_track                = tracked_cells(k).y_pos                     ;

    displacement_x         = diff(x_track)                              ;
    displacement_y         = diff(y_track)                              ;
    displacement           = sqrt(displacement_x.^2 + displacement_y.^2);
    displacement_magn(:,k) = displacement                               ;
    velocity_magn(:,k)     = displacement/delta_time                    ;

end

% Assuming all tracks belong to the same frame time (that's one of the
% reasons we have filtered the tracks as suggested), compute mean
% displacement and velocity

mean_displacement  =  mean(displacement_magn,2);
mean_velocity      =  mean(velocity_magn,2)    ;

plot(time_hou(2:end),mean_displacement)
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Mean cell displacement (um)','Fontsize',20,'Color','k');
title('Cell displacement from cellpose');
saveas(gcf, 'displacement.tif');

plot(time_hou(2:end),mean_velocity)
xlabel('Time (h)','Fontsize',20,'Color','k');
ylabel('Mean cell velocity (um/h)','Fontsize',20,'Color','k');
title('Cell velocity from cellpose');
saveas(gcf, 'velocity.tif');
close all



