%# THIS CODE WORKS WITH CELLPROFILER CSV OUTPUT FILE.
%# CELLPROFILER VERSION 4.3.1
%# WE COLLECT DATA AND CLASSIFY CELL CYCLE BASED ON NEGATIVE CONTROL OF EACH SAMPLE...
%# BY USING 'cellcyclestages' FUNCTION.
%# Ref. Kriengkrai Phongkitkarun. 2021. All rights reserved.

%# INITIAL PROGRAM BEFORE RUN
close all; clear; clc;
warning('off','all');

%# IMPORT CELLPROFILER OUTPUT FILE
fpath = '/Users/heng/Dropbox/Macbook Pro/Projects/Tooba project/CellCycle Example';
tbl = readtable(fullfile(fpath,'MyExpt_Nuclei_edited.csv'));

%# FILTER SELECTED COLUMNS AND CHANGE NAMES
columns = {'Metadata_Row','Metadata_Column','Metadata_Well','AreaShape_Area',...
    'Intensity_IntegratedIntensity_CorrSubDAPI','Intensity_IntegratedIntensity_CorrSubEDU',...
    'Intensity_IntegratedIntensity_CorrSubHoechst','Intensity_MeanIntensity_CorrSubCherry',...
    'Intensity_MeanIntensity_CorrSubVenus'};
tbl = tbl(:,columns);
tbl.Properties.VariableNames = {'rowid','columnid','well','nuclei_area',...
    'nuclei_integrated_dapi','nuclei_integrated_edu','nuclei_integrated_hoechst',...
    'nuclei_mean_cherry','nuclei_mean_venus'};

%# DEFINE PLATEMAP
%#   Write in matrix ; M{S,T} = [R,C];
%#   >> S:sample
%#   >> T:treatment
%#   >> R:row_id
%#   >> C:column_id
%#   Note! If there is replicate well, write each well in each row...
%#         as following example;
%#         platemap{1,1} = [ 1 1 ; 1 2 ; 1 3 ];
platemap{1,1} = [ 2 2 ]; % no TMP
platemap{1,2} = [ 2 3 ]; % conc 1
platemap{1,3} = [ 2 4 ]; % conc 2
platemap{1,4} = [ 2 5 ]; % conc 3

%# TRANSFORM TABLE DATA TO MATRIX DATA
%#   Write in matrix array ; M{S,T} = D;
%#   >> S:sample
%#   >> T:treatment
%#   >> D:data array without index column
Data = cell(size(platemap));
for i = 1:size(platemap,1)
    for j = 1:size(platemap,2)
        Dtemp = [];
        for k = 1:size(platemap{i,j},1) % Loop in case there is a replicate well.
            idx = tbl.rowid == platemap{i,j}(k,1) & ...
                tbl.columnid == platemap{i,j}(k,2);
            Dtemp = [ Dtemp ; table2array(tbl(idx,4:end)) ]; % The first 3 columns are an index. Start with column 4.
        end
        Data{i,j} = Dtemp;
    end
end
clear i j idx Dtemp;
Data_colname = tbl.Properties.VariableNames(4:end)'; % Collect the column name.

%# THIS CODE IS FOR SHOW ALL TREATMENT (DAPI&EDU PLOT)
%{
fig = figure('units','normalized','position',[ 0.1 0.1 0.8 0.4]);
ti = tiledlayout(size(platemap));

for i = 1:size(platemap,1)
    for j = 1:size(platemap,2)
        x = Data{i,j}(:,3);
        y = Data{i,j}(:,4);
        
        nexttile;
        sc = scatter(x,y,1,'filled');
        sc.MarkerFaceColor = 'k';
        
        axis([ 0 60 0 80 ]);
    end
end
%}

%# Next, we need to classify cellcycle using negative control of each sample...
%# and apply that threshold to all treatment.

%# DEFINE ALL NECESSARY VARIABLES
S = []; % New data matrix variable.
G = []; % Decoration setting variable.
params = [ 2 3 ]; % Use Nuclei Integrated Intensity of DAPI and EDU column.
opts = struct();
% opts.is_decor = true; % Use this in case there need to decorated data.

%# RUN CONTROL GROUP
for i = 1:size(platemap,1)
    [S{i,1},G{i,1}] = cellcyclestages('f',Data{i,1},params,opts);
end
clear i;

%# APPLY 'G' TO TREATMENT
for i = 1:size(platemap,1)
    for j = 2:size(platemap,2)
        [S{i,j},G{i,j}] = cellcyclestages('a',Data{i,j},params,G{i,1});
    end
end
clear i j;

%# SAVE ALL VARIABLES TO .MAT FILE
S_colname = [ Data_colname ; 'cellcycle_4stage' ];
save(fullfile(fpath,'mydata.mat'));
fprintf('Done.\n');
