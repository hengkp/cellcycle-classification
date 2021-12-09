%# THIS CODE SHOWS HOW TO PLOT CELL CYCLE GRAPH.
%# 
%# WE COLLECT DATA AND CLASSIFY CELL CYCLE BASED ON NEGATIVE CONTROL OF EACH SAMPLE...
%# BY USING 'cellcyclestages' FUNCTION.
%# Ref. Kriengkrai Phongkitkarun. 2021. All rights reserved.

%# INITIAL PROGRAM BEFORE RUN
close all; clear; clc;
warning('off','all');

%# IMPORT CLUSTERING DATA .MAT FILE
fpath = '/Users/heng/Dropbox/Macbook Pro/Projects/Tooba project/CellCycle Example';
fname = 'mydata.mat';
load(fullfile(fpath,fname));

%# DEFINE FOR SAVING ALL GRAPH
is_save = true; % true or false
if is_save==1; mkdir(fullfile(fpath,'Figure')); end

%# DEFINE COLOR FOR EACH STAGE OF CELL CYCLE
%#   Start with G1, ES, LS, and G2M respectively.
color = [ 30 94 255 ; 43 252 90 ; 251 200 21 ; 238 49 0 ]./255;

%# DEFINE AXIS PARAMETERS
param_x = 'nuclei_integrated_dapi';
param_y = 'nuclei_integrated_edu';

%# DEFINE CONC (FOR LABELLING THE TITLE)
conc = [ 0 0.4 4 40 ];

%# LOOP PLOT BAR OBJECT COUNT EACH SAMPLE AND TREATMENT
for k = 1:size(platemap,1)
    for i = 1:size(platemap,2)
        x = S{k,i}(:,ismember(S_colname,param_x)); % Select column for "X" axis
        y = S{k,i}(:,ismember(S_colname,param_y)); % Select column for "Y" axis
        z = []; % Select column for "Size"
        c = S{k,i}(:,end); % Select column 'cellcycle_4stage' for "Color"
        % Define optinal parameters
        opt.title_name = char(64+platemap{k,i}(1)) + string(platemap{k,i}(2)); % Define as 'wellname'
        opt.title_name = sprintf('Clone %d , Conc %s',k,string(conc(i)));
        opt.colormap = color;
        fig = plotcc(x,y,x,c,opt);
        if is_save==1; saveas(fig,fullfile(fpath,'Figure',sprintf('CellCycle Plot - Clone %d - Conc %s',k, strrep(string(conc(i)),'.',','))),'png'); close all; end
        clear fig;
    end
    clear i x y z;
    
end
clear k;

function fig = plotcc(x,y,z,c,varargin)

g = struct();
g.colormap = [];
g.xlim = [];
g.ylim = [];
g.linewidth = 2;
g.fontsize = 32;
g.sc_linewidth = [];
g.sc_facecolor = [];
g.sc_facealpha = [];
g.sc_edgecolor = [];
g.sc_edgealpha = [];
g.title_name = '';
g.title_fontsize = 36;

if nargin > 4
    opt = varargin{1};
    if ~isstruct(opt)
        warning("""opt"" must be a structural variable.");
    end
    fieldname = fieldnames(opt);
    for i = 1:numel(fieldname)
        if ismember(fieldname{i},fieldnames(g))
            g.(fieldname{i}) = opt.(fieldname{i});
        end
    end
end

fig = figure('units','normalized','position',[0.2 0.25 0.4 0.5]); hold on;
sc = scatter(x,y,z,c,"filled");
if ~isempty(g.sc_facecolor); sc.MarkerFaceColor = g.sc_facecolor; end
if ~isempty(g.sc_facealpha); sc.MarkerFaceAlpha = g.sc_facealpha; end
if ~isempty(g.sc_edgecolor); sc.MarkerEdgeColor = g.sc_edgecolor; end
if ~isempty(g.sc_edgealpha); sc.MarkerEdgeAlpha = g.sc_edgealpha; end
if ~isempty(g.colormap); colormap(g.colormap); end
if ~isempty(g.xlim); xlim(g.xlim); end
if ~isempty(g.ylim); ylim(g.ylim); end
box on;
set(gca,'fontsize',g.fontsize,'linewidth',g.linewidth);
if ~isempty(g.title_name); title(g.title_name,'fontsize',g.title_fontsize); end
end