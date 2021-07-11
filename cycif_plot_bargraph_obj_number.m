%# THIS CODE SHOWS HOW TO PLOT CELL CYCLE BAR GRAPH.
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

%# LOOP PLOT BAR OBJECT COUNT EACH SAMPLE AND TREATMENT
for k = 1:size(platemap,1)
    ccount = zeros(size(platemap,2),4); dcount = zeros(size(platemap,2),4);
    for i = 1:size(platemap,2)
        obj = S{k,i}(:,end); % Select column 'cellcycle_4stage' which is the last column
        ccount(i,:) = [ sum(obj==1) sum(obj==2) sum(obj==3) sum(obj==4) ]; % Raw object number
        dcount(i,:) = ccount(i,:)./numel(obj).*100; % Percent object number
    end
    clear i obj;
    x = [ 0.8 2 3.2 4.4 ];
    fig = figure('units','normalized','position',[0.2 0.25 0.3 0.5]); hold on;
    for i = 1:4
        br = bar(x(i),dcount(i,:),'stacked');
        for j = 1:4
            br(j).FaceColor = color(j,:);
            br(j).LineWidth = 3;
            br(j).BarWidth = 0.8;
        end
    end
    clear i j;
    xlim([ 0 5.2 ]);
    ylim([ 0 100 ]);
    xticks(x);
    xticklabels([0,0.4,4,40]);
    yticks(0:50:100);
    box on;
    set(gca,'fontsize',36,'linewidth',2);
    title(sprintf('Clone %d',k),'fontsize',42);
    if is_save==1; saveas(fig,fullfile(fpath,'Figure',sprintf('Bar Plot - Clone %d',k)),'png'); close all; end
    clear fig br ccount dcount;
end
clear k;