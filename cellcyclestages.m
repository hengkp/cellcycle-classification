function [B,varargout] = cellcyclestages(mode,A,params,varargin)
%cellcyclestages(mode,A,params,opts) is a function for calculate stage of
%cell cycle in each cell.
%
% INPUT PARAMETERS:
%      mode   : character string
%          >>> 'find'  or 'f' for find cutpoint
%          >>> 'apply' or 'a' for apply parameter
%      A      : matrix of data
%      params : column id for classify stage of cell cycle e.g. [ 1 2 ]
%          >>> Nuclei_Integrated_DAPI
%          >>> Nuclei_Integrated_EDU
%          >>> Nuclei_Area (Optional)
%      opt    : optional structure
%          >>> outliermode : 'percentile' (default) or 'manual' or '' for
%          pass this process
%          >>> percentile  : [ 5 95 5 95 ] (default)
%          >>> lowerbound  : [ min_x min_y max_x min_y ] (default is [])
%          >>> is_decor    : false (default) or true for decoration data in
%          find cutpoint mode
%          >>> decor       : [x,y] decoration bound datapoint
%          >>> cpmode      : 'auto4(a4)' or 'auto3(a3)' or 'auto2s(a2s)' or 'manual(m)'
%          (default is 'manual')
%          >>> sharedcov   : false (default) or true
%          >>> gmfit4      : gmdistribution for 'auto4' cpmode
%          >>> gmfit3      : gmdistribution for 'auto3' cpmode
%          >>> cp          : cutpoint for 'semi-auto' cpmode
%          matrix of [ x-value , y-value ] 
%          >>> cpm         : manual cutpoint for 'manual' cpmode
%          matrix of [ x1 , y1 ; x2 , y2 ; x3 , y3 ]
%
% OUTPUT PARAMETERS:
%      B : matrix (A) with column stage of cell cycle
%      g : structure of 'opt' parameter contains cutpoint information.
%
%Ref. Kriengkrai Phongkitkarun. 2021. All rights reserved.


g = struct();
g.percentile = [ 1 99 1 99 ]; %# percentile of [ xmin xmax ymin ymax ]
g.outliermode = 'percentile';
g.outlier = [];
g.lowerbound = []; %# two data points [ x1, y1, x2, y2 ]
g.is_decor = false;
g.decor = [];
g.cpmode = 'manual'; %# cut point mode : '' , 'partial' , 'manual'
g.sharedcov = false;
g.gmfit4 = [];
g.gmfit3 = [];
g.cp = [];
g.cpm = [];

if nargin == 4
    opts = varargin{1};
    fns = fieldnames(opts);
    for i = 1:numel(fns)
        if ismember(fns{i},fieldnames(g))
            g.(fns{i}) = opts.(fns{i});
        end
    end
    clear fns i opts;
end

%# RUNNING CODE ...
switch mode
    case {'find','f'} %# Find optimize parameter for cut point
        
        %# Filter Outlier
        if ~isempty(g.outliermode)
            while true
                A1 = A(:,params(1));
                A2 = A(:,params(2));
                switch g.outliermode
                    case 'percentile'
                        [~,idx1] = rmoutliers(A1,'percentiles',g.percentile(1:2));
                        [~,idx2] = rmoutliers(A2,'percentiles',g.percentile(3:4));
                        idx = ~(idx1|idx2);
                    case 'manual'
                        if ~isempty(g.outlier)
                            idx = A1 > g.outlier(1) & A1 <= g.outlier(2) & ...
                                A2 > g.outlier(3) & A2 <= g.outlier(4);
                        else
                            idx = zeros(numel(A1),1);
                        end
                end
                
                close all;
                fig = figure(); hold on;
                sc = scatter(A1,A2,35,idx,'filled');
                sc.MarkerFaceAlpha = 0.5;
                colormap([ 188 185 173 ; 38 71 255 ]./255);
                switch g.outliermode
                    case 'percentile'
                        title(sprintf('Current percentile is %.2f , %.2f , %.2f , %.2f', ...
                            g.percentile(1),g.percentile(2),g.percentile(3),g.percentile(4)));
                    case 'manual'
                        title('Current outlier mode is ''manual''.');
                end
                
                satisfy = input('Satisfy? (y/[n]): ','s');
                if ismember(satisfy,{'y','Y'})
                    A = A(idx,:);
                    g.outlier = [ min(A(:,params(1))) ...
                        max(A(:,params(1))) ...
                        min(A(:,params(2))) ...
                        max(A(:,params(2))) ];
                    break;
                else
                    while true
                        fprintf('Your current mode is ''%s''. Press Enter for continue on same mode or \nspecify new mode by types: percentile(p), manual(m).\n',g.outliermode);
                        changemode = input('Change CP mode? (p/m): ','s');
                        if ismember(changemode,{'percentile','p'})
                            g.outliermode = 'percentile';
                        elseif ismember(changemode,{'manual','m'})
                            g.outliermode = 'manual';
                        end
                        clear changemode;
                        switch g.outliermode
                            case 'percentile'
                                fprintf('Specify percentile in x-axis and y-axis by using ''space'' between value.\n');
                                str = input('New value: ','s');
                                newval = str2double(split(str));
                                if size(newval,1) == 4
                                    g.percentile = newval';
                                    break;
                                end
                            case 'manual'
                                xo = get(gca,'xlim');
                                yo = get(gca,'ylim');
                                sc.MarkerFaceColor = [ 38 71 255 ]./255;
                                title('Specify first data points for min outlier in x-axis.');
                                [x1,~] = ginput(1);
                                li = line([x1 x1],yo,'color','k','linewidth',2);
                                title('Specify second data points for max outler in x-axis.');
                                [x2,~] = ginput(1);
                                li = line([x2 x2],yo,'color','k','linewidth',2);
                                title('Specify third data points for min outlier in y-axis.');
                                [~,y1] = ginput(1);
                                li = line(xo,[y1 y1],'color','k','linewidth',2);
                                title('Specify fourth data points for max outler in y-axis.');
                                [~,y2] = ginput(1);
                                li = line(xo,[y2 y2],'color','k','linewidth',2);
                                g.outlier = [ x1 x2 y1 y2 ];
                                clear li x1 x2 y1 y2 xo yo;
                                break;
                        end
                    end
                end
            end
            clear A1 A2 idx1 idx2 idx fig sc satisfy str newval;
            close all;
        end
        
        %# Filter lower bound of G1 and G2M
        if ~isempty(g.outlier)
            while true
                
                A1 = A(:,params(1));
                A2 = A(:,params(2));
                min1 = min(A1); max1 = max(A1);
                min2 = min(A2);
                lbflag = false;
                if isempty(g.lowerbound)
                    lbflag = true;
                    g.lowerbound = [ min1 min2 max1 min2 ];
                end
                ymb = (g.lowerbound(4)-g.lowerbound(2))/(g.lowerbound(3)-g.lowerbound(1));
                ycb = g.lowerbound(2) - ymb*g.lowerbound(1);
                ybp.xmin = min1;
                ybp.ymb = ymb;
                ybp.yminb = g.lowerbound(2);
                ybfunc = @(x,y,p) y > p.ymb*(x-p.xmin)+p.yminb;
                idx = ybfunc(A1,A2,ybp);
                
                close all;
                fig = figure(); hold on;
                sc = scatter(A1,A2,35,idx,'filled');
                sc.MarkerFaceAlpha = 0.5;
                colormap([ 188 185 173 ; 38 71 255 ]./255);
                pl = plot([min1,max1],ymb.*[min1,max1]+ycb,'k--','linewidth',2);
                title(sprintf('Current lower bound is %.2f , %.2f , %.2f , %.2f', ...
                    g.lowerbound(1),g.lowerbound(2),g.lowerbound(3),g.lowerbound(4)));
                
                if lbflag == 1
                    g.lowerbound = [];
                end
                
                satisfy = input('Satisfy? (y/[n]): ','s');
                if ismember(satisfy,{'y','Y'})
                    A = A(idx,:);
                    break;
                else
                    delete(pl);
                    sc.MarkerFaceColor = [ 38 71 255 ]./255;
                    title('Specify first data points for lower bound in y-axis.');
                    [x1,y1] = ginput(1);
                    pl = plot(x1,y1,'k.','markersize',40);
                    title('Specify second data points for lower bound in y-axis.');
                    [x2,y2] = ginput(1);
                    pl = plot(x2,y2,'k.','markersize',40);
                    
                    g.lowerbound = [ x1 y1 x2 y2 ];
                    clear pl x1 y1 x2 y2;
                end
            end
            clear A1 A2 min1 max1 min2 ymb ycb ybp ybfunc idx fig sc pl satisfy pl1 pl2 lbflag; 
            close all;
        end
        
        %# Manual Decoration
        if g.is_decor == 1
            
            if ~isempty(g.decor)
                cont = input('The decoration setting is found!\nDo you want to re-decorated data? (y/[n]): ','s');
            else
                cont = 'y';
            end
            
            if ismember(cont,{'y','Y','yes'})
                fig = figure(); hold on;
                sc = scatter(A(:,params(1)),A(:,params(2)),50,'filled');
                sc.MarkerFaceColor = [ 117 156 255 ]./255;
                sc.MarkerEdgeColor = [ 30 94 255 ]./255;
                title('Specify decoration point.');
                
                x = []; y = []; j = 0; pl = plot([]);
                while true
                    j = j + 1;
                    [xg,yg,button] = ginput(1);
                    if isempty(button)
                        break;
                    end
                    x(j,1) = xg;
                    y(j,1) = yg;
                    clear xg yg button;
                    delete(pl);
                    pl = plot(x,y,'k.-');
                    pl.LineWidth = 3;
                end
                if ~isempty(x)
                    x(end) = x(1); y(end) = y(1);
                    g.decor = [ x y ];
                    clear x y;
                end
                close all;
                clear fig sc x y j pl;
                
                idx = inpolygon(A(:,params(1)),A(:,params(2)),g.decor(:,1),g.decor(:,2));
                A = A(idx,:);
            else
                idx = inpolygon(A(:,params(1)),A(:,params(2)),g.decor(:,1),g.decor(:,2));
                A = A(idx,:);
            end
            clear cont;
        end
        
        %# Find Cut Point
        while true
            
            A1 = A(:,params(1));
            A2 = A(:,params(2));
            
            switch g.cpmode
                case {'auto4','a4'}
                    
                    gkmean = kmeans(A(:,params),4);
                    gmean = [ mean(A(gkmean==1,params(1))) mean(A(gkmean==2,params(1)))  mean(A(gkmean==3,params(1)))  mean(A(gkmean==4,params(1))) ]';
                    [gmean,idg] = sortrows(gmean);
                    gkmean = 1.*(gkmean==idg(1)) + 2.*(gkmean==idg(2)) + 3.*(gkmean==idg(3)) + 4.*(gkmean==idg(4));
                    g.gmfit4 = fitgmdist(A(:,params),4,'CovarianceType','full',...
                        'SharedCovariance',g.sharedcov,'start',gkmean);
                    gcluster = cluster(g.gmfit4,A(:,params));
                    gmean = [ mean(A(gcluster==1,1)) mean(A(gcluster==2,1)) mean(A(gcluster==3,1)) mean(A(gcluster==4,1)) ]';
                    [gmean,idg] = sortrows(gmean);
                    Acp = 1.*(gcluster==idg(1)) + 2.*(gcluster==idg(2)) + 3.*(gcluster==idg(3)) + 4.*(gcluster==idg(4));
                    
                    clear gkmean gmean idg gcluster;
                    
                case {'auto3','a3'}
                    
                    gkmean = kmeans(A(:,params),3);
                    gmean = [ mean(A(gkmean==1,1)) mean(A(gkmean==2,1)) mean(A(gkmean==3,1)) ]';
                    [gmean,idg] = sortrows(gmean);
                    gkmean = 1.*(gkmean==idg(1)) + 2.*(gkmean==idg(2)) + 3.*(gkmean==idg(3));
                    g.gmfit3 = fitgmdist(A(:,1:2),3,'CovarianceType','full',...
                        'SharedCovariance',g.sharedcov,'start',gkmean);
                    gcluster = cluster(g.gmfit3,A(:,params));
                    gmean = [ mean(A(gcluster==1,1)) mean(A(gcluster==2,1)) mean(A(gcluster==3,1)) ]';
                    [gmean,idg] = sortrows(gmean);
                    Acp = 1.*(gcluster==idg(1)) + 2.*(gcluster==idg(2)) + 3.*(gcluster==idg(3));
                    
                    clear gkmean gmean idg gcluster;
                    
                case {'auto2s','a2s'}
                    
                    gmfit1 = fitgmdist(A1,2,'CovarianceType','full',...
                        'SharedCovariance',g.sharedcov);
                    gmfit2 = fitgmdist(A2,2,'CovarianceType','full',...
                        'SharedCovariance',g.sharedcov);
                    gcluster1 = cluster(gmfit1,A1);
                    gcluster2 = cluster(gmfit2,A2);
                    cp1 = mean([ max([ min(A1(gcluster1==1)) min(A1(gcluster1==2)) ]) ...
                        min([ max(A1(gcluster1==1)) max(A1(gcluster1==2)) ]) ]);
                    
                    cp2 = mean([ max([ min(A2(gcluster2==1)) min(A2(gcluster2==2)) ]) ...
                        min([ max(A2(gcluster2==1)) max(A2(gcluster2==2)) ]) ]);
                    
                    g.cp = [ cp1 cp2 ];
                    
                    Acp = 1 .* (A1 <= cp1 & A2 <= cp2) + ...
                          2 .* (A1 <= cp1 & A2 >  cp2) + ...
                          3 .* (A1 >  cp1 & A2 >  cp2) + ...
                          4 .* (A1 >  cp1 & A2 <= cp2) ;
                    
                    clear gmfit1 gmfit2 gcluster1 gcluster2;
                    
                case {'manual','m'}
                    
                    close all;
                    fig = figure(); hold on;
                    sc = scatter(A1,A2,35,'filled');
                    sc.MarkerFaceAlpha = 0.5;
                    sc.MarkerFaceColor = [ 38 71 255 ]./255;
                    
                    title('Specify first data points for CP in X-axis.');
                    [x1,y1] = ginput(1);
                    pl = plot(x1,y1,'k.','markersize',40);
                    pl = plot([x1,x1],get(gca,'ylim'),'k--','linewidth',2);
                    
                    title('Specify first data points for CP in Y-axis.');
                    [x2,y2] = ginput(1);
                    pl = plot(x2,y2,'k.','markersize',40);
                    
                    title('Specify second data points for CP in Y-axis.');
                    [x3,y3] = ginput(1);
                    pl = plot(x3,y3,'k.','markersize',40);
                    
                    min1 = min(A1);
                    max1 = max(A1);
                    ym = (y3-y2)/(x3-x2);
                    yc = y2 - ym*x2;
                    ypr.xmin = min1;
                    ypr.ym = ym;
                    ypr.ycp1 = y2;
                    ycpfunc = @(x,y,p) y > p.ym*(x-p.xmin)+p.ycp1;
                    
                    g.cpm = [ x1 y1 ; x2 y2 ; x3 y3 ];
                    
                    pl = plot([min1,max1],ym.*[min1,max1]+yc,'k--','linewidth',2);
                    pause(2);
                    
                    Acp = 1.*(A1 <= x1 & ~ycpfunc(A1,A2,ypr)) + ...
                          2.*(A1 <= x1 &  ycpfunc(A1,A2,ypr)) + ...
                          3.*(A1 >  x1 &  ycpfunc(A1,A2,ypr)) + ...
                          4.*(A1 >  x1 & ~ycpfunc(A1,A2,ypr));
            end
            
            close all;
            fig = figure(); hold on;
            sc = scatter(A1,A2,35,Acp,'filled');
            sc.MarkerFaceAlpha = 0.5;
            colormap([ 30 94 255 ; 43 252 90 ; 251 200 21 ; 238 49 0 ]./255);
            
            satisfy = input('Satisfy? (y/[n]): ','s');
            if ismember(satisfy,{'y','Y'})
                A = [ A Acp ];
                break;
            else
                fprintf('Your current mode is ''%s''. Press Enter for continue on same mode or \nspecify new mode by types: auto(a), semi-auto(s), manual(m).\n',g.cpmode);
                changemode = input('Change CP mode? (a4/a3/s/m): ','s');
                if ismember(changemode,{'auto4','a4'})
                    g.cpmode = 'auto4';
                elseif ismember(changemode,{'auto2s','a2s'})
                    g.cpmode = 'auto2s';
                elseif ismember(changemode,{'manual','m'})
                    g.cpmode = 'manual';
                end
                clear changemode;
            end
            clear A1 A2 Acp fig sc;
        end
        
    case {'apply','a'} %# Apply parameter to data
        
        %# Apply Outlier
        if ~isempty(g.outlier) && numel(g.outlier) == 4
            A1 = A(:,params(1));
            A2 = A(:,params(2));
            idx = A1 >  g.outlier(1) & ...
                A1 <= g.outlier(2) & ...
                A2 >  g.outlier(3) & ...
                A2 <= g.outlier(4) ;
            A = A(idx,:); clear idx A1 A2;
        end
        
        %# Apply Lower bound
        if ~isempty(g.lowerbound) && numel(g.lowerbound) == 4
            A1 = A(:,params(1));
            A2 = A(:,params(2));
            min1 = min(A1); max1 = max(A1);
            ymb = (g.lowerbound(4)-g.lowerbound(2))/(g.lowerbound(3)-g.lowerbound(1));
            ycb = g.lowerbound(2) - ymb*g.lowerbound(1);
            ybp.xmin = min1;
            ybp.ymb = ymb;
            ybp.yminb = g.lowerbound(2);
            ybfunc = @(x,y,p) y > p.ymb*(x-p.xmin)+p.yminb;
            idx = ybfunc(A1,A2,ybp);
            A = A(idx,:); clear idx A1 A2 min1 max1 ymb ycb ybp ybfunc;
        end
        
        %# Apply decoration
        if ~isempty(g.decor)
            idx = inpolygon(A(:,params(1)),A(:,params(2)),g.decor(:,1),g.decor(:,2));
            A = A(idx,:); clear idx;
        end
        
        %# Apply Cut Point
        Acp = [];
        switch g.cpmode
            case {'auto4','a4'}
                if ~isempty(g.gmfit4)
                    try
                        gcluster = cluster(g.gmfit4,A(:,params));
                    catch
                        error('''g.gmfit4'' is not in gmdistribution format.\nPlease re-calculate it or try other method.\n');
                    end
                    gmean = [ mean(A(gcluster==1,1)) mean(A(gcluster==2,1)) mean(A(gcluster==3,1)) mean(A(gcluster==4,1)) ]';
                    [gmean,idg] = sortrows(gmean);
                    Acp = 1.*(gcluster==idg(1)) + 2.*(gcluster==idg(2)) + 3.*(gcluster==idg(3)) + 4.*(gcluster==idg(4));
                    clear gcluster gmean idg;
                end
                
            case {'auto3','a3'}
                if ~isempty(g.gmfit3)
                    try
                        gcluster = cluster(g.gmfit3,A(:,params));
                    catch
                        error('''g.gmfit3'' is not in gmdistribution format.\nPlease re-calculate it or try other method.\n');
                    end
                    gmean = [ mean(A(gcluster==1,1)) mean(A(gcluster==2,1)) mean(A(gcluster==3,1)) ]';
                    [gmean,idg] = sortrows(gmean);
                    Acp = 1.*(gcluster==idg(1)) + 2.*(gcluster==idg(2)) + 3.*(gcluster==idg(3));
                    clear gcluster gmean idg;
                end
                
            case {'auto2s','a2s'}
                if ~isempty(g.cp)
                    cp1 = g.cp(1);
                    cp2 = g.cp(2);
                    Acp = 1 .* (A1 <= cp1 & A2 <= cp2) + ...
                          2 .* (A1 <= cp1 & A2 >  cp2) + ...
                          3 .* (A1 >  cp1 & A2 >  cp2) + ...
                          4 .* (A1 >  cp1 & A2 <= cp2) ;
                    clear cp1 cp2;
                end
                
            case {'manual','m'}
                if ~isempty(g.cpm) && size(g.cpm,1) == 3 && size(g.cpm,2) == 2
                    A1 = A(:,params(1));
                    A2 = A(:,params(2));
                    x1 = g.cpm(1,1); y1 = g.cpm(1,2);
                    x2 = g.cpm(2,1); y2 = g.cpm(2,2);
                    x3 = g.cpm(3,1); y3 = g.cpm(3,2);
                    min1 = min(A1);
                    max1 = max(A1);
                    ym = (y3-y2)/(x3-x2);
                    yc = y2 - ym*x2;
                    ypr.xmin = min1;
                    ypr.ym = ym;
                    ypr.ycp1 = y2;
                    ycpfunc = @(x,y,p) y > p.ym*(x-p.xmin)+p.ycp1;
                    Acp = 1.*(A1 <= x1 & ~ycpfunc(A1,A2,ypr)) + ...
                          2.*(A1 <= x1 &  ycpfunc(A1,A2,ypr)) + ...
                          3.*(A1 >  x1 &  ycpfunc(A1,A2,ypr)) + ...
                          4.*(A1 >  x1 & ~ycpfunc(A1,A2,ypr));
                    clear A1 A2 x1 y1 x2 y2 x3 y3 min1 max1 ym yc ypr ycpfunc;
                end
        end
        A = [ A Acp ];
end

%# Return OUTPUT parameters
B = A;
if nargout == 2
    varargout{1} = g;
end

end

