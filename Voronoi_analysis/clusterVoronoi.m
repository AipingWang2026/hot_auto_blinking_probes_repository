
% Inputs
%   convFact = 1- or 2-element array with the first value being pixels to 
%       nm conversion, for converting [x,y] positions from pixels to
%       nanometers; the second element (if desired) being the pixel size of 
%       the final mask.
%
% Outputs
%   cluster - structure having the following fields
%       .nLocs = number of localizations in a cluster.
%       .area = summed area of all voronoi polygons comprising a cluster.
%       .center = [x,y] position of the cluster's center of mass weighted 
%           by the density of each voronoi polygon in a cluster.
%           Center = <rho*X> of particles where rho = 1/area.
%           Units are the same as input xy array.
%       .Locs = cell matrix containing the original localizations contained
%           within the cluster in columns [x,y] in columns 1,2, zero-rank
%           area in column 3, first-rank area in column 4, index values
%           column 5.
%       .ellipticity = ellipticity ratio for the cluster
%       .HullArea = Area contained in the convex hull drawn around the
%           localizations comprising a cluster; this area is larger than
%           the summed area in the .areas field
%       .NND = Nearest neighbor distance (column 1) in units of input xy
%           array, and the cluster index of the closest neighbor (column 2), 
%           mean distance to all neighbors (column 3), median distance to
%           all neighbors (column 4); neighbors found via triangulation.
%   fig_h - handle to the figure generated if shoIM == true
%   xy - input localization coordinates with zero & first order Voronoi
%       areas in columns 3 & 4, respectively
%   thresh - if input ~isempty(thresh), then the output thresh corresponds
%       to the manually selected Voronoi area threshold

function [cluster, fig_h, xy, thresh] = clusterVoronoi(xy,thresh,neighborList,minNLoc,convFact,showIM,mask,figTitle)

if size(xy,2) < 2 && size(xy,1) < 2
    error('unknown [x,y] localization list format')
elseif size(xy,2) == 2
    if ~exist('neighborList','var') || isempty(neighborList)
        [xy,~,~,~,neighborList] = VoronoiAreas(xy,true);
    else
        xy = VoronoiAreas(xy,false);
    end
end
% extract the Voronoi areas
Vareas = xy(:,3);
% extract number of localizations
nPts = length(xy);
% check for threshold input
if ~exist('thresh','var') || isempty(thresh)
    disp('No area threshold input, select now on plot...')
    VareasPlot = Vareas(~isnan(Vareas));
    % set the histogram binning limits
%     Nbins = round(2*nPts^(1/3));
%     areaLim = 4*median(VareasPlot);
%     [counts, centers] = hist(VareasPlot(VareasPlot<areaLim),Nbins);

areaLim = interpercentilerange(VareasPlot,[ .99 .995]);
areaLim = areaLim(1);
[counts, centers] = histcounts(VareasPlot(VareasPlot<areaLim), 'BinMethod','fd');
stp = centers(2)-centers(1);
centers = stp*0.5+centers(1:end-1);
    
    base = [min(centers)*0.05 min(counts)*0.05];
    ylimTop = ceil(1.05*max(counts));
    accept = false;
    while ~accept
        f1 = figure;
        plot(centers,counts,'b-')
        xlabel('Voronoi Polygon Area')
        ylabel('Count')
        title('Click once to set the threshold area')
        axis([0 centers(end), 0 ylimTop])
        [x_ui,~] = ginput(1);
        thresh = x_ui(1);    
        hold on
        box on
        area([base(1),thresh],[1,1]*ylimTop-base(2),base(2),...intersection(2),...
            'LineStyle','none',...
            'FaceColor',[175 238 238]/255)
        plot(centers,counts,'b-')
        commandwindow
        accept = input(['Area threshold selected to be: ' num2str(thresh,'%.4g') ', do you accept? (y/n) ==> '],'s');
        if ~strcmpi(accept,'y')
            accept = false;
        end
        close(f1)
    end
end
    
% check for input conversion factors to go from pixels to nm
if ~exist('convFact','var') || isempty(convFact)
    pix2nm = 1;
    p = 1;
elseif length(convFact)==2
    pix2nm = convFact(1);
    p = convFact(2);
else
    pix2nm = convFact(1);
    p = 1;
end
% [pix2nm,p]
% check for input neighbor list
if ~exist('neighborList','var') || isempty(neighborList)
    disp('Neighbor list for Voronoi Polygons not included, calculating now...')
    [~,~,~,~,neighborList] = VoronoiAreas(xy,true);
end

% check for the minimum # localizations/cluster
if ~exist('minNLoc','var') || isempty(minNLoc)
    minNLoc = 5;
    disp(['MINIMUM NUMBER OF LOCALIZATIONS PER CLUSTER SET TO DEFAULT = ' num2str(minNLoc)])
end

%% check to see if user wants to plot the mask with localizations
if ~exist('showIM','var') || isempty(showIM)
    showIM = false;
end
if ~exist('figTitle','var') || isempty(figTitle)
    figTitle = '';
else
    if strcmp(figTitle(end-3),'.') % then input is full filename
        figTitle = figTitle(1:end-4); % remove extension
    end
    % remove folder separators if present
    idx = strfind(figTitle,filesep);
    if ~isempty(idx)
        figTitle = figTitle(idx(end)+1:end);
    end
    % now replace underscores with spaces
    figTitle(figTitle=='_') = ' ';
    figTitle = [figTitle '; Threshold = ' num2str(thresh,'%.3g')];
end
%% there had better be a mask if they want to plot it
% [isempty(mask) showIM ((~exist('mask','var') || isempty(mask)) && showIM)]
if (~exist('mask','var') || isempty(mask)) && showIM
    disp(['No mask input, calculating now using conversions: pix2nm = ' ...
        num2str(pix2nm) ', final pixel size = ' num2str(p) ' nm'])
    [mask,~] = Locs2Mask( xy(:,1:2));
    mask = imresize(mask,pix2nm/p);
end

%% Begin Algorithm

% define the localizations above the threshold that are to be kept
kppt = Vareas <= thresh; % X(:,3)<=thresh; % nanometers

% these are the [x,y] values that will be plotted
xypt = xy(:,1:2).*pix2nm./p; % X(:,1:2)./p;
% check to make sure the mask is big enough for plotting
if showIM && ((max(xypt(:,1)) > size(mask,2)) || (max(xypt(:,2)) > size(mask,1)))
    % then the mask is too small, so re-size it
    mask = imresize(mask,pix2nm/p);
    disp('mask resized for input conversion factors')
end

% %% testing
% minNLoc = 5;
% showIM = true;
fig_h = [];
if showIM
    fig_h = figure('units','pixels','Name','Voronoi Clusters');
    disp('Figure will resize when calculation completes, please be patient...')
    imshow(mask)
    fig_ax = gca;
    hold on
    % plot voronoi diagram
    hv = voronoi(xypt(:,1),xypt(:,2),'g');
    % change color of voronoi diag
    for i=1:2
        set(hv(i),'Color',[0.9290    0.6940    0.1250])%[0 0.6 0])
    end
    % blue - plot all localizations having high density
    plot(xypt(kppt,1),xypt(kppt,2),'bs','MarkerSize',3)
    [y,x] = find(mask);
    lm = [min(x) max(x) min(y) max(y)];
    if lm(1) > 3, lm(1) = lm(1)-2; end
    if lm(2) < size(mask,2)-2, lm(2) = lm(2)+2; end
    if lm(3) > 3, lm(3) = lm(3)-2; end
    if lm(4) < size(mask,1)-2, lm(4) = lm(4)+2; end
    set(fig_h,'position',[501 125 918 841])
    axis(fig_ax,lm)
    set(fig_ax,'units','normalized')
    set(fig_ax,'position',[0.0928 0.0777 0.8126 0.8870])
end

% %%

c = 0;
usedPt.nClusters = zeros(nPts,1);
usedPt.cluster = cell(nPts,1);
% cluster = struct('areas',[],'center',[],'ellipticity',[],'Locs',{},'nLocs',[],'NND',[]);
strCR = '';
for pt = 1:nPts;
    strout = sprintf('Grouping localizations into clusters: %d / %d \n',pt,nPts);
    fprintf([strCR strout]);
    strCR = repmat('\b',1,length(strout));
    if kppt(pt) && ~usedPt.nClusters(pt)
        % get the neighbors surrounding the seed-point
        idxNebHi = neighborList{pt,1}(kppt(neighborList{pt,1}));
        if ~isempty(idxNebHi)
            nNbHi = 0;
            sz = size(idxNebHi,1);
            %% find all connected neighbors above thresh
            while sz ~= nNbHi
                %%
                nNbHi = length(idxNebHi);
                tmp = cell(nNbHi,1);
                for n = 1:nNbHi
                    tmp{n} = neighborList{idxNebHi(n),1};
                end
                idxNebHi2 = unique(cell2mat(tmp));
                if ~any(ismember(idxNebHi2,idxNebHi))
                    idxNebHi = sort([idxNebHi2;idxNebHi]);
                else
                    idxNebHi = idxNebHi2;
                end
                idxNebHi = idxNebHi(kppt(idxNebHi));
                sz = size(idxNebHi,1);
            end
        else
            idxNebHi = pt;
        end
        
        %%        
        nLocCluster = length(idxNebHi);
        if nLocCluster >= minNLoc
            % extract cluster data
            c = c+1;
            % mark all of the indicies idxNebHi as used in a cluster
            for i = 1:nLocCluster
                %             usedPt.nClusters(idxNebHi(i)) = usedPt.nClusters(idxNebHi(i))+1;
                usedPt.cluster{idxNebHi(i)} = [c;usedPt.cluster{idxNebHi(i)}];
                usedPt.nClusters(idxNebHi(i)) = length( usedPt.cluster(idxNebHi(i)) );
            end
            % 1) number of localizations in the cluster
            cluster.nLocs(c,1) = nLocCluster;
            % 2) total area
            nebAreas = Vareas(idxNebHi);
            cluster.areas(c,1) = sum(nebAreas);
            % 3) center of cluster
            % % method 1: <X> of particle
            % Xbar = mean(xypt(idxNebHi,:));
            % plot(Xbar(1),Xbar(2),'k+','MarkerSize',12)
            %     % reject for now in favor of using the voronoi area for centroid
            %     % localization as well
            % % method 2: <X> of verticies
            % % find the verticies of the cluster
            % idxV = cell2mat( ...
            %     cellfun(@transpose,...
            %     VorDat{2}(idxNebHi),...
            %     'UniformOutput',false)...
            %     );
            % Vbar = mean(clustV,1);
            % plot(Vbar(1),Vbar(2),'kx','MarkerSize',12)
            %     % No, gives weight to the outer-most verticies
            % method 3: <rho*X> of particles where rho = 1/area
            rho = 1./nebAreas;
            rXbar = [sum(rho.*xypt(idxNebHi,1))/sum(rho),sum(rho.*xypt(idxNebHi,2))/sum(rho)];
            % cluster centroid found by weighting each localization by it's voronoi
            % density
            cluster.center(c,:) = rXbar*p/pix2nm;
            % 4) listing of which localizations are in the cluster
            cluster.Locs{c,1} = [xy(idxNebHi,:), idxNebHi];
%             if length(idxNebHi) > 1
            % 5) send localizations for calculation of cluster ellipticity
                cluster.ellipticity(c,1) = ...
                    clusterEllipticity(cluster.Locs{c,1}(:,1:2),cluster.center(c,1:2));
%                 % ellipticity of the cluster (make subfunction)
%                 xsz = max(xypt(idxNebHi,1))-min(xypt(idxNebHi,1));
%                 ysz = max(xypt(idxNebHi,2))-min(xypt(idxNebHi,2));
%                 if xsz < ysz, b = xsz; a=ysz;
%                 else b = ysz; a=xsz; end
%                 cluster.ellipticity(c,1) = sqrt( 1-b^2/a^2 );
%             else
%                 cluster.ellipticity(c,1) = 1;
%             end
            
            % find edges of the identified cluster
%                 % 1) using the verticies of the included Voronoi Polygons
%                 clustV = VorDat{1}(idxV,:)./p;
%                 clustV = unique(clustV,'rows');
%                 % indicate the verticies
%                 % plot(clustV(:,1),clustV(:,2),'m*','MarkerSize',4)
%                 % show the area delimiting each cluster
%                 dtclust = delaunayTriangulation( clustV ); % using Voronoi verticies for convHull
%               > POOR, GETS BIGGER
                % 2) Using the localizations comprising the cluster
                clustLoc = xypt(idxNebHi,:); % using localizations for convHull
                dtclust = delaunayTriangulation( clustLoc );
                % Convex Hull
                [ConHull, clustVol] = convexHull(dtclust);
                 cluster.HullArea(c,1) = clustVol;
                 % NOT IDEAL, BUT EASIER THAN TRYING TO TRACE THE POINTS
                 % FOR NOW
            
            
            % optionally plot Voronoi Diagram 
            if showIM
                %     % show the seed-point for cluster identification
                %     plot(xypt(pt,1),xypt(pt,2),'s','MarkerSize',7,'Color',[1 0 1].*0.85)
                plot(rXbar(1),rXbar(2),'o','MarkerSize',8,'MarkerFaceColor','k')
                %             plot(clustV(ConHull,1),clustV(ConHull,2),'r-') % using Voronoi verticies for convHull
                plot(clustLoc(ConHull,1),clustLoc(ConHull,2),'r-') % using localizations for convHull
            end
        end
    end
end

% future expansion: perform some type of water shedding on the polynomials
%   to separate clusters that are in contact with each other, but form a
%   type of dumb-bell shape that would be an unusual cluster shape.

%% find NND for each cluster
nClust = length(cluster.nLocs);
DTctr = delaunayTriangulation(cluster.center);
NL = findVoronoiNeighbors(DTctr,false);
for c = 1:nClust
%     % set the NND as the nearest centroid-to-centroid distance
%     r = sqrt( ...
%         (cluster.center(c,1)-cluster.center(:,1)).^2 + ...
%         (cluster.center(c,2)-cluster.center(:,2)).^2 ...
%         );
%     [minNND,idx_r] = min( r(r>0) );
%     if idx_r > c
%         idx_r = idx_r+1;
%     end
%     cluster.NND(c) = [minNND, idx_r];
% end

% alternative: performa  delaunay Triangulation, then find the neighbors
% and take an average NND over them
    
    r_neighbors = sqrt( ...
        (cluster.center(c,1)-cluster.center(NL{c,1},1)).^2 + ...
        (cluster.center(c,2)-cluster.center(NL{c,1},2)).^2 ...
        );
    [minNND,idx_r] = min( r_neighbors );
    % note, this minimum result gives a double-counting since each pair of 
    % closest clusters will have the same NND
    
    % [minimum NND, index of closest cluster, mean & median distance  to neighbors]
    cluster.NND(c,:) = [minNND, NL{c,1}(idx_r), mean(r_neighbors), median(r_neighbors)];
    
end

% update the title if necessary
if showIM
    title(fig_ax,figTitle)
end


% future expansion idea: use the external-most verticies of a cluster to
% determine the nearest cluster, then calculate the NND between centers


end % of function