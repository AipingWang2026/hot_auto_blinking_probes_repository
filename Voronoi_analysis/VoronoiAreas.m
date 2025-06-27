% Function for performing 2D Voronoi Tesellation and calculating the zero 
% and first rank area for each [x,y] input point.
%   J.Otterstrom MATLAB 2016a
%
% INPUTS
%   x (required) - either a vector of 'x' positions or a matrix of '[x,y]' 
%       positions
%   y (optional) - a vector of 'y' positions; include this input ONLY if x 
%       is a vector
%   rankAreas (optional) - a logical true (1) or false (0) input to either
%       calculate only the zero rank areas (false input) or also the first
%       rank areas (true input). Default value is false (0).
%   showText (optional) - logical true (1) or false (0) input to determine
%       whether the user is updated regarding the status of the Voronoi
%       area calculations
%
% OUTPUT
%   X - either a 3-column (rankAreas=false) or 4-column(rankAreas=true)
%       matrix with columns corresponding to:
%       column 1: input 'x' vector, or first column of x matrix
%       column 2: input 'y' vector, or second column of x matrix
%       column 3: zero rank area per point
%       column 4: first rank area per point iff rankAreas=true
%   DT - Delaunay Triangulation object output by MATLAB function:
%       delaunayTriangulation()
%   neighborList - a 2-column cell matrix where each row corresponds to one
%       input position in [x,y]. The first column contains a list of
%       neighboring positions for each input as determined by
%       triangulation. The second column is a count for the number of
%       neighbors.
%   VorDat - 1x2 cell array {V, C} containing the Voronoi verticies as the 
%       first element and the Voronoi regions are the second element, both 
%       being returned from the voronoiDiagram method for 
%       delaunayTriangulation objects
%   repidx - mx4 cell array reporting the indicies of m number of repeated
%       localizations in the input xy localization list. 
%       column 1 = indicies of repeated localizations
%       column 2 = [x,y] values of the repeated localizations (for checking)
%       column 3 = difference between localizations reported in column 2, 
%           all elements should be [0,0], otherwise there's a problem
%       column 4 = difference between the indicies of repeated localizations
%           there should not be any values of zero
%
% to call function use commands:
% X = VoronoiAreas(x,y,rankAreas) where x and y are vectors and rankAreas 
%   is logical true (1) or false (0) 
% X = VoronoiAreas(x,rankAreas) where x is a matrix and rankArea is logical
%   1 or 0
%   equivalently: X = VoronoiAreas(x,[],rankAreas)
% X = VoronoiAreas(x,y) where x and y are vectors
% X = VoronoiAreas(x) where x is a matrix
%
% Finished: 160913
% update: 160915, included neighborList output
% update: 161010, use delaunayTriangulation function output for Voronoi
%       calculation
% update: 161017, included showText variable
% update: 161021, include check for unique localizations and return the
%       repeated indicies in repidx variable



function [X,DT,VorDat,repidx,neighborList] = VoronoiAreas(x,y,rankAreas,showText)

if nargin < 2 || isempty(y) || length(x)~=length(y)% (length(y)==1 && (~exist('rankAreas','var') || isempty(rankAreas)) )
    if exist('y','var') && ~isempty(y)
        if length(y)==1
            if exist('rankAreas','var')
                showText = rankAreas;
            end
            rankAreas = y;
            clear y
        else
            error('unrecognized second input')
        end
    end
    % input assumed to be a 2-column matrix of [x,y] values
    if size(x,2)>size(x,1)%~iscolumn(x) % then the data is along columns, so flip it
        rot = 1;
        X = transpose(x);
%         disp('Sec 1 yes transpose')
    else
        rot = 0;
        X = x;
%         disp('Sec 1 no transpose')
    end
else
    if ~iscolumn(x), 
        x = transpose(x); 
%         disp('Sec 2 yes transpose')
    end
    if ~iscolumn(y), 
        y = transpose(y); 
    end
    X = [x,y];
    rot = 0;
%     disp('Sec 2')
end

% if nargin == 3 && ~isempty(origInput3forText)
%     showText = origInput3forText;
% else
if ~exist('showText','var') || isempty(showText)
    showText = true;
end

% if ~exist('rankAreas','var') || isempty(rankAreas)
%     rankAreas = false;
% end

if nargout > 4 && (~exist('rankAreas','var') || isempty(rankAreas) || rankAreas~=1) % then they want the neighbor list
    rankAreas = true;
    disp('You request the neighbor list, so the 1st rank areas will also be calculated')
elseif ~exist('rankAreas','var') || isempty(rankAreas)
    rankAreas = false;
end

nPoints = size(X,1);

tic

% Estimate the amount of time needed for the calculation
if rankAreas
    X = [X,nan(nPoints,2)];
%     estT = polyval([0.97945 -3.7775],log10(nPoints));
    estT = polyval([1.0015 -3.5355],log10(nPoints));
    % initialize output
    neighborList = cell(nPoints,2);
else
    X = [X,nan(nPoints,1)];
%     estT = polyval([1.055 -4.8924],log10(nPoints));
    estT = polyval([0.97777 -3.9502],log10(nPoints));
end
estT = 10^estT;
if showText
    [estT,estt_units] = detTunits(estT);
    disp(['Estimated calculation time is ' num2str(estT,'%.4g') estt_units])
end

% calculate the Voroni verticies "V" and their connectivity matrix "C"
if showText
    strout = 'Calculating Voronoi polygons...';
    fprintf([strout '\n'])
end
% need to ensure the positions are unique
[xu,ixy,ixu] = unique(X(:,1:2),'rows');
% xu = X(ixy,:);  X = xu(ixu,:);
if length(xu) ~= nPoints
%     Xorig = X;
    % find the indicies of repeated localizations
    bincts = histc(ixu,unique(ixu)); %size(find(bincts==2))
    k = find(bincts==2);
    repidx = cell(length(k),4);
    for j = 1:length(k)
        repidx{j,1} = intersect(find(X(:,1)==xu(k(j),1)),find(X(:,2)==xu(k(j),2)));
        repidx{j,2} = X(repidx{j},1:2); % original Locs
        repidx{j,3} = diff(repidx{j,2}); % difference between [x,y] of repeats, should be [0,0]
        repidx{j,4} = diff(repidx{j,1}); % index separation
    end
    X = X(ixy,:);
    nPoints = size(X,1);
else
    repidx = {};
%     Xorig = [];
end
% [V,C]=voronoin(X(:,1:2));
DT = delaunayTriangulation(X(:,1:2));
[V,C] = voronoiDiagram(DT);
VorDat = {V, C};
if showText
    fprintf( [repmat('\b',1,length(strout)+1), strout, 'DONE!\n'] )
end

% Calculate the area for each Voronoi polygon
strCR = '';
for pt=1:nPoints
%     X(pt,3) = polyarea(V(C{pt},1),V(C{pt},2));
    X(pt,3) = polyarea(V(C{pt},1),V(C{pt},2));
    if showText
        strout = sprintf('Calculating zero rank area for each polygon: %d / %d \n',pt,nPoints);
        fprintf([strCR strout]);
        strCR = repmat('\b',1,length(strout));
    end
end

if rankAreas
    % use DT to find neighboring polynomials and calculate summed 1st rank
    % area
    neighborList = findVoronoiNeighbors(DT, [], showText);
    % using the sortrows function to identify neighbors is much faster than
    % using the 'find' function
%     [ConSort,idx] = sortrows(DT.ConnectivityList(:));
%     nCon = length(ConSort);
%     szConMtx = size(DT.ConnectivityList);
%     prevIdx = 1;
    strCR = '';
    for pt=1:nPoints 
%         % find 1st rank neighbors using delaunay triangulation
%         nextIdx = prevIdx+1;
%         while nextIdx<=nCon && pt == ConSort(nextIdx)
%             nextIdx = nextIdx+1;
%         end
%         nextIdx = nextIdx-1; % remove increment since test was actually for ConSort(nextIdx-1)
%         [r,~]=ind2sub(szConMtx,idx(prevIdx:1:nextIdx));
%         % determine the indicies for the point in question & neighbors
%         pt_neighbor = unique(DT.ConnectivityList(r,:));
%         % remove the point in question and save a listing of neighbors for each point
%         neighborList{pt,1} = pt_neighbor(~ismember(pt_neighbor,pt));
%         % save the number of neighbors
%         neighborList{pt,2} = nextIdx-prevIdx+1;
%         prevIdx = nextIdx+1; % start next search from subsequent index
        
        % if the polygon is closed, calculate first-rank area
        if ~isnan(X(pt,3))
%             neighborAreas = X(pt_neighbor,3);
            neighbPts = [pt;neighborList{pt,1}];
            neighborAreas = X(neighbPts,3);
            neighborAreas = neighborAreas(~isnan(neighborAreas));
            X(pt,4) = sum(neighborAreas);
        end
        
        if showText
            strout = sprintf('Calculating 1st rank area for each polygon: %d / %d \n',pt,nPoints);
            fprintf([strCR strout]);
            strCR = repmat('\b',1,length(strout));
        end
    end
end
% % now return the list as before removal of duplicate localizaionts
%>> can't return the original list because it fucks up cluster extraction
%   since the neighborlist is for the 'unique' set of localizations
% if ~isempty(Xorig)
%     if rankAreas
%         Xorig(:,3:4) = X(ixu,3:4);
%     else
%         Xorig(:,3) = X(ixu,3);
%     end
%     X = Xorig;
% end

% find how long the calculation took and report to user
t = toc;
if showText
    [t,t_units] = detTunits(t);
    disp(['For ' num2str(nPoints) ' localizations, total Voronoi area calc time = ' num2str(t,'%.3g') t_units])
    disp(' ')
end

if rot, X = transpose(X); end


% subfunction
    function [tin,tmunits] = detTunits(tin)
        if tin<60
            tmunits = ' seconds';
        elseif tin<60*60
            tin=tin/60;
            tmunits = ' minutes';
        elseif tin<60*60*24
            tin=tin/(60*60);
            tmunits = ' hours';
        else
            tin=tin/(60*60*24);
            tmunits = ' days';
        end
    end
end