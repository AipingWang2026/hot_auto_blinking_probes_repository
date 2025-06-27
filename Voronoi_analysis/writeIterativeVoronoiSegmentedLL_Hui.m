% Function writes a new localization list for Insight3 where localizations
% are assigned to a channel depending on the Voronoi area they occupy,
% where segmentation is obtained using the function
%   iterativeVoronoiSegmentation.m
% 
% Inputs
%   LL - Insight3 object 
%   xy - array of [x,y] coordinates in columns [1,2], the zero rank Voronoi 
%       areas in column 3 and, optionally, first rank Voronoi areas in
%       column 4
%   thresh - vector of Voronoi area threshold values for segmenting the 
%       localization list
%   spath - (optional) path where the localization lists should be saved


function writeIterativeVoronoiSegmentedLL_Hui(LL,xy,thresh,spath)

% thresh = Allthresholds;
[fpath,fname,ext] = fileparts( LL.filename );
if ~exist('spath','var') || exist('spath','dir')
    spath = fpath;
end
saveName = fullfile( spath, [fname '_iterVorSeg' ext] );

dat = LL.getData();
ch = LL.getColumnIndex('channel');
% delete all channel information
dat(:,ch) = nan;

% find localizations on the edge with infinite Voroni area
idxSet = isnan(xy(:,3)); % channel 0
dat(idxSet,ch) = 0;
% %%
for t = 1:length(thresh)
    %
%     t = 1;
    % find localizations above the threshold
    idx = xy(:,3) > thresh(t);
    
    % avoid previously chosen points
    idx(idxSet) = false;
    
    % reset the channel assignment 
    dat(idx,ch) = t-1;
    
    % update the list of utilized localizations
    idxSet(idx) = 1;
%     t
end

dat(~idxSet,ch) = t;
%%update the missing two point
%dat(length(dat)-1:length(dat),ch) = t;
%update the nan position number
dat(isnan(dat(:,ch)),ch) = t;


if ~isempty( find(isnan(dat(:,ch)), 1) )
    error('  Uh oh, some localizations didn''t get their channel reset')
end


% save the thresholded .bin list
LL.setFilename( saveName );
LL.forceFileOverwrite(true);
LL.setData( dat );
LL.write;

end % of function