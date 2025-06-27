% Function that iteratively identifies thresholds for Voronoi Areas given
% an input of Voronoi areas from data and from a Monte Carlo simulation
%
% input is intended to be the output from the associated function:
%     VoronoiMonteCarlo_JO.m
%
% Input
%   xy = array containing [x,y] coordinates in columns 1 & 2, the zero-rank
%       Voronoi area in column 3 and (optionally) the first-rank Voronoi
%       area in column 4
%   Varea_rnd = cell vector with each element containing the Voronoi areas
%       obtained from randomized localizations and that is an output from
%       VoronoiMonteCarlo_JO.m
%   maxloop = maximum number of thresholds to identify.  Default is 8 since
%       Insight3 - presumed destination program for image rendering - only
%       has 9 channels
%   signif = confidence bounds for the Monte Carlo simulation.  If no input
%       is provided, a default value of 95% will be assigned
%   startThresh = Optional argument that can take different forms:
%       logical 
%           true - (default) indicates that a threshold has been applied to 
%           BOTH the xy(:,3) and Varea_rnd inputs
%           false - indicates that a threshold has been applied to  xy(:,3) 
%           before input, but not to Varea_rnd
%       double 
%           The value of the threshold will be applied to both xy(:,3) and
%           Varea_rnd before data is binned for identification of the first
%           threshold value
%     Note: The program attempts to check if the user has applied a
%     threshold only to xy(:,3) data and not to Varea_rnd, but did not
%     indicate this fact. This would be equivalent to an input value of 
%     'false'.
%     If a user inputs value of 'true', but in reality applied it only 
%     xy(:,3), then erroneous threshold values will be returned.
%   showplot = logical true/false value for plotting the distributions of 
%       Voronoi areas for the data & randomized localizations together for
%       each iteratively identified threshold
%
% Output
%   thresholds = array with first column containing thresholds identified 
%       by iteratively applying a threshold to the voronoi area of the 
%       input data and finding where its area distribution crosses the that 
%       of the random data.
%       The second column are the p-values associated with the ttest
%       performed to identify when the area distribution for xy(:,3) data
%       is statistically below the upper confidence bound of the randomized
%       data, thereby invalidating any thresholds identified

function thresholds = iterativeVoronoiSegmentation(xy,Varea_rnd,maxloop,signif,startThresh,showplot)

% % testing
% startThresh = unithresh; showplot=true

if ~exist('maxloop','var') || isempty(maxloop)
    maxloop = 8;
end
if (~exist('signif','var') || isempty(signif))
    signif = 95;
    disp(['Applying default value of ' num2str(signif) '% confidence bounds'])
end

% z = norminv(1-(100-signif)*0.01/2);
thresholds = nan(maxloop,2);
iter = size(Varea_rnd,1);
loop = 0;

Vareas = xy(~isnan(xy(:,3)),3);
areaLim = interpercentilerange(Vareas,[ .99 .999]);
areaLim = areaLim(2);

if ~exist('startThresh','var') || isempty(startThresh)
    startThresh = true;
    startThreshInput = false;
else
    if ~islogical(startThresh) && sum(startThresh==[0 1])
        % input was 0 or 1, so convert to logical for later testing
        startThresh = logical(startThresh);
    else
        %make space for the starting threshold value
        thresholds = thresholds(1:end-1,:);
    end
    startThreshInput = true;
end
if ~exist('showplot','var') || isempty(showplot)
    showplot = false;
elseif ~islogical(showplot)
    showplot = logical(showplot);
end


stoploop = false;
%%
while loop <= maxloop && ~stoploop
    %%
    loop = loop+1;
    if loop > 1
        % apply first threshold to Voronoi Areas of data
        xy = xy(xy(:,3) < thresholds(loop-1,1),:);
    elseif ~islogical(startThresh)
        % apply starting threshold unless it is a logical value
        xy = xy(xy(:,3) < startThresh,:);
%         disp('applied startThresh to data')
%     else
%         disp('No thresholds applied to first loop')
    end
    % Select out Varea of data
    Varea_dat = xy(~isnan(xy(:,3)),3);
    % set the histogram binning limits
%     Nbins = 2*round(2*(size(xy,1))^(1/3));
    % select Voronoi areas within the limit for binning
    Varea_dat = Varea_dat(Varea_dat<areaLim);
    nptsAnalyzed = size(Varea_dat,1);
    % obtain data's Voronoi Area distribution
%     [counts, centers] = hist(Varea_dat,Nbins);
    [counts, centers] = histcounts(Varea_dat, 'BinMethod','fd');
    % try to determine if a previous threshold has been applied
    if loop==1 && ~startThreshInput % user didn't specify if a previous threshold applied
        if 100*counts(end)/size(xy,1) > 0.1% of input data, likley no previous threshold
            startThresh = false;
            disp('Previous application of Voronoi area threshold detected')
        end
    end
    
    % obtain randomized points' Voronoi Area distribution
    counts_rnd = nan(iter,length(counts));
%     counts_rnd = nan(iter,Nbins);
%     centers_rnd = nan(iter,Nbins);
    for j = 1:size(Varea_rnd,1)
        % apply thresholds as appropriate
        Varea_rnd_used = Varea_rnd{j};%(1:nRndPts);
%         Varea_rnd_used = Varea_rnd_used(Varea_rnd_used<areaLim);
        if loop > 1
            % apply the threshold to randomized data also
            Varea_rnd_used = Varea_rnd_used( Varea_rnd_used<thresholds(loop-1,1) );
        elseif ~islogical(startThresh)
            Varea_rnd_used = Varea_rnd_used( Varea_rnd_used<startThresh );
%             if j == 1
%             disp('applied startThresh to randomized data')
%             end
%         else
%             if j ==1
%             disp('No thresholds applied to randomized data in first loop')
%             end
        end
        % select randomized voronoi areas from the same number of points as 
        % there are number of data points below the threshold(i)
        nRndPts = min([nptsAnalyzed,size(Varea_rnd_used,1)]);
        Varea_rnd_used = Varea_rnd_used(1:nRndPts);
        [counts_rnd(j,:), ~] = histcounts(Varea_rnd_used,centers);
        if islogical(startThresh) && ~startThresh && loop == 1
            % a threshold has been applied to data, but not random Locs,
            % Set final bin to zero, since 'centers' will span a majority
            % of the distribution
            counts_rnd(j,end) = 0;
            disp('Attempted to correct for previously applied threshold')
%             % so control the bin centers
%             [counts_rnd(j,:), centers_rnd(j,:)] = ...
%                 hist(Varea_rnd_used,centers);
%             if j == 1
%             disp('Randomized data is binned equally to actual data for first threshold identification')
%             end
%         else
%             % let bin centers vary
%             [counts_rnd(j,:), centers_rnd(j,:)] = ...
%                 hist(Varea_rnd_used,Nbins);
        end
        % Adjust the counts to have: sum(counts_rnd(j,:))==sum(counts)
        counts_rnd(j,:) = counts_rnd(j,:)*sum(counts)/sum(counts_rnd(j,:));
    end
    
    stp = centers(2)-centers(1);
    centers = stp*0.5+centers(1:end-1);
    
    MeanCounts = mean(counts_rnd,1);
    StdCounts = std(counts_rnd,1,1);
    z = norminv(1-(100-signif)*0.01/2);
    upperConf = MeanCounts + z * StdCounts;
    lowerConf = MeanCounts - z * StdCounts;
%     
%     centers_rnd = mean(centers_rnd,1);
%     MeanCounts = mean(counts_rnd,1);
%     StdCounts = std(counts_rnd,1,1);
%     upperConf = MeanCounts + z * StdCounts;

    % make a spline fit for the data & random data, avoiding the
    % untrustworthy final bin
    splineData = fit(centers(1:end-1)', counts(1:end-1)','smoothingspline');
    splineRnd = fit(centers(1:end-1)', MeanCounts(1:end-1)','smoothingspline');
    xtst = linspace(centers(1),centers(end),length(centers)*5);
    ySpDat = feval( splineData, xtst );
    ySpRnd = feval( splineRnd, xtst );
    [xi, yi] = polyxpoly(xtst, ySpDat, xtst, ySpRnd);
    
%     % the last bin is untrustworthy for determining intersections
%     [xi, yi] = polyxpoly(centers(1:end-1), counts(1:end-1),...
%         centers(1:end-1), MeanCounts(1:end-1));
    intersection = [xi, yi];
%     [xUC, ~] = polyxpoly(centers(1:end-1), counts(1:end-1),...
%         centers(1:end-1), upperConf(1:end-1));
    % I expect at most two intersections; if there are more, then I presume
    % it's from a noisy signal and so the data must first be smoothed
    % to find intersections that may or may not exist
    if size(intersection,1) > 2 && length(counts)>15 % limit in slideFilter.m
        counts_sm = slideFilter('avg',counts,4,2);
        MeanCounts_sm = slideFilter('avg',MeanCounts,4,2);
%         upperConf_sm = slideFilter('avg',upperConf,4,2);
        splineData = fit(centers(1:end-1)', counts_sm(1:end-1)','smoothingspline');
        splineRnd = fit(centers(1:end-1)', MeanCounts_sm(1:end-1)','smoothingspline');
        xtst = linspace(centers(1),centers(end),length(centers)*5);
        ySpDat = feval( splineData, xtst );
        ySpRnd = feval( splineRnd, xtst );
        [xi, yi] = polyxpoly(xtst, ySpDat, xtst, ySpRnd);
%         [xi, yi] = polyxpoly(centers(1:end-4), counts_sm(1:end-4), ...
%             centers(1:end-4), MeanCounts_sm(1:end-4));
        intersection = [xi, yi];
%         [xUC, ~] = polyxpoly(centers(1:end-1), counts_sm(1:end-1),...
%             centers(1:end-1), upperConf_sm(1:end-1));
    end
    
    if ~isempty(intersection) %&& size(intersection,1) == 2
        % See if data has fallen within confidence bounds, so unreliable intersection
        idx = find(centers < intersection(1,1),1,'last');
%         idx = find(centers < min([intersection(1) max(xUC)]),1,'last');
%         idx = find(centers < min(xUC),1,'last');
%         idx = find(centers < min(intersection(:,1)),1,'last');
%         y1 =  counts(1:idx)-upperConf(1:idx) ;
%         y1 = abs( counts(1:idx)-MeanCounts(1:idx) );
%         y2 = abs( counts(1:idx)-(MeanCounts(1:idx) + z * StdCounts(1:idx)) );
%         [hpth,pval] = ttest(y1-y2);
        [hpth,pval] = ttest(counts(1:idx),upperConf(1:idx),...
            'Tail','right','Alpha',(100-signif)/100);    
%         figure,plot(centers(1:idx),counts(1:idx)-upperConf(1:idx)),title(['p-val = ' num2str(pval)]),xlabel('Voronoi Area'),ylabel('Data - UpperConf')
        if ~isnan(hpth)&&hpth % mean(y1) > mean(y2) --------Di
            if loop > 1
                % ensure leftward-moving thresholds
                if min(intersection(:,1)) < thresholds(loop-1) 
                    thresholds(loop,1) = min(intersection(:,1));
                    thresholds(loop,2) = pval;
                else stoploop = true;
                end
            else
                thresholds(loop,1) = min(intersection(:,1));
                thresholds(loop,2) = pval;
            end
        else
            stoploop = true;
            if loop == 1
                % inform user that no intersections were found with
                % p-value below their specified criteria
                disp('NO THRESHOLDS IDENTIFIED FROM AREA DISTRIBUTION INTERSECTION ANALYSIS')
                disp(['Intersection at Voronoi Area = ' num2str(xi(1),'%.3g') ' has p-value = ' num2str(pval,'%.3g') ])
                disp(['     This value is larger than the user-defined cutoff of p < ' num2str((100-signif)/100)])
                disp(['To have this intersection selected, re-analyze data with confidence bounds of ' num2str(floor(100*(1-pval))) '%'])
            end
        end
    else stoploop = true;
    end
    
    if showplot
        plotVoronoiMCdat([centers', counts', MeanCounts', ...
            (lowerConf)', (upperConf)'],...
            thresholds(loop,1), signif);
        
%         idx = find(centers < min(intersection(:,1)),1,'last');
%         y1 = counts(1:idx)-MeanCounts(1:idx);
%         y2 = counts(1:idx)-(MeanCounts(1:idx) + z * StdCounts(1:idx));
%         y3 = y1-y2;
%         [hpth,pval] = kstest2(y1,y2);
%         [hpth,pval] = ttest2(y1,y2);
%         [hpth,pval] = ttest(y1,0,'Tail','right');
%         figure
%         hold on
%         plot(centers(1:idx),y1,'-','Color',[0 0.7 0]) %green, dist to mean
%         plot(centers(1:idx),y2,'-','Color',[0.8 0.3 0]) % orange, dist to upp conf
%         plot(centers(1:idx),y3,'k-')% black, average separation
%         title({['to Mean = ' num2str(mean(y1),'%.3g') ' ± ' num2str(std(y1,1),'%.3g')] ...            ', to upper Conf = ' num2str(mean(y2),'%.3g') ' ± ' num2str(std(y2,1),'%.3g')];...
%             ['Ttest = [' num2str(hpth) ' ' num2str(pval,'%.3g') ']' ]})
    end
end
% remove nan values
thresholds = thresholds( ~isnan(thresholds(:,1)),: );
% include the starting threshold, if input
if isa(startThresh,'double')
    thresholds = [startThresh 0; thresholds];
end

end % of function
