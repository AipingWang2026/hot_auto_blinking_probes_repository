        clear
% clc

[pixVals,pix2nm,finalpix,clusterID,iter,signif,minNLoc,maxloop,showIM,files] = ... 
        setVoronoiClusterParams();
nfiles = size(files.data,1);

% If doing manual thresholding, ask user if they have an a priori value
% they want to use
switch clusterID
    case 'manual'
        options.WindowStyle = 'normal';
        manualThresh = inputdlg([...
            'Set the threshold [pix^2] you desire for cluster identification.  '...
            'This threshold will be used for all files analyzed.               '...
            'If no value is input, then you will be prompted to choose the threshold from '...
            'the distribution of Voronoi areas for each file.'],...
            'Manual Voronoi Threshold',[1 68],{''},options);
        clear options
        if ~isempty( manualThresh )
            if strcmp(manualThresh{1},''), manualThresh = [];
            else manualThresh = str2double( manualThresh{1} );
            end
        end
        
end

%%%%%% begin algorithm %%%%%% 
for f = 1:nfiles

% load localization list

LL = Insight3( fullfile( files.data{f,2},files.data{f,1} ) );

% get coordinates in units of pixels
xy = LL.getXYcorr;

switch clusterID
    case 'automatic'
        
        timer1 = tic;
        [Histograms, ~, xy, neighborList, mask, DT, VorDat, repidx, Varea_rnd] = ...
            VoronoiMonteCarlo_JO(xy,iter,signif,pixVals);
        t1 = toc(timer1);
        
        % save the output
        savefile = [LL.filename(1:end-4) '_iterVorSegData.mat'];
        save(savefile,...
            'pix2nm','finalpix','pixVals','clusterID','iter','signif','minNLoc','maxloop','showIM',... inputs
            'Histograms', 'xy', 'neighborList', 'mask', 'DT', 'VorDat', 'repidx', 'Varea_rnd',... outputs
            'savefile')

        %% Determine the threshold for area/localization assuming a uniform
        % distribution, a.la. SR-Tesseler
        unithresh = 2*((finalpix/pix2nm)^2)*mask.Area/size(xy,1);
        
        % send localizations with small Voronoi areas for seq seg
        % Apply uniform threshold to establish the base object for sequential
        % cluster segmentation
        Allthresholds = iterativeVoronoiSegmentation_Wenhui(xy,Varea_rnd,maxloop,signif,unithresh,showIM);
        
        % save threshold info
        save(savefile,'Allthresholds','unithresh','-append')

        % plot Voronoi distributions with obtained thresholds
        plotVoronoiMCdat(Histograms, Allthresholds(:,1), signif);
        set(gca,'XLim',[0 unithresh*1.6])
        % save cluster area vs Loc/cluster image as .png file
        saveas(gcf,[savefile(1:end-8) '_thresholds.png'],'png')

        timer2 = tic;
        % write a new molecule list assigning localizations to channels 0-9
        % according to their Voronoi area        
        writeIterativeVoronoiSegmentedLL_Hui(LL,xy,Allthresholds)
        
        % send localizations for clustering
        cluster = cell(size(Allthresholds,1),1);
        for th = 1:size(Allthresholds,1)
            if showIM
                [cluster{th}, fig_h] = ...
                    clusterVoronoi(xy,Allthresholds(th,1),neighborList,minNLoc,...
                    [pix2nm,finalpix],showIM,mask.BW,LL.filename);
                % save cluster stats for application of threshold 'th'
                saveas(gcf,[savefile(1:end-8) '_VoronoiThresh#' num2str(th) '.png'],'png')
            else
                cluster{th} = ...
                    clusterVoronoi(xy,Allthresholds(th,1),neighborList,minNLoc);
            end
            % plot the cluster statistics
            figAnnotation = [files.data{f,1}(1:end-4) ', thresh = ' num2str(Allthresholds(th,1),'%.3g') '    '];
            VClustFig = plotVoronoiClusterStats2(...
                cluster{th}.nLocs, ...
                cluster{th}.areas*(pix2nm^2), ...
                cluster{th}.NND(:,1)*pix2nm, ...
                figAnnotation);
            % save cluster stats for application of threshold 'th'
            saveas(gcf,[savefile(1:end-8) '_ClustThresh#' num2str(th) '.png'],'png')
            
            % save the clusters in a localization list for visualization
            writeVoronoiClusteredLL(LL,cluster{th},Allthresholds(th,1))
        end
        % save cluster info
        save(savefile,'cluster','-append')

        if size(Allthresholds,1) > 1
            plotVoronoiArea_nLocs(cluster,Allthresholds,pix2nm,LL.filename( find(LL.filename==filesep,1,'last')+1:end-4))
            % save cluster area vs Loc/cluster image as .png file
            saveas(gcf,[savefile(1:end-8) '_AreaVLoc.png'],'png')
        end
        t2 = toc(timer2);
        
        % inform user about total calculation time
        disp(' ')
        disp(['Analysis complete: ' num2str(LL.numMolecules) ' localizations '])
        disp(['Total analysis time: ' num2str( minutes( seconds(t1+t2)),'%.3g') ' minutes' ])
        disp(['Monte Carlo simulation required ' num2str( minutes( seconds(t1)),'%.3g') ' minutes' ])
        disp(['Localization clustering required ' num2str( minutes( seconds(t2)),'%.3g') ' minutes' ])
        disp(' ')
       
        
    case 'manual'
           
        [cluster, fig_h, xy, Allthresholds] = clusterVoronoi(xy,manualThresh,[],minNLoc,[pix2nm,finalpix],showIM);
        % plot the cluster statistics
        figAnnotation = [files.data{f,1}(1:end-4) ', thresh = ' num2str(Allthresholds,'%.3g') '    '];
        VClustFig = plotVoronoiClusterStats2(...
            cluster.nLocs, ...
            cluster.areas*(pix2nm^2), ...
            cluster.NND(:,1)*pix2nm, ...
            figAnnotation);
        % save cluster info
        savefile = [LL.filename(1:end-4) '_manualVorSegData.mat'];
        save(savefile,'xy','cluster','Allthresholds','pix2nm','finalpix','minNLoc')
        disp(['Analysis complete: ' num2str(LL.numMolecules) ' localizations '])
        
        % write a new molecule list assigning localizations to channels 0-9
        % according to their Voronoi area        
        writeIterativeVoronoiSegmentedLL_Hui(LL,xy,Allthresholds)
        
        % save the clusters in a localization list for visualization
        writeVoronoiClusteredLL(LL,cluster,Allthresholds)

end % cluster switch


end % file loop 