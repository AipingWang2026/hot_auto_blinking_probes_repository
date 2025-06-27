%This function takes several conditions and creates
% 1 Plot of precision distribution of localizations
% Numerical data of loc precision for all files, divided by name
%1 excel with mean and median loc precision, number of localizations and
%precision
%Alvaro Castells Garcia 2024-04



clear all
close all




%% Colors
% Here i am defining some colors (self-explanatory) so when plotting we can
% just put the names of the color variable. To be used with the variable
% "colors" in next part
lightblue = [0, 0.4470, 0.7410];
orange=     [0.8500, 0.3250, 0.0980];
mustard=    [0.9290, 0.6940, 0.1250];
purple=     [0.4940, 0.1840, 0.5560];
lightgreen= [0.4660, 0.6740, 0.1880];
skyblue =   [0.3010, 0.7450, 0.9330];
granate=    [0.6350, 0.0780, 0.1840];
blue =      [0, 0, 1];
darkgreen = [0, 0.5, 0];
red =       [1, 0, 0];
magenta =   [0.75, 0, 0.75];
grey =      [0.25, 0.25, 0.25];

%% Variables
% Here we are setting the variables for the computation 
categ = 1;  % number of conditions (for example DMSO and ActD)
% catName = [{'HSVMock'},{'HSV2hpi'}, {'HSV3hpi'}, {'HSV8hpi'}]; %modify names of categories. Due to how matlab writes tables, please avoid spaces
catName = [{'SelectedCategory'},]; %modify names of categories. Due to how matlab writes tables, please avoid spaces
colors = [lightblue; orange; purple; lightgreen];
pix2nm = 117; %Pixels to nm value. Oni 117, Nstorm 160.  
Directory = ''; %Folder to save the data in 

%% Variables for histogram
Precisionhistlimnm =50; %limit histogram of loc precision, in nm.
Precisionhistlimstep = 1; %steps histogram of loc precision, in nm.

%Locperframehistlimnm = 10000; %limit histogram of loc per frame, in number of frames.
%Locperframehistlimstep = 10; %steps histogram of loc precision, in number of frames.
locperframebinwidth = 50; %Bin widtch for localizations per frame. Now it is set at 10 frames. Important that all histograms should have the same number of frames, otherwise it will break


%% Script variables

%One matrix in which each line is the Histogram results of localizations
%per frame. Increments of 50, abolute numbers

Histallcellslocpercell = [];
%edgeslocpercell = 10:Locperframehistlimstep:Locperframehistlimnm; % frames

%One matrix in which each line is the Histogram results of loc precision in
%increments of 5 nm from 0 to 300 nm, distributed as probability
%per frame
Histallcellslocprecisionpercell = [];
edgeslocprecision = 0:Precisionhistlimstep:Precisionhistlimnm; % nm 

Meanlocpres = [];
Medianlocpres = [];
Nlocperframe= [];
EDGESlocperframe = [];
Nlocprecision =[];
EDGESlocprecision =[];



%% Load data
% load .bin data of localizations 
% Note: You can use the selected ROIs or the raw ROIs, it should not be
% that different. But it should be better if there was only one file per
% cell, as otherwise the number of loc per frame will vary wildly

for m = 1:categ
    % loop over the categories
    % load .bin

    Dataname{m} = uipickfiles;
    
end

for m = 1:categ

    data = Dataname{1,m};
    
    for k = 1:length(data)
            
        DD = Insight3(data{1,k});
        % take the Loc Precision Column
        %%%%%%VERYIMPORTANTNOTE: SEE LOC PRECISION COLUMN, SI NO ESTA AHI HAY QUE VOLVER A CONVERTIR A .BIN
        LocPrecSingleCell = (DD.data(:,(6)))./20; %Column 6 is the results of the loc precision in nm of x and y multiplied by 10. If we divide it by 20 we eliminate the multiplication and we have the average of x and y. 
        FrameSingleCell = DD.data(:, 14);
        
        %Hist loc per frame
        [Nlocperframe(k,:),EDGESlocperframe(k,:)] = histcounts(FrameSingleCell,'BinWidth',locperframebinwidth);
        
        %Hist locprecision
        [Nlocprecision(k,:),EDGESlocprecision(k,:)] = histcounts(LocPrecSingleCell,edgeslocprecision,'normalization','probability');
        %[Nlocprecision{1,m}{k,1},EDGESlocprecision{1,m}{k,1}] = histcounts(D32{k},edgeslocprecision,'normalization','probability');

        %Mean loc precision
        Medianlocpres(k,1) = median(LocPrecSingleCell);
        Meanlocpres (k,1) = mean(LocPrecSingleCell);
        %Median loc precision
        
        
        %Make an histogram of Loc Pres single cell
        %Save the histogram as a variable
        %Make an histogram for Frames in single cell
        
        
    end
end

%% Histogram average
    
   %Make an average of several histograms of Loc Pres single cell
   
   Averagelocperframe = mean(Nlocperframe, 1);
   Averagelocprecision = mean(Nlocprecision, 1);
   
   StdevLocperframe = std(Nlocperframe, 1);
   StdevLocprecision = std(Nlocperframe, 1);
   
   %% Create table and save table
   
Cellnames = transpose(Dataname{1,1});
Data_summary = table (Cellnames, Medianlocpres, Meanlocpres);
Data_histlocprec = table (Cellnames, Nlocprecision);
Data_histperframe = table (Cellnames, Nlocperframe);
writetable(Data_summary,strcat(Directory,'\Loc_precission_summary.xlsx'));
writetable(Data_histlocprec,strcat(Directory,'\Loc_precission_histogram.xlsx'));
writetable(Data_histperframe,strcat(Directory,'\Locperframe_summary.xlsx'));
   
   
   
   
   
   
   
   
   
   
   
   
   
% %    
% %    
% %    
% %    
% %    
% %    
% %    
% %    
% %    %Make an average of several histogrames for frame precision in single cell
% %         
% %         
% %         
% %         
% %         
% %         % generate density map (20 refers to size of SR pixel)
% %         [Dens] = QuickDensity(D(:,1),D(:,2),30,minX,maxX,minY,maxY);
% %         % smooth density map with gaussian filter ("2" refers to gaussian sigma) 
% %         h = fspecial('gaussian',[10 10], 5);
% %         Dens = filter2(h, Dens);
% %         % binarize DNA (you can adjust threshold (ie. 0.001) depending on
% %         % your dataset,scan different threshold values to find appropriate) 
% %         MaskDNA = imbinarize(Dens,'adaptive','Sensitivity',1);
% % %         figure, imshow(MaskDNA)
% % %       dilate the mask
% %         se = strel('disk',18);
% %         Filled = imdilate(MaskDNA, se);
% %         % fill the holes        
% %         Filled = imfill(Filled,'holes');
% %         % erode the mask
% %         Filled = imerode(Filled, se);
% %         
% % % %         COMMENT part below to skip visualization of density maps and masks generated,
% % %           figure, subplot(1,3,1),imagesc(Dens), axis equal
% % %                      subplot(1,3,2),imagesc(MaskDNA),axis equal
% % %                      subplot(1,3,3),imagesc(Filled),axis equal
% %                                         
% %                    % m 
% %                    % k
% %                     
% %         percentage_black{m,k} = 100-((sum(sum(MaskDNA))*100)/sum(sum(Filled)));
% % 
% %     end 
% %     end
% % 
% % 
% %     
% %            DD = Insight3(data{1,k});
% %         % take the Loc Precision Column
% %         %%%%%%VERYIMPORTANTNOTE: SEE LOC PRECISION COLUMN, SI NO ESTA AHI HAY QUE VOLVER A CONVERTIR A .BIN
% %         LocPrecSingleCell = DD.data(:,COLUMNAAAA);
% %         FrameSingleCell = DD.data(:, columnaframe);
% %         
% %         %Hist loc per frame
% %         [Nlocperframe{1,m}{k,1},EDGESlocperframe{1,m}{k,1}] = histcounts(FrameSingleCell,edgeslocpercell);
% %         
% %         
% %         %Hist locprecision
% %         [Nlocprecision{1,m}{k,1},EDGESlocprecision{1,m}{k,1}] = histcounts(LocPrecSingleCell,edgeslocprecision,'normalization','probability');
% %         %[Nlocprecision{1,m}{k,1},EDGESlocprecision{1,m}{k,1}] = histcounts(D32{k},edgeslocprecision,'normalization','probability');
% % 
% %         %Mean loc precision
% %         Medianlocpres(1,k) = median(LocPrecSingleCell);
% %         Meanlocpres (1,k) = mean(LocPrecSingleCell);
% %         %Median loc precision
