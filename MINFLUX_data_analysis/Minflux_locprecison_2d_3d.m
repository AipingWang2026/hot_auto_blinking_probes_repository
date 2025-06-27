% Script Description:
% This script processes all .mat files in the specified folder by:
% 1. Calculating localization precision using MINFLUX data 
% 2. Generating histograms with median precision indicators
% 3. Saving results to 'localization_precisions.csv'
% 4. Saving figures in PDF, TIFF, and FIG formats
% 5. Using sky blue colors for histograms with darker edges

% Initialize variables
resultsTable = table([], [], [], [],...
    'VariableNames', {'file_name', 'loc_precision_X', 'loc_precision_Y', 'loc_precision_Z'});

% Configure paths
dataFolder = '';
matFiles = dir(fullfile(dataFolder, '*.mat'));

% Define colors
faceColor = [0.35 0.71 0.87];    % Sky blue (RGB)
edgeColor = [0.2 0.5 0.7];       % Darker sky blue (RGB)

% Main processing loop
for k = 1:length(matFiles)
    fileName = matFiles(k).name;
    fprintf('Processing file: %s\n', fileName);
    data = load(fullfile(dataFolder, fileName));
    
    % Extract and validate localizations
    locations = squeeze(data.itr.loc(:,end,:));
    validIdx = data.vld == 1;
    locations = locations(validIdx, :);
    traceIDs = data.tid(validIdx);
    
    % Filter short traces
    [~, ~, groupIdx] = unique(traceIDs);
    counts = accumarray(groupIdx, 1);
    validTraces = counts >= 2;
    locations = locations(ismember(groupIdx, find(validTraces)), :);
    traceIDs = traceIDs(ismember(groupIdx, find(validTraces)));
    
    % Calculate statistics
    [~, ~, newGroupIdx] = unique(traceIDs);
    traceSTD = splitapply(@(x) std(x,0,1), locations, newGroupIdx);
    
    % Compute median precision (convert to nm)
    precisionX = median(traceSTD(:,1)) * 1e9;
    precisionY = median(traceSTD(:,2)) * 1e9;
    precisionZ = median(traceSTD(:,3)) * 1e9;
    
    % Update results table
    newRow = table(string(fileName), precisionX, precisionY, precisionZ,...
        'VariableNames', resultsTable.Properties.VariableNames);
    resultsTable = [resultsTable; newRow];
    
    % Save incremental results
    outputPath = fullfile(dataFolder, 'localization_precisions.csv');
    if exist(outputPath, 'file')
        existingData = readtable(outputPath);
        writetable([existingData; resultsTable], outputPath);
    else
        writetable(resultsTable, outputPath);
    end
    
    % Generate visualization
    binSize = 0.5;  % in nm
    fig = figure('Name', fileName, 'NumberTitle', 'off',...
        'Color', 'w', 'Visible', 'off',...
        'PaperPositionMode', 'auto',...
        'InvertHardcopy', 'off');
    
    if precisionZ == 0
        % 2D data visualization
        subplot(1,2,1);
        createHistogram(traceSTD(:,1)*1e9, binSize, 'X', precisionX, faceColor, edgeColor);
        
        subplot(1,2,2);
        createHistogram(traceSTD(:,2)*1e9, binSize, 'Y', precisionY, faceColor, edgeColor);
    else
        % 3D data visualization
        subplot(1,3,1);
        createHistogram(traceSTD(:,1)*1e9, binSize, 'X', precisionX, faceColor, edgeColor);
        
        subplot(1,3,2);
        createHistogram(traceSTD(:,2)*1e9, binSize, 'Y', precisionY, faceColor, edgeColor);
        
        subplot(1,3,3);
        createHistogram(traceSTD(:,3)*1e9, binSize, 'Z', precisionZ, faceColor, edgeColor);
    end
    
    % Save figures in multiple formats
    [~, baseName] = fileparts(fileName);
    savePath = fullfile(dataFolder, baseName);
    
    % Save as FIG
    saveas(fig, [savePath '.fig']);
    
    % Save as PDF (vector format)
    print(fig, [savePath '.pdf'], '-dpdf', '-painters', '-r300')
    
    % Save as TIFF (300 dpi)
    print(fig, [savePath '.tiff'], '-dtiff', '-r300')
    
    close(fig);  % Close figure after saving
end

% Helper function for histogram creation
function createHistogram(data, binWidth, axisLabel, medianValue, faceColor, edgeColor)
    histogram(data,...
        'BinWidth', binWidth,...
        'FaceColor', faceColor,...
        'EdgeColor', edgeColor,...
        'LineWidth', 0.5);
    
    hold on;
    xline(medianValue,...
        'LineStyle', '--',...
        'Color', 'r',...
        'LineWidth', 2);
    
    title(sprintf('%s Axis Precision: %.2f nm', axisLabel, medianValue));
    xlabel('Standard Deviation (nm)');
    ylabel('Frequency');
    xlim([0, max(20, prctile(data, 95))]);
    grid on;
    set(gca, 'FontSize', 10, 'LineWidth', 1);
end