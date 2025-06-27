% Define file path and name
filePath = '';
fileName = '.csv';
fullPath = fullfile(filePath, fileName);

% Read CSV data
data = readtable(fullPath);

% Extract required columns
photons = data{:, 14};         % Column 14: psf-photon-count
xPrecision = data{:, 29};      % Column 29: X precision (nm)
yPrecision = data{:, 30};      % Column 30: Y precision (nm)
zPrecision = data{:, 31};      % Column 31: Z precision (nm)

% Calculate mean values
meanPhotons = mean(photons, 'omitnan');
meanX = mean(xPrecision, 'omitnan');
meanY = mean(yPrecision, 'omitnan');
meanZ = mean(zPrecision, 'omitnan');

% Calculate median values
medianPhotons = median(photons, 'omitnan');
medianX = median(xPrecision, 'omitnan');
medianY = median(yPrecision, 'omitnan');
medianZ = median(zPrecision, 'omitnan');

% Display calculated means
fprintf('Mean Photons = %.2f\n', meanPhotons);
fprintf('Mean X precision = %.2f nm\n', meanX);
fprintf('Mean Y precision = %.2f nm\n', meanY);
fprintf('Mean Z precision = %.2f nm\n', meanZ);

% Generate output filenames base
[path, name, ~] = fileparts(fullPath);

% Automatically calculate optimal bin width for histograms
calculateOptimalBinWidth = @(data) range(data(~isnan(data))) / (1 + 3.322 * log10(sum(~isnan(data))));

% Adjustment factors (smaller than 1 makes bins thinner, larger than 1 makes bins thicker)
photonAdjustment = 0.8;      % Adjustment factor for photon count histogram
xyPrecisionAdjustment = 0.7; % Adjustment factor for X/Y precision histograms
zPrecisionAdjustment = 0.8;  % Adjustment factor for Z precision histogram (thinner by 20%)

% Calculate optimal bin widths with adjustments
photonBinWidth = calculateOptimalBinWidth(photons) * photonAdjustment;
xPrecisionBinWidth = calculateOptimalBinWidth(xPrecision) * xyPrecisionAdjustment;
yPrecisionBinWidth = calculateOptimalBinWidth(yPrecision) * xyPrecisionAdjustment;
zPrecisionBinWidth = calculateOptimalBinWidth(zPrecision) * zPrecisionAdjustment;

% 1. Create and save Photons histogram
fig = figure('Visible', 'off');
histogram(photons, 'BinWidth', photonBinWidth);
hold on;

yLimits = ylim;
line([medianPhotons, medianPhotons], yLimits, ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
xlabel('Photons');
ylabel('Frequency');
title('Distribution of Photons');
dim = [0.7 0.8 0.1 0.1];
str = sprintf('Mean = %.2f', meanPhotons);
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
           'BackgroundColor', 'white', 'EdgeColor', 'none', ...
           'FontSize', 12, 'FontWeight', 'bold');
set(gcf, 'Position', [100, 100, 800, 600]);
outputBase = fullfile(path, [name ' Photons']);
print([outputBase '.pdf'], '-dpdf', '-bestfit');
saveas(gcf, [outputBase '.fig']);
print([outputBase '.tiff'], '-dtiff', '-r300');
close(fig);

% 2. Create and save X precision histogram
fig = figure('Visible', 'off');
histogram(xPrecision, 'BinWidth', xPrecisionBinWidth);
hold on;

yLimits = ylim;
line([medianX, medianX], yLimits, ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
xlabel('X precision (nm)');
ylabel('Frequency');
title('Distribution of X precision');
dim = [0.7 0.8 0.1 0.1];
str = sprintf('Mean = %.2f', meanX);
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
           'BackgroundColor', 'white', 'EdgeColor', 'none', ...
           'FontSize', 12, 'FontWeight', 'bold');
set(gcf, 'Position', [100, 100, 800, 600]);
outputBase = fullfile(path, [name ' X precision']);
print([outputBase '.pdf'], '-dpdf', '-bestfit');
saveas(gcf, [outputBase '.fig']);
print([outputBase '.tiff'], '-dtiff', '-r300');
close(fig);

% 3. Create and save Y precision histogram
fig = figure('Visible', 'off');
histogram(yPrecision, 'BinWidth', yPrecisionBinWidth);
hold on;

yLimits = ylim;
line([medianY, medianY], yLimits, ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
xlabel('Y precision (nm)');
ylabel('Frequency');
title('Distribution of Y precision');
dim = [0.7 0.8 0.1 0.1];
str = sprintf('Mean = %.2f', meanY);
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
           'BackgroundColor', 'white', 'EdgeColor', 'none', ...
           'FontSize', 12, 'FontWeight', 'bold');
set(gcf, 'Position', [100, 100, 800, 600]);
outputBase = fullfile(path, [name ' Y precision']);
print([outputBase '.pdf'], '-dpdf', '-bestfit');
saveas(gcf, [outputBase '.fig']);
print([outputBase '.tiff'], '-dtiff', '-r300');
close(fig);

% 4. Create and save Z precision histogram
fig = figure('Visible', 'off');
histogram(zPrecision, 'BinWidth', zPrecisionBinWidth);
hold on;

yLimits = ylim;
line([medianZ, medianZ], yLimits, ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
xlabel('Z precision (nm)');
ylabel('Frequency');
title('Distribution of Z precision');
dim = [0.7 0.8 0.1 0.1];
str = sprintf('Mean = %.2f', meanZ);
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
           'BackgroundColor', 'white', 'EdgeColor', 'none', ...
           'FontSize', 12, 'FontWeight', 'bold');
set(gcf, 'Position', [100, 100, 800, 600]);
outputBase = fullfile(path, [name ' Z precision']);
print([outputBase '.pdf'], '-dpdf', '-bestfit');
saveas(gcf, [outputBase '.fig']);
print([outputBase '.tiff'], '-dtiff', '-r300');
close(fig);

disp('All histograms created and saved successfully!');