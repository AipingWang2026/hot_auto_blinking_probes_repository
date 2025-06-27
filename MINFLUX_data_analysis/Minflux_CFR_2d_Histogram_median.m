% Enhanced script with sky blue bars 
% Load data file
data_path = '';
load(data_path); 

% Extract 5th column of cfr
cfr_column5 = itr.cfr(:,5); 
median_value = median(cfr_column5); 

% Create histogram with customized bars
figure;
h = histogram(cfr_column5,...
    'EdgeColor', [0.1 0.4 0.7],...      
    'FaceColor', [0.35 0.75 1],...      
    'LineWidth', 1.2,...                
    'BinWidth', 0.8*range(cfr_column5)/30); 

title('Distribution of CFR Column 5 Data');
xlabel('CFR Column 5 Values');
ylabel('Occurrence (Counts)');
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
box off;
grid off;

% Add median features
hold on;
line([median_value median_value], ylim, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
text(0.8, 0.8, sprintf('Median = %.2f', median_value),...  
    'Units', 'normalized',...
    'Color', 'r',...
    'FontSize', 12,...
    'HorizontalAlignment', 'left');

% Save configurations
[filepath, name, ~] = fileparts(data_path);
output_filename = fullfile(filepath, [name '_cfr']);

set(gcf, 'Renderer', 'painters');
saveas(gcf, [output_filename '.fig']);
print(gcf, [output_filename '.pdf'], '-dpdf');
print(gcf, [output_filename '.tif'], '-dtiff', '-r300');

disp(['Files saved to: ' filepath]);