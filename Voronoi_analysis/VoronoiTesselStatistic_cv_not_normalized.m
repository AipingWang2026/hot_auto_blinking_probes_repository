%% Voronoi tessels statistic

% adjsut the paramter in the first part and then run the code
% output: figure 2x2 with histograms (left) and cdf (rigth) of all the
% single experiments (top) and the mean values (bottom)

% mean_hist contains m cell (one for each category) with mean values of
% tessels area and std 

% cvicario 20/06/2018
% update on 13/07/2018 (normalization factor density of localizations)
% update on 25/10/2018 (add standard deviation on plots)
% update on 28/11/2018 (add shade error bar plot)
% update on 25/02/2019 (double plot, tessel area & 1/tessel area)

clear all
close all

%% parameter to add:

% number of categories (each category may include many experiments)
categ =2  ;
% change jet colormap or add new colors as you wish
% Colors =jet(categ);
lightblue = [0, 0.4470, 0.7410];
orange=     [0.8500, 0.3250, 0.0980];
% mustard=    [0.9290, 0.6940, 0.1250];
% purple=     [0.4940, 0.1840, 0.5560];


% DMSO/ACTD/WAPLko/WAPLkoActD are magenta/green/blue/red
% Colors below are: blue, orange, yellow, red, green and aquamarine
%Colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880; 0 0.75 0.75];
%Colors =  [lightblue; orange; purple;lightgreen;red;grey ];
Colors =  [lightblue; orange;];
%Colors =  [0.6314    0.0039    0.7569; 0.5176    0.8902    0.0039; 0 0 1;1 0 0];
%Colors =  [0.5176    0.8902    0.0039; 0.6314    0.0039    0.7569; 0 0 1;1 0 0];
% Colors = [0.49 0.18 0.56; 0.93 0.69 0.13; 0.85 0.32 0.09 ;0 0.48 0.74];
% Colors = [0.49 0.18 0.56];

% size of the line for the single experiment
sizeline_exp = 0.5;
% size of the line for the mean
sizeline_mean = 2;

% load data .mat from voronoi analysis
for m = 1:categ
    % loop over the categories
    % load .mat file from Voronoi analysis of all the experiment in the
    % category
    DataNames{m} = uipickfiles;
    
end

% change the limits and the number of bins as you wish (try to keep high number of bins for edges_cdf )
% edges_histogram = linspace(0,0.06,100);
% edges_cdf = linspace(0,10,5000);
edges_histogram = logspace(0,5,500);
edges_cdf = logspace(0,5,500);



for m = 1:categ

    data = DataNames{1,m};
    
    for k = 1:length(data)
        
% %         USE THIS PART IF YOU LOADED .bin
%         DD = Insight3(data{1,k});
%         DD = DD.data(:,3:4);
%         [V,C] = voronoin(DD);
%         VoronoiAreas{m,1}{k} = VorArea(V,C);
%         VoronoiAreas{m,1}{k}  = VoronoiAreas{m,1}{k}*160*160;  
        

% %        % USE THIS PART IF YOU LOADED .mat 
        D = data{1,k};
        DD = importdata(D);
%       extract voronoi areas from voronoi.mat analysis file 

        VoronoiAreas{m,1}{k}  = DD.xy(:,3) ; % areas in px 
        VoronoiAreas{m,1}{k}  = VoronoiAreas{m,1}{k}*160*160;  % convert px to nm

        
        %% FIGURE 
        
        % % % % % % with tessel area and cumulative
        
%         % figure 2x2, top left
%         % overaly of all histograms of the experiments in the specified
%         % category (thin line, color of the category)
        figure(7)
        subplot(2,2,1)
        hold all
%         [N_hist{k},EDGES_hist{k}] = HistStepPlot(VoronoiAreas{m,1}{k},edges_histogram,Colors(m,:),sizeline_exp,'-');
        [N_hist{k},~] = histcounts(VoronoiAreas{m,1}{k},edges_histogram,'normalization','probability');
        % x axis is tessel area in nm2!!!
         plot(edges_histogram(1:end-1),N_hist{k},...
             '-', 'LineWidth',sizeline_exp,...
            'Color',Colors(m,:));
        set (gca, 'Fontsize',14);
        grid on;
        set(gcf,'color','w');
        title('Tessels areas histogram')
        xlabel('nm^{2}')
        ylabel('Probability')  
        set(gca, 'XScale', 'log') 
        
        % figure 2x2, top right
        % overlay of all cfd of the experiments in the specified 
        % category (thin line, color of the category)
        figure(7)
        subplot(2,2,2)
        hold all
        [N{k},~] = histcounts(VoronoiAreas{m,1}{k},edges_cdf,'normalization','cdf');
         % x axis is tessel area in nm2!!!
        plot(edges_cdf(:,1:end-1),N{k},'Color',Colors(m,:),'linewidth',sizeline_exp);
        set(gca, 'XScale', 'log') 
        set (gca, 'Fontsize',14);
        grid on;
        set(gcf,'color','w');  
        title('Tessels areas CDF')
        xlabel('nm^{2}')
        ylabel('Probability')
%     
        % % % % % % with 1/tessel area and anti - cumulative
        
%         % figure 2x2, top left
%         % overaly of all histograms of the experiments in the specified
%         % category (thin line, color of the category)
        figure(17)
        subplot(2,2,1)
        hold all
%        % x axis is 1/tessel area in nm-2!!!
         plot(1./edges_histogram(1:end-1),N_hist{k},...
             '-', 'LineWidth',sizeline_exp,...
            'Color',Colors(m,:));
        set (gca, 'Fontsize',14);
        grid on;
        set(gcf,'color','w');
        title('Inverse Tessels areas histogram')
        xlabel('nm^{-2}')
        ylabel('Probability')  
        set(gca, 'XScale', 'log') 
        
        % figure 2x2, top right
        % overlay of all cfd of the experiments in the specified 
        % category (thin line, color of the category)
        figure(17)
        subplot(2,2,2)
        hold all
      %        % x axis is 1/tessel area in nm-2!!!
        plot(1./edges_cdf(:,1:end-1),1-N{k},'Color',Colors(m,:),'linewidth',sizeline_exp);
        set(gca, 'XScale', 'log') 
        set (gca, 'Fontsize',14);
        grid on;
        set(gcf,'color','w');  
        title('Inverse Tessels areas CDF')
        xlabel('nm^{-2}')
        ylabel('Probability')


    end

    general_hist{m} = vertcat(N_hist{:});
    general_cum{m} = vertcat(N{:});
    
    % mean values and errors 
    mean_hist{m}(1,:) = mean(general_hist{m},1);
    mean_hist{m}(2,:) = std(general_hist{m},1);
    mean_hist{m}(3,:) = std(general_hist{m},1) / sqrt(size(general_hist{m},1));
    mean_hist{m}(4,:) = median(general_hist{m},1);
    mean_hist{m}(5,:) = iqr(general_hist{m},1);
    
    mean_cum{m}(1,:) = mean(general_cum{m},1);
    mean_cum{m}(2,:) = std(general_cum{m},1);
    mean_cum{m}(3,:) = std(general_cum{m},1) / sqrt(size(general_cum{m},1));
    mean_cum{m}(4,:) = median(general_cum{m},1);
    mean_cum{m}(5,:) = iqr(general_cum{m},1);

    
    % % % % % % Figure with TESSEL AREA and CUMULATIVE on x-axis
    
    % figure 2x2, bottom left
    % mean histogram of all the experiments in the specified category
    figure(7)
    subplot(2,2,3)
    hold all
    p = plot(edges_histogram(:,1:end-1),mean_hist{m}(1,:),...
        'linewidth',sizeline_mean,...
        'Color',Colors(m,:));
%      p.Color(4) = 0.7; 
%     % with std
%         p = plot(1./EDGES_hist{1}(:,1:end-1),(mean_hist{m}(1,:)+mean_hist{m}(2,:)),...
%         'linewidth',0.5,...
%         'linestyle','--',...
%         'Color',Colors(m,:));
%      p.Color(4) = 0.7;
%             p = plot(1./EDGES_hist{1}(:,1:end-1),(mean_hist{m}(1,:)-mean_hist{m}(2,:)),...
%         'linewidth',0.5,...
%         'linestyle','--',...
%         'Color',Colors(m,:));
%      p.Color(4) = 0.7;
    set (gca, 'Fontsize',14);
    grid on;
    set(gcf,'color','w');
    xlabel('nm^{2}')
    ylabel('Probability')
    set(gca, 'XScale', 'log') 
    title('Tessels areas histogram average')
    
    
    % figure 2x2, bottom right
    % mean cdf of all the experiments in the specified category
    figure(7)
    subplot(2,2,4)
    hold all
    p = plot(edges_cdf(:,1:end-1),mean_cum{m}(1,:),...
        'linewidth',sizeline_mean,...
        'Color',Colors(m,:));
%     p.Color(4) = 0.7;   
%     % with std
%         p = plot(1./EDGES{1}(:,1:end-1),1-mean_cum{m}(1,:)+mean_cum{m}(2,:),...
%         'linewidth',0.5,...
%         'linestyle','--',...
%         'Color',Colors(m,:));
%      p.Color(4) = 0.7;
%             p = plot(1./EDGES{1}(:,1:end-1),1-mean_cum{m}(1,:)-mean_cum{m}(2,:),...
%         'linewidth',0.5,...
%         'linestyle','--',...
%         'Color',Colors(m,:));
%      p.Color(4) = 0.7;
% %     
    set(gca, 'XScale', 'log') 
    set (gca, 'Fontsize',14);
    grid on;
    set(gcf,'color','w');
    xlabel('nm^{2}')
    ylabel('Probability')
    title('Tessels areas CDF average')
    
    
    
    % % % % % % Figure with 1/TESSEL AREA and ANTICUMULATIVE on x-axis
    
    % figure 2x2, bottom left
    % mean histogram of all the experiments in the specified category
    figure(17)
    subplot(2,2,3)
    hold all
    p = plot(1./edges_histogram(:,1:end-1),mean_hist{m}(1,:),...
        'linewidth',sizeline_mean,...
        'Color',Colors(m,:));
%      p.Color(4) = 0.7; 
%     % with std
%         p = plot(1./EDGES_hist{1}(:,1:end-1),(mean_hist{m}(1,:)+mean_hist{m}(2,:)),...
%         'linewidth',0.5,...
%         'linestyle','--',...
%         'Color',Colors(m,:));
%      p.Color(4) = 0.7;
%             p = plot(1./EDGES_hist{1}(:,1:end-1),(mean_hist{m}(1,:)-mean_hist{m}(2,:)),...
%         'linewidth',0.5,...
%         'linestyle','--',...
%         'Color',Colors(m,:));
%      p.Color(4) = 0.7;
    set (gca, 'Fontsize',14);
    grid on;
    set(gcf,'color','w');
    xlabel('nm^{-2}')
    ylabel('Probability')
    set(gca, 'XScale', 'log') 
    title('Inverse Tessels areas histogram average')
    
    
    % figure 2x2, bottom right
    % mean cdf of all the experiments in the specified category
    figure(17)
    subplot(2,2,4)
    hold all
    p = plot(1./edges_cdf(:,1:end-1),1-mean_cum{m}(1,:),...
        'linewidth',sizeline_mean,...
        'Color',Colors(m,:));
%     p.Color(4) = 0.7;   
%     % with std
%         p = plot(1./EDGES{1}(:,1:end-1),1-mean_cum{m}(1,:)+mean_cum{m}(2,:),...
%         'linewidth',0.5,...
%         'linestyle','--',...
%         'Color',Colors(m,:));
%      p.Color(4) = 0.7;
%             p = plot(1./EDGES{1}(:,1:end-1),1-mean_cum{m}(1,:)-mean_cum{m}(2,:),...
%         'linewidth',0.5,...
%         'linestyle','--',...
%         'Color',Colors(m,:));
%      p.Color(4) = 0.7;
% %     
    set(gca, 'XScale', 'log') 
    set (gca, 'Fontsize',14);
    grid on;
    set(gcf,'color','w');
    xlabel('nm^{-2}')
    ylabel('Probability')
    title('Inverse Tessels areas CDF average')
    
    
end


%% plot the median anticumulative with shaded error bars of standard error


for m = 1:categ

    figure(70)
    lineProps.width = 2;
    lineProps.edgestyle = '-';
    lineProps.col = {Colors(m,:)};
    
    %% CHECK THE COMBINATION YOU PREFER TO PLOT
    
    % Plot median anticumulative vs 1/tessel are with shade error SE
    mseb(1./edges_cdf(:,1:end-1),1-mean_cum{m}(4,:),(mean_cum{m}(3,:)),lineProps,0.9); 
    
        
%     % Plot median anticumulative vs 1/tessel are with shade error STD
%     mseb(1./edges_cdf(:,1:end-1),1-mean_cum{m}(4,:),(mean_cum{m}(2,:)),lineProps,0.9); 
%     
%         
%     % Plot mean anticumulative vs 1/tessel are with shade error SE
%     mseb(1./edges_cdf(:,1:end-1),1-mean_cum{m}(1,:),(mean_cum{m}(3,:)),lineProps,0.9); 
%     
%         
%     % Plot mean anticumulative vs 1/tessel are with shade error STD
%     mseb(1./edges_cdf(:,1:end-1),1-mean_cum{m}(1,:),(mean_cum{m}(2,:)),lineProps,0.9); 
%     
%     
    hold all
    set (gca, 'Fontsize',14);
    grid on;
    set(gcf,'color','w');
    xlabel('nm^{-2}')
    ylabel('Probability')
    set(gca, 'XScale', 'log') 
%     xlim([0.0005 1])
%     ylim([0 1])
    title('Inverse Tessels areas CDF median')
    
end






















