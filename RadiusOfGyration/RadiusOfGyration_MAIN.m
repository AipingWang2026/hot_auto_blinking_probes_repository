
%% RADIUS OF GYRATION MAIN SCRIPT. Computes radius of gyration (in um) of localizations in an image or selected ROI
%Script by Chiara VICARIO 2019-02 and Laura MARTIN 2023-10
%Revised and annotated by Victoria NEGUEMBOR 2021-06
%Added create results table and save in .txt and .xlsx by Laura MARTIN 2023-10
%Added measurement of G_std (=root square of previous G) and Gyr_Radius by Laura MARTIN 2023-10
%Changed to analyse files from 1 Condition only per Run. 
%It creates a subfolder with the desired name where xlsx is saved. 
%%by Laura MARTIN 2023-11
%Modified Alvaro Castells added variable nmperpixel in it and the dependent scripts so it can be used with different microscopes
%Di modified subfolder and add parameter nmperpixel in each functions
close all
clear all 
%% parameter to add:

% % Change 'CONDITION_1' with the desired name in both code lines.
% mkdir('CONDITION_1');                                 % % Creates subfolder. 

SubFolder = 'C:\Users\HP\Desktop\code WAP';

categ = 1;
nmperpixel= 160;
% STORM  = 160 nm per pixel; ONI = 117 nm per pixel

% load data .bin of the coordinates 

for m = 1:categ
    % loop over the categories
    % load .bin file
    DataNames{m} = uipickfiles;
    %List = transpose(DataNames);
end

for m = 1:categ

    data = DataNames{1,m};
    
    for k = 1:length(data)
        D = Insight3(data{1,k});
        D = D.data(:,3:4);
        G_old_px2{m,1}(1,k) = RadiusOfGyration(D);
        G_var_nm2{m,1}(1,k) = RadiusOfGyration_var(D,nmperpixel);
        G_std_nm{m,1}(1,k) = RadiusOfGyration_std(D,nmperpixel);
        Gyr_Radius_nm{m,1}(1,k) = RadiusOfGyration_LM(D,nmperpixel);
        
    end
end

        Rg_old_px2 = transpose(G_old_px2{m,1});
        Rg_var_nm2 = transpose(G_var_nm2{m,1});
        Rg_std_nm = transpose(G_std_nm{m,1});
        R_GYR_nm = transpose(Gyr_Radius_nm{m,1});

Datanames = transpose(DataNames{m});
Radius_of_Gyration = table (Datanames, Rg_old_px2, Rg_var_nm2, Rg_std_nm, R_GYR_nm);


writetable(Radius_of_Gyration,strcat(SubFolder,'\Radius_of_Gyration.csv'));
%writetable(Radius_of_Gyration,strcat(SubFolder,'\CONDITION1\_Radius_of_Gyration.txt'),'Delimiter','tab');

%% 
