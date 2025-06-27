



% loadONIcsv

% script to transform an ONI csv file into a BIN file for Insight3
%Modified by Alvaro on 17/07/20
%% Data location
% A: This was hard coded but we added so it so you can select the file. Select the .csv.
% dataLoc = 'E:\Code\MATLAB_NonCode\ONI';
% dataName = 'ptfLC3_Mcherry tubulin_076.csv';
%210304_Mdified for the sigma values so it can be corrected for filtering
%Paint localizations
%211021_Areamodified to be a multiplication of the precision x10, as
%otherwise the value is too small
%       Intensity modified?as it derives from the area
[dataName, dataLoc, ~] = uigetfile('*.csv','select the data file','multiselect','on');
%3d modified

for ii = 1:length(dataName)
data(ii) = fullfile(dataLoc,dataName(ii));

nmperpixel = 117; % A: Change this to adjust the nm per pixel. if the csv is already in pixels comment this part

M = csvread(data{ii},1,0);
X = M(:,3)/nmperpixel;
Y = M(:,4)/nmperpixel;
Z = M(:,5)/nmperpixel;
Channel = M(:,1);
Xc = M(:,3)/nmperpixel;
Yc = M(:,4)/nmperpixel;
Zc = M(:,5)/nmperpixel;
Height = M(:,13);
Area = (M(:,6)+ M(:,7)).*10;
%Area = (ones([(size(M,1)) 1])).*10000;
%Area = M(:,15)*nmperpixel;
%Width = (M(:,15).*M(:,16)) *nmperpixel;
%Width = (ones([(size(M,1)) 1])).*300;
%Width = (ones([(size(M,1)) 1])).*300;
Width = ((M(:,15))).*nmperpixel;
Phi = zeros([(size(M,1)) 1]);
Ax = ones([(size(M,1)) 1]);
BG = M(:,16);
I = Area;
Frame = (M(:,2))+1;
Length = ones([(size(M,1)) 1]);
Link = (ones([(size(M,1)) 1])).*-1;
Valid = ones([(size(M,1)) 1]);

% validInd = logical(M(:,31));
% X = X(validInd);
% Y = Y(validInd);
% Z = Z(validInd);
i3 = Insight3();
%i3.setData( [Channel,X,Y,Xc,Yc,Height,Area,Width,Phi,Ax,BG,I,Frame,Length,Link,Valid,Z,Zc] );
i3.setData( [X,Y,Xc,Yc,Height,Area,Width,Phi,Ax,BG,I,Channel,Valid,Frame,Length,Link,Z,Zc] );
i3.write(fullfile(dataLoc,[dataName{ii}(1:end-4) '.bin']));
end
 