
%% Radius of Gyration (Standard Deviation) FUNCTION Script.
%Script by Laura MARTIN 2023-11
% modified by Alvaro Castells added nmperpixel variable
%% This function measure the max radius of dispersion of coordinates around a Centroid
%% as the 

function Gyr_Radius_nm = RadiusOfGyration_LM(D,nmperpixel)

            Cx = sum(D(:,1))./length(D);    % x centroid
            Cy = sum(D(:,2))./length(D);    % y centroid  
            coorX = D(:,1);                 % x list coordinates
            coorY = D(:,2);                 % y list coordinates
            xstd = sqrt((1./length(D)).*sum((coorX-Cx).^2)); % std x 
            ystd = sqrt((1./length(D)).*sum((coorY-Cy).^2)); % std y
            
            mstd= (xstd+ystd)/2;            %mean of the std x and y

            Gyr_Radius_px = sqrt((2*((mstd)^2)));  %approximated radius of dispersion (in px)
            Gyr_Radius_nm = Gyr_Radius_px * nmperpixel;   %approximated radius of dispersion (in nm)


            
end
