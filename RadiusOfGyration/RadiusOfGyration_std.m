
%% Radius of Gyration (Standard Deviation) FUNCTION Script.
%Script by Laura MARTIN 2023-10
%Added nmperpixel so it can be used by different microscopes Alvaro Castells 2023-12

function Gyr_std_nm = RadiusOfGyration_std(D,nmperpixel)

            Cx = sum(D(:,1))./length(D);    % x centroid
            Cy = sum(D(:,2))./length(D);    % y centroid  
            r = sqrt(Cx.^2+Cy.^2);          % r mean
            
            ri = sqrt((D(:,1).^2)+(D(:,2).^2)); % ri 
            
            %%Std (in px)
            Gyr_std_px = sqrt((1./length(D)).*sum((ri-r).^2));  
            %%Std (in nm)
            Gyr_std_nm = Gyr_std_px * nmperpixel;
           


            
end
