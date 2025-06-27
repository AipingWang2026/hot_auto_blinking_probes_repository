
%% Radius of Gyration (variance) FUNCTION Script.
%Script by Chiara VICARIO 2019-02
%Revised and annotated by Victoria NEGUEMBOR 2021-06 and by Laura MARTIN 2023-10
%Added nmperpixel so it can be used with different microscopes Alvaro Castells 2023-12
function Gyr_nm2 = RadiusOfGyration_var(D,nmperpixel)

            Cx = sum(D(:,1))./length(D);    % x centroid
            Cy = sum(D(:,2))./length(D);    % y centroid  
            r = sqrt(Cx.^2+Cy.^2);          % r mean
            
            ri = sqrt((D(:,1).^2)+(D(:,2).^2)); % ri 
            
            %%Variance (in px^2)
            Gyr_px = (1./length(D)).*sum((ri-r).^2);
            %%Variance (in nm^2)
            Gyr_nm2 = Gyr_px * nmperpixel * nmperpixel;
  
          
            
end
            
