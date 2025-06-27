
%% Radius of Gyration (variance) FUNCTION Script.
%Script by Chiara VICARIO 2019-02
%Revised and annotated by Victoria NEGUEMBOR 2021-06 and by Laura MARTIN 2023-10
function Gyr_px2 = RadiusOfGyration(D)

            Cx = sum(D(:,1))./length(D);    % x centroid
            Cy = sum(D(:,2))./length(D);    % y centroid  
            r = sqrt(Cx.^2+Cy.^2);          % r mean
            
            ri = sqrt((D(:,1).^2)+(D(:,2).^2)); % ri 
            
            %%Variance (in px^2)
            Gyr_px2 = (1./length(D)).*sum((ri-r).^2);
  
          
            
end
            