for m=1:length(VoronoiAreas)
    for k=1:length (VoronoiAreas{m,1})
        median (k,m)=nanmedian(VoronoiAreas{m,1} {1,k} (:,1));
    end
end
