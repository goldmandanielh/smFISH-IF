function [distMat] = pairwiseDistanceMat(TS_summary, pixSize)
%Given TS_summary data structure, calculates pairwise distance between all
%translation sites
%points
%   Detailed explanation goes here

distMat=NaN(size(TS_summary,1), size(TS_summary,1)); %initate distance matrix

for i = 1:size(TS_summary,1)
    for j = i:size(TS_summary,1)
        pos1=[TS_summary.x_pos(i), TS_summary.y_pos(i)];
        pos2=[TS_summary.x_pos(j), TS_summary.y_pos(j)];
        distMat(i,j)=pix2nm(euclidDistance(pos1, pos2), pixSize); %calculates distance matrix in nm
        if i==j
           distMat(i,j)=NaN;
        end
    end
end


end

