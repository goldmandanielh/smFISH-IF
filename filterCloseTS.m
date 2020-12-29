function [TS_summary_distFilter] = filterCloseTS(TS_summary, distMat, thresh)
%Using previously calculated distance matrix, eliminate all translation
%sites that are closer than threshold
%   Only filters translation sites that are in the SAME image


indxcloseTS=find(distMat<thresh); %returns linear index values for all distances below thresh

coordcloseTS=zeros(size(indxcloseTS,1), 2); %intiates matrix with row numbers of close sites

for i = 1:size(indxcloseTS,1)
    [pos_1, pos_2]=ind2sub(size(distMat), indxcloseTS(i, 1));
    coordcloseTS(i,:)=[pos_1, pos_2];  
end

indxNotSameCell=[];
for j = 1:size(coordcloseTS, 1) %removes pairs that are not in the same image
    if ~(TS_summary.FILE(coordcloseTS(j,1)) == TS_summary.FILE(coordcloseTS(j,2)))
       indxNotSameCell=vertcat(indxNotSameCell, j); 
    end
end

coordcloseTS(indxNotSameCell,:)=[];

rowcloseTS=vertcat(coordcloseTS(:,1), coordcloseTS(:,2)); %vert cat all row indexes from close points

TS_summary(rowcloseTS, :)=[];

TS_summary_distFilter=TS_summary;

end

