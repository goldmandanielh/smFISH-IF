resultsDir='Z:\users\nliving5\2020\E4.44 FISH-IF 5UTR Reporters\Results\ST\';
tsData='Translation_Sites_Filtered_Final.txt';
TS_summary=FISH_TS_Summary_Import(fullfile(resultsDir,tsData));


fracTranslating=[];

for image=[2:16];
    fileName=['C3-ST_' num2str(image) '_TS_outline.txt']; %Defines file name
    rowsImage=find(TS_summary.FILE==fileName); %Finds rows that have the file name
    cellID=unique(TS_summary.CELL(min(rowsImage):max(rowsImage))); %Finds unique cell IDs only in image rows
    imageTS=TS_summary(min(rowsImage):max(rowsImage),:); %Takes out only rows from table that are in image
    
    for i = 1:numel(cellID) %Finds unique number of cells within file
        fracTranslating_temp=[];
        indxCell=find(imageTS.CELL==cellID(i)); %Finds indexes of rows in sub table that are in each cell
        numRNA=size(indxCell,1); %Defines number of RNA in cell
        translating=(imageTS.N_IntInt(indxCell)>0);
        fracTranslating_temp=sum(translating)/numRNA;
        
        if numRNA>5
        fracTranslating=vertcat(fracTranslating, fracTranslating_temp);
        end
        
    end
    
end

plotSpread(fracTranslating)