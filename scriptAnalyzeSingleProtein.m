%Script to analyze single protein intensities from samples

% path='Z:\users\nliving5\2020\E4.44 FISH-IF 5UTR Reporters\';
% conditions={'ST'};
% 
% scMeanProteinInt=[];
% imageMeanProteinInt=[];
% 
% for i=1:numel(conditions)
%     resultDir=fullfile(path, 'Results\', char(conditions(i)), '\');
%     %scMeanProteinInt_temp=[];
%     
% for image=[2:16] %Image indexes
%     proteinResult=fullfile(resultDir,['C3-ST_' num2str(image) '_TS_outline_spots_200902.txt']);
%     [cell_prop_protein, par_microscope_protein, file_names, flag_file, version,size_img,comment]=FQ_load_results_WRAPPER_v2(proteinResult,[]);
%     scMeanProteinInt_temp=[];
%     imageProtein=[];
%     
%     for cellID=[1:size(cell_prop_protein,2)]
%         cellProtein=cell_prop_protein(cellID).spots_fit(:,4);
%         cellMeanProtein=mean(cellProtein); 
%         imageProtein=vertcat(imageProtein, cellProtein);
%         scMeanProteinInt=vertcat(scMeanProteinInt, cellMeanProtein);  
%     end
%     
%     imageMeanProteinInt=vertcat(imageMeanProteinInt, mean(imageProtein));
%     
% end
% end
% 
% %boxplot(scMeanProteinInt)
% %hold on
% plotSpread(scMeanProteinInt);
% title('Mean Single Protein Intensities (Single Cell)');
% ylabel('Amplitude');
% ylim([0 3600])
% 
% figure
% plotSpread(imageMeanProteinInt);
% %xlabel('Image Index');
% ylabel('Mean Image Single Protein Amplitude');
% title('Means Single Protein Intensities (Image)')
% 
% figure
% scatter(1:length(scMeanProteinInt), scMeanProteinInt)
% xlabel('Cell Index')
% ylabel('Mean Intensity')
% title('Mean Cell Intensity Indexed by Cell')

%% Histogram of all single proteins amplitude

path='Z:\users\nliving5\2020\E4.50 Repeat FISH-IF 5UTR Reporters\';
conditions={'ST_dNluc'};

totalProtein=[];

for i=1:numel(conditions)
    resultDir=fullfile(path, 'Results\', char(conditions(i)), '\');
    
for image=[2:16]
    proteinResult=fullfile(resultDir,['C3-ST_dNluc_' num2str(image) '_TS_outline_spots_200902.txt']);
    [cell_prop_protein, par_microscope_protein, file_names, flag_file, version,size_img,comment]=FQ_load_results_WRAPPER_v2(proteinResult,[]);
    
    for cellID=[1:size(cell_prop_protein,2)]
        cellProtein=cell_prop_protein(cellID).spots_fit(:,4); 
        totalProtein=vertcat(totalProtein, cellProtein); 
    end
end

end


histogram(totalProtein, 'BinWidth', 100, 'normalization','probability')
%hold on
%histogram(ST_ZNF598_single, 'BinWidth', 100, 'normalization','probability')
% histogram(A60_ZNF598_single, 'BinWidth', 100, 'normalization','probability')
% title('Single Protein Amplitude (All Cells)')
% legend
% hold off
%% Relationship betwen nubmer of single proteins and intensity

path='Z:\users\nliving5\2020\E4.1 FISH-IF ATF4, ATF4-A60, U2TF2\200106_RQC\';
conditions={'ST_NT'};

scMeanProteinInt=[];
    
for image=[1:20]
    proteinResult=fullfile(resultDir,['C2-' char(conditions) '_' num2str(image) '_TS_outline_spots_200324.txt']);
    [cell_prop_protein, par_microscope_protein, file_names, flag_file, version,size_img,comment]=FQ_load_results_WRAPPER_v2(proteinResult,[]);
    
    scMeanProteinInt_temp=[];
    for cellID=[1:size(cell_prop_protein,2)]
        cellMeanProtein=mean(cell_prop_protein(cellID).spots_fit(:,4));
        cellNumProtein=length(cell_prop_protein(cellID).spots_fit(:,4));
        scMeanProteinInt_temp=vertcat(scMeanProteinInt_temp, [cellMeanProtein cellNumProtein]); 
    end
  
    scMeanProteinInt=vertcat(scMeanProteinInt, scMeanProteinInt_temp);
    
end

% scatter(ST(:,2), ST(:,1))
% hold on
% scatter(A60(:,2), A60(:,1))
% scatter(A60_ZNF598(:,2), A60_ZNF598(:,1))
% ylabel('Mean Single Protein Intensity')
% xlabel('# Single Proteins')
% ylim([0 2500])
%% Comparing mean protein intensity to number of mRNA

%conditions={'ST_NT'};
path='Z:\users\nliving5\2020\E4.1 FISH-IF ATF4, ATF4-A60, U2TF2\200106_RQC\';
resultDir=fullfile(path, 'Results_Updated\', char(conditions), '\');
tsData='Translation_Sites_Filtered_Final.txt';
TS_summary=FISH_TS_Summary_Import(fullfile(resultDir, tsData));
meanTSInt=[];

for image=[1:20];
    fileName=['C2-ST_' num2str(image) '_TS_outline.txt']; %Defines file name
    rowsImage=find(TS_summary.FILE==fileName); %Finds rows that have the file name
    cellID=unique(TS_summary.CELL(min(rowsImage):max(rowsImage))); %Finds unique cell IDs only in image rows
    
    for i = 1:numel(cellID) %Finds unique number of cells within file
        %meanTSInt_temp=[]; %defines mean translation site intensity matrix
        imageTS=TS_summary(min(rowsImage):max(rowsImage),:); %Takes out only rows from table that are in image
        indxNonZero=find(imageTS.CELL==cellID(i)); %Finds indexes of rows in sub table that are in each cell
        %meanTSInt_temp=mean(imageTS.N_IntInt(indxNonZero)); %Finds mean oftranslation sites
        numRNA=size(imageTS.N_IntInt(imageTS.CELL==cellID(i)),1); %Defines number of RNA
        %meanTSInt=vertcat(meanTSInt, [meanTSInt_temp numRNA]); %concatinates mean protein intensity for each cell 
    end
    
end

% scatter(ST(:,2),ST(:,1))
% hold on
% scatter(A60(:,2),A60(:,1))
% scatter(A60_ZNF598(:,2),A60_ZNF598(:,1))
% ylabel('Mean Translation Site Itensity')
% xlabel('#mRNA')
%%