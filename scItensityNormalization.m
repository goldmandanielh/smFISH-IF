% Analysis script to do normalization using single cell intensities
%% Initiate paths and load summary data structure into matlab

homeDir='Z:\common\nliving5\For Nate\puro\';
tsData='Translation_Sites_Filtered_Final.txt';
resultsDir=fullfile(homeDir, 'Results_Updated\untreated\');

TS_summary=FISH_TS_Summary_Import(fullfile(resultsDir, tsData));
TS_summary_global=TS_summary;
%% Executes script to replace normalized data points with single cell normalized data points

for image=[1:10]
    proteinResult=fullfile(resultsDir,['C2-190731_untreated_2_' num2str(image) '_TS_outline_spots_200403.txt']);
    [cell_prop_protein, par_microscope_protein, file_names, flag_file, version,size_img,comment]=FQ_load_results_WRAPPER_v2(proteinResult,[]);
    
    outlineName=['C2-190731_untreated_2_' num2str(image) '_TS_outline.txt']; %Defines outline file name
    rowsImage=find(TS_summary.FILE==outlineName); %Finds rows that have the file name
    cells=unique(TS_summary.CELL(min(rowsImage):max(rowsImage))); %Finds unique cell IDs only in image rows
    imageTS=TS_summary(min(rowsImage):max(rowsImage),:);
    
    for cellID=[1:numel(cells)]
        cellMeanProteinAMP=mean(cell_prop_protein(cellID).spots_fit(:,4)); %mean amplitude
        cellMeanProteinXY=mean(cell_prop_protein(cellID).spots_fit(:,7));
        cellMeanProteinZ=mean(cell_prop_protein(cellID).spots_fit(:,9));
        scMeanIntInt=AMP2IntInt(cellMeanProteinAMP, cellMeanProteinXY, cellMeanProteinZ);
        
        for i=1:size(TS_summary,1)
            if TS_summary.FILE(i)==outlineName && TS_summary.CELL(i)==cells(cellID) && ...
                    TS_summary.N_IntInt(i)>0
            TS_AMP=TS_summary.AMP(i);
            TS_sigmaXY=TS_summary.sigma_xy(i);
            TS_sigmaZ=TS_summary.sigma_z(i);
            TS_summary.N_IntInt(i)=AMP2IntInt(TS_AMP,TS_sigmaXY,TS_sigmaZ)/scMeanIntInt;
            end    
        end
    end    
end
%% Plot Comparison of Single Cell and Global Mean

histogram(TS_summary_global.N_IntInt(TS_summary_global.N_IntInt>0), 'BinWidth', 1, 'normalization', 'probability')
hold on
histogram(TS_summary.N_IntInt(TS_summary.N_IntInt>0), 'BinWidth', 1, 'normalization', 'probability')
hold off

mean(TS_summary_global.N_IntInt(TS_summary_global.N_IntInt>0))
mean(TS_summary.N_IntInt(TS_summary.N_IntInt>0))
   