%%Script Analyze mRNA FISH intensities
%We were asked by a reviewer to show that mRNA intensity did not change
%with the translation site intensity, this script will call the FISH data
%and do this analysis.
%% Load IF Data Results
homeDir='Z:\users\nliving5\2019\E3.128 FISH-IF Repeat ST, A36, A60 and ATF4 and ATF4-A60\191212_RQC\Results_Updated\ST_A60\';
tsData='A60_2.txt';

resultsDir=fullfile(homeDir, tsData);
TS_summary=FISH_TS_Summary_Import(resultsDir);

%% Begin with loading mRNA experiment I will use the data from E3.95

mRNAItensities=[];
for image=[1:15]
mRNA_Result=fullfile(homeDir, ['C1-A60_' num2str(image) '__outline_spots_200407.txt']);
[cell_prop_RNA, par_microscope_RNA, file_names, flag_file, version,size_img,comment]=FQ_load_results_WRAPPER_v2(mRNA_Result,'');
nCell=size(cell_prop_RNA, 2);

for cellID=[1:nCell]
    mRNAIntensity_temp=[];
    mRNAIntensity_temp=cell_prop_RNA(cellID).spots_fit(:,4);
    mRNAItensities=vertcat(mRNAItensities, mRNAIntensity_temp);
end

end

A60_mRNA_2=mRNAItensities;
%% Plot mRNA Intensity Histograms

histogram(mRNAItensities,'BinWidth',100,'Normalization','Probability')

