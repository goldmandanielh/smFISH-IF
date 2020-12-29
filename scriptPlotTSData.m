%Plotting of fitting from translation site data site data

homeDir='Z:\users\nliving5\2020\E4.50 Repeat FISH-IF 5UTR Reporters\';
tsData='Translation_Sites_Filtered_Final.txt';

%NormFactor=1.301; %ST_nLuc_BFP_AID
NormFactor=1.544;%ST_AID
resultsDir=fullfile(homeDir, 'Results\ST_dNluc\', tsData);

TS_summary=FISH_TS_Summary_Import(resultsDir);

tsInt=TS_summary.N_IntInt;
tsInt_translating=tsInt(tsInt>0);
%tsInt_2=TS_summary.N_Amp;

figure
histogram(tsInt(tsInt>0)*NormFactor, 'BinWidth', 1, 'normalization', 'probability')
hold on
%histogram(tsInt_2(tsInt_2>0), 'BinWidth', 1, 'normalization', 'probability')

mean(tsInt(tsInt>0))*NormFactor
%mean(tsInt_2(tsInt_2>0))

fracTranslating=sum(tsInt>0)/size(tsInt,1)
