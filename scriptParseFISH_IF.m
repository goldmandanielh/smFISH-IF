% Analysis script for new FISH-IF Interpretation
%% Initiate paths and load summary data structure into matlab

path='Z:\users\nliving5\2020\E4.50 Repeat FISH-IF 5UTR Reporters\Results';
resultsFile=fullfile(path, 'ST_dNluc\');
summary_file_name='C3-FQ_batch_summary_NASCENT_200902.txt';

TS_summary_file=fullfile(resultsFile, summary_file_name);

TS_summary=importTS_Summary(TS_summary_file);
%% Add mRNA positions to the table and change coordinates to proper system

mRNAx=[];
mRNAy=[];

for image=[2:8 10:16]
mRNA_Result=fullfile(resultsFile, ['C1-ST_dNluc_' num2str(image) '__outline_spots_200902.txt']);
[cell_prop_RNA, par_microscope_RNA, file_names, flag_file, version,size_img,comment]=FQ_load_results_WRAPPER_v2(mRNA_Result,'');
nCell=size(cell_prop_RNA, 2);

for cell=[1:nCell]
    mRNAx_temp=[];
    mRNAy_temp=[];
    
    if size(cell_prop_RNA(cell).spots_fit,1)
    mRNAx_temp=nm2pix(cell_prop_RNA(cell).spots_fit(:,2), 107.5); %There is an issue with this coordinate system
    mRNAy_temp=nm2pix(cell_prop_RNA(cell).spots_fit(:,1), 107.5); %Y is inverted, to plot Y it is 2048-coordinate
    mRNAx=vertcat(mRNAx, mRNAx_temp);
    mRNAy=vertcat(mRNAy, mRNAy_temp);
    end
    
end
    
end

TS_summary.mRNAx_pos=mRNAx;
TS_summary.mRNAy_pos=mRNAy;
%TS_summary.y_pos=2048-TS_summary.y_pos; %Changes y positions to correct coordinate system
%% Filter positions that are far and dim
lowdistThresh=200; %distance threshold in nm for dim spots
highdistThresh=400;%distance threshold in nm for all spots
lowintensityThresh=2; %intensity threshold, will already be normalized
highintensityThresh=20; %Ignores high intensity threshold if it is really bright
absDis=400; %absolute furthest distance allowed, even for very bright spots
sigmaXY_thresh=200;
distTot=[];

for i = [1:size(TS_summary, 1)]
    mRNA_pos=[TS_summary.mRNAx_pos(i), TS_summary.mRNAy_pos(i)];
    TS_pos=[TS_summary.x_pos(i), TS_summary.y_pos(i)];
    sigmaXY=[TS_summary.sigma_xy(i)];
    
    crit_1 = pix2nm(euclidDistance(mRNA_pos, TS_pos), 107.5) > lowdistThresh;
    crit_2 = TS_summary.N_IntInt(i)<lowintensityThresh;
    crit_3 = pix2nm(euclidDistance(mRNA_pos, TS_pos), 107.5) > highdistThresh;
    crit_4 = TS_summary.N_IntInt(i)<highintensityThresh;
    crit_5 = pix2nm(euclidDistance(mRNA_pos, TS_pos), 107.5) > absDis;
    crit_6 = sigmaXY>sigmaXY_thresh;
    
    if crit_1 && crit_2  
       TS_summary.N_IntInt(i)=0;
       TS_summary.N_Amp(i)=0;
    end 
    
    if crit_3 && crit_4
       TS_summary.N_IntInt(i)=0;
       TS_summary.N_Amp(i)=0;
    end 
    
    if crit_5
       TS_summary.N_IntInt(i)=0;
       TS_summary.N_Amp(i)=0;
    end 
    
    if crit_6
       TS_summary.N_IntInt(i)=0;
       TS_summary.N_Amp(i)=0;
    end 
    
    dist=pix2nm(euclidDistance(mRNA_pos, TS_pos), 107.5);
    distTot=[distTot dist];
end

TS_distThresh=500; %using fitted positions of translation sites, eliminates translation sites that are close
TS_summary(:,24)=num2cell(distTot');
distMat=pairwiseDistanceMat(TS_summary, 107.5); %measure pairwise distance between all spots (uses translations site coordinates)
TS_summary_distanceFiltered=filterCloseTS(TS_summary, distMat, TS_distThresh);
%% Filter cells based on the number of RNA (optional)
% 
% limRNA=35;
% rowDelete=[];
% for imageID=[2:16]
%     fileName=['C3-ST_' num2str(imageID) '_TS_outline.txt']; %Defines file name
%     rowsImage=find(TS_summary_distanceFiltered.FILE==fileName); %Finds rows that have the file name
%     cellID=unique(TS_summary_distanceFiltered.CELL(min(rowsImage):max(rowsImage)));
%     
%     if imageID==1
%         rowDelete=vertcat(rowDelete, rowsImage);
%     else
%         
%     for i=1:numel(cellID)
%         rowsCell=find(TS_summary_distanceFiltered.FILE==fileName & ...
%             TS_summary_distanceFiltered.CELL==cellID(i));
%         if length(rowsCell)>limRNA
%             rowDelete=vertcat(rowDelete, rowsCell);
%         end  
%     end  
%     end
% end
% 
% TS_summary_distanceFiltered(rowDelete,:)=[];
%% Save filtered translation sites

file_name_full=fullfile(resultsFile,'Translation_Sites_Filtered_Final.txt');

fid = fopen(file_name_full,'w');  

         fprintf(fid,'FILE\tCELL\tTS\tN_IntInt\tN_PSFsup\tN_PSFsup_std\tN_SumPix\tN_Amp\tN_MaxInt\tsigma_xy\tsigma_z\tAMP\tBGD\tSize_mean[nm]\tSize_std[nm]\tBGD_cell\tTS_PixSum\tPSF_BgdSum\tTS_MaxInt\tx_pos\ty_pos\tmRNAx_pos\tmRNAy_pos\tDistance_nm\n'); 
         string_write = ['%s', repmat('\t%g', 1, 23), '\n'];
         
         for j=1:size(TS_summary_distanceFiltered,1)
            fprintf(fid, string_write, TS_summary_distanceFiltered.FILE(j), table2array(TS_summary_distanceFiltered(j,2:24)));
         end
         
fclose(fid);
