%%%%%%Scripts for generating fPowertional connectivity matrices using Power atlas and
%%%%%%calculating SW/clustering coefficient/path length


%%STNAMIC MATRICES GENERATION 
nsubj = [176]; %number of participants
timeseries_nos = [871];% number of timepoints 
nroi = 218 %%number of regions of interest

% Loop for all participants
for subjind = 1:nsubj
    subjind
    data = avTC{1,subjind}; %%avTC is time courses from regions of interest 
    STN_Power_Mat{1,subjind}= corrcoef(data);
    %% Power Atlases
    %number of nodes/regions 
    for no = 1:nroi
        STN_Power_Mat{1,subjind}(no,no) = 0;
    end
    STN_Power_Mat{1,subjind}(STN_Power_Mat{1,subjind}<0) = 0;
end



%%%%To calculate small worldness,etc
%%%% threshold the matrices
Thresholds = [0.1:0.01:0.3]
numbofthr = numel(Thresholds)

%Change according to no of networks 
% participants 
for  subjind = 1:nsubj%%%%%%%%%
   %number of graphs
       Data_Adults =(STN_Power_Mat{1,subjind}); 
       matrix_Adults = weight_conversion(Data_Adults, 'autofix');
       %FC
       % Calculate mean functional connectivity
       [A_FC_Adults(subjind,1), mean_neg(subjind,1), mean_abs(subjind,1)] = calculate_mean_FC(matrix_Adults);
       % for each threshold
       for  y = 1:numbofthr
           %%%%Apply thresholds to the FC matrices
           thr_Adults{subjind,y} = threshold_proportional(matrix_Adults, Thresholds(1,y)); 
                     
           %%%%Calculate characteristic path length
           len_Adults{subjind,y} = weight_conversion(thr_Adults{subjind,y}, 'lengths');
           Dist_Adults{y,subjind} = distance_wei(len_Adults{subjind,y});
           CharPath_Adults{y,subjind} = charpath(Dist_Adults{y,subjind},0,0);
           
           %%%%Calculate clustering coefficient
           norm_Adults{subjind,y} = weight_conversion(thr_Adults{subjind,y}, 'normalize');
           CC_Adults{y,subjind} = clustering_coef_wu(norm_Adults{subjind,y});
           
           %%%%Calculate small-world propensity, normalized clustering
           %%%%coefficient and normalized characteristic path length
           [A_PHI_Adults(y,subjind), delta_C_Adults(y,subjind), delta_L_Adults(y,subjind)] = small_world_propensity(thr_Adults{subjind,y});
       end
       av_A_PHI_Adults(subjind,1) = nanmean(A_PHI_Adults(:,subjind));
       av_delta_C_Adults(subjind,1) = nanmean(delta_C_Adults(:,subjind));
       av_delta_L_Adults(subjind,1) = nanmean(delta_L_Adults(:,subjind)); 
       %av_CC_Adults(subjind,1) = nanmean(CON_CC{y,n}(:,:));    
end

%%%%
for subjind = 1:nsubj
    for thr = 1:21
        AutoFix_Dist_Adults{thr,subjind} = weight_conversion(Dist_Adults{thr,subjind}, 'autofix');
    end
end

for subjind = 1:nsubj
    for roi1 = 1:nroi-1
        for roi2 = roi1+1:nroi
            for thr = 1:numbofthr
                tempDist_Adults(thr,1) = AutoFix_Dist_Adults{thr,subjind}(roi1,roi2);
            end
            avDist_Adults{subjind}(roi1,roi2) = mean(nonzeros(tempDist_Adults));
            avDist_Adults{subjind}(roi2,roi1) = avDist_Adults{subjind}(roi1,roi2);
        end
    end
end

for subjind = 1:nsubj
    AutoFix_avDist_Adults{1,subjind} = weight_conversion(avDist_Adults{subjind}, 'autofix');
end

for subjind = 1:nsubj
    for roi = 1:nroi
        tempDist = nonzeros(AutoFix_avDist_Adults{subjind}(roi,:));
        invertDist = 1./tempDist;
        numbofnodes = numel(tempDist);
        LocEff{subjind}(roi,1) = sum(invertDist)/numbofnodes;
    end
    LocEff{subjind}(isnan(LocEff{subjind}))=0;
    clear tempDist invertDist
end
