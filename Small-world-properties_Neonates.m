%%%%%%Scripts for generating funtional connectivity matrices using Power atlas and
%%%%%%calculating small-world propensity/clustering coefficient/path length/nodal efficiency
 
nsubj = [420]; %number of participants
timeseries_nos = [1600];% number of timepoints 
nroi = [217] %%number of regions of interest

% Loop for all participants
for subjind = 1:nsubj
    subjind
    data = TimeCourses{1,subjind};
    STN_Power_Mat{1,subjind}= corrcoef(data);
    %% Power Atlases
    for no = 1:nroi
        STN_Power_Mat{1,subjind}(no,no) = 0;
    end
    STN_Power_Mat{1,subjind}(STN_Power_Mat{1,subjind}<0) = 0;
end


%%%%To calculate small-world propensity etc.
% threshold the matrices
Thresholds = [0.1:0.01:0.3];
numbofthr = numel(Thresholds);
 
for  subjind = 1:nsubj
   %%%number of graphs
       Data_Neonates =(STN_Power_Mat{1,subjind}); 
       matrix_Neonates = weight_conversion(Data_Neonates, 'autofix');
       %FC
       % this is our own script
       [A_FC_Neonates(subjind,1), mean_neg(subjind,1), mean_abs(subjind,1)] = calculate_mean_FC(matrix_Neonates);
       
       % for each threshold
       for  y = 1:numbofthr
           thr_Neonates{subjind,y} = threshold_proportional(matrix_Neonates, Thresholds(1,y));   
                      
           %%%%Calculate characteristic path length and clustering
           %%%%coefficient according to Rubinov and Sporns (2010; 2011)
           %%%%(https://sites.google.com/site/bctnet/)
           %%%%Characteristic path length
           len_Neonates{subjind,y} = weight_conversion(thr_Neonates{subjind,y}, 'lengths');
           Dist_Neonates{y,subjind} = distance_wei(len_Neonates{subjind,y});
           CharPath_Neonates{y,subjind} = charpath(Dist_Neonates{y,subjind},0,0);
           
           %%%%Clustering coefficient
           norm_Neonates{subjind,y} = weight_conversion(thr_Neonates{subjind,y}, 'normalize');
           CC_Neonates{y,subjind} = clustering_coef_wu(norm_Neonates{subjind,y});
           
           %%%%Calculate small-world propensity, normalizied clustering
           %%%%coefficient and normalized characteristic path length
           %%%%according to Muldoon et al.(2016)(http://www.seas.upenn.edu/~dsb/)
           [A_PHI_Neonates(y,subjind), delta_C_Neonates(y,subjind), delta_L_Neonates(y,subjind)] = small_world_propensity(thr_Neonates{subjind,y});
       end
       av_A_PHI_Neonates(subjind,1) = nanmean(nonzeros(A_PHI_Neonates(:,subjind)));
       av_delta_C_Neonates(subjind,1) = nanmean(nonzeros(delta_C_Neonates(:,subjind)));
       av_delta_L_Neonates(subjind,1) = nanmean(nonzeros(delta_L_Neonates(:,subjind))); 
end

%%%%To calculate nodal efficiency
for subjind = 1:nsubj
    for thr = 1:numbofthr
        AutoFix_Dist_Neonates{thr,subjind} = weight_conversion(Dist_Neonates{thr,subjind}, 'autofix');
    end
end

for subjind = 1:nsubj
    for roi1 = 1:nroi-1
        for roi2 = roi1+1:nroi
            for thr = 1:numbofthr
                tempDist_Neonates(thr,1) = AutoFix_Dist_Neonates{thr,subjind}(roi1,roi2);
            end
            avDist_Neonates{subjind}(roi1,roi2) = mean(nonzeros(tempDist_Neonates));
            avDist_Neonates{subjind}(roi2,roi1) = avDist_Neonates{subjind}(roi1,roi2);
        end
    end
end
%%%%Calculate nodal efficiency according to Achard & Bullmore(2007)
for subjind = 1:nsubj
    for roi = 1:nroi
        tempDist = nonzeros(avDist_Neonates{subjind}(roi,:));
        invertDist = 1./tempDist;
        numbofnodes = numel(tempDist);
        LocEff{subjind}(roi,1) = sum(invertDist)/numbofnodes;
    end
    LocEff{subjind}(isnan(LocEff{subjind}))=0;
    clear tempDist invertDist
end


