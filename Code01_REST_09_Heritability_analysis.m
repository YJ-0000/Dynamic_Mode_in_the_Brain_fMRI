clear; clc; close all;
current_path = pwd;
%%
load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm
load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_indiv_10_B

[sub_ids,sorted_idx] = sort(sub_ids);
D = D(:,sorted_idx);
B = B(:,sorted_idx);

D(1,:) = [];
D(:,tau<2000) = [];
B(:,tau<2000) = [];
sub_ids(tau<2000) = [];
% tau(tau<2000) = [];

%%
gene_data_table = readtable('data/RESTRICTED_tbk5452_10_15_2024_18_49_21.csv');
behav_data_table = readtable('data/Behavior_Data_8_18_2024_5_18_26.csv');
freesurfer_data_table = readtable('data/HCP_freesurfer_ys1j13_10_16_2024_0_33_16.csv');
for nrow = size(gene_data_table,1):-1:1
    if ~any(sub_ids==gene_data_table(nrow,'Subject').Variables)
        gene_data_table(nrow,:) = [];
    end
end
for nrow = size(behav_data_table,1):-1:1
    if ~any(sub_ids==behav_data_table(nrow,'Subject').Variables)
        behav_data_table(nrow,:) = [];
    end
end

for nrow = size(freesurfer_data_table,1):-1:1
    if ~any(sub_ids==freesurfer_data_table(nrow,'Subject').Variables)
        freesurfer_data_table(nrow,:) = [];
    end
end
gene_data_table = sortrows(gene_data_table, 'Subject');
behav_data_table = sortrows(behav_data_table, 'Subject');
freesurfer_data_table = sortrows(freesurfer_data_table, 'Subject');

nonEmptyIdx = ~cellfun(@isempty, gene_data_table.ZygositySR) & ...
              ~ismissing(gene_data_table.ZygositySR);
gene_data_table = gene_data_table(nonEmptyIdx, :);
behav_data_table = behav_data_table(nonEmptyIdx, :);
freesurfer_data_table = freesurfer_data_table(nonEmptyIdx, :);
D = D(:,nonEmptyIdx);
B = B(:,nonEmptyIdx);

age = gene_data_table.Age_in_Yrs;
sex = strcmp(behav_data_table.Gender,'F');
bmi = gene_data_table.BMI;
race1 = strcmp(gene_data_table.Race,'Am. Indian/Alaskan Nat.');
race2 = strcmp(gene_data_table.Race,'Asian/Nat. Hawaiian/Othr Pacific Is.');
race3 = strcmp(gene_data_table.Race,'Black or African Am.');
race4 = strcmp(gene_data_table.Race,'More than one');
race5 = strcmp(gene_data_table.Race,'Unknown or Not Reported');
race6 = strcmp(gene_data_table.Race,'White');
ICV = freesurfer_data_table.FS_IntraCranial_Vol;
TGMV = freesurfer_data_table.FS_Total_GM_Vol;
handedness = gene_data_table.Handedness;
edu_year = gene_data_table.SSAGA_Educ;
try 
    load results/mean_FD
catch
    load results/HCP_timeseries_noise_regressors
    mean_FD = nan(size(noise_regressors));
    for n_sub = 1:size(mean_FD,1)
        for n_ses = 1:size(mean_FD,2)
            if ~isempty(noise_regressors{n_sub,n_ses})
                mean_FD(n_sub,n_ses) = cal_mean_FD(noise_regressors{n_sub,n_ses});
            end
        end
    end
    save results/mean_FD mean_FD
end
mean_FD(tau<2000,:) = [];
mean_FD(~nonEmptyIdx,:) = [];


X = [ones(size(age,1),1),age,sex,age.*sex,age.^2,(age.^2).*sex,bmi,ICV.^(1/3),TGMV.^(1/3),race1,race2,race3,race4,race5,race6,handedness,edu_year,mean(mean_FD,2)];

nan_sub = sum(isnan(X),2)>0;
X(nan_sub,:) = [];
gene_data_table(nan_sub,:) = [];
D(:,nan_sub) = [];
B(:,nan_sub) = [];

new_table = gene_data_table(:,[1,7,8,4]);
new_table = renamevars(new_table, 'ZygositySR', 'Zygosity');


Y = [];
abs_D = abs(D(1:2:end,:));
angle_D = angle(D(1:2:end,:));
BB = B(1:2:end,:);
for ii = 1:size(abs_D,1)
%     dd = abs_D(ii,:)';
%     dd = angle_D(ii,:)';
    dd = BB(ii,:)';
    new_dd = dd - X * pinv(X) * dd;
    Y = [new_dd,Y];
end

writetable(new_table,'data/zygosity.csv');

%%
ACEfit_Par.Model     = 'ACE'; 
ACEfit_Par.P_nm      =  Y';
ACEfit_Par.InfMx     = [current_path,'/data/zygosity.csv']; 
% cd(current_path); mkdir('ACE_model_results_mag');
% ACEfit_Par.ResDir    = './ACE_model_results_mag'; 
% cd(current_path); mkdir('ACE_model_results_angle');
% ACEfit_Par.ResDir    = './ACE_model_results_angle'; 
cd(current_path); mkdir('ACE_model_results_B');
ACEfit_Par.ResDir    = './ACE_model_results_B'; 
ACEfit_Par.Subset    = [];
ACEfit_Par.Pmask     = '';                % Brain mask image (default: 
                                          % whole volume)

ACEfit_Par.Dsnmtx    = '';                % Design matrix (default: 
                                          % all-ones vector)

ACEfit_Par.Nlz       = 1;                 % Inverse Gaussian normalisation 
ACEfit_Par.AggNlz    = 0;                 % Aggregate heritability 
ACEfit_Par.NoImg     = 0;                 % If 1, suppress image-wise 
                                          % inference, and only compute
                                          % summaries.
                                        
ACEfit_Par.alpha_CFT = [];
ACEfit_Par.nPerm     = 5000;
ACEfit_Par.nBoot     = 1000; 
ACEfit_Par.nParallel = 1; 

%%% 2) Update 'ACEfit_Par' with the input data information
ACEfit_Par = PrepData(ACEfit_Par);

%%% 3) Run the original data once
ACEfit_Par = ACEfit(ACEfit_Par);

%%% 4) Add permutation and bootstrapping information, and save
%%%    "ACEfit_Par.mat"
PrepParallel(ACEfit_Par);


%
% Permutation inference for computing FWE/FDR-corrected p-values
%

if ACEfit_Par.nPerm>0
    
    %%% 1)
    %%%%% The following code can be scripted or parallelized as you wish.  
    %%%%% Please refer to "README_APACE_intro.pdf" for the example snippets
    %%%%% for parallelization.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    RunID = 1;
    ACEfit_Perm_Parallel(ACEfit_Par,RunID);
    % Inside ResDir, it will create result sets, ACEfit_Parallel_XXXX.mat,
    % one for each RunID.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 2) This will merge together all the results from the specified
    %%%    RunID's
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    ACEfit_Perm_Parallel_Results(ACEfit_Par);
    
    %%% 3)
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    ACEfit_Results(ACEfit_Par);
    
end


%
% Bootstrapping inference for constructing CIs
%

if ACEfit_Par.nBoot>0
    
    %%% 1)
    %%%%% The following code can be scripted or parallelized as you wish.
    %%%%% Please refer to "README_APACE_intro.pdf" for the example snippets
    %%%%% for parallelization.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    RunID = 1;
    ACEfit_Boot_Parallel(ACEfit_Par,RunID);
    % Inside ResDir, it will create result sets, BootCI_Parallel_XXXX.mat,
    % one for each RunID.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 2) This will merge together all the results from the specified RunID's
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    ACEfit_Boot_Parallel_Results(ACEfit_Par);
    
    %%% 3) Construct CIs   
    load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
    Boot_CIs(ACEfit_Par);
    
end


%
% Aggregate heritability (aka "Steve's method") for multiple phenotypes, 
% with P-values via permutation and CI's via boostrapping.
%
% Note that permuation (steps 1-3) and bootstrapping (steps 1-3) can be
% skipped by setting ACEfit_Par.nPerm=0 and ACEfit_Par.nBoot=0 separately; 
% the following code can be run immediately after "PrepParallel". 
%

load(fullfile(ACEfit_Par.ResDir,'ACEfit_Par.mat'));
AgHe_Method(ACEfit_Par);

% % Once with no variance normalisation
% ACEfit_Par.AggNlz = 0; % de-meaning only
% AgHe_Method(ACEfit_Par,'_NoNorm');
% 
% % Now, again, but with variance normalisation
% ACEfit_Par.AggNlz = 1; % de-meaning and scaling to have stdev of 1.0
% AgHe_Method(ACEfit_Par,'_Norm');


%
% Generate summary file 
%

APACEsummary(ACEfit_Par,'ResultSummary'); % Save results to "ResultsSummary.csv"
% APACEsummary(ACEfit_Par); % Print results to the screen

%% Violin plot

% Validate Dimensions
[num_items, num_subjects_D] = size(abs_D);
[num_subjects_table, ~] = size(new_table);

if num_subjects_D ~= num_subjects_table
    error('Number of subjects in D (%d) does not match number of rows in subjectsTable (%d).', num_subjects_D, num_subjects_table);
end

% -------------------------------------------------------------------------
% Step 2: Filter Out Subjects Without Siblings or Twins
% -------------------------------------------------------------------------

% Identify subjects sharing the same mother
[~, ~, mother_idx] = unique(new_table.Mother_ID);
mother_counts = histcounts(mother_idx, 1:max(mother_idx)+1);
sameMother = mother_counts(mother_idx) > 1;

% Identify subjects sharing the same father
[~, ~, father_idx] = unique(new_table.Father_ID);
father_counts = histcounts(father_idx, 1:max(father_idx)+1);
sameFather = father_counts(father_idx) > 1;

% Subjects with at least one shared parent
hasSiblingOrTwin = sameMother | sameFather;

% Display number of subjects before filtering
fprintf('Number of subjects before filtering: %d\n', num_subjects_D);

% Filter D and subjectsTable
% abs_D_filtered = abs_D(:, hasSiblingOrTwin);
abs_D_filtered = Y(hasSiblingOrTwin,:)';
new_table_filtered = new_table(hasSiblingOrTwin, :);

% Display number of subjects after filtering
[num_items_filtered, num_subjects_filtered] = size(abs_D_filtered);
fprintf('Number of subjects after filtering: %d\n', num_subjects_filtered);

% Update num_subjects_D to reflect the filtered data
num_subjects_D = num_subjects_filtered;

% -------------------------------------------------------------------------
% Step 2: Compute Pairwise Cosine Distances
% -------------------------------------------------------------------------

% Transpose D to have subjects as rows
D_transposed = abs_D_filtered'; % [num_subjects x num_items]

% Compute pairwise cosine distances using pdist
% 'cosine' distance in pdist is defined as 1 - cosine similarity
cosine_distances = pdist(D_transposed, 'cosine'); % [nchoosek(num_subjects,2) x 1]

% -------------------------------------------------------------------------
% Step 3: Generate All Unique Subject Pairs
% -------------------------------------------------------------------------

% Generate all unique subject pairs using nchoosek
pairs = nchoosek(1:num_subjects_D, 2); % [num_pairs x 2]

% Extract zygosity and parent IDs for each pair
zygosity1 = new_table_filtered.Zygosity(pairs(:,1));
zygosity2 = new_table_filtered.Zygosity(pairs(:,2));

mother1 = new_table_filtered.Mother_ID(pairs(:,1));
mother2 = new_table_filtered.Mother_ID(pairs(:,2));

father1 = new_table_filtered.Father_ID(pairs(:,1));
father2 = new_table_filtered.Father_ID(pairs(:,2));

% -------------------------------------------------------------------------
% Step 4: Assign Groups Based on Zygosity and Parent IDs
% -------------------------------------------------------------------------

% Initialize group labels as 'Unrelated'
groupLabels = strings(size(pairs,1),1);
groupLabels(:) = "Unrelated";

% Vectorized Assignment for Efficiency

% Identify MZ pairs: Both 'MZ' and share both parents
isMZ_pair = strcmp(zygosity1, 'MZ') & strcmp(zygosity2, 'MZ') & ...
            (mother1 == mother2) & (father1 == father2);
groupLabels(isMZ_pair) = "MZ";

% Identify DZ pairs: Both 'notMZ' and share both parents
isDZ_pair = strcmp(zygosity1, 'NotMZ') & strcmp(zygosity2, 'NotMZ') & ...
            (mother1 == mother2) & (father1 == father2);
groupLabels(isDZ_pair) = "DZ";

% Identify Sibling pairs: Not MZ or DZ, but share at least one parent
isSibling_pair = ~isMZ_pair & ~isDZ_pair & (mother1 == mother2 | father1 == father2);
groupLabels(isSibling_pair) = "Sibling";

% Remaining pairs are 'Unrelated' (already set)

% -------------------------------------------------------------------------
% Step 5: Create a Table with Group and Cosine Distance
% -------------------------------------------------------------------------

% Create a table for plotting
distanceTable = table(groupLabels, cosine_distances', 'VariableNames', {'Group','CosineDistance'});

% -------------------------------------------------------------------------
% Step 6: Generate Violin Plot
% -------------------------------------------------------------------------

groupNames = {'MZ', 'DZ', 'Sibling', 'Unrelated'};
[~, group_inx] = ismember(distanceTable.Group, groupNames);

% an alternative color scheme for some plots
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30]; 

figure('Name','Cosine Distance Box Plot','NumberTitle','off');
h = daviolinplot(distanceTable.CosineDistance,'groups',group_inx,'color',c,'xtlabels', groupNames,'violin','full');
% boxplot(distanceTable.CosineDistance, distanceTable.Group);
xlabel('Group');
ylabel('Cosine Distance');
title('Cosine Distance Distribution by Group');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);




%% functions
function mean_FD = cal_mean_FD(X)
    r = 50; % mm
    delta_trans = X(2:end,1:3)-X(1:end-1,1:3);
    delta_rot = deg2rad(X(2:end,4:6)-X(1:end-1,4:6));
    delta_rot_mm = delta_rot * r;
    FD = sum(abs(delta_trans),2) + sum(abs(delta_rot_mm),2);
    mean_FD = mean(FD);
end