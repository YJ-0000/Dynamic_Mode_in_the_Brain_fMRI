% clear;
load DMs/DM_ABIDE_cortical_subcortical_noROInorm_indiv_10

D(1,:) = [];

% ASD<->TD
% Age->ASD,TD
% ADOS
% site, sex: covariates
current_path = pwd;
fullTable=readtable('data/Phenotypic_V1_0b_preprocessed1.csv');
cd(current_path);

fileNameOriginal=fullTable.FILE_ID;
isASDarrayOriginal=fullTable.DX_GROUP;
isASDarray=zeros(length(image_file_list),1);
ageArrayOriginal=fullTable.AGE_AT_SCAN;
ageArray=zeros(length(image_file_list),1);
ADOSarrayOriginal=fullTable.ADOS_TOTAL;
ADOSarray=zeros(length(image_file_list),1);
siteArrayOriginal=fullTable.SITE_ID;
siteArray=cell(length(image_file_list),1);
sexArrayOriginal=fullTable.SEX;
sexArray=zeros(length(image_file_list),1);
for i=1:length(image_file_list)
    for j=1:length(fileNameOriginal)
        if contains(image_file_list(i).name,fileNameOriginal{j})
            isASDarray(i)=isASDarrayOriginal(j);
            ageArray(i)=ageArrayOriginal(j);
            ADOSarray(i)=ADOSarrayOriginal(j);
            siteArray(i)=siteArrayOriginal(j);
            sexArray(i)=sexArrayOriginal(j);
        end
    end
end

siteList = unique(siteArray)';
x_group=isASDarray;
x_group(x_group==2)=0; % 1:ASD, 2->0:TD
x_site=zeros(size(D,2),length(siteList));
for s=1:size(siteList,2)
    for i=1:size(D,2)
        if strcmp(siteArray{i},siteList{s})==1
            x_site(i,s)=1;
        end
    end
end

x_site(:,sum(x_site)==0) = [];

x_sex=sexArray;
x_sex(x_sex==2)=0; % 1:M, 2->0:F
%X=[x_group,x_site,x_sex];
X_prime=[x_site,x_sex,ageArray];
b=glmfit(X_prime,x_group,'normal','constant','off');
x_group_residual=x_group-X_prime*b;
X=[x_group_residual,X_prime];

%%
% permutation test on the first t-value
Y = [abs(D)',angle(D)'];
observed_tlist=zeros(size(Y,2),1);
observed_plist=zeros(size(Y,2),1);
for eig=1:size(Y,2)
    y = Y(:,eig);
    [~,~,stats]=glmfit([x_group,X_prime],y,'normal','constant','off');
    observed_tlist(eig)=stats.t(1); % only the 1st t-value is meaningful
    observed_plist(eig)=stats.p(1); % only the 1st t-value is meaningful
end
%%
nperm=10000;
perm_tlist=zeros(nperm,1);
parfor p=1:nperm
    perm_temp=zeros(size(Y,2),1);
    Xperm=[X(randperm(size(X,1)),1),X_prime];
    for eig=1:size(Y,2)
        y = Y(:,eig);
        [~,~,stats]=glmfit(Xperm,y,'normal','constant','off');
        perm_temp(eig)=stats.t(1);
    end
    perm_tlist(p)=max(abs(perm_temp)); % 2-sided
end
p_tlist=zeros(size(D,1)*2,1);
for eig=1:size(Y,2)
    p_tlist(eig)=sum( abs(observed_tlist(eig))<=perm_tlist )/nperm;
end
