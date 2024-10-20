clear;clc;
current_path = pwd;
%%

task_names = {'WM','EMOTION','MOTOR','LANGUAGE','GAMBLING','SOCIAL','RELATIONAL'};
len_time_list = [405,176,284,316,253,274,232];

for n_task = 1:length(task_names)
    current_task = task_names{n_task};
    disp(current_task);
    load(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted_meta.mat']);
    task_input = cell(length(sub_ids),2);
    for nsub = 1:length(sub_ids)
        disp(nsub);
        if strcmp(current_task,'LANGUAGE') && nsub == 615
            continue
        end
        for nses = 1:2
            % % % create input
            cd([folder_preproc_list(nsub).folder, filesep, folder_preproc_list(nsub).name]);
            cd(num2str(sub_ids(nsub)));
            cd('MNINonLinear\Results');
            try
                if nses == 1
                    cd(['tfMRI_',current_task,'_LR']);
                else
                    cd(['tfMRI_',current_task,'_RL']);
                end
                cd('EVs');
            
                switch current_task
                    case 'WM'
%                         WM_items = {'body','faces','places','tools'};
%                         onset_0bk = zeros(length(WM_items),1);
%                         duration_0bk = zeros(length(WM_items),1);
%                         for n_item = 1:length(WM_items)
%                             task_info = readtable(['0bk_',WM_items{n_item},'.txt'], "FileType","text",'Delimiter', '\t');
%                             onset_0bk(n_item) = task_info(1,1).Variables;
%                             duration_0bk(n_item) = task_info(1,2).Variables;
%                         end
%                         [onset_0bk,idx_sort] = sort(onset_0bk);
%                         duration_0bk = duration_0bk(idx_sort);
%                         task_info_0bk = table(onset_0bk,duration_0bk);
% 
%                         onset_2bk = zeros(length(WM_items),1);
%                         duration_2bk = zeros(length(WM_items),1);
%                         for n_item = 1:length(WM_items)
%                             task_info = readtable(['2bk_',WM_items{n_item},'.txt'], "FileType","text",'Delimiter', '\t');
%                             onset_2bk(n_item) = task_info(1,1).Variables;
%                             duration_2bk(n_item) = task_info(1,2).Variables;
%                         end
%                         [onset_2bk,idx_sort] = sort(onset_2bk);
%                         duration_2bk = duration_2bk(idx_sort);
%                         task_info_2bk = table(onset_2bk,duration_2bk);
% 
%                         trial_types = {'0bk','2bk'};
%                         task_input{nsub,nses} = {task_info_0bk,task_info_2bk};
                        
                        trial_types = {'0bk_body','0bk_faces','0bk_places','0bk_tools','2bk_body','2bk_faces','2bk_places','2bk_tools'};
                        temp_cell = cell(1,length(trial_types));
                        for n_trial = 1:length(trial_types)
                            task_info = readtable([trial_types{n_trial},'.txt'], "FileType","text",'Delimiter', '\t');
                            temp_cell{n_trial} = task_info;
                        end
                        task_input{nsub,nses} = temp_cell;

                    case 'EMOTION'
                        trial_types = {'fear','neut'};
                        task_info_fear = readtable('fear.txt', "FileType","text",'Delimiter', '\t');
                        task_info_neut = readtable('neut.txt', "FileType","text",'Delimiter', '\t');
                        task_input{nsub,nses} = {task_info_fear,task_info_neut};
                    case 'MOTOR'
                        trial_types = {'cue','lh','rh','lf','rf','t'};
                        temp_cell = cell(1,length(trial_types));
                        for n_trial = 1:length(trial_types)
                            task_info = readtable([trial_types{n_trial},'.txt'], "FileType","text",'Delimiter', '\t');
                            temp_cell{n_trial} = task_info;
                        end
                        task_input{nsub,nses} = temp_cell;
                    case 'LANGUAGE'
                        trial_types = {'cue','present_math','present_story','question_math','question_story','response_math','response_story'};
                        temp_cell = cell(1,length(trial_types));
                        for n_trial = 1:length(trial_types)
                            task_info = readtable([trial_types{n_trial},'.txt'], "FileType","text",'Delimiter', '\t');
                            temp_cell{n_trial} = task_info;
                        end
                        task_input{nsub,nses} = temp_cell;
                    case 'GAMBLING'
                        trial_types = {'win_event','loss_event','neut_event'};
                        temp_cell = cell(1,length(trial_types));
                        for n_trial = 1:length(trial_types)
                            task_info = readtable([trial_types{n_trial},'.txt'], "FileType","text",'Delimiter', '\t');
                            temp_cell{n_trial} = task_info;
                        end
                        task_input{nsub,nses} = temp_cell;
                    case 'SOCIAL'
                        trial_types = {'mental','rnd','mental_resp','other_resp'};
                        temp_cell = cell(1,length(trial_types));
                        for n_trial = 1:length(trial_types)
                            task_info = readtable([trial_types{n_trial},'.txt'], "FileType","text",'Delimiter', '\t');
                            temp_cell{n_trial} = task_info;
                        end
                        task_input{nsub,nses} = temp_cell;
                    case 'RELATIONAL'
                        trial_types = {'relation','match','error'};
                        temp_cell = cell(1,length(trial_types));
                        for n_trial = 1:length(trial_types)
                            task_info = readtable([trial_types{n_trial},'.txt'], "FileType","text",'Delimiter', '\t');
                            temp_cell{n_trial} = task_info;
                        end
                        task_input{nsub,nses} = temp_cell;
                end
            
            catch
                warning('Something wrong!')
                continue
            end

            
        end

    end
    cd(current_path);
    save(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted_task_input'], 'task_input','trial_types');
end