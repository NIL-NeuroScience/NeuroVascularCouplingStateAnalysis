%% load data
% data root folder
root = '/projectnb/devorlab/bcraus/HRF/1P/23-10-25/Thy1_215/';
% run number
Run = 1;

% load data
root_processed = [root filesep 'processed'];
rfp_norm=h5read(sprintf('%s%srun%04i.h5',root_processed ,filesep,Run),'/rfp/norm');
gfp_normHD=h5read(sprintf('%s%srun%04i.h5',root_processed ,filesep,Run),'/gfp/normHD');
HbO=h5read(sprintf('%s%srun%04i.h5',root_processed ,filesep,Run),'/hemodynamics/HbO');
HbR=h5read(sprintf('%s%srun%04i.h5',root_processed ,filesep,Run),'/hemodynamics/Hb');
HbT=HbR+HbO;
load([root filesep 'Triggers' filesep sprintf('Run%03i.mat',Run)]);
load([root 'dataIn.mat']); 

% load signals
load([root '/DataAnalysis/Run0' num2str(Run) '/Signals.mat'])
jRGECO= rfp_norm;
GRAB = gfp_normHD;

% Save results in this folder
folder_to_save = '/projectnb/devorlab/skura/HRF/cross_corr_results/run_1_231025_Thy1_215_xcorrfft';
if ~isfolder(folder_to_save)
    mkdir(folder_to_save); 
end


%% Load brain_mask if it is in the results_save folder. Otherwise make a brain mask.

% load brain_mask
if exist([folder_to_save filesep 'brain_mask.mat'],'file')
    load([folder_to_save filesep 'brain_mask.mat']);
else
    img = dataIn(Run).template;
    img = img(:,:,1);
    brain_mask = zeros(size(img));
    h = msgbox(' Please select one side of the brain region');
    uiwait(h)
    figure; colormap('gray');
    while(1)
        imagesc(img); axis image;
        [brain_mask1,Xi1,Yi1] = roipoly;
        hold on;
        plot(Xi1,Yi1,'color','k');
        hold off;
        button = questdlg('Are you statisfied with ROI');
        if strcmp(button, 'Yes')
          break;
        end
    end
    h = msgbox(' Please select other side of the brain region');
    uiwait(h)
    while(1)
        imagesc(img); axis image; hold on;
        plot(Xi1,Yi1,'color','k');
        hold off;
        [brain_mask2,Xi2,Yi2] = roipoly;
        hold on;
        plot(Xi2,Yi2,'color','k');
        hold off;
        button = questdlg('Are you statisfied with ROI');
        if strcmp(button, 'Yes')
          break;
        end
    end
    imagesc(img); hold on;
    plot(Xi1,Yi1,'color','k');
    plot(Xi2,Yi2,'color','k');
    hold off;
    brain_mask(brain_mask1==1) = 1;
    brain_mask(brain_mask2==1) = 1;

    save([folder_to_save filesep 'brain_mask.mat'],'brain_mask1','Xi1','Yi1','brain_mask2','Xi2','Yi2','brain_mask')
end

%% Get lag and cross corrlearion between Ca and HbT and cluster lag and correation to define states

% Spatial smoothing of the data
smooth_scale = 3;
HbT_smoothed = smooth2d(HbT,smooth_scale);
jRGECO_smoothed = smooth2d(jRGECO ,smooth_scale);
GRAB_smoothed = smooth2d(GRAB ,smooth_scale);

[sX, sY, sZ] = size(HbT);
% apply filtering
brain_idx = find(brain_mask == 1);
HbT_smoothed = reshape(HbT_smoothed,[sX*sY sZ]);
HbT_smoothed_brain = HbT_smoothed(brain_idx,:);

fc= 0.01;
fs = dataIn(Run).settings.fs;
[b,a] = butter(6,fc/(fs/2));
HbT_lowpass_brain =  filtfilt(b,a,HbT_smoothed_brain');
HbT_higpass_brain = HbT_smoothed_brain'-HbT_lowpass_brain;


jRGECO_smoothed = reshape(jRGECO_smoothed,[sX*sY sZ]);
jRGECO_smoothed_brain = jRGECO_smoothed(brain_idx,:);
jRGECO_lowpass_brain =  filtfilt(b,a,jRGECO_smoothed_brain');
jRGECO_higpass_brain = jRGECO_smoothed_brain'-jRGECO_lowpass_brain;

GRAB_smoothed = reshape(GRAB_smoothed,[sX*sY sZ]);
GRAB_smoothed_brain = GRAB_smoothed(brain_idx,:);
GRAB_lowpass_brain =  filtfilt(b,a,GRAB_smoothed_brain');
GRAB_higpass_brain = GRAB_smoothed_brain'-GRAB_lowpass_brain;

analyse = 'highpass';
% analyse = 'lowpass';
% analyse = 'allpass';
sr = dataIn(Run).settings.fs;
time_shift = 2;
corr_time_length = 30;
idx_shift = round(time_shift*sr);
n_steps = floor(size(jRGECO,3)/idx_shift);

max_corr_oneside_brain = zeros(length(brain_idx),n_steps);
max_corr_lag_oneside_brain = zeros(length(brain_idx),n_steps);


for w = 1:n_steps
    idx_start = floor(1+idx_shift*(w-1));
    idx_end = min(floor(idx_start+sr*30),sZ);
    if strcmp(analyse, 'allpass')
        Ca_act_30sec = smoothdata(jRGECO_smoothed_brain(idx_start:idx_end,:),1,"movmean",6);
        HbT_act_30sec = smoothdata(HbT_smoothed_brain(idx_start:idx_end,:),1,"movmean",6);
    elseif strcmp(analyse, 'highpass')
        Ca_act_30sec = smoothdata(jRGECO_higpass_brain(idx_start:idx_end,:),1,"movmean",6);
        HbT_act_30sec = smoothdata(HbT_higpass_brain(idx_start:idx_end,:),1,"movmean",6);
    else
        Ca_act_30sec = smoothdata(jRGECO_lowpass_brain(idx_start:idx_end,:),1,"movmean",6);
        HbT_act_30sec = smoothdata(HbT_lowpass_brain(idx_start:idx_end,:),1,"movmean",6);
    end
   
    HbT_act_30sec = HbT_act_30sec-mean(HbT_act_30sec,1);
    Ca_act_30sec = Ca_act_30sec-mean(Ca_act_30sec,1);
   
    R = f_xcorr(Ca_act_30sec,HbT_act_30sec,3*sr);
    
    id = ceil(size(R,1)/2);
    [max_R, max_idx] = max(abs(R(1:id,:)),[],1);
    ind = sub2ind(size(R),max_idx, 1:size(R,2));
    max_corr_oneside_brain(:,w) = R(ind);
    max_corr_lag_oneside_brain(:,w) = (max_idx-(3*sr+1));
end
max_corr_oneside = zeros(sX*sY,n_steps);
max_corr_lag_oneside = zeros(sX*sY,n_steps);
for u = 1:length(brain_idx)
    max_corr_oneside(brain_idx(u),:) = max_corr_oneside_brain(u,:);
    max_corr_lag_oneside(brain_idx(u),:) = max_corr_lag_oneside_brain(u,:);
end

max_corr_oneside= reshape(max_corr_oneside,[sX sY n_steps]);
max_corr_lag_oneside = reshape(max_corr_lag_oneside,[sX sY n_steps]);

%%
% Get clusters at each seed
cluster_idx = zeros(sX,sY,n_steps);
C = [-1 0.5*m_factor; -0.4 0.5*m_factor; -0.5 -0.4*m_factor; -2 -0.4*m_factor];
m_factor = 5;
for u = 1:sX
    u
    for v = 1:sY
        if brain_mask(u,v) == 1
            lag_arr_time = squeeze(max_corr_lag_oneside(u,v,:))/(fs);
            data = [lag_arr_time squeeze(max_corr_oneside(u,v,:))*m_factor];
            cluster_idx(u,v,:) = kmeans(data,4, 'MaxIter',10000,'Start',C);
        end
    end
end

% Get centroids of each cluster at each seed
centroids = zeros(sX,sY,4,2);
for u = 1:sX
    u
    for v = 1:sY
        if brain_mask(u,v) == 1
            for cidx = 1:4
                idx1 = find(squeeze(cluster_idx(u,v,:) == cidx));
                centroids(u,v,cidx,:) = [mean(max_corr_lag_oneside(u,v,idx1)) mean(max_corr_oneside(u,v,idx1))];
            end
        end
    end
end

colors =[1 0 0;
     0 1 0;
     0 0 1;
     0 1 1];
 pcolors ={'r.','b.','g.','c.'};
centroids = reshape(centroids,[sX*sY,4,2]);
figure
hold on
for cidx = 1:size(colors,1)
    cidx
%     scatter(centroids(brain_idx,cidx,1),centroids(brain_idx,cidx,2),5,'MarkerFaceColor',colors(cidx,:));
    plot(centroids(brain_idx,cidx,1)/fs,centroids(brain_idx,cidx,2),pcolors{cidx});
end
xlabel('lag'); ylabel('correlation coeff');xlim([-3 3]); ylim([-1 1])
title('center of mass for each cluster at each pixel')
hold off

% Get clusters for cluster centroids
m_factor = 20;
C = [-1 0.5*m_factor; -2.7 0.3*m_factor; -0.45 -0.4*m_factor; -4 -0.5*m_factor];
% C = [-0.88 0.42*m_factor; -4 0.1*m_factor; -0.5 -0.4*m_factor; -2.5 -0.3*m_factor];
com = centroids(brain_idx,:,:);
com = reshape(com,[length(brain_idx)*4,2]);
figure; plot(com(:,1)/fs,com(:,2),'k.'); xlabel('lag'); ylabel('correlation coeff');xlim([-3 3]); ylim([-1 1])
data_com = [com(:,1) com(:,2)*m_factor];
cluster_com = kmeans(data_com,4, 'MaxIter',100000,'Start',C);
transparent_value = [0.01 1];
for u = 1:2
    figure;
    scatter(data_com(cluster_com==1,1)/fs,data_com(cluster_com==1,2)/m_factor,2,'ro', 'MarkerEdgeAlpha', transparent_value(u))
    hold on
    scatter(data_com(cluster_com==2,1)/fs,data_com(cluster_com==2,2)/m_factor,2,'bo','MarkerEdgeAlpha', transparent_value(u))
    scatter(data_com(cluster_com==3,1)/fs,data_com(cluster_com==3,2)/m_factor,2,'go', 'MarkerEdgeAlpha', transparent_value(u))
    scatter(data_com(cluster_com==4,1)/fs,data_com(cluster_com==4,2)/m_factor,2,'co', 'MarkerEdgeAlpha', transparent_value(u))
    hold off
    xlabel('lag'); ylabel('correlation coeff');xlim([-3 3]); ylim([-1 1])
    title('center of mass for each cluster at each pixel')
    saveas(gcf,[folder_to_save filesep 'cluster_com_plot' '_' num2str(transparent_value(u)) '.png']);
end


cluster_com = reshape(cluster_com,[length(brain_idx),4]);
cluster_idx_com = zeros(sX*sY,4);
cluster_idx_com(brain_idx,:) = cluster_com;
cluster_idx_com = reshape(cluster_idx_com,[sX, sY, 4]);

% update clusters
cluster_idx_updated = cluster_idx;
for u = 1:sX
    u
    for v = 1:sY
        if brain_mask(u,v) == 1
            for cidx = 1:4
                idx1 = find(squeeze(cluster_idx(u,v,:) == cidx));      
                cluster_idx_updated(u,v,idx1) = cluster_idx_com(u,v,cidx);
            end
        end
    end
end

save([folder_to_save filesep 'cross_corr_results_' analyse '.mat'],'max_corr_oneside','max_corr_lag_oneside','cluster_idx_updated','cluster_idx');

%%

NE = gfp_normHD(:,:,fs*corr_time_length/2:fs*2:end);
NE_size = size(NE,3);
pupil_diameter = Signals.pupil(fs*corr_time_length/2:fs*2:end);
pupil_diam_size = length(pupil_diameter);

pts = [250 125;
       170 125;
       150 105;
       50 125];
%   
% colors =[1 0 0;
%      0 0 1;
%      0 1 0;
%      0 1 1];
% for pt = 1:size(pts,1) 
%     figure;
%     hold on
%     for u = 1:size(colors,1)
%         idx1 = find(squeeze(cluster_idx_updated(pts(pt,1),pts(pt,2),1:pupil_diam_size))==u);
%         scatter(squeeze(NE(seed_pts(pt,1),seed_pts(pt,2),idx1)),pupil_diameter(idx1),40,'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',colors(u,:),...
%               'LineWidth',1.5);
%     end
%     xlabel('Norepinephrine'); ylabel('Pupil dimater'); 
%     ylim([0 1]); 
%     xlim([-0.08 0.08])
% % xlim([0 1])
%     saveas(gcf,[folder_to_save filesep 'updated_classification_scatter_plot' '_' num2str(pt) '.png']);
%     hold off;
% end

% %%
% for pt = 1:size(pts,1)
%     figure;
%     hold on
%     for u = 1:size(colors,1)
%         idx1 = find(squeeze(cluster_idx_updated(pts(pt,1),pts(pt,2),1:pupil_diam_size))==u);
%         scatter(squeeze(max_corr_lag_oneside(seed_pts(pt,1),seed_pts(pt,2),idx1))/fs,max_corr_oneside(seed_pts(pt,1),seed_pts(pt,2),idx1),40,'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',colors(u,:),...
%               'LineWidth',1.5);
%     end
%     xlim([-3 3]); ylim([-1 1]); xticks([-3:1:3]);
%     xlabel('lag (seconds)'); ylabel('correlation coefficient'); 
%     hold off;
% end

% plot to select mouse alertness levels
for pt = 4
    figure;
    hold on
    for u = 1:length(pupil_diameter)
        h(u) = scatter(squeeze(NE(pts(pt,1),pts(pt,2),u)),pupil_diameter(u),40,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',colors(cluster_idx_updated(pts(pt,1),pts(pt,2),u),:),...
              'LineWidth',1.5);
    end
    xlabel('Norepinephrine'); ylabel('Pupil dimater'); 
    ylim([0 1]); 
    xlim([-0.08 0.08])
%     xlim([0 1])
    hold off;
end

select_point = true;
line_pts = [];
while select_point
    [xpt, ypt]  = ginput(2);
    hold on;
    plot(xpt, ypt, 'k--','LineWidth',2);
    line_pts = [line_pts; xpt ypt];
    answer = questdlg('Do you want to select another line'); 
    if ~strcmp(answer,'Yes')
        select_point = false;
    end
end

n_splits = size(line_pts,1)/2+1; 
folder_to_save_imges = [folder_to_save filesep 'state_images_' num2str(n_splits)];
if ~isfolder(folder_to_save_imges)
    mkdir(folder_to_save_imges); 
end
saveas(gcf,[folder_to_save_imges filesep 'selected_lines_plot.png']);


line_vectors = line_pts(1:2:end,:)-line_pts(2:2:end,:);
img_state1 = zeros(sX,sY,size(line_vectors,1)+2);
img_state2 = zeros(sX,sY,size(line_vectors,1)+2);
img_state3 = zeros(sX,sY,size(line_vectors,1)+2);
img_state4 = zeros(sX,sY,size(line_vectors,1)+2);
img_state_all = zeros(sX,sY,size(line_vectors,1)+2);
for u = 1:sX
    u
    for v = 1:sY
        if brain_mask(u,v) == 1
            all_pts = [squeeze(NE(u,v,:)) squeeze(pupil_diameter')];
            for w = 1:size(line_vectors,1)
                all_pts_vectors = all_pts-line_pts(2*w,:);
                vec_cross_product = all_pts_vectors(:,1).*line_vectors(w,2)-all_pts_vectors(:,2).*line_vectors(w,1);
                idx1 = find(vec_cross_product >= 0);
                if w == 1
                    idx = idx1;
                else
                    idx = intersect(idx1,idx2);
                end
                img_state1(u,v,w) = length(find(cluster_idx_updated(u,v,idx) == 1))/length(idx);
                img_state2(u,v,w) = length(find(cluster_idx_updated(u,v,idx) == 2))/length(idx);
                img_state3(u,v,w) = length(find(cluster_idx_updated(u,v,idx) == 3))/length(idx);
                img_state4(u,v,w) = length(find(cluster_idx_updated(u,v,idx) == 4))/length(idx);
                img_state_all(u,v,w) = length(idx)/size(NE,3);
                idx2 = find(vec_cross_product < 0);
                if w == size(line_vectors,1)
                    idx = idx2;
                    img_state1(u,v,w+1) = length(find(cluster_idx_updated(u,v,idx) == 1))/length(idx);
                    img_state2(u,v,w+1) = length(find(cluster_idx_updated(u,v,idx) == 2))/length(idx);
                    img_state3(u,v,w+1) = length(find(cluster_idx_updated(u,v,idx) == 3))/length(idx);
                    img_state4(u,v,w+1) = length(find(cluster_idx_updated(u,v,idx) == 4))/length(idx);
                    img_state_all(u,v,w+1) = length(idx)/size(NE,3);
                    
                    img_state1(u,v,w+2) = length(find(cluster_idx_updated(u,v,:) == 1))/size(NE,3);
                    img_state2(u,v,w+2) = length(find(cluster_idx_updated(u,v,:) == 2))/size(NE,3);
                    img_state3(u,v,w+2) = length(find(cluster_idx_updated(u,v,:) == 3))/size(NE,3);
                    img_state4(u,v,w+2) = length(find(cluster_idx_updated(u,v,:) == 4))/size(NE,3);
                    img_state_all(u,v,w+2) = 1;
                end
               
            end
            
        end
    end
end


for u = 1:size(img_state1,3)
    figure; colormap('jet'); imagesc(img_state1(:,:,u),[0 1]); axis image; colorbar;
    saveas(gcf,[folder_to_save_imges filesep 'state1_' num2str(u) '.png']);
    figure; colormap('jet'); imagesc(img_state2(:,:,u),[0 1]); axis image; colorbar;
    saveas(gcf,[folder_to_save_imges filesep 'state2_' num2str(u) '.png']);
    figure; colormap('jet'); imagesc(img_state3(:,:,u),[0 1]); axis image; colorbar;
    saveas(gcf,[folder_to_save_imges filesep 'state3_' num2str(u) '.png']);
    figure; colormap('jet'); imagesc(img_state4(:,:,u),[0 1]); axis image; colorbar;
    saveas(gcf,[folder_to_save_imges filesep 'state4_' num2str(u) '.png']);
    figure; colormap('jet'); imagesc(img_state_all(:,:,u),[0 1]); axis image; colorbar;
    saveas(gcf,[folder_to_save_imges filesep 'state_all_' num2str(u) '.png']);
%     close all
end
%%
% Thy1_215_231023_run5
% seed_pts = [169 130;
%             269 154;
%             131 19;
%             294 293];

% Thy1_215_231023_run1
% seed_pts = [196 118;
%     54 283;
%     207 415;
%     26 130];

% Thy1_215_231023_run5
% seed_pts = [196 118;
%     57 320;
%     203 421;
%     38 129];

% Thy1_215_231025_run1
seed_pts = [195 310;
    48 319;
    203 430;
    59 149];

pupil_diam_size = length(pupil_diameter);
 transperancy = [0.8 0.6 0.4 0.2];
for pt = 1:size(seed_pts,1)
    figure;
    hold on
    for u = 1:size(colors,1)
        idx1 = find(squeeze(cluster_idx_updated(seed_pts(pt,1),seed_pts(pt,2),1:pupil_diam_size))==u);
        scatter(squeeze(NE(seed_pts(pt,1),seed_pts(pt,2),idx1)),pupil_diameter(idx1),40,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',colors(u,:),...
              'LineWidth',1.5);
    end
    for u=1:size(line_pts,1)/2
        ii = 2*(u-1)+1;
        plot([line_pts(ii,1) line_pts(ii+1,1)] , [line_pts(ii,2) line_pts(ii+1,2)] , 'k--','LineWidth',2);
    end
    
    xlabel('Norepinephrine'); ylabel('Pupil dimater'); 
    ylim([0 1]); 
    xlim([-0.08 0.08])
%      xlim([0 1])
    saveas(gcf,[folder_to_save_imges filesep 'updated_classification_scatter_plot' '_' num2str(pt) '.png']);
    hold off;
end

for pt = 1:size(seed_pts,1)
    figure;
    hold on
    for u = 1:size(colors,1)
        idx1 = find(squeeze(cluster_idx_updated(seed_pts(pt,1),seed_pts(pt,2),1:pupil_diam_size))==u);
        scatter(squeeze(max_corr_lag_oneside(seed_pts(pt,1),seed_pts(pt,2),idx1))/fs,max_corr_oneside(seed_pts(pt,1),seed_pts(pt,2),idx1),40,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',colors(u,:),...
              'LineWidth',1.5);
    end

    xlim([-3 3]); ylim([-1 1]); xticks([-3:1:3]);
    xlabel('lag (seconds)'); ylabel('correlation coefficient'); 
    saveas(gcf,[folder_to_save_imges filesep 'updated_cluster_scatter_plot' '_' num2str(pt) '.png']);
    hold off;
end

