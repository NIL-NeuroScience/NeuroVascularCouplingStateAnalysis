%% load data
% data root folder
root = '/projectnb/devorlab/bcraus/HRF/1P/23-10-20/Thy1_215/';
% run number
Run = 4;

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
folder_to_save = '/projectnb/devorlab/skura/HRF/cross_corr_results/run_4_231020_Thy1_215_newcorr';
if ~isfolder(folder_to_save)
    mkdir(folder_to_save); 
end


%% Run this block if you don't have brain mask in the folder where results are saved, otherwise skip this block
% Select brain mask
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

%% Get lag and cross corrlearion between Ca and HbT and cluster lag and correation to define states

% initial points for cluster
m_factor = 5;
C = [-1 0.5*m_factor; -3 0.5*m_factor; -0.5 -0.4*m_factor; -2 -0.4*m_factor];

% Spatial smoothing of the data
smooth_scale = 3;
HbT_smoothed = smooth2d(HbT,smooth_scale);
jRGECO_smoothed = smooth2d(jRGECO ,smooth_scale);
% GRAB_smoothed = smooth2d(GRAB ,smooth_scale);

% load brain_mask
if exist([folder_to_save filesep 'brain_mask.mat'],'file')
    load([folder_to_save filesep 'brain_mask.mat']);
end

[sX, sY, sZ] = size(HbT);
% apply filtering
brain_idx = find(brain_mask == 1);
HbT_smoothed = reshape(HbT_smoothed,[sX*sY sZ]);
HbT_smoothed_brain = HbT_smoothed(brain_idx,:);

fc= 0.01;
fs = dataIn(Run).settings.fs;
[b,a] = butter(6,fc/(fs/2));
HbT_lowpass_brain =  filtfilt(b,a,HbT_smoothed_brain');
HbT_lowpass = zeros(size(HbT_smoothed ));
HbT_lowpass(brain_idx,:) = HbT_lowpass_brain';
HbT_lowpass = reshape(HbT_lowpass,[sX sY sZ]);

HbT_smoothed = reshape(HbT_smoothed,[sX sY sZ]);
HbT_higpass = HbT_smoothed-HbT_lowpass;


jRGECO_smoothed = reshape(jRGECO_smoothed,[sX*sY sZ]);
jRGECO_smoothed_brain = jRGECO_smoothed(brain_idx,:);

[b,a] = butter(6,fc/(fs/2));
jRGECO_lowpass_brain =  filtfilt(b,a,jRGECO_smoothed_brain');
jRGECO_lowpass = zeros(size(jRGECO_smoothed ));
jRGECO_lowpass(brain_idx,:) = jRGECO_lowpass_brain';
jRGECO_lowpass = reshape(jRGECO_lowpass,[sX sY sZ]);

jRGECO_smoothed = reshape(jRGECO_smoothed,[sX sY sZ]);
jRGECO_higpass = jRGECO_smoothed-jRGECO_lowpass;

% %%
GRAB_smoothed = smooth2d(GRAB ,smooth_scale);
GRAB_smoothed = reshape(GRAB_smoothed,[sX*sY sZ]);
GRAB_smoothed_brain = GRAB_smoothed(brain_idx,:);

[b,a] = butter(6,fc/(fs/2));
GRAB_lowpass_brain =  filtfilt(b,a,GRAB_smoothed_brain');
GRAB_lowpass = zeros(size(GRAB_smoothed ));
GRAB_lowpass(brain_idx,:) = GRAB_lowpass_brain';
GRAB_lowpass = reshape(GRAB_lowpass,[sX sY sZ]);

GRAB_smoothed = reshape(GRAB_smoothed,[sX sY sZ]);
GRAB_higpass = GRAB_smoothed-GRAB_lowpass;

analyse = 'highpass';
% analyse = 'lowpass';
% analyse = 'allpass';
sr = dataIn(Run).settings.fs;

% sr = dataIn.settings.fs;
time_shift = 2;
corr_time_length = 30;
idx_shift = round(time_shift*sr);
n_steps = floor(size(jRGECO,3)/idx_shift);

max_corr = zeros(sX,sY,n_steps);
max_corr_lag= zeros(sX,sY,n_steps);

max_corr_oneside = zeros(sX,sY,n_steps);
max_corr_lag_oneside = zeros(sX,sY,n_steps);

cluster_idx = zeros(sX,sY,n_steps);

for u = 1:sX
    u
    for v = 1:sY
        if brain_mask(u,v) == 1
            for w = 1:n_steps
                idx_start = floor(1+idx_shift*(w-1));
                idx_end = min(floor(idx_start+sr*30),sZ);
                if strcmp(analyse, 'allpass')
                    Ca_act_30sec = smooth(squeeze(jRGECO_smoothed(u,v,idx_start:idx_end)));
                    HbT_act_30sec = smooth(squeeze(HbT_smoothed(u,v,idx_start:idx_end)));
                elseif strcmp(analyse, 'highpass')
                    Ca_act_30sec = smooth(squeeze(jRGECO_higpass(u,v,idx_start:idx_end)));
                    HbT_act_30sec = smooth(squeeze(HbT_higpass(u,v,idx_start:idx_end)));
                else
                    Ca_act_30sec = smooth(squeeze(jRGECO_lowpass(u,v,idx_start:idx_end)));
                    HbT_act_30sec = smooth(squeeze(HbT_lowpass(u,v,idx_start:idx_end)));
                end
                                
                Ca_act_30sec = Ca_act_30sec-mean(Ca_act_30sec);
                HbT_act_30sec = HbT_act_30sec-mean(HbT_act_30sec);
                [r,lags] = xcorr(Ca_act_30sec,HbT_act_30sec,3*sr,'normalized');
                
                [max_r, lag_idx] = max(abs(r));
                max_corr(u,v,w) = r(lag_idx);
                max_corr_lag(u,v,w) = lags(lag_idx);
                
                % one sided cross corr
                id = ceil(length(lags)/2);
                [max_r, lag_idx] = max(abs(r(1:id)));
                max_corr_oneside(u,v,w) = r(lag_idx);
                max_corr_lag_oneside(u,v,w) = lags(lag_idx);
                
            end
            lag_arr_time = squeeze(max_corr_lag_oneside(u,v,:))/(fs);
            data = [lag_arr_time squeeze(max_corr_oneside(u,v,:))*m_factor];
            cluster_idx(u,v,:) = kmeans(data,4, 'MaxIter',10000,'Start',C);
        end
    end
end
save([folder_to_save filesep 'cross_corr_results_' analyse '.mat'],'max_corr','max_corr_lag','max_corr_oneside','max_corr_lag_oneside','cluster_idx');