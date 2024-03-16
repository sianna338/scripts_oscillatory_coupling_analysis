%% fieldtrip
addpath('C:\Users\siann\Downloads\fieldtrip-20231113\fieldtrip-20231113')
ft_defaults
%% import so and spindle onset times
spindles = readtable("df_spindles.csv"); 
Slow_osc = readtable("df_sos.csv"); 

%% define trial-matrix based on SO onsets
sf = 5000; 
% pre = 0.4; % time before neg SO peak
% post = 0.8; % time after neg SO peak
trialmatrix = zeros(height(Slow_osc), 3); 
trialmatrix(:, 1) = Slow_osc.NegPeak*sf-3*sf; 
trialmatrix(:,2) = Slow_osc.NegPeak*sf+3*sf; 
trialmatrix(:,3) = -3*sf; 
trialmatrix = array2table(trialmatrix)
%% read-in data
cfg = [];
cfg.datafile = 'C:\Users\siann\Documents\SleepMem5a_S01_SpindleEEGTMS_20230705215233.eeg'
cfg.headerfile = 'C:\Users\siann\Documents\SleepMem5a_S01_SpindleEEGTMS_20230705215233.vhdr'
cfg.channel = {'C3'}; % indicate the channels we would like to read and/or exclude.
cfg.trl = trialmatrix;
cfg.continous = 'yes'
data_raw = ft_preprocessing(cfg);

%% ERPs
cfg = []; 
cfg.channels = {'C3'} 
cfg.keeptrials = 'no'; 
data_timelock = ft_timelockanalysis(cfg, data_raw)
% find peaks of slow oscillation
[Max, Idx_max]  = max(data_timelock.avg)
[Min, Idx_min]  = min(data_timelock.avg)
%% freq analysis
cfg = [];
cfg.method    = 'mtmconvol';
cfg.channels = {'C3'} 
cfg.foi       = 8:1:20;
% cfg.t_ftimwin = 7./cfg.foi; % time window depends on frequency, here: 7 cycles per second
% cfg.toi       = 0:0.01:max(Slow_osc.Duration)
% cfg.toi       = 0:0.1:mean(cellfun(@numel, data_raw.time)/data_raw.fsample)/2 % half of the trial length
% cfg.toi       = 0:0.05:(mean(cellfun(@numel,data_raw.time))-1)/data_raw.fsample % time vector
%cfg.toi       = 4.5:0.05:7.5 % time vector
cfg.toi = -1.5:0.05:1.5
cfg.t_ftimwin = ones(length(cfg.foi),1).*0.5; % fixed window length, 500 ms
% cfg.t_ftimwin = ones(length(cfg.foi),1).*max(Slow_osc.Duration)/1000 
cfg.taper     = 'hanning';
% cfg.pad = ceil(max(cellfun(@numel, data_raw.time)/data_raw.fsample))
freq_sos  = ft_freqanalysis(cfg, data_raw);

%% plot
cfg = [];
cfg.baselinetype = 'relative'; 
% cfg.baselinetype = 'db'
% cfg.maskstyle    = 'saturation';
figure;
% yyaxis left; 
ft_singleplotTFR(cfg, freq_sos);
xlim([(Idx_min-2000)/data_raw.fsample-trialmatrix(1,3) (Idx_min+5000)/data_raw.fsample-trialmatrix(1,3)])
yyaxis right; 
% plot(data_timelock.time(1, find(data_timelock.time==round(4.5)):find(data_timelock.time==round(7.5))), data_timelock.avg(1, find(data_timelock.time==round(4.5)):find(data_timelock.time==round(7.5))), 'k');
plot(data_timelock.time(1, Idx_min-2000:Idx_min+5000), data_timelock.avg(1, Idx_min-2000:Idx_min+5000), 'k');