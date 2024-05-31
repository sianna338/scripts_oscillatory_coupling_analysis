%%
datapath = 'E:\spindle_ppTMS\EEG'
subject = 'sub-06'
session = 'ses-adapt'
sf = 5000
sw_times = readtable([datapath filesep subject filesep session filesep subject '_sw_times'])
So_troughs = sw_times.NegPeak

% fieldtrip
path_ft   = 'C:\Users\siann\Downloads\fieldtrip-20231113\fieldtrip-20231113';
addpath(path_ft);
ft_defaults;

%% 
filter_01.order = 5000
filter_01.range = [0.5 2]
filter_02.order = 1660
filter_02.range = [12 15]
%% segment data around the detected SO troughs

cfg = [];
cfg.trl = [So_troughs*sf-(2.5*sf) So_troughs*sf+(2.5*sf) zeros(length(So_troughs),1)]

cfg.datafile = [datapath, filesep, subject, filesep, session, filesep, 'spindle-ppTMS_', subject, '_', session, '.eeg']
cfg.headerfile = [datapath, filesep, subject, filesep, session, filesep, 'spindle-ppTMS_', subject, '_', session, '.vhdr']
cfg.continous = 'yes';

cfg = ft_definetrial(cfg);
trial_matrix = cfg.trl


% read-in data
cfg.channel = {'Cz'};

data_raw_cz = ft_preprocessing(cfg);

%% convert data to trial x channel x timepoints matrix
num_channels = size(data_raw_cz.trial{1},1)
num_trials = size(data_raw_cz.trial,2)
num_timepoints = size(data_raw_cz.trial{1},2)
input_dat = zeros(num_trials, num_channels, num_timepoints)
for itrial = 1:num_trials
    input_dat(itrial,:,:) = data_raw_cz.trial{itrial}(:,:)
end 
%% filter the data in the spindle and SO frequency range 
first_filter = fir1(filter_01.order, filter_01.range/(sf/2))
second_filter = fir1(filter_02.order, filter_02.range/(sf/2))
filter1_dat = filter(first_filter, 1, input_dat, [], 3)
filter2_dat = filter(second_filter, 1, input_dat, [], 3)
 
%% extract phase and amplitude time series

SO_phase_timeseries = zeros(num_trials, num_channels, num_timepoints);
sp_amp_timeseries = zeros(num_trials, num_channels, num_timepoints);
for itrial=1:num_trials
    for ichan = 1:num_channels
        SO_phase_timeseries(itrial,ichan,:) = angle(hilbert(squeeze(filter1_dat(itrial,ichan,:))))
        sp_amp_timeseries(itrial,ichan,:) = abs(hilbert(squeeze(filter2_dat(itrial,ichan,:)))) 
    end 
end 

%% make input window smaller (-2, +2 seconds) around spindle troughs to avoid filter artifacts
SO_phase_timeseries = SO_phase_timeseries(:,:,0.5*sf:end-0.5*sf)
sp_amp_timeseries = sp_amp_timeseries(:,:,0.5*sf:end-0.5*sf)

num_timepoints = size(sp_amp_timeseries,3)

%% find for each trial, SO phase for which there is spindle amplitude peak 
sp_amp_max = zeros(num_trials, num_channels, 1);
idx_max_sp = zeros(num_trials, num_channels, 1);
so_phase_sp_max = zeros(num_trials, num_channels, 1);

for itrial=1:num_trials
    for ichan = 1:num_channels
        [sp_amp_max(itrial, ichan), idx_max_sp(itrial, ichan)] = max(sp_amp_timeseries(itrial, ichan,:))
        so_phase_sp_max(itrial, ichan) = SO_phase_timeseries(itrial, ichan, idx_max_sp(itrial, ichan))
    end 
end


%% extract mean vector direction and vector length with circ_stat toolbox (SO phases with max sp amp as input)

[mu ul ll] = circ_mean(so_phase_sp_max)
circ_plot(so_phase_sp_max,'pretty','ro',true,'linewidth',2,'color','r')