function [plv] = compute_plv(input_dat, sf, filter_01, filter_02)
% This computes the phase locking value for two different EEG signals 
%
% Input parameters:
%   input_dat is the raw data, in num trials x num channels x time point
%   format 
%   sf is the sampling rate
%   filter_01.range is the lower and higher cutoff for the FIR filter
%   filter_01.order is the order of the FIR filter 
%   the same needs to be specified for the second frequency band of
%   interest
% Output parameters:
%   plv is

num_trials = size(input_dat, 1)
num_channels = size(input_dat,2)

% discard the first and last couple of samples to avoid artifcats
input_dat = input_dat(:, :, filter_01.order:end-filter_01.order)
% input_dat = input_dat(:, :, 10000:20000)

num_timepoints = size(input_dat,3)

% create filters and filter data in the two frequency bands of interest
first_filter = fir1(filter_01.order, filter_01.range/(sf/2))
second_filter = fir1(filter_02.order, filter_02.range/(sf/2))
filter1_dat = filter(first_filter, 1, input_dat, [], 3)
filter2_dat = filter(second_filter, 1, input_dat, [], 3)
 
% hilbert transform and extract phase time series
first_phase_timeseries = zeros(num_trials, num_channels, num_timepoints);
second_phase_timeseries = zeros(num_trials, num_channels, num_timepoints);
for itrial=1:num_trials
    for ichan = 1:num_channels
        first_phase_timeseries(itrial,ichan,:) = angle(hilbert(squeeze(filter1_dat(itrial,ichan,:))))
        second_phase_timeseries(itrial,ichan,:) = angle(hilbert(squeeze(filter2_dat(itrial,ichan,:)))) 
    end 
end 

% compute PLV
plv = zeros(num_trials, num_channels);
for itrial=1:num_trials
    for ichan=1:num_channels
        plv(itrial, ichan) = 1/num_timepoints*abs(sum(exp(1i*(first_phase_timeseries(itrial,ichan,:)-second_phase_timeseries(itrial,ichan,:)))))
    end 
end 

plv = squeeze(plv);
return;