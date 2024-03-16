% run phase-amplitude coupling anaylsis (SO-spindle), extract phase
% timeseries from the slower oscillation (SO) and amplitude timeseries from
% the faster oscillation (spindle), then identify SO-phase where power in
% the spindle band is maximal
%
% Fieldtrip and circular statistics toolbox need to be installed
%
% Input:
% rootidr   data directory containing the subfolder 'data' with seperate
% folders for each subject 
function [mean_direction, vector_length, pval, z] = SO_spindle_PAC_analysis(rootdir)

%% settings
datadir   = fullfile(rootdir,'\data\');
addpath(datadir)
cd (datadir)

subjects = {'P01' 'P02' 'P03' 'P04' 'P05' 'P06' 'P07' 'P08'};


for isubject=1:numel(subjects)

    %% load-in the data 
    cfg = [];
    cfg.demean  = 'yes';
    cfg.detrend = 'yes';
    cfg.dataset = fullfile(datadir, subjects{isubject}, 'SN00.edf')
    % cfg.dataset = fullfile(datadir,subjects{isubject}, 'data')
    data_sleep   = ft_preprocessing(cfg);
    assignin('base', 'data_sleep', data_sleep)
    % filter the data in the SO frequency range and extract the phase timeseries  
    cfg                     = [];
    cfg.bpfilter            = 'yes';              
    cfg.bpfreq              = [0.5 2];
    cfg.bpinstabilityfix    = 'reduce';
    cfg.bpfiltdir           = 'twopass';
    data_phase          = ft_preprocessing(cfg, data_sleep); % filtering step

    cfg             = []; 
    cfg.hilbert     = 'angle'; 
    data_phase  = ft_preprocessing(cfg, data_phase); % extract phase timeseries for the signal filtered in the SO frequency band 

   % bring data into better shape
    cfg = [];
    cfg.keeptrials = 'yes';
    cfg.latency    = [-1.5 1.5];
    data_phase = ft_timelockanalysis (cfg, data_phase);
    
    % filter the data in the spindle frequency range and extract the amplitude time series      
    cfg                     = [];
    cfg.bpfilter            = 'yes';              
    cfg.bpfreq              = [12 15];
    cfg.bpinstabilityfix    = 'reduce';
    cfg.bpfiltdir           = 'twopass'; 
    data_amp            = ft_preprocessing(cfg, data_sleep); % filtering step

    cfg             = []; 
    cfg.hilbert     = 'abs'; 
    data_amp    = ft_preprocessing(cfg, data_amp); % extract the amplitude timeseries for the signal filtered in the spindle frequency band 

    % bring data into better shape
    cfg            = [];
    cfg.keeptrials = 'yes';
    cfg.latency    = [-1.5 1.5];
    data_amp   = ft_timelockanalysis (cfg, data_amp);    

    % extract maximum amplitude values for the spindle frequency band 
    val = [];trialnum = [];
    ind = [];ichan    = [];
    for trialnum      = 1:size(data_amp.trial,1); 
        for ichan     = 1:size(data_amp.trial,2);
            [val(trialnum ,ichan) ind(trialnum ,ichan)] = max (data_amp.trial(trialnum,ichan,:));
        end
    end

    % extract corresponding phase of the SO (at the times of maximum
    % spindle activity)
    val_phase         = [];
    for trialnum      = 1:size(data_amp.trial,1); 
        for ichan     = 1:size(data_amp.trial,2);
            [val_phase(trialnum ,ichan)] = data_phase.trial(trialnum,ichan, ind(trialnum ,ichan));
        end
    end
    val_phase_all{(isubject)} = val_phase;
end
    assignin('base', 'val_phase', val_phase)
    assignin('base', 'val', val)

    assignin('base', 'val_phase_all', val_phase_all)
% circ_mean gives you average direction of vector, so average phase angle
% where spindle activity is maximal

for iSub  = 1:numel(subjects) 
        %[mean_direction(iSub), ul(iSub), ll(iSub)] = circ_mean(val_phase_all{iSub}(:,2)) ;
        % mean_direction(iSub) = squeeze(circ_mean(val_phase_all{iSub}(:,2))) ;
        mean_direction(iSub) = circ_mean(val_phase_all{iSub}(:,2)) % select channel C4
        mean_direction_squeezed = mean(mean_direction)
        
    end    
        assignin('base', 'mean_direction', mean_direction)
        assignin('base', 'mean_direction_squeezed', mean_direction_squeezed)
    % get strength of PAC
    %  Computes mean resultant vector length for circular data.

    % signficance test for preferred phase 
    [pval, z] = circ_rtest(mean_direction); % Computes Rayleigh test for non-uniformity of circular data.
    vector_length = circ_r(mean_direction);
    assignin('base', 'vector_length', vector_length)
    % plot circular histogram 
    figure;
    for iSub  = 1:numel(subjects) 
        h = circ_plot(mean_direction(iSub), 'pretty', 'ro', true, 'color','r')% handle to the individual dots
        hold on;
    end 
    h_2 = circ_plot(mean(mean_direction, 2),'pretty','go',true,'linewidth',1,'color','g')

    % h = circ_plot(z, 'pretty', 'ro', true)

    % save all the output variables    

    %% step 2: correlate spindle-phase specific MEP amplitude modulation with spindle-SO phase coupling 

    % [RHO,PVAL] = circ_corrcl(alpha_all, brain);
