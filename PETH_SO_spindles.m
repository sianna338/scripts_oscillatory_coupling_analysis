
function [rst,detectwndw]=PETH_SO_spindles(slow_oscillations,spindles,pre,pos); 

spindles = 6000*rand(300, 1)
slow_oscillations = 6000*rand(300, 1)
pre = 1.2;
pos = 1.2;
% slow_oscillations is vector of detected SO troughs
% spindles is vector of detected fast spindles centers
detectwndw=[slow_oscillations(:,1)-pre slow_oscillations(:,1)+pos]; % create for each detected SO center timewindow for PETH analysis 
rst=[]; 
% for i = 1:length(spindles)
    for j = 1:length(detectwndw)
            count = sum(spindles >= detectwndw(j,1) & spindles <= detectwndw(j,2)); % Count how many spindles fall within the time window 
% for i=1:length(detectwndw)
%     count=sum(spindles<=detectwndw(i,2) & spindles>=detectwndw(i,1)) %count how many times spindle beginning time is bigger or equal to lower limit and smaller or equal to upper limit
        if any (count~=0) % when spindles fall into the time window of interest 
        t= find(spindles >= detectwndw(j,1) & spindles <= detectwndw(j,2)) % identify indices for desired spindle events 
             for k=1:length(t)
                 rst(end+1,1)=spindles(t(k))-slow_oscillations(j,1); % time difference between spindle beginning and SO center
                 rst(end,2)=j; % which slow oscillation is spindle coupled to
                 rst(end,3)=slow_oscillations(j,1);  % time of spindle-coupled SO-trough 
            end
        else % when none of the spindles fall into the timewindow of SO-center - pre + post
             rst(end+1,1)=NaN;
             rst(end,2)=j;
             rst(end,3)=slow_oscillations(j,1);
        end
    end
end 