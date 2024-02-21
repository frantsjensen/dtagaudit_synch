function [cues] = detectEvents(tag,tcue,NS)

% Humpback whale song note detector
CH=1;           % Analysis channel
AFS_RES = 8000; % Desired sample rate after downsampling
HP_filt = 40;   % High pass filter, Hz
mindur  = 0.10  % Minimum duration (in s)
gap     = 0.25; % Minimum gap (in s) between subsequent detections


% Main function

% Default if detection fails
cues=[];

% Load audio data
[x,afs_org] = dtagwavread(tag,tcue,NS);  
if isempty(x), return, end
   
% Resample to conserve memory
x=resample(x(:,CH)-mean(x(:,CH)),AFS_RES,afs_org);

% Filter to get rid of zero offset, then calculate envelope
[B,A]=butter(4,HP_filt/(AFS_RES/2),'high');
xf=filtfilt(B,A,x);
env = abs(hilbert(xf));

% Smooth envelope with 10 ms running average, then convert to dB
NN = AFS_RES/100;
env_smooth = 20*log10(conv(env,ones(NN,1)/NN,'same'));

% Calculate some statistics for finding a good threshold
env_dist = prctile(env_smooth,[2.5 10 50 90 97.5]);

% Take mean of 2.5th percentile and 97.5th percentile as starting threshold
thr1 = mean(env_dist([1,5]);


% find events exceeding threshold
env_smooth>thr1;
k=find(env>5);

if ~isempty(k)
    detect_samples = [k(1) k(end)]; % Default is one long detection
else
    return
end


% Segment into different events using critical gap length
dk = diff(k);
crit = AFS_RES*gap;
sep = find(dk>crit);

% Segmentation loop
for j=1:length(sep)
    % Move last sample to end of last detection
    detect_samples(j+1,2) = detect_samples(j,2);

    % Update end of first block
    detect_samples(j,2) = k(sep(j));

    % Update start of next block
    detect_samples(j+1,1) = k(sep(j)+1);
end

% Convert to start time vs duration
cues(:,1) = tcue+detect_samples(:,1)/AFS_RES;
cues(:,2) = (detect_samples(:,2)-detect_samples(:,1))/AFS_RES;

% Remove detections that are too short
cues(find(cues(:,2)<mindur),:) = [];

figure(90), clf, hold on
plot([1:length(env_smooth)]/AFS_RES,env_smooth,'k')
for j=1:size(cues,1)
    plot(detect_samples(j,:)/AFS_RES),max(env_smooth(detect_samples(j,1):detect_samples(j,2)),'r','linewidth',3)
end

return



