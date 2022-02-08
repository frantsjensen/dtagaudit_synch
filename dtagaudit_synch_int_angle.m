function [AXS] = dtagaudit_synch_int_angle(tag,tcue,NS,AOA_FH)
% QUICK SCRIPT TO PLOT ANGLE OF ARRIVAL
% COMPATIBLE W DTAG2 AND DTAG3

global CH AFS_RES SPEC_BL SPEC_OL SPEC_CLIM SPEC_FLIM AOA_SCF

tagver = dtagtype(tag);

% Manually choose cutoff filter
% make it possible to set lowpass criterion, high-pass criterion, and
% amplitude threshold
settings = {'High-pass filter (kHz)','Low-pass filter (kHz)','Amplitude threshold (dB)'}; 
set_values = {'5','20','-74'}; % Default
if nargin>3
    set_values(1) = num2str(AOA_FH/1000);
end
set_values = inputdlg(settings,'Change LP/HP/Threshold settings',[1,60],set_values);

if isempty(set_values{1})
    AOA_FH = 4000;
else
    AOA_FH = 1000*[str2num(set_values{1}) str2num(set_values{2})];
end

if isempty(set_values{3})
    THRESH = -74;
else
    THRESH = str2num(set_values{3});
end

% Default settings
AOA_BL = 0.020 ;    % Block length (ms) of running angle of arrival estimate
AOA_OL = 0.75 ; % Overlap in percent % Was 0.5

% Load data
[x,afs_org] = dtagwavread(tag,tcue,NS);
if isempty(x), return, end

if size(x,2)<2
    disp(' Angle of arrival not possible with single-channel tags')
    return
end

% Prepare figure first
figure(97),clf,
set(gcf,'name','Angle')

% Create figure axes
AXs = axes('position',[0.11,0.53,0.8,0.40]) ;
AXa = axes('position',[0.11,0.11,0.8,0.40]) ; box on, hold on

% Resample to conserve memory
for j=1:size(x,2)
    x(:,j)=x(:,j)-mean(x(:,j));
end
x=resample(x,AFS_RES,afs_org);
afs = AFS_RES;

% Create and plot spetrogram
[B,F,T] = spectrogram(x(:,CH),hamming(SPEC_BL),round(SPEC_BL*SPEC_OL),SPEC_BL,afs, 'yaxis') ;
BB = adjust2Axis(20*log10(abs(B))) ;
axes(AXs), imagesc(tcue+T,F/1000,BB,SPEC_CLIM) ; axis xy; grid ;

% Restrict y axis in spectrogram
set(AXs,'XTickLabel',[],'YLim',SPEC_FLIM,'XLim',tcue+[0 NS]) ;
ylabel('Frequency, kHz')
set(gca,'XTickLabel','')

% Filter before calculating angle of arrival
if length(AOA_FH)==1
    [baoa aaoa] = butter(4,AOA_FH/(afs/2),'high') ;
elseif length(AOA_FH)==2
    [baoa aaoa] = butter(4,AOA_FH/(afs/2)) ;
end
[xf]=filter(baoa,aaoa,x);

% Set up analysis
L = round(AOA_BL*afs);
N = 1;
i = 0;

% Work through data sequence
while N+L<=length(xf)
    k=N+[1:L] ;

    % Only analyse angle of arrival data if rms level above threshold
    if 20*log10(std(xf(k,1)))>THRESH
        [aa,qq] = xc_tdoa(xf(k,1),xf(k,2),afs/AOA_SCF*1.02) ;
        if qq>0.7
            i=i+1 ; 
            time(i) = (N+round(0.5*L))/afs ;
            RL(i)   = 20*log10(std(xf(k,1))) ;
            aoa(i)  = real(asin(aa*AOA_SCF/afs)*180/pi) ;
        end
    end
    N=N+round(L*(1-AOA_OL));
end

% Now plot aoa data
axes(AXa), scatter(tcue+time,aoa,[],RL); 
axis xy, grid ; box on
set(gca,'YLim',[-90 90])
 
% Adjust time limits
set([AXs AXa],'XLim',tcue+[0 NS]) 
AXS=[AXs AXa];