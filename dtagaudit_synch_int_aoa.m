function [] = dtagaudit_synch_int_aoa(tag,tcue,RES,rtag,time_shift)
% R=dtagaudit_synch_int_aoa(tag,tcue,RES,rtag,NS,time_shift)
% where tag (tagged whale) and rtag (receiving tag) 
% are string names for the two tags to compare,
% tcue is a start time cue for the audit plot, and
% R is an audit structure from the tagged whale.
% 
% Optional input arguments include
% time_shift (extra time shift, default 0)

global AFS_RES NS

if nargin<5
    time_shift = 0 ; % default is no time shift
elseif isempty(NS)
    time_shift = 0 ;
end

% Check if multi-channel tags
[x,~] = dtagwavread(tag,tcue,NS);
if size(x,2)<2
    disp(' Angle of arrival not possible with single-channel tags')
    return
end

% Make sure that there is a time shift for each rtag
if length(time_shift)==1
    time_shift=ones(length(rtag))*time_shift;
end

% Manually choose cutoff filter
% make it possible to set lowpass criterion, high-pass criterion, and
% amplitude threshold
settings = {'High-pass filter (kHz)','Low-pass filter (kHz)','Amplitude threshold (dB)'}; 
set_values = {'5','20','-74'}; % Default
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

% Prepare figure first
figure(98),clf

% Make corresponding axes for angle-of-arrival plots
switch length(rtag)
    case 1
        AXa(1) = axes('position',[0.11,0.42,0.78,0.30]) ;
        AXa(2) = axes('position',[0.11,0.11,0.78,0.30]) ;
        
    case 2
        AXa(1) = axes('position',[0.11,0.52,0.78,0.20]) ;
        AXa(2) = axes('position',[0.11,0.31,0.78,0.20]) ;
        AXa(3) = axes('position',[0.11,0.10,0.78,0.20]) ;
       
    case 3
        AXa(1) = axes('position',[0.11,0.57,0.78,0.15]) ;
        AXa(2) = axes('position',[0.11,0.41,0.78,0.15]) ;
        AXa(3) = axes('position',[0.11,0.25,0.78,0.15]) ;
        AXa(4) = axes('position',[0.11,0.09,0.78,0.15]) ;
        
    case 4
        AXa(1) = axes('position',[0.11,0.60,0.78,0.12]) ;
        AXa(2) = axes('position',[0.11,0.47,0.78,0.12]) ;
        AXa(3) = axes('position',[0.11,0.34,0.78,0.12]) ;
        AXa(4) = axes('position',[0.11,0.21,0.78,0.12]) ;
        AXa(5) = axes('position',[0.11,0.08,0.78,0.12]) ;

    case 5
        AXa(1) = axes('position',[0.11,0.62,0.78,0.10]) ;
        AXa(2) = axes('position',[0.11,0.51,0.78,0.10]) ;
        AXa(3) = axes('position',[0.11,0.40,0.78,0.10]) ;
        AXa(4) = axes('position',[0.11,0.29,0.78,0.10]) ;
        AXa(5) = axes('position',[0.11,0.18,0.78,0.10]) ;
        AXa(6) = axes('position',[0.11,0.07,0.78,0.10]) ;
        
    otherwise
        disp('AOA script only supports up to 6 channels')
end

% Make axes for first tag cues+spectrogram
AXc = axes('position',[0.11,0.94,0.78,0.02]) ;
AXs = axes('position',[0.11,0.73,0.78,0.20]) ;
        
% Change settings for cue axes
bc = get(gcf,'Color') ;
set([AXc],'XLim',[0 1],'YLim',[0 1]) ;
set([AXc],'Box','off','XTick',[],'YTick',[],'XColor',bc,'YColor',bc,'Color',bc) ;
set([AXs AXa],'Box','on')

% Find time difference between tags
cuediff = dtagtimediff(tag,rtag,time_shift);

% Plot the info from tag 1
plotspec(tag,tcue,RES,[AXc AXs]);
[RES]=plotaoa(tag,tcue,RES,[AXa(1)],NS,AOA_FH,THRESH);

% Add info from other tags
for i=1:length(rtag)
    newtag = char(rtag(i)) ;
    try
        [RES2]=plotaoa(newtag,tcue+cuediff(i),[],[AXa(i+1)],NS,AOA_FH,THRESH);
    catch
        disp(['Error plotting angle-of-arrival for tag ' char(rtag(i)) ' at cue ' num2str(tcue+cuediff(i))])
    end
end



function     [RES] = plotaoa(tag,tcue,RES,AXS,NS,AOA_FH,THRESH)
% AXS a one-element vector of axis handles

global AOA_SCF AFS_RES

if nargin<7
    THRESH  = -74  ; % Amplitude threshold in dB
end

if nargin<6
    AOA_FH  = 6e3 ;     % high-pass filter for angle-of-arrival measurement
end

if nargin<5
    NS = 4 ;           % number of seconds to display
end

% Settings for AOA measurements
BL      = 0.040 ; % Block length for AOA measurements, s % Was 0.01
OL      = 0.75 ; % Overlap in percent % Was 0.5

% Get axes handles from input variables
AXs=AXS(1);

if isempty(RES)
    RES=loadaudit(tag);
end

if nargin<3 | isempty(RES),
   RES.cue = [] ;
   RES.comment = [] ;
end

% Read in data
[x,afs_org] = dtagwavread(tag,tcue,NS);
if isempty(x), return, end

% Resample
x   = resample(x,AFS_RES,afs_org);
afs = AFS_RES ;

% Filter data
if length(AOA_FH)==1
    [baoa aaoa] = butter(4,AOA_FH/(afs/2),'high') ;
elseif length(AOA_FH)==2
    [baoa aaoa] = butter(4,AOA_FH/(afs/2)) ;
end
xf      = filter(baoa,aaoa,x) ;

% Set up analysis
L = round(BL*afs);
N = 1;
i = 0;

% Work through data sequence
while N+L<=length(xf)
    k=N+[1:L] ;

    % Only analyse angle of arrival data if rms level above threshold
    if 20*log10(std(xf(k,1)))>THRESH
        i=i+1 ; 
        time(i) = (N+round(0.5*L))/afs ;
        [aa,qq] = xc_tdoa(xf(k,1),xf(k,2),afs/AOA_SCF*1.02) ;
        aoa(i)  = real(asin(aa*AOA_SCF/afs)*180/pi) ;
        RL(i)   = 20*log10(std(xf(k,1))) ;
    end
    N=N+round(L*(1-OL));
end

% Now plot aoa data
axes(AXs), scatter(time,aoa,[],RL); 
text(0.1,80,[tag(1:4) '\_' tag(6:9)],'VerticalAlignment','Top','HorizontalAlignment','Left');
axis xy, grid ; box on

if max(real(aoa))>90
    disp([' Error: Max angle is ' num2str(max(abs(aoa)))])
end

% Restrict y axis in spectrogram
set(AXs,'XTickLabel',[],'YLim',[-90 90],'XLim',[0 NS]) ;
ylabel('AOA, deg')

return


function [x,afs] = plotspec(tag,tcue,RES,AXS,col)

global CH NS AFS_RES SPEC_BL SPEC_OL SPEC_CLIM SPEC_FLIM

if nargin<5
    col = [1 1 1];
end

% Get axes handles from input variables
AXc=AXS(1); AXs=AXS(2);

% Read audio and resample to conserve memory
[x,afs_org] = dtagwavread(tag,tcue,NS) ;
if isempty(x), return, end
x=resample(x(:,CH)-mean(x(:,CH)),AFS_RES,afs_org);
afs = AFS_RES;

% Calculate spectrogram
[B,F,T] = spectrogram(x(:,CH),hamming(SPEC_BL),round(SPEC_BL*SPEC_OL),SPEC_BL,afs, 'yaxis') ;

% Format and plot spectrogram   
BB = adjust2Axis(20*log10(abs(B))) ;
axes(AXs), imagesc(tcue+T,F/1000,BB,SPEC_CLIM) ;
text(tcue+0.1,0.9*SPEC_FLIM(2),[tag(1:4) '\_' tag(6:9)],...
    'VerticalAlignment','Top','HorizontalAlignment','Left',...
    'Backgroundcolor',[0.2 0.2 0.2],'Color',col,'FontWeight','bold');
axis xy, grid ;

% Restrict y axis in spectrogram
set(AXs,'XTickLabel',[],'YLim',SPEC_FLIM) ;
ylabel('Frequency, kHz')

% Check data for labels
if isempty(RES), RES=loadaudit(tag); end

% Add labels
if ~(nargin<3 || isempty(RES))
    dtagaudit_plotRES(AXc,RES,[tcue tcue+NS]) ;
end

return