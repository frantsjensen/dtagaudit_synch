function [AXS] = dtagaudit_synch_int(tag,tcue,RES,rtag,time_shift)
% R=tagaudit_synch_int(tag,tcue,R,rtag,NS,time_shift)
% where tag (tagged whale) and rtag (receiving tag) 
% are string names for the two tags to compare,
% tcue is a start time cue for the audit plot, and
% R is an audit structure from the tagged whale,
% rtag is a 9-character tag deployment string OR a cell array of 
% tag deployment strings to compare with.
% 
% Optional input arguments include
% time_shift (extra time shift for each tag, default 0)

global NS ENV_LIM

% First, check input parameters
if nargin<5
    time_shift = 0 ; % default is no time shift
elseif isempty(time_shift)
    time_shift = 0 ;
end

% Make sure that there is a time shift for each rtag
if length(time_shift)~=length(rtag)
    time_shift=zeros(length(rtag));
end

% Make rtag string array into cell array
if ischar(rtag)
    rtag={rtag};
end
        
for i=1:length(rtag)
    if ~strcmp(tag(1:8),rtag{i}(1:8))
        disp('Error: Tag deployments must be from same day')
    end
end

% Prepare figure first
figure(99),clf,
set(gcf,'name','Simultaneous Spectrogram Plot')

% Create axes for tag comparisons
switch length(rtag)
    case 1
        % Make axes for first tag cues+spectrogram
        AXc = axes('position',[0.11,0.93,0.78,0.03]) ;
        AXs = axes('position',[0.11,0.64,0.78,0.28]) ;

        % Make corresponding axes for second tag cues+spectrogram
        AXc2(1) = axes('position',[0.11,0.60,0.78,0.03]) ;
        AXs2(1) = axes('position',[0.11,0.31,0.78,0.28]) ;

        % Make axis for amplitude envelopes for all whales
        AXa = axes('position',[0.11,0.10,0.78,0.18]) ;

    case 2
        % Make axes for first tag cues+spectrogram
        AXc = axes('position',[0.11,0.94,0.78,0.02]) ;
        AXs = axes('position',[0.11,0.71,0.78,0.20]) ;

        % Make corresponding axes for second tag cues+spectrogram
        AXc2(1) = axes('position',[0.11,0.68,0.78,0.02]) ;
        AXs2(1) = axes('position',[0.11,0.47,0.78,0.20]) ;
        AXc2(2) = axes('position',[0.11,0.44,0.78,0.02]) ;
        AXs2(2) = axes('position',[0.11,0.23,0.78,0.20]) ;
        
        % Make axis for amplitude envelopes for all whales
        AXa = axes('position',[0.11,0.09,0.78,0.12]) ;
        
    case 3
        % Make axes for first tag cues+spectrogram
        AXc = axes('position',[0.11,0.94,0.78,0.02]) ;
        AXs = axes('position',[0.11,0.78,0.78,0.15]) ;

        % Make corresponding axes for other tag cues+spectrogram
        AXc2(1) = axes('position',[0.11,0.75,0.78,0.02]) ;
        AXs2(1) = axes('position',[0.11,0.59,0.78,0.15]) ;
        AXc2(2) = axes('position',[0.11,0.56,0.78,0.02]) ;
        AXs2(2) = axes('position',[0.11,0.40,0.78,0.15]) ;
        AXc2(3) = axes('position',[0.11,0.37,0.78,0.02]) ;
        AXs2(3) = axes('position',[0.11,0.21,0.78,0.15]) ;

        % Make axis for amplitude envelopes for all whales
        AXa = axes('position',[0.11,0.09,0.78,0.10]) ;
        
    case 4
        % Make axes for first tag cues+spectrogram
        AXc = axes('position',[0.11,0.935,0.78,0.015]) ;
        AXs = axes('position',[0.11,0.810,0.78,0.120]) ;

        % Make corresponding axes for other tag cues+spectrogram
        AXc2(1) = axes('position',[0.11,0.785,0.78,0.015]) ;
        AXs2(1) = axes('position',[0.11,0.660,0.78,0.120]) ;
        AXc2(2) = axes('position',[0.11,0.635,0.78,0.015]) ;
        AXs2(2) = axes('position',[0.11,0.510,0.78,0.120]) ;
        AXc2(3) = axes('position',[0.11,0.485,0.78,0.015]) ;
        AXs2(3) = axes('position',[0.11,0.360,0.78,0.120]) ;
        AXc2(4) = axes('position',[0.11,0.335,0.78,0.015]) ;
        AXs2(4) = axes('position',[0.11,0.210,0.78,0.120]) ;

        % Make axis for amplitude envelopes for all whales
        AXa = axes('position',[0.11,0.09,0.78,0.10]) ;

    otherwise
        disp('Simultaneous tag auditing up to 5 tags only')
end   

% Balanced color scheme, categorical
CMAT = [...
    242 089 094 ;... % Red
    089 154 210 ;... % Blue
    249 166 090 ;... % Orange    
    158 103 171 ;... % Purple
    204 112 087 ;... % Rust
    121 195 106 ;... % Green    
    215 127 178 ;... % Pink/purple
    114 114 114 ]/255;% Grey    
        
% Change settings for cue axes
bc = get(gcf,'Color') ;
set([AXc AXc2],'XLim',[0 1],'YLim',[0 1]) ;
set([AXc AXc2],'Box','off','XTick',[],'YTick',[],'XColor',bc,'YColor',bc,'Color',bc) ;
set([AXs AXs2 AXa],'Box','on')

% Find time difference between tags
cuediff = dtagtimediff(tag,rtag,time_shift);

% Plot the info from tag 1
[t,rl] = plotspec(tag,tcue,RES,[AXc AXs]);
axes(AXa), hold off, plot(tcue+t,rl,'k','LineWidth',1) ;
set(gca,'XAxisLocation','top','XLim',[tcue, tcue+NS],'YLim',ENV_LIM,'XTickLabel','') ;
ylabel('RL (dB)'), grid on

for i=1:length(rtag)
    % Get variables
    newtag = char(rtag(i)) ;
    newtime = tcue+cuediff(i);

    % Add spectrogram information
    [t,rl] = plotspec(newtag,newtime,[],[AXc2(i),AXs2(i)],CMAT(i,:));
    axes(AXa), hold on, plot(tcue+t,rl,'color',CMAT(i,:),'LineWidth',1) ; grid
end

% Define spectrogram handles for output to mother function
AXS=[AXs AXc ; AXs2' AXc2'];



function [t,rl] = plotspec(tag,tcue,RES,AXS,col)

global CH NS AFS_RES SPEC_BL SPEC_OL SPEC_CLIM SPEC_FLIM ENV_HP

if nargin<5
    col = [1 1 1];
end

% Get axes handles from input variables
AXc=AXS(1); AXs=AXS(2);

% Restrict y axis in spectrogram
axes(AXs), axis xy, grid ;
set(AXs,'XTickLabel',[],'YLim',SPEC_FLIM,'XLim',[tcue tcue+NS]) ;
ylabel('F, kHz')

% Read audio
[x,afs_org] = dtagwavread(tag,tcue,NS) ;

% Check if there's data
if isempty(x), 
    t=[]; rl = []; 

    text(tcue+0.02*NS,0.9*SPEC_FLIM(2),[tag(1:4) '\_' tag(6:9) ' no data'],...
    'VerticalAlignment','Top','HorizontalAlignment','Left',...
    'Backgroundcolor',[0.2 0.2 0.2],'Color',col,'FontWeight','bold');

    return, 
end

% Resample data
x=resample(x(:,CH)-mean(x(:,CH)),AFS_RES,afs_org);
afs = AFS_RES;

% Calculate spectrogram
[B,F,T] = spectrogram(x(:,CH),hamming(SPEC_BL),round(SPEC_BL*SPEC_OL),SPEC_BL,afs, 'yaxis') ;

% Format and plot spectrogram   
BB = adjust2Axis(20*log10(abs(B))) ;
axes(AXs), imagesc(tcue+T,F/1000,BB,SPEC_CLIM) ; axis xy, grid ;
text(tcue+0.02*NS,0.9*SPEC_FLIM(2),[tag(1:4) '\_' tag(6:9)],...
    'VerticalAlignment','Top','HorizontalAlignment','Left',...
    'Backgroundcolor',[0.2 0.2 0.2],'Color',col,'FontWeight','bold');
set(AXs,'XTickLabel',[],'YLim',SPEC_FLIM,'XLim',[tcue tcue+NS]) ;
ylabel('F, kHz')

% Check data for labels
if isempty(RES), RES=loadaudit(tag); end

% Add labels
if ~ (nargin<3 || isempty(RES))
    dtagaudit_plotRES(AXc,RES,[tcue tcue+NS]) ;
end

% Add RL
[BBB,AAA]=butter(4,ENV_HP/(afs/2),'high');
xf=filtfilt(BBB,AAA,x);  
N_env = AFS_RES*NS/3000;
xx = buffer(xf,N_env,N_env/2,'nodelay');
rl = 20*log10(std(xx));
t  = [1:length(rl)]*N_env/2/afs;

return
