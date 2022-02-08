function [AXS] = dtagaudit_synch_int_revise(tag,tcue,rtag,NS,time_shift)
% QUICK SCRIPT TO SEE AMPLITUDE

% INTEGRATE GLOBAL HP FILTER


if nargin<5,
    time_shift = 0 ; % default is no time shift
elseif isempty(time_shift)
    time_shift = 0 ;
end

% Make sure that there is a time shift for each rtag
if length(time_shift)==1,
    time_shift=ones(length(rtag))*time_shift;
end

if nargin<5
    NS=15;
elseif isempty(NS)
    NS=15;
end

% Make rtag string array into cell array
if ischar(rtag)
    rtag={rtag};
end
        
for i=1:length(rtag),
    if ~strcmp(tag(1:8),rtag{i}(1:8))
        disp('Error: Tag deployments must be from same day')
    end
end

% Prepare figure first
figure(97),clf,
set(gcf,'name','Revise cues')

% Create figure axes
        
% Make axis for amplitude envelopes for all whales
AXs = axes('position',[0.11,0.73,0.8,0.20]) ;
AXa = axes('position',[0.11,0.11,0.8,0.60]) ; box on, hold on

% Colors
CMAT = [238 045 046 ;... % Red
        024 089 169 ;... % Blue
        000 139 071 ;... % Green
        243 125 035 ;... % Orange
        105 042 147 ;... % Purple
        162 029 032 ;... % Rust
        179 056 147 ;... % Pink/purple
        001 001 001]/255;% Black

% Find time difference between tags
cuediff = dtagtimediff(tag,rtag,time_shift);

% Find random tag for spectrogram
k = round(1+length(rtag)*rand(1,1)) ; 
if k == 1,
    thistag = tag;
    thiscue = tcue ;
else
    thistag = char(rtag(k-1)) ;
    thiscue = tcue + cuediff(k-1) ;
end

% Construct spectrogram
plotspec_mod(thistag,thiscue,AXs,NS) ;
set(gca,'XLim',thiscue+[0 NS])

% Construct amplitude plot
[x,afs] = dtagwavread(tag,tcue,NS);
[cl xx] = findclicks(x(:,1),0,afs,1000) ;
kk = 1:5:length(xx) ;
axes(AXa), set(AXa,'XLim',tcue+[0 NS],'YLim',[-80 -10]), 
ylabel('RL (dB)'), grid on, hold on
h(1) = plot(tcue+kk/afs,20*log10(xx(kk)),'k','LineWidth',2,'Color',CMAT(1,:));
legg{1} = char(tag) ;

for i=1:length(rtag)
    thistag = char(rtag(i)) ;
    thiscue = tcue+cuediff(i) ;
    [x,afs] = dtagwavread(thistag,thiscue,NS);
    [cl xx] = findclicks(x(:,1),0,afs,1000) ;
    kk = 1:5:length(xx) ;
    h(i+1)=plot(tcue+kk/afs,20*log10(xx(kk)),'k','LineWidth',2);
    set(h(i+1),'Color',CMAT(i+1,:));
    legg{i+1} = char(rtag(i)) ;
end

for i=1:length(legg),
    legg{i}(findstr(legg{i},'_'))='-' ;
end

legend(h,legg,'Location','Northwest')



function   [cc,xx] = findclicks(x,thresh,fs,fh)
%
%     clicks = findclicks(x,thresh,fs,fh)
%     Return the time cue to each click in x, in seconds
%     x is a signal vector
%     thresh is the click detection threshold - try a value of 0.01 to
%     begin with
%     fs is the sampling rate of x [default = 32000]
%     fh is the high pass filter frequency [default = 10000]
%
%     mark johnson, WHOI
%     August, 2000
%     modified June 2003

% for sperm whales use:
tc = 2.5e-3 ;           % power averaging time constant in seconds

% for ziphius use:
tc = 0.5e-3 ;           % power averaging time constant in seconds

blanking = 20e-3 ;      % blanking time after a click is detected before another

[b,a] = cheby1(6,0.5,fh/fs*2,'high') ;
pp = 1/tc/fs ;

xf = filter(b,a,x);
xx = filter(pp,[1 -(1-pp)],abs(xf)) ;
cc = [] ;

if thresh==0,
   return
end

cc = find(diff(xx>thresh)>0)/fs ;
done = 0 ;

if isempty(cc),
   return ;
end

while ~done,
   kg = find(diff(cc)>blanking) ;
   done = length(kg) == (length(cc)-1) ;
   cc = cc([1;kg+1]) ;
end
return

function     [] = plotspec_mod(tag,tcue,AXs,NS)
% AXS a two-element vector of axis handles
% AXc=AXS(1); AXs=AXS(2);

global CH AFS_RES SPEC_BL SPEC_OL SPEC_CLIM SPEC_FLIM

if nargin<4
    NS = 4 ;           % number of seconds to display
end


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
    'Backgroundcolor',[0.2 0.2 0.2],'Color',[1 1 1],'FontWeight','bold');
axis xy, grid ;

% Restrict y axis in spectrogram
set(AXs,'XTickLabel','','YLim',SPEC_FLIM) ;
ylabel('Frequency, kHz')

return