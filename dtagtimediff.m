% cuediff = dtagtimediff(tag,rtag,time_shift)
% 
% Estimate the difference in time frames between two tags deployed on the same day
% tag is the primary tag, given as a character string
% rtag are the secondary tags, entered either as a character string, or as
%   a cell array, of tag names.
% cuediff will be the number of seconds added to the time cue from primary
%   tag in order to get the time cue for the secondary tag
% Optional input argument time_shift is a vector of time offsets to add to
%   nominal time offset as reported by tag clocks
%
% Made compatible with dtag2 and dtag3/2015
%
% F. H. Jensen, 2013

function cuediff = dtagtimediff(tag,rtag,time_shift)

if nargin < 3
    time_shift = zeros(length(rtag));
end

% Convert character string to cell array
if isstr(rtag), rtag = {rtag}; end

% Preallocate cue difference vector
cuediff = zeros(length(rtag),1);

% INSERT HERE EXACT TIME DIFFERENCE METHOD


% Find approximate time difference between tags

% Find tag version
tagver = dtagtype(tag);

% Find time of tag on for D2
if tagver==2,
    [c t s ttype fs id] = tagcue(10,tag) ;
    if isempty(c), return, end
    t1=datevec(datenum(t-[zeros(1,5) 10])) ; % subtract the 10s from tagcue call
    for i=1:length(rtag),
        [c t s ttype fs id] = tagcue(10,char(rtag(i))) ;
        if isempty(c), return, end
        t2=datevec(datenum(t-[zeros(1,5) 10])) ; % subtract the 10s from tagcue call
        cuediff(i) = etime(t1,t2) + time_shift(i) ; % Number of seconds added to tcue to get t2cue
    end
    
% Find time of tag on for D3    
elseif tagver==3,
    [CAL,DEPLOY] = d3loadcal(tag);
    t1 = DEPLOY.SCUES.TIME(1,:) ;
    for i=1:length(rtag)
        [CAL,DEPLOY] = d3loadcal(char(rtag(i)));
        t2 = DEPLOY.SCUES.TIME(1,:) ;
        cuediff(i) = etime(t1,t2) + time_shift(i) ; % Number of seconds added to tcue to get t2cue
    end
end

% Check that cuediff is within the same half-day
if max(cuediff)>12*3600,
    disp('Warning: More than 12 hour difference in tag start time. Check CAL files')
end