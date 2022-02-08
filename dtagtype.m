function [DTAG] = dtagtype(tag)
% Find if DTAG is version 2 or 3 and switch between them
% Looks for type of cal file, so only works if you have cal files in cal
% directory

shortname = tag([1:2 6:9]) ;
subdir = tag(1:4) ;

global TAG_PATHS
if isempty(TAG_PATHS) | ~isfield(TAG_PATHS,'CAL'),
   if isempty(SILENT),
      fprintf(' No %s file path - use settagpath\n', 'CAL') ;
   end
   return
end

% If DTAG2, CAL file is stored as .MAT
d2suffix = strcat(tag,'cal.mat') ;
d2cal = sprintf('%s/%s',getfield(TAG_PATHS,'CAL'),d2suffix) ;

% If DTAG3, CAL file is stored as .XML
d3suffix = strcat(tag,'cal.xml') ;
d3cal = sprintf('%s/%s',getfield(TAG_PATHS,'CAL'),d3suffix) ;

% Now check which one is true
if exist(d2cal)
    DTAG=2;
elseif exist(d3cal)
    DTAG=3;
else
    DTAG=[];
    disp('No CAL file found when evaluating function dtagtype')
end
    