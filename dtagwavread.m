function [x,afs] = dtagwavread(tag,tcue,tdur);
% function [x,afs] = dtagwavread(tag,tcue,tdur);
% Switch between dtag2 and dtag3 wavread functions
% FHJ

tagver = dtagtype(tag);

if tagver==3,
    recdir = d3makefname(tag,'RECDIR'); % get audio path
    prefix = [tag(1:2) tag(6:9)];
    
    D = dir([recdir prefix '*.wav']);
    if isempty(D)
        prefix = tag;
        D = dir([recdir prefix '*.wav']);
    end
    [x,afs] = d3wavread([tcue tcue+tdur],recdir,prefix , 'wav' ) ;
elseif tagver==2,
    [x,afs] = tagwavread(tag,tcue,tdur);
end
    
