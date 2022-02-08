function dtagaudit_settings(audit_type)

% Configure DTAGAUDIT settings

if nargin<1
    audit_type = 'default';
end

% All variables will be made global
global CH NS AFS_RES SPEC_BL SPEC_OL SPEC_CLIM SPEC_FLIM ENV_LIM ENV_HP

global foc_labels nonfoc_labels
foc_labels = {'focal','foc'};
nonfoc_labels = {'nonfoc','non','nf'};

switch lower(audit_type)

    % Settings for Bottlenose dolphins
    case 'tt' % Tursiops truncatus

        CH = 1 ;            % which channel to display if multichannel audio
        NS = 10;            % number of seconds to display
        AFS_RES = 60e3 ;    % Resample sound to limit data and speed up specgram
        SPEC_BL = 512 ;     % spectrogram (fft) block size
        SPEC_OL = 0.75 ;    % spectrogram overlap
        SPEC_CLIM = [-90,-10] ; % color axis limits in dB for specgram
        SPEC_FLIM = [0,30]; % Spectrogram frequency limit in kHz (but max half of AFS_RES)   
        ENV_LIM = [-70,-20];% Envelope limits
        ENV_HP  = 1000;     % Envelope high-pass filter

        foc_labels = {'sw','swc','w','wM','rs','bp'};
        nonfoc_labels = {'nf'};
        
    % Settings for pilot whales
	case {'pw','gm'} 

        CH = 1 ;            % which channel to display if multichannel audio
        NS = 15;            % number of seconds to display
        AFS_RES = 60e3 ;    % Resample sound to limit data and speed up specgram
        SPEC_BL = 512 ;     % spectrogram (fft) block size
        SPEC_OL = 0.75 ;    % spectrogram overlap
        SPEC_CLIM = [-90,-10] ; % color axis limits in dB for specgram
        SPEC_FLIM = [0,30]; % Spectrogram frequency limit in kHz (but max half of AFS_RES)   
        ENV_LIM = [-70,-20];% Envelope limits
        ENV_HP  = 2000;     % Envelope high-pass filter

        foc_labels = {'sk','skp','bp','rs','w'};
        nonfoc_labels = {'nf'};        
        
    % Settings for baleen whales (may need to be differentiated)
    case {'mn','bp','ea'}
        
        CH = 1 ;            % which channel to display if multichannel audio
        NS = 30;            % number of seconds to display
        AFS_RES = 24e3 ;    % Resample sound to limit data and speed up specgram
        SPEC_BL = 1024 ;    % spectrogram (fft) block size
        SPEC_OL = 0.75 ;    % spectrogram overlap
        SPEC_CLIM = [-60,10];% color axis limits in dB for specgram
        SPEC_FLIM = [0,8];  % Spectrogram frequency limit in kHz (but max half of AFS_RES)   
        ENV_LIM = [-60,0];% Envelope limits
        ENV_HP  = 100;       % Envelope high-pass filter    

    % Default settings
    otherwise
        
        CH = 1 ;            % which channel to display if multichannel audio
        NS = 10;            % number of seconds to display
        AFS_RES = 60e3 ;    % Resample sound to limit data and speed up specgram
        SPEC_BL = 512 ;     % spectrogram (fft) block size
        SPEC_OL = 0.75 ;    % spectrogram overlap
        SPEC_CLIM = [-90,-10] ; % color axis limits in dB for specgram
        SPEC_FLIM = [0,30]; % Spectrogram frequency limit in kHz (but max half of AFS_RES)   
        ENV_LIM = [-70,-20];% Envelope limits
        ENV_HP  = 2000;     % Envelope high-pass filter  
        
end

return