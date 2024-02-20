    function     [tcue] = dtagaudit_main(tag,tcue)
%
%     tcue = dtagaudit_main(tag,tcue)
%
%     DTAG audit tool compatible with DTAG versions 1 through 3
%     [Support for simultaneous tags plus multiple smaller hacks]
%
%     Input:
%     - tag:    deployment ID of tag data 
%     - tcue:   Start time of audio block to load
%     Output:
%     - tcue:   Start time of last audioblock
%
%     Audit files will be saved in configured audit folder;
%       file name will be given as deploymentIDaud.txt
%
%     OPERATION
%     Type or click on the display for the following functions:
%     - type 'f' to go to the next block
%     - type 'b' to go to the previous block
%     - click on the graph to get the time cue, depth, time-to-lastf
%       and frequency of an event. Time-to-last is the elapsed time 
%       between the current click point and the point last clicked. 
%       Results display in the matlab command window.
%     - type 's' to select the current segment and add it to the audit.
%       You will be prompted to enter a sound type on the matlab command
%       window. Enter a single word and type return when complete.
%     - type 'l' to select the currect cursor position and add it to the 
%       audit as a 0-length event. You will be prompted to enter a sound 
%       type on the matlab command window. Enter a single word and type 
%       return when complete.
%     - type 'x' to change the audit entry at the cursor position.
%       If there is no audit entry at the cursor, nothing happens.
%       Otherwise, a prompt will appear with the current cue type
%       for all cues around cursor position. Delete text to remove that
%       cue, or change it to update it
%     - type 'i' to play the current segment and simultaneously save to
%       current folder as tempsound.wav. Rename to keep sound file.
%     - type 'p' to play the displayed sound segment 
%       through the computer speaker/headphone jack.
%     - type 'q' or press the right hand mouse button to finish auditing.
%     - type 'a' to report the angle of arrival of the selected segment
%
%     Multiple tag comparison package (separate files required):
%     - type 'o' to set up multiple tag comparison.
%     - type 't' to compare spectrograms and received levels across tags
%                (will open separate window)
%     - type 'r' to compare intensity across tags
%     - type 'T' to compare angle-of-arrival across tags (separate window)
%     - type 'S' to start marking sequences on non-focal tags (requires 't' first):
%           -NOTE crosshairs should indicate which window is active!
%                Press twice in same axis to activate sequence
%                Press 's' to save active sequence as cue
%                Press 'l' to mark instant point
%                Press 'x' to edit any cue types close to cursor
%                   -change text to update cue type
%                   -remove text to delete cue
%                Press any other button to go back to auditing window
%
%     OBS: Change focal and non-focal labels in dtagaudit_plotRES.m to
%     colorcode focal and non-focal sounds
%
%   Original code by Mark Johnson, www.soundtags.org
%   Modified by Frants Jensen, frants.jensen@gmail.com

% SHARED GLOBAL VARIABLES
global CH NS AFS_RES SPEC_BL SPEC_OL SPEC_CLIM SPEC_FLIM ENV_LIM ENV_HP AOA_SCF
dtagaudit_settings(tag(1:2));

% UNSHARED VARIABLES
SOUND_FH = 0 ;   % high-pass filter for sound playback - 0 for no filter
SOUND_DF = 1 ;     % decimation factor for playing sound
volume = 8 ;       % amplification factor for audio output - often needed to
                   % hear weak signals (if volume>1, loud transients will
                   % be clipped when playing the sound cut
QUERY_LABEL = 0;   % Ask for extra label when pressing "i" and saving sound

% Warning if using other than CH1
if ~CH==1
    disp([ ' Dtagaudit_synch warning: using ch ' num2str(CH)]);
end

% Default Settings for multiple tag comparisons
possibletags = '' ; othertags = '' ; AXS_synch = [] ; time_shift = 0  ;

% Find tag version
tagver = dtagtype(tag);

% Decide on hydrophone separation
if tagver==3
    AOA_SCF = 1500/0.045 ; % sound speed / hydrophone separation (adjusted for D3)
elseif tagver==2
    AOA_SCF = 1500/0.0228 ; % v/h % This should be /0.0228 (0.9 inch)
end

% Check if there's a recent synch audit saved settings:
synch_fname = [tag(1:9) '_synch.mat'] ;
if exist(synch_fname)
    load(synch_fname);
    disp(['Loading settings from ' synch_fname])
end

% Check if audit structure exists
% Load pre-existing audit information
RES = loadaudit(tag);
if isempty(RES)
	disp(' WARNING - NO PRE-EXISTING AUDIT FOUND!')
	RES.cue = [] ;
	RES.stype = [];
end

% Check frequency limits
SPEC_FLIM(2) = min( [SPEC_FLIM(2) AFS_RES/2000]);

% Load prh
k = loadprh(tag,'p','fs') ;           % read p and fs from the sensor file
if k==0
   fprintf('Unable to find a PRH file - continuing without\n') ;
   p = [] ; fs = [] ;
end


%%%%%%%%%%%%%%%%%%% Prepare filters %%%%%%%%%%%%%%%%%%%%%

% check sampling rate and reset to AFS_RESAMPLE if too high
[~,afs] = dtagwavread(tag,tcue,0.1);
if afs>AFS_RES, afs=AFS_RES; end

% filters for sound playback
if SOUND_FH > 0
   [bs,as] = butter(6,SOUND_FH/(afs/2),'high') ;
else
   bs = [] ;
end

%%%%%%%%%%%%%%%%% Prepare figure plot %%%%%%%%%%%%%%%%%%%

current = [0 0] ;
figure(1),clf
if ~isempty(p)
   kb = 1:floor(NS*fs) ;
   AXm = axes('position',[0.11,0.76,0.78,0.18]) ;
   AXc = axes('position',[0.11,0.70,0.78,0.05]) ;
   AXs = axes('position',[0.11,0.34,0.78,0.35]) ;
   AXp = axes('position',[0.11,0.11,0.78,0.2]) ;
else
   AXm = axes('position',[0.11,0.60,0.78,0.34]) ;
   AXc = axes('position',[0.11,0.52,0.78,0.07]) ;
   AXs = axes('position',[0.11,0.11,0.78,0.38]) ;
end

bc = get(gcf,'Color') ;
set(AXc,'XLim',[0 1],'YLim',[0 1]) ;
set(AXc,'Box','off','XTick',[],'YTick',[],'XColor',bc,'YColor',bc,'Color',bc) ;

while 1
   [x,afs_org] = dtagwavread(tag,tcue,NS);
   if isempty(x), return, end
   
   % Resample to conserve memory
   x=resample(x(:,CH)-mean(x(:,CH)),AFS_RES,afs_org);

   % Filter to get rid of zero offset
   [BBB,AAA]=butter(4,ENV_HP/(afs/2),'high');
   xf=filtfilt(BBB,AAA,x);      

   % Create top plot
   N_env = AFS_RES*NS/3000;
   xx = buffer(xf,N_env,N_env/2,'nodelay');
   rl = 20*log10(std(xx));
   t  = [1:length(rl)]*N_env/2/afs;
   axes(AXm), plot(tcue+t,rl,'k') ; grid
   set(AXm,'XAxisLocation','top','XLim',[tcue, tcue+NS],'YLim',ENV_LIM) ;
   ylabel('Intensity, dB')
   
   % Add audit cues
   dtagaudit_plotRES(AXc,RES,[tcue tcue+NS]) ;

   % Create bottom plot
   if ~isempty(p)
      ks = kb + round(tcue*fs) ;
      axes(AXp),plot(ks/fs,p(ks)), grid
   	  set(gca,'YDir','reverse','XLim',[tcue tcue+NS]);
      xlabel('Time, s'), ylabel('Depth, m')
   end
   
   % Construct spectrogram plot
   [B,F,T] = spectrogram(x,hamming(SPEC_BL),round(SPEC_BL*SPEC_OL),SPEC_BL,afs, 'yaxis') ;
   BB = adjust2Axis(20*log10(abs(B))) ;
   axes(AXs), imagesc(tcue+T,F/1000,BB,SPEC_CLIM) ; axis xy; grid ;
   if ~isempty(p)
      set(AXs,'XTickLabel',[]) ;
   else
      xlabel('Time, s')
   end
   ylabel('Frequency, kHz')
   
   % Set axes
   axis([tcue, tcue+NS, SPEC_FLIM])
   
   % Add current selection to spectrogram plot   
   hold on
   hhh = plot([0 0],0.8*min([afs/2000 SPEC_FLIM(2)])*[1 1],'k*-') ;    % plot cursor
   hold off
   

   done = 0 ;
   while done == 0
      axes(AXs) ; pause(0) ;
      [gx, gy, button] = ginput(1) ;
      if isempty(button)
         % This counts for volume up/down buttons
      

% QUIT AUDITING PROGRAM AND RETURN CUE STRUCTURE
      elseif button==3 || button=='q'
         saveaudit(tag,RES);
         disp(['Audit saved to audit folder: ' gettagpath('AUDIT')])
         return

% INSERT SEQUENCE CUE
      elseif button=='s'
         ss = input(' Enter comment... ','s') ;
         cc = sort(current) ;
         RES.cue = [RES.cue;[cc(1) diff(cc)]] ;
         RES.stype{size(RES.cue,1)} = ss ;
         saveaudit(tag,RES);
         dtagaudit_plotRES(AXc,RES,[tcue tcue+NS]) ;

% CHECK ANGLE OF ARRIVAL OF SEQUENCE         
      elseif button=='a'
          try
              length_temp=0.1*ceil(abs(diff(current)*10));             
              dtagaudit_synch_int_angle(tag,min(current),length_temp);
          catch
              disp(' An error occurred during dtagaudit_synch_int_angle and action was aborted ')
          end
          
% INSERT INSTANTANEOUS CUE AT CURSOR POSITION
      elseif button=='l'
         ss = input(' Enter comment... ','s') ;
         RES.cue = [RES.cue;[gx 0]] ;
         RES.stype{size(RES.cue,1)} = ss ;
         saveaudit(tag,RES);
         dtagaudit_plotRES(AXc,RES,[tcue tcue+NS]) ;

% DELETE OR CHANGE CUE AT CURSOR POSITION
      elseif button=='x'
        % Find all cues around synch cursor
        kres =(find(gx>=RES.cue(:,1)-0.1 & gx<sum(RES.cue')'+0.1)) ;
        if ~isempty(kres)
            % Find the type of cues and allow users to edit them
            ktype = RES.stype(kres) ;
            ktype=inputdlg(ktype,'Edit/remove cues',[1,40],ktype);
            % If cancel button is pressed, ktype will be empty cell
            if ~isempty(ktype)
                % Otherwise, go through cell array of new sound types
                for i=1:length(ktype)
                    RES.stype(kres(i))=ktype(i);
                    kempty(i)=isempty(ktype{i});
                end

                % If there are empty sound types, remove them
                kempty=kres(find(kempty));
                kkeep = setxor(1:size(RES.cue,1),kempty) ;
                RES.cue = RES.cue(kkeep,:) ;
                RES.stype = {RES.stype{kkeep}} ;
            end
            saveaudit(tag,RES);
            dtagaudit_plotRES(AXc,RES,[tcue tcue+NS]) ;
            clear kempty kres

         else
            fprintf(' No saved cue at cursor\n') ;
         end          

% GO FORWARD ONE SCREEN
      elseif button=='f'
            tcue = tcue+floor(NS)-0.5 ;
            done = 1 ;

% GO BACK ONE SCREEN            
      elseif button=='b'
            tcue = max([0 tcue-NS+0.5]) ;
            done = 1 ;

% PLAY ALL SOUND            
      elseif button=='p'
            if ~isempty(bs)
               xf = filter(bs,as,x) ; 
               sound(volume*xf,afs/SOUND_DF,16) ;
            else
               sound(volume*x,afs/SOUND_DF,16) ;
            end

% PLAY SHORT SELECTION AND SAVE TEMP FILE            
      elseif button=='i'
            if all(current>0)
                length_temp=0.1*ceil(abs(diff(current)*10));
                [x_temp,afs_org] = dtagwavread(tag,min(current),length_temp);
                if isempty(x), return, end
   
                % First, query for sound label
                if QUERY_LABEL
                    label=inputdlg({'Sound label'},'Add label (skA, skB, etc) to save',[1,40],{''});
                else
                    label=[];
                end
                if ~isempty(label)
                    %outputname = ['call types/' tag '_' num2str(round(min(current))) '_' char(label) '.wav'];
                    outputname = [tag '_' num2str(round(min(current))) '_' label '.wav'];
                else
                    outputname = 'tempsound.wav';
                end
                
                % Then try to save
                try
                    audiowrite(outputname,x_temp(:,CH),afs_org);
                catch
                    disp([' Error writing to ' outputname ', close other programs and try again'])
                end
                
                % Finally, filter and play sound
                xf = resample(x_temp,AFS_RES,afs_org);                    
                if ~isempty(bs)
                    xf = filter(bs,as,xf) ;
                end
                xf = xf.*tukeywin(length(xf),0.2./(length(xf)/AFS_RES));
                sound(volume*xf,AFS_RES/SOUND_DF,16) ;
            else
                fprintf('Invalid click: Need to mark interval before pressing i to play sound from this interval\n')                
            end

% ADJUST SYNCH AUDITING PARAMETERS            
      elseif button=='o'
          % Find out which tags are possible
          if isempty(possibletags)
              letters='abcdefghij'; n=0;
              % Check for other tags from same day
              for i=1:length(letters)
                  % Rewrite this to find all possible tags, figure out when
                  % they were on and off the animal, and find estimated
                  % time difference - then keep track of this in structure
                  if tagver==2
                      othertagfile = makefname([tag(1:8) letters(i)],'AUDIO',1);
                  elseif tagver==3
                      othertagfile = d3makefname([tag(1:8) letters(i)],'AUDIO',1);
                  end
                  slashes = [union(findstr(othertagfile,'/'),findstr(othertagfile,'\'))];
                  othertagfile = othertagfile(1:slashes(end));
                  if exist(othertagfile) % compatible w d2
                      if ~strcmp(tag,[tag(1:8) letters(i)]) % Exclude focal tag
                          n=n+1;                            % Add tag to possible tag list
                          possibletags{n}=[tag(1:8) letters(i)]; 
                      end
                  end
              end
          end

          % Choose tags for comparing
          InitValue=[];
          if ~isempty(othertags)
              for i=1:length(othertags)
                  if strmatch(char(othertags(i)),possibletags,'exact')
                      InitValue=[InitValue strmatch(char(othertags(i)),possibletags,'exact')];
                  end
              end
          end
          [kd,ok] = listdlg('ListString',possibletags,'SelectionMode','multiple',...
          'InitialValue',InitValue,'OKString','Accept','Name','Tags',...
          'PromptString',['Select tags to compare with focal ' tag],...
          'ListSize',[300 100]);
        
          if ~ok
              continue
          end
          othertags = possibletags(kd) ;

          % Check time shifts available for each tag
          if ~(length(othertags)==length(time_shift))
              time_shift=zeros(1,length(othertags));
          end

          % Adjust initialvalue
          clear InitValue
          for i=1:length(othertags)
              InitValue{i}=num2str(time_shift(i));
          end

          % Adjust time shifts for each tag
          kd=inputdlg(othertags,'Time offset',[1,30],InitValue);
          time_shift = zeros(length(othertags),1);
          for i=1:length(othertags)
              if ~isnan(str2double(kd{i}))
                  time_shift(i)=str2double(kd{i});
              else
                  time_shift(i)=0;
              end
          end
          
          % Save synch audit settings for quick retrieval
          synch_fname = [tag(1:9) '_synch.mat'] ;
          save (synch_fname,'othertags','time_shift')

      elseif button=='r'
          if isempty(othertags)
              disp(' First adjust synch parameters using "o" ')
              continue
          end
          
          try
              length_temp=0.1*ceil(abs(diff(current)*10));
              dtagaudit_synch_int_revise(tag,min(current),othertags,length_temp,time_shift);
          catch
              disp(' An error occurred during tagaudit_synch_int_revise and action was aborted ')
              disp([' Trying to compare ' tag ' with the following tags:' ])
              disp(char(othertags))
          end

          
% SYNCH SOUND COMPARISON WITH OTHER TAGS          
      elseif button=='t'
          if isempty(othertags)
              disp(' First adjust synch parameters using "o" ')
              continue
          end
          
          if length(othertags)<6
              try
                  [AXS_synch]=dtagaudit_synch_int(tag,tcue,RES,othertags,time_shift);
              catch
                  disp(' An error occurred during dtagaudit_synch_int and action was aborted ')
                  disp([' Trying to compare ' tag ' with the following tags:' ])
                  disp(char(othertags))
              end
          else
              disp('dtagaudit_synch currently supports comparing focal tag with up to 4 other tags')
          end

          
% SYNCH ANALYSIS OF ANGLE-OF-ARRIVAL ON MULTIPLE TAGS
      elseif button=='T'
          if isempty(othertags)
              disp(' First adjust synch parameters using "o" ')
              continue
          end
          
          try
              dtagaudit_synch_int_aoa(tag,tcue,RES,othertags,time_shift);
          catch
              disp(' An error occurred during dtagaudit_synch_aoa_int and action was aborted ')
              disp([' Trying to compare ' tag ' with the following tags:' ])
              disp(char(othertags))
          end
          
% MARK SEQUENCE ON NON-FOCAL TAGS
      elseif button=='S'
          if isempty(othertags)
              disp(' Cannot mark multiple sequences. First adjust synch parameters using "o" ')
              continue
          end
          
          % Check if simultaneous spectrogram plot and handles exist
          if isempty(findobj('type','figure','name','Simultaneous Spectrogram Plot'))
              disp(' Cannot mark multiple sequences. Need to display simultaneous spectrograms ')
              continue
          elseif isempty(AXS_synch)
              disp(' Cannot mark multiple sequences. Need to display simultaneous spectrograms first')
              continue
          end
          

          % Call that figure
          figure(findobj('type','figure','name','Simultaneous Spectrogram Plot'))

          % Set up artificial cursors in each synch plot, make them
          % invisible until changed to visible
          for i=1:length(AXS_synch(:,1))
              axes(AXS_synch(i,1));
              hold on, hhh_synch(i) = plot([0 0],0.8*min([afs/2000 SPEC_FLIM(2)])*[1 1],'k*-','visible','off') ; hold off
          end              
          
          MARK=1;
          current_synch = [0 0] ;
          axis_selected = [0 0] ;

          % Activate relevant axis to bring out crosshairs from ginput
          axes(AXS_synch(1,1));
              
          while MARK
            
              % Find two points
              [gx_synch, gy_synch, button_synch] = ginput(1) ;
          
              % Verify that user clicked in right figure
              if gcf~=findobj('type','figure','name','Simultaneous Spectrogram Plot') 
                  continue,
              end
          
              if button_synch==1
                  if gy_synch<0 || gx_synch<0 || isempty(find(AXS_synch(:,1)==gca,1))
                      disp('Invalid click')
                  else
                      % Update selection
                      current_synch=[current_synch(2) gx_synch]; 
                      axis_selected=[axis_selected(2) gca];
                      % Make all synch cursors invisible
                      set(hhh_synch,'visible','off')
                      % Then activate the cursor that is active and update position
                      thistag=find(AXS_synch(:,1)==axis_selected(2));
                      set(hhh_synch(thistag),'visible','on','XData',current_synch(axis_selected==gca),'YData',0.8*min([afs/2000 SPEC_FLIM(2)])*ones(length(find(axis_selected==gca)),1));
                  end

                  
              elseif button_synch=='s'
                  
                  if ~(axis_selected(1)==axis_selected(2)) || all(axis_selected==0)
                      disp('Click twice in same axis before marking sequence')
                  else
                  
                      % Find axis that was clicked in and corresponding tag
                      thistag=find(AXS_synch(:,1)==axis_selected(2));

                      % Decide cue type
                      if thistag==1
                          ss=char(inputdlg(tag,'Cue type',[1,40]));
                          if ~isempty(char(ss))
                              cc = sort([current_synch]) ;
                              RES.cue = [RES.cue;[cc(1) diff(cc)]] ;
                              RES.stype{size(RES.cue,1)} = ss ;
                              saveaudit(tag,RES);
                              dtagaudit_plotRES(AXc,RES,[tcue tcue+NS]); % Update auditing plot
                              dtagaudit_plotRES(AXS_synch(thistag,2),RES,[tcue tcue+NS]); % Update synch plot
                          end
                      else
                          % Load audits from other tags, 
                          % add cue to cue list, and save again
                          ss=char(inputdlg(othertags(thistag-1),'Cue type',[1,40]));
                          if ~isempty(char(ss))
                              cc = sort([current_synch]) ;
                              RESTEMP=loadaudit(char(othertags(thistag-1)));
                              RESTEMP.cue = [RESTEMP.cue;[cc(1) diff(cc)]] ;
                              RESTEMP.stype{size(RESTEMP.cue,1)} = ss ;
                              saveaudit(char(othertags(thistag-1)),RESTEMP) ;
                              dtagaudit_plotRES(AXS_synch(thistag,2),RESTEMP,get(AXS_synch(thistag,2),'XLIM')) ;
                              clear RESTEMP
                              fclose('all') % Force matlab to close any files that are still open
                          end
                      end
                  end


              elseif button_synch=='i'
                  
                  if ~eq(axis_selected(1),axis_selected(2)) || axis_selected(1)==0 || diff(current_synch)>NS
                      disp('Click twice in same axis before playing sequence')
                  else
                  
                      % Find axis that was clicked in and corresponding tag
                      thistag=find(AXS_synch(:,1)==axis_selected(2));

                      % Define what tag recording is selected
                      if thistag==1
                          play_tag = tag ;
                      else
                          play_tag = char(othertags(thistag-1)) ;
                      end
                      
                      if all(current_synch>0)
                            length_temp=0.1*ceil(abs(diff(current_synch)*10));
                            [x_temp,afs_org] = dtagwavread(play_tag,min(current_synch),length_temp) ;
                            if isempty(x), return, end

                            % First, query for sound label
                            if QUERY_LABEL
                                label=inputdlg({'Sound label'},'Add label (skA, skB, etc) to save',[1,40],{''});
                            else
                                label=[];
                            end
                            if ~isempty(label)
                                outputname = ['call types/' play_tag '_' num2str(round(min(current_synch))) '_' char(label) '.wav'];
                                %outputname = [play_tag '_' num2str(round(min(current_synch))) '.wav'];
                            else
                                outputname = 'tempsound.wav';
                            end

                            % Then, write original sound extract to disk
                            try
                                audiowrite(outputname,x_temp(:,CH),afs_org);
                            catch
                                disp([ ' Error writing to ' outputname ', close other programs and try again'])
                            end

                            % Finally, play back sound
                            if ~isempty(bs)
                                x_temp = filter(bs,as,x_temp) ;
                            end
                            xf = resample(x_temp,48e3,afs_org);
                            sound(volume*xf,48e3/SOUND_DF,16) ;
                      else
                            fprintf('Invalid click: Need to mark interval before pressing i to play sound from this interval\n')                
                      end

                  end
                  
              elseif button_synch=='l'
                  
                  % Cancel if selection is not within a valid axis
                  if isempty(find(AXS_synch(:,1)==gca,1))
                      disp('Click in valid spectrogram axis to activate that axis first')
                      continue
                  end
                  
                  % Buttons other than left clicking does not update selected
                  % axis - instead, check that click is within selected axis
                  axislimits=get(axis_selected(2),'YLim');
                  if gy_synch<axislimits(1) || gy_synch>axislimits(2)
                      disp('To change audit cue, press x while hovering in active spectrogram')
                      continue
                  end
                
                  % Find axis that was clicked in and corresponding tag
                  thistag=find(AXS_synch(:,1)==axis_selected(2));

                  if thistag==1
                      ss=char(inputdlg(tag,'Cue type',[1,40]));
                      if ~isempty(char(ss))
                          RES.cue = [RES.cue;[gx_synch 0]] ;
                          RES.stype{size(RES.cue,1)} = ss ;
                          saveaudit(tag,RES);
                          dtagaudit_plotRES(AXc,RES,[tcue tcue+NS]) ;
                          dtagaudit_plotRES(AXS_synch(thistag,2),RES,[tcue tcue+NS]); % Update synch plot
                      end
                  else
                      % Load audits from other tags, 
                      % add cue to cue list, and save again
                      ss=char(inputdlg(othertags(thistag-1),'Cue type',[1,40]));
                      if ~isempty(char(ss))
                          RESTEMP=loadaudit(char(othertags(thistag-1)));
                          RESTEMP.cue = [RESTEMP.cue;[gx_synch 0]] ;
                          RESTEMP.stype{size(RESTEMP.cue,1)} = ss ;
                          saveaudit(char(othertags(thistag-1)),RESTEMP) ;
                          dtagaudit_plotRES(AXS_synch(thistag,2),RESTEMP,get(AXS_synch(thistag,2),'XLIM')) ;
                          clear RESTEMP
                      end
                  end
                 
              % Delete cue at cursor
              elseif button_synch=='x'

                    % Buttons other than left clicking does not update selected
                    % axis - instead, check that click is within selected axis
                    axislimits=get(axis_selected(2),'YLim');
                    if gy_synch<axislimits(1) | gy_synch>axislimits(2)
                        disp('To change audit cue, press x while hovering in active spectrogram')
                        continue
                    end

                    % Find axis that was clicked in and corresponding tag
                    thistag=find(AXS_synch(:,1)==gca);

                    % If tag is focal, RES is current audit structure
                    if thistag==1
                        % Find all cues around synch cursor
                        kres =(find(gx_synch>=RES.cue(:,1)-0.1 & gx_synch<sum(RES.cue')'+0.1)) ;
                        if ~isempty(kres)
                            % Find the type of cues and allow users to edit them
                            ktype = RES.stype(kres) ;
                            ktype=inputdlg(ktype,'Edit/remove cues',[1,40],ktype);
                            % If cancel button is pressed, ktype will be empty cell
                            if ~isempty(ktype)
                                % Otherwise, go through cell array of new sound types
                                for i=1:length(ktype)
                                    RES.stype(kres(i))=ktype(i);
                                    kempty(i)=isempty(ktype{i});
                                end

                                % If there are empty sound types, remove them
                                kempty=kres(find(kempty));
                                kkeep = setxor(1:size(RES.cue,1),kempty) ;
                                RES.cue = RES.cue(kkeep,:) ;
                                RES.stype = {RES.stype{kkeep}} ;

                                % Update information in single audit and 
                                % simultaneous audit screen
                                saveaudit(tag,RES);                                
                                dtagaudit_plotRES(AXc,RES,[tcue tcue+NS]) ;
                                dtagaudit_plotRES(AXS_synch(thistag,2),RES,[tcue tcue+NS]); % Update synch plot
                            end
                            clear kempty kres

                        else
                            fprintf(' No saved cue at cursor\n') ;
                        end

                    % If tag is non-focal, load audit as RESTEMP    
                    else
                        % WARNING: BACK UP AUDIT TXT FILES AND TEST 
                        % THOROUGHLY ANY CHANGES MADE TO THIS SECTION  

                        % Load audits from other tag, 
                        RESTEMP=loadaudit(char(othertags(thistag-1)));
                        % Find all cues around cursor
                        kres =(find(gx_synch>=RESTEMP.cue(:,1)-0.1 & gx_synch<sum(RESTEMP.cue')'+0.1)) ;
                        if ~isempty(kres)
                            % If cancel button is pressed, ktype will be empty cell
                            ktype = RESTEMP.stype(kres) ;
                            ktype=inputdlg(ktype,'Edit cues (empty to remove)',[1,30],ktype);
                            if ~isempty(ktype)
                                for i=1:length(ktype)
                                    % Otherwise, go through cell array of new sound types
                                    RESTEMP.stype(kres(i))=ktype(i);
                                    kempty(i)=isempty(ktype{i});
                                end
                                % Remove empty sound types
                                kempty=kres(find(kempty));
                                kkeep = setxor(1:size(RESTEMP.cue,1),kempty) ;
                                RESTEMP.cue = RESTEMP.cue(kkeep,:) ;
                                RESTEMP.stype = {RESTEMP.stype{kkeep}} ;
                                % Now save audit again
                                saveaudit(char(othertags(thistag-1)),RESTEMP) ;
                                % Update synch plot
                                dtagaudit_plotRES(AXS_synch(thistag,2),RESTEMP,get(AXS_synch(thistag,2),'XLIM')); 
                            end
                            clear kempty kres
                        else
                            fprintf(' No saved cue at cursor\n') ;
                        end
                        clear RESTEMP
                    end

              % If other button is pressed, including q, go back to single auditing screen   
              else
                  disp(' Aborting nonfocal call marking')
                  MARK=0; continue
              end
          end
                


% ADJUST AUDIT SELECTION          
      elseif button==1
         if gy<0 || gx<tcue || gx>tcue+NS
            fprintf('Invalid click: commands are f b s l p x q\n')

         else
            current = [current(2) gx] ;
            set(hhh,'XData',current) ;
            if ~isempty(p)
               fprintf(' -> %6.1f\t\tdiff to last = %6.2f\t\tp = %6.1f\t\tfreq. = %4.2f kHz\n', ...
                 gx,diff(current),p(round(gx*fs)),gy) ;
			else
               fprintf(' -> %6.1f\t\tdiff to last = %6.1f\t\tfreq. = %4.2f kHz\n', ...
                 gx,diff(current),gy) ;
	        end
         end
      end
   end
end

return
