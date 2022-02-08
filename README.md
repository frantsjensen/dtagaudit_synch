# dtagaudit_synch toolbox

This is a Matlab Toolbox for synchronizing and auditing multiple concurrently deployed DTAGs, and with a variety of functions to enable concurrent data labelling and facilitate identification of source of sounds in audio recordings. 
Code has been tested with Matlab 2021a.

Code is derived from auditing scripts developed by Dr. Mark Johnson (www.soundtags.org). Both the dtag2 and dtag3 toolboxes are required for these tools to function, and they can be downloaded here: https://www.soundtags.org/dtags/dtag-toolbox/

Additionally, tools require that tag paths have been configured. In particular, the following tag paths have to be configured, either through a startup.m file or manually at the start of each session:
settagpath('AUDIO',audiopath) - path to parent directory in which deployment year folders (mn09, tt12, pw17, depending on dataset) are located. 
settagpath('CAL',calpath) - full path where all cal files are located
settagpath('PRH',prhpath) - full path where all prh files are located
settagpath('AUDIT',audpath) - full path where all audits are located

Main function for the toolbox is dtagaudit_synch. To use:
-install matlab from www.mathworks.com (tested with Matlab 2021a but fairly backwards compatible to at least Matlab 2016a)
-download dtag2 and dtag3 matlab tools from www.soundtags.org
-start up a new matlab sesssion and make sure you have the dtag2 and dtag3 matlab tools and the dtagaudit_synch folder added to your matlab path
-define tag deployment ID, load existing audit information, and start audit at an arbitrary time
tag = 'tt12_131a';
R = loadaudit(R);
R = dtagaudit_synch(tag,sttime,R);

Press help dtagaudit_synch for a description of different options while running audit script.

-at the end of the audit session, save your updated audit.

Adjust settings such as spectrogram parameters and other variables in the dtagaudit_settings file. 
Settings are based on the first two letters (species ID) in the tag deployment id, so define a new option if needed.



