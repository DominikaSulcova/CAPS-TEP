function  filename = initialize_logfile(session_info)
% ----------------------------------------------------------------------
% Author: Dominika
% Fcn: creates a text logfile that serves as documentation for a CAPS-TEP
%      experimental session
% Input: session_info = a cell array with 5 cells
%           --> project name, subject number, session number, condition,
%           treated side
% ----------------------------------------------------------------------
% get the date
c = fix(clock);
date = [num2str(c(3)) '/' num2str(c(2)) '/' num2str(c(1))];
clear c

% define the filename
filename = [session_info{1} '_' session_info{2} '_S' session_info{3} '.txt'];

% get the experimental condition
switch session_info{4}
    case '1'
        condition = 'capsaicin';
    case '0'
        condition = 'vehicle';
end

% get the side of application
switch session_info{5}
    case '1'
        side = 'right';
        hemisphere = 'left';
    case '0'
        side = 'left';
        hemisphere = 'right';
end

% write the file 
fileID = fopen(filename,'w');
fprintf(fileID, '******************************************************************************************************\r\n');
fprintf(fileID, ['study: ', session_info{1}, '\r\n']); 
fprintf(fileID, ['subject: ', session_info{2}, '\r\n']);
fprintf(fileID, ['session: S', session_info{3}, '\r\n']);
fprintf(fileID, ['condition: ', condition, '\r\n']);
fprintf(fileID, ['recorded: ', date, '\r\n']);
fprintf(fileID, '******************************************************************************************************\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, 'DATA ACQUISITION\r\n');
fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '- TMS-EEG-EMG is recorded in 10 mins blocks, first at baseline and then at 8 timepoints (onset each 15 mins)\r\n');
fprintf(fileID, '  following the application of a patch (either capsaicin or vehicle alone) on the hand dorsum.\r\n');
fprintf(fileID, '  Pain intensity ratings are recorded before and after each stimulation block.\r\n');
fprintf(fileID, '- TMS:\r\n');
fprintf(fileID, '     over M1, 120 %%rMT, 80 stimuli per timepoint\r\n');
fprintf(fileID, '     MagVenture MagPro X100 (biphasic sin pulse) + Visor2 neuronavigation (generic brain model)\r\n');
fprintf(fileID, '- EEG: \r\n');
fprintf(fileID, '     Bittium system - TESLA amplifiers, NeurOne recording software, 20kHz sampling rate, 3500Hz LP device filter\r\n');
fprintf(fileID, '     setup with 32 electrodes (10-20), from which M1(29) + M2(30) were used as references, ground at the AFz position\r\n');
fprintf(fileID, '- EMG:\r\n');
fprintf(fileID, '     recorded from contralateral FDI muscle\r\n');
fprintf(fileID, '     MOBI system connected to Visor2, 1024Hz sampling rate\r\n');
fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, ['side of patch application: ', side, ' --> ', hemisphere, ' hemisphere stimulated\r\n']); 
fclose(fileID);
end