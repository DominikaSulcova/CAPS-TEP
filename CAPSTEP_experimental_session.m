%% CAPS-TEP experimental session script
% Author: Dominika
% Description: The script runs along the CAPS-TEP experimental session
% It allows to: 1) initialize a logfile for current participant and fill in
%                  relevant information
%               2) time individual steps of the sessiom (starting at T0) 
%               3) write the intensity ratings into the appropriate line of 
%                  an existing MATLAB table
%% clean out, set parameters
clear all; clc;
timepoints = [0 10 15 25 30 40 45 55 60 70 75 85 90 100];
%% enter subject and experiment information
% enter the info
prompt = {'Study:', 'Subject number:', 'Session:'};
dlgtitle = 'Input information';
dims = [1 35];
definput = {'CAPS-TEP', '', '1'};
session_info = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput
% condition
answer = questdlg('Which solution is applied?', 'Experimental condition',...
    'capsaicin', 'vehicle', 'none');
switch answer
    case 'capsaicin'
        session_info{4} = '1';
    case 'vehicle'
        session_info{4} = '0';
end
% side of patch application
answer = questdlg('To which side is the patch applied?', 'Treated side',...
    'right', 'left', 'none');
switch answer
    case 'right'
        session_info{5} = '1';
    case 'left'
        session_info{5} = '0';
end
clear answer
%% initialize output files
% logfile
filename = initialize_logfile(session_info);
% load the table of intensity ratings
load('intensity_ratings.mat');
% fill in the condition
index_row = find([intensity_ratings{2:end,1}] == 0);
index_row = index_row(1)+1;
intensity_ratings{index_row, 1} = str2num(session_info{2});
intensity_ratings{index_row, 2} = str2num(session_info{3});
intensity_ratings{index_row, 3} = str2num(session_info{4});
%% baseline rMT and closest electrodes
% rMT
session_info{6} = inputdlg('What is the rMT at the beginning of the sesion?', 'rMT');
% register electrodes 
prompt = {'Electrode 1:', 'Electrode 2:', 'Electrode 3:'};
dlgtitle = 'Hotspot - closest electrodes';
dims = [1 50];
definput = {'', '', ''};
electrodes = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput
% associate numbers with electrode names
load('labels.mat')
labels{31} = 'C1'; labels{32} = 'C2'; 
for e = 1:length(electrodes)
    electrodes{e} = [electrodes{e} ' (' labels{str2double(electrodes{e})} ')'];
end
%% create the timer
% parameters
blocks = length(timepoints) / 2;
event_time = 60;
block_time = 60*15;
total_time = blocks * block_time;
timer_count = 1;

% create a timer
block_timer = timer;
block_timer.UserData = {timer_count};
block_timer.StartFcn = @start_function;
block_timer.TimerFcn = @timer_function;
block_timer.StopFcn = @timer_cleanup;
block_timer.Period = event_time;
block_timer.StartDelay = event_time;
block_timer.TasksToExecute = (blocks * block_time)/event_time - 1;
block_timer.ExecutionMode = 'fixedSpacing';
%% experimental session
% press SPACE to start timing the session 
disp('Press SPACE to start timing the session.');
ok=0; 
while ok==0;    
    WaitSecs(1); 
    [onset_time,keycode,deltasecs]=KbWait([]);   
    [a,b]=max(keycode); 
    if b==32; 
        start_session(filename); 
        ok=1;     
    end; 
end;    
clear a b keycode ok onset_time deltasecs 

% timed loop 
counter = 1; 
start(block_timer);
for i = 1: length(timepoints)
    intensity_ratings = rate(intensity_ratings, index_row, counter, timepoints);
    counter = counter + 1;
end

% stop(block_timer);
% get(block_timer, 'Running');
% delete(block_timer);
%% verify rMT 
session_info{7} = inputdlg('What is the rMT at the end of the sesion?', 'rMT');
fileID = fopen(filename,'a');
fprintf(fileID, '\r\n');
fprintf(fileID, ['rMT (beginning): ' char(session_info{6}) ' %%MSO\r\n']);
fprintf(fileID, ['rMT (end): ' char(session_info{7}) ' %%MSO\r\n']);
fprintf(fileID, '\r\n');
fprintf(fileID, ['M1 closest electrodes: ' electrodes{1} ' - ' electrodes{2} ' - ' electrodes{3} '\r\n']);
fprintf(fileID, '\r\n');
fprintf(fileID, 'intensity ratings saved to MATLAB cell array --> intensity_ratings.mat\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, 'notes:\r\n');
fclose(fileID);

% save the ratings
save('intensity_ratings.mat','intensity_ratings');







