%% CAPS-TEP: merge blocs into timepoints 
clear all
clc

% parameters
participant = [1, 2];
condition = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
block = {'b1' 'b2'};
prefix = 'ds art-sup dc ep reref chan-select EEG';
channel_n = 30;

%% merge block 1 and 2
% loop through participants
for p = 1:length(participant)
    for c = 1:length(condition)
        % dataset core name
        if participant(p) < 10
            name = [prefix ' 0' num2str(participant(p)) ' ' char(condition(c))];  
        else
            name = [prefix ' ' num2str(participant(p)) ' ' char(condition(c))];  
        end
        
        % loop through timepoints
        for t = 1:length(time)
           % load both datasets
           for b = 1:length(block)
               load([name ' ' char(time(t)) ' ' char(block(b)) '.mat']);
               load([name ' ' char(time(t)) ' ' char(block(b)) '.lw6'], '-mat');    
               statement = ['data_' num2str(b) ' = data; header_' num2str(b) ' = header;'];
               eval(statement);
           end

           % merge the data together
           data = cat(1, data_1, data_2); 

           % merge events and update in header
           T_1 = struct2table(header_1.events); 
           T_2 = struct2table(header_2.events); 
           T = [T_1; T_2];
           n = size(T, 1);
           T.epoch(1:end) = [1:n];
           header.events = table2struct(T);                  
           header.events = header.events';  

           % update the header info
           header.name = [name ' ' char(time(t))];
           if header_1.datasize(2) == channel_n & header_2.datasize(2) == channel_n
               header.datasize = [n, channel_n, 1, 1, 1, 6000];  
           else
               disp('ERROR: uneaven number of channels!')
               return
           end

           % save the merged dataset
           save([name ' ' char(time(t)) '.mat'], 'data');
           save([name ' ' char(time(t)) '.lw6'], 'header');
        end
    end
end
clear p c t b name statement data_1 data_2 data header_1 header_2 header T T_1 T_2 n 