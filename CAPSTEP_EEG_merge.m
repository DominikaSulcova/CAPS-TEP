%% CAPS-TEP: merge blocs into timepoints 
clear all
clc

% parameters
participant = [12];
treatment = {'caps'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
blocks = {'b1' 'b2'};
prefix = 'reref ds art-sup ep dc EEG';
%% merge 
% loop through participants
for p = 1:length(participant)
    for a = 1:length(treatment)

        name = [prefix ' ' num2str(participant(p)) ' ' char(treatment(a))];  

        for b = 1:length(time)       
           % load both datasets
           for c = 1: length(blocks)
               load([name ' ' char(time(b)) ' ' char(blocks(c)) '.mat']);
               load([name ' ' char(time(b)) ' ' char(blocks(c)) '.lw6'], '-mat');    
               command = ['data_' num2str(c) ' = data; header_' num2str(c) ' = header;'];
               eval(command);
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
           header.name = [name ' ' char(time(b))];
           header.datasize = [n,32,1,1,1,6000];   

           % save the merged dataset
           save([name ' ' char(time(b)) '.mat'], 'data');
           save([name ' ' char(time(b)) '.lw6'], 'header');
        end
    end
end