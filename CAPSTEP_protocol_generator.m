%% CAPS-TEP stimulation protocol creator
% Author: Dominika
% Description: The following script creates text files compatible with 
% MagPro X100 TMS stimulator. According to the experimental design, each 
% protocol consists of 40 single-pulse biphasic TMS stimuli with random
% ISI <6, 8>s.
%% parameters
target_cortex = 'M1';   
n_protocols = 20;
amplitude = 120;                            % amplitude in percentage of the rMT
ISI = [6000, 7000, 8000; 15, 15, 9];        % matrix of times disponible for ISI
                                            % r1: times (ms); r2: repetitions
max_rep = 3;                                % maximum repetition of an ISI
%% write the experiment file
repetitions = sum(ISI(2,:))+2;
for n = 1:n_protocols
    waiting_times = randomize_isi(ISI, max_rep);
    filename = [target_cortex '_' num2str(n) '.CG3'];
    initializeMagProFilebi(filename,repetitions);
    for i= 1:repetitions
        if i==1
            delay = 10000;
            writeMagPro_singlepulsebi(filename,i,delay,5);
        elseif i == 2
            delay = 10000;
            writeMagPro_singlepulsebi(filename,i,delay,amplitude);
        else
            delay = waiting_times(i-2);
            writeMagPro_singlepulsebi(filename,i,delay,amplitude);
        end        
    end
end



