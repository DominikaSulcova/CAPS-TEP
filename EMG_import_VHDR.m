function [out_header,out_data]=EMG_import_VHDR(filename);
%RLW_import_VHDR
%
%Import Brainvision VHDR
%filename : name of VHDR file 
%
% Author : 
% Andre Mouraux
% Institute of Neurosciences (IONS)
% Universite catholique de louvain (UCL)
% Belgium
% 
% modified by Dominika (13/08/2020)
% 
% Contact : andre.mouraux@uclouvain.be
% This function is part of Letswave 6
% See http://nocions.webnode.com/letswave for additional information
%

out_header=[];
out_data=[];

%load the VHDR file
%load data
dat=ft_read_data(filename);
%load header
hdr=ft_read_header(filename);
%load events
trg=ft_read_event(filename);

%set header
out_header.filetype='time_amplitude';
out_header.name=filename;
out_header.tags='';
out_header.history(1).configuration=[];
out_header.datasize=[hdr.nTrials hdr.nChans 1 1 1 hdr.nSamples];
out_header.xstart=(hdr.nSamplesPre/hdr.Fs)*-1;
out_header.ystart=0;
out_header.zstart=0;
out_header.xstep=1/hdr.Fs;
out_header.ystep=1;
out_header.zstep=1;

%dummy chanloc
chanloc.labels='';
chanloc.topo_enabled=0;
chanloc.SEEG_enabled=0;
%set chanlocs
for chanpos=1:hdr.nChans;
    chanloc.labels=hdr.label{chanpos};
    out_header.chanlocs(chanpos)=chanloc;
end;

%set events
numevents=size(trg,2);
%set events
if numevents==0;
    out_header.events=[];
else
    for eventpos=1:numevents;
        event.code='unknown';
        if isempty(trg(eventpos).value);
            event.code=trg(eventpos).type;
        else
            event.code=trg(eventpos).value;
        end;
        if isnumeric(event.code);
            event.code=num2str(event.code);
        end;
        event.latency=(trg(eventpos).sample*out_header.xstep)+out_header.xstart;
        event.epoch=1;
        out_header.events(eventpos)=event;
    end;
end;

%set data
out_data=zeros(out_header.datasize);
for chanpos=1:out_header.datasize(2);
    for epochpos=1:out_header.datasize(1);
        out_data(epochpos,chanpos,1,1,1,:)=squeeze(dat(chanpos,:,epochpos));
    end;
end;

