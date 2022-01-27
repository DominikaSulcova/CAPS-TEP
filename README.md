# CAPS-TEP
Project for fun :)

------------------------------------------------------------
Experimental protocol
------------------------------------------------------------
15 (16) young healthy volunteers
Datasets for each subject: corresponding EEG and EMG epochs
- 2 sessions --> topical capsaicin (2% solution) OR vehicle (ethanol) applied in a patch over FDI muscle
- 7 timepoint per session: baseline + 6 recordings every 15 mins following the patch application

TMS:
- MagVenture MagPro X100, biphasic sin pulse
- 80 stimuli per timepoint (split in 2 blocks per 40 stims)
- over M1 - FDI hotspot
- single pulse TMS of 120 %rMT 
- ISI randomly 6 - 8 s

EEG:
- NeurOne EEG system (MEGA/Bittium)
- 32 (30) electrodes – referenced to mastoids (ground an extra electrode on the front)
- SR 20 kHz

EMG:
- Visor2 MobiEMG amplifier
- SR 1024 Hz

------------------------------------------------------------
AVAILABLE SCRIPTS & FUNCTIONS
------------------------------------------------------------
Experimental session
------------------------------------------------------------
- CAPSTEP_protocol_generator.m
- CAPSTEP_experimental_session.m
- CAPSTEP_initialize_logfile.m


EEG --> TEPs 
------------------------------------------------------------
**EEG data preprocessing:**
- CAPSTEP_EEG_import.m
- EEG_import_MEGA.m
- EEG_history_import.mat
- CAPSTEP_EEG_preprocess.lwscript
- CAPSTEP_EEG_merge.m
- CAPSTEP_EEG_ica_timecourse.lwscript
- CAPSTEP_EEG_ica_FFT.lwscript
- CAPSTEP_EEG_filter.lwscript

**TEP analysis:**
- CAPSTEP_TEP_process.m


EMG -->  MEPs
------------------------------------------------------------
**EMG data preprocessing:**
- CAPSTEP_EMG_import.m
- EMG_import_VHDR.m
- EMG_history_import.mat

**MEP analysis:**
- CAPSTEP_MEP_process.m
