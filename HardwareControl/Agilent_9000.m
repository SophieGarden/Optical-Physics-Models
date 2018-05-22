% Agilent_9000.m
% Stand-alone program template for taking data with Agilent 9000 Oscilloscope
% Tom Galvin - University of Illinois at Urbana-Champaign - 2016
% v1.0


%% Meter and Scan Parameters
% If you wanted to turn this into a function, make these the inputs
scope_address = 'USB0::0x0957::0x900A::MY52370145::INSTR'; % VISA address of oscilloscope
NUMTRIGS = 20;                                      % Number of data points desired
REPRATE  = 1;                                       % Repetition rate of trigger


%% Open Connection to Oscilloscope
% Create visa object of oscilloscope
visaObj = visa('ni',scope_address);
visaObj.InputBufferSize = 50e6;     % Set the buffer size
visaObj.Timeout = 2;                % Set the timeout value
visaObj.ByteOrder = 'littleEndian'; % Set the Byte order
fopen(visaObj);                     % Open the connection


%% Configure channels on oscilloscope
% Refer to Keysight Infiniium Programming Guide for command reference

% fprintf(visaObj,'*RST');                       % Reset the instrument
fprintf(visaObj,':STOP');                        % Stop any acquisition which may be in progress
fprintf(visaObj,'*CLS');                         % Clear the output buffer

%% Set up Trigger
% Set trigger source to AUX and rising edge triggering
fprintf(visaObj,':TRIGger:EDGE:SOUR AUX');
fprintf(visaObj,':TRIGger:EDGE:SLOPe POS');
fprintf(visaObj,':TRIGger:LEVel AUX, 2');
fprintf(visaObj,':TRIGger:SWEep TRIGgered');


%% Set up channels
% Change input impedance to 1 MOhm and coupling to DC for all four channels
fprintf(visaObj,':CHANnel1:INPut DC');
fprintf(visaObj,':CHANnel2:INPut DC');
fprintf(visaObj,':CHANnel3:INPut DC');
fprintf(visaObj,':CHANnel4:INPut DC');

% Not Using Probes, so set attenuation to 1
fprintf(visaObj,':CHANnel1:PROBe 1,RAT');
fprintf(visaObj,':CHANnel2:PROBe 1,RAT');
fprintf(visaObj,':CHANnel3:PROBe 1,RAT');
fprintf(visaObj,':CHANnel4:PROBe 1,RAT');

% Turn off bandwidth limit
fprintf(visaObj,':CHANnel1:BWLimit OFF');
fprintf(visaObj,':CHANnel2:BWLimit OFF');
fprintf(visaObj,':CHANnel3:BWLimit OFF');
fprintf(visaObj,':CHANnel4:BWLimit OFF');

% Turn on relevant channels
fprintf(visaObj,':CHANnel1:DISPlay OFF');
fprintf(visaObj,':CHANnel2:DISPlay ON');
fprintf(visaObj,':CHANnel3:DISPlay ON');
fprintf(visaObj,':CHANnel4:DISPlay OFF');

% COMMENTED - user setss these parameters with buttons on oscilloscope front panel
% % Set Channel offsets                        
% fprintf(visaObj,':CHANnel1:OFFSet 1');
% fprintf(visaObj,':CHANnel2:OFFSet  1');
% fprintf(visaObj,':CHANnel3:OFFSet  .1');
% fprintf(visaObj,':CHANnel4:OFFSet  1');

% % Set Channel Scales
% fprintf(visaObj,':CHANnel1:SCALe 1');
% fprintf(visaObj,':CHANnel2:SCALe .5');
% fprintf(visaObj,':CHANnel3:SCALe .05');
% fprintf(visaObj,':CHANnel4:SCALe 1');



%% Configure Acquisition Subsytem and Waveform
fprintf(visaObj,':WAVEFORM:VIEW MAIN');                % Only record points visible in window
fprintf(visaObj,':SYSTem:HEADer OFF');                 % Don't get header data when reading from oscilloscope
fprintf('\n\n\n Reading in as ascii...\n')
fprintf(visaObj,':ACQuire:POINts AUTO');               % Automatically match the number of points in the window
fprintf(visaObj,':ACQuire:INTerpolate OFF');           % Turn off interpolation... Reduced number of data points
fprintf(visaObj,':WAVeform:FORMat ASCii');             % Return data in ASCII format

% Configure segment mode -- In this mode the part of the waveform visible
% in the window is recorded for the specified number of triggers
% (NUMTRIGS), then data can be read.
fprintf(visaObj,sprintf(':ACQuire:SEGMented:COUNt %i',NUMTRIGS)); % NUMTRIGS number of waveforms recorded
fprintf(visaObj,':WAVeform:SEGMented:ALL ON');                    % Download all segments at once
fprintf(visaObj,':ACQuire:MODE SEGMented');                       % Turn on SEGMENTED mode

% Acquire the data
fprintf(visaObj,':DIGitize CHANnel2,CHANnel3');   % Capture data from channeles 2 and 3
visaObj.Timeout = NUMTRIGS/REPRATE+5;             % Don't time out until scope has had chance to acquire data


% Read in data from Channel 2
fprintf(visaObj,':WAVeform:SOURce CHANnel2');
fprintf(visaObj,':WAVeform:DATA?');
raw = str2num(fscanf(visaObj,'%s'));
chan2 = reshape(raw,length(raw)/NUMTRIGS,NUMTRIGS);  % Reshape data

% Read in Data from Channel 3
fprintf(visaObj,':WAVeform:SOURce CHANnel3');
fprintf(visaObj,':WAVeform:DATA?');
raw = str2num(fscanf(visaObj,'%s'));
chan3 = reshape(raw,length(raw)/NUMTRIGS,NUMTRIGS);   % Reshape data

% Create timebase
t_window = str2double(query(visaObj,':TIMebase:WINDow:RANGe?'));    % Time range in visible window on scope
Dt = t_window/(size(chan3,1)-1);                                    % Spacing between data samples
t_center = str2double(query(visaObj,':TIMebase:WINDow:POSition?')); % Time at center 
t_base = Dt*((1:size(chan3,1))-(size(chan3,1)-1)/2) + t_center;     % Time at every sample relative to trigger


%% Reset display so oscilloscope is running again
% Set trigger source to AUX and rising edge triggering
fprintf(visaObj,':ACQuire:MODE RTIME');       % Acquire all data on one trigger
fprintf(visaObj,':TRIGger:SWEep TRIGgered');
fprintf(visaObj,':RUN');

% Turn on relevant channels
fprintf(visaObj,':CHANnel1:DISPlay OFF');
fprintf(visaObj,':CHANnel2:DISPlay ON');
fprintf(visaObj,':CHANnel3:DISPlay ON');
fprintf(visaObj,':CHANnel4:DISPlay OFF');



%% Close Connection and delete visa object
% query(visaObj,':SYSTEM:ERR? STR')
fclose(visaObj);
delete(visaObj);
clear visaObj;


%% Plot Captured Waveforms
plot(1e9*t_base,chan2)
hold on
plot(1e9*t_base,chan3)
hold off