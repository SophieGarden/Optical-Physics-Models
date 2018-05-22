% EnergyMAX.m
% Stand-alone program template for taking data with Coherent EnergyMAX meters
% Tom Galvin - University of Illinois at Urbana-Champaign - 2016
% v1.0


%% Meter and Scan Parameters
% If you wanted to turn this into a function, make these the inputs
METERS.address = {'ASRL24::INSTR','ASRL23::INSTR'};  % VISA address of meter
METERS.wavelength = [532 660];                       % nm -- wavelength of meters
METERS.range = {'MIN','MIN'};                        % range of meters 'MAX' or 'MIN'
NUMTRIGS = 200;                                      % Number of data points desired
REPRATE  = 10;                                       % Repetition rate of trigger




%% Open connections to each power meter and configure settings
for ind_meter = 1:length(METERS.address)
    
    % Set up connection for power meter
    pm(ind_meter) = visa('ni',METERS.address{ind_meter});
    set(pm(ind_meter),'InputBufferSize',max(1e5,NUMTRIGS*20+1e4));
    set(pm(ind_meter),'Timeout',2);
    set(pm(ind_meter),'BaudRate',9600);
    set(pm(ind_meter),'DataBits',8);
    set(pm(ind_meter),'StopBits',1);
    set(pm(ind_meter),'Parity','none');
    set(pm(ind_meter),'Terminator',13);               % USE THE LINE FEED as terminator
    set(pm(ind_meter),'FlowControl','none');
    
    % Open the connection
    fopen(pm(ind_meter));
    
    % Make sure meter is not running and clear its buffer
    fprintf(pm(ind_meter),'ABORt');                    % Stop streaming data collection
    if pm(ind_meter).BytesAvailable > 0                % Clear buffer, if necessary
        fread(pm(ind_meter),pm(ind_meter).BytesAvailable,'char');
    end
    fprintf(pm(ind_meter),'SYSTem:ERRor:CLEar');      % Clear any existing errors
    
    % Configure Power Meter
    fprintf(pm(ind_meter),'CONFigure:MEASure J');      % Set measurement to Joules
    fprintf(pm(ind_meter),sprintf('CONFigure:WAVElength %i'  ,METERS.wavelength(ind_meter))); % Set wavelength of this power meter to 532 nm
    fprintf(pm(ind_meter),sprintf('CONFigure:RANGe:SELect %s',METERS.range{ind_meter}));
    fprintf(pm(ind_meter),'TRIGger:SOURce EXT');       % Set
    fprintf(pm(ind_meter),'TRIGger:LEVel 2');
    fprintf(pm(ind_meter),'TRIGger:SLOPe POSitive');
    fprintf(pm(ind_meter),'CONFigure:ITEMselect PULS,SEQ');
    
    % Print Model Numbers
    fprintf('Meter %i: %s',ind_meter,query(pm(ind_meter),'SYSTem:INFormation:MODel?'));  
end

%% Start power meters taking data
for ind_meter = 1:length(METERS.address)
    fprintf(pm(ind_meter),'TRIGger:SEQuence 0');       % Set the sequence ID to zero for both channels
    fprintf(pm(ind_meter),'INIT');                     % Begin streaming data collection
end

%% Wait while data is acquired
pause(NUMTRIGS/REPRATE+1);
% Stop streaming data collection
for ind_meter = 1:length(METERS.address)
    fprintf(pm(ind_meter),'ABORt');  
end

%% Read data from power meters and close the connections
power_meter_data = zeros(NUMTRIGS,2,length(METERS.wavelength));
for ind_meter = 1:length(METERS.address)
    % Read in entire buffer
    raw_data = fread(pm(ind_meter),pm(ind_meter).BytesAvailable,'char');
    raw_data = raw_data - 128;                              % Subtract parity bit added by Coherent
    data1 = str2num(char(raw_data'));                       % Convert data to array
    power_meter_data(:,:,ind_meter) = data1(1:NUMTRIGS,:);  % Save data in output array
    
    % Close connections and delete the references to them
    fclose(pm(ind_meter))
    delete(pm(ind_meter))
    clear pm(ind_meter)
end

%% Plot results
plot( power_meter_data(:,2,1),power_meter_data(:,1,1), power_meter_data(:,2,2),power_meter_data(:,1,2));