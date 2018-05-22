clear all;
clc;

addpath('C:\Users\tgalvin2\Documents\MATLAB\lpxcontrol');

%% Open connection to specrometer controller
% Make sure no connection is open already
out = instrfind();
for ind = 1:length(out)
    if strcmp(out(ind).name,'Serial-COM13')
        fclose(out(ind));
        delete(out(ind));
    end
end

%% Set up connection to spectrometer and ADC board
S = serial('COM13');
fopen(S);
pause(3)

%% initialize   
% Parameters of scan
freq = 20; % unit:Hz; external triggered by the DG535
shots = 3; % # of shots per data point
wvlStart = 4220; % Ang %don't scan in one dot, it will stops when trying to make the axis
wvlEnd = 4390; % Ang
wvlStep = 0.2; % Ang
intTime = 250; % ms
rpm=60; %revolution per minute,speed set for the motor
spa=50; %50 steps per Angstrom
Scantime=(abs(wvlEnd-wvlStart)/wvlStep*(shots/freq*3+(abs(wvlStep*spa)/200/rpm*60)*1.2))/60; % unit: miniute
fprintf('scan time: %f minutes= %f hours\n',Scantime,Scantime/60);

% send the right deriction to move during scan
if wvlStart > wvlEnd
    wvlStep = -wvlStep;
end

% Recall Previous position
load('Wendy_Position');

% Prompt user for starting wavelength
    fprintf('Does the wavelength shown on Spectrometer Wheel match Matlab memory for current position?\n 1) Matlab Value: %f Angstrom is correct\n  2) Enter New Position\n',Matlab_Position);
    the_answer = input('Response: ');
    if the_answer == 1     % Matlab memory is correct
        Stopped_Position = Matlab_Position;
    elseif the_answer == 2
        Stopped_Position = input('Please enter position on spectrometer: ');
        Matlab_Position = Stopped_Position;
    end
    
save('Wendy_Position','Matlab_Position');
wvlCurrent = Stopped_Position; %Stopped_Position: current wavelength position read from wheel of Spectrometer

% wavelength check in case manually set some wrong numbers
wvlDiff=abs(wvlStart-wvlCurrent);
if (wvlDiff > 510)  %unit:Angstrom
	fclose(S);
    delete(S);
    clear S
    error('Wavelength difference > 50 nm ');
end

if (wvlCurrent < 1700 ) || (wvlCurrent > 12000 ) || (wvlStart < 1700 ) || (wvlStart > 12000 ) || (wvlEnd < 1700 ) || (wvlEnd > 12000 ) 
	fclose(S);
    delete(S);
    clear S
	error('Wavelength out of Range: wavlength start set as %f Angstrom; wavelength End set as %f Angstrom',wvlStart,wvlEnd',wvlStart,wvlEnd);
end
    
if (wvlStart == wvlEnd)
	fclose(S);
    delete(S);
    clear S
	error('Start Wavelenth equal to End Wavelength: wavlength start set as %f Angstrom; wavelength End set as %f Angstrom',wvlStart,wvlEnd);
    end

%% Set up Digitizer and move spectrometer to starting line

% Move the spectrometer to the starting position
% in Arduino: Steps are long :Long variables are extended size variables for number storage, and store
%               32 bits (4 bytes), from -2,147,483,648 to 2,147,483,647. 
steps1=(wvlStart-wvlCurrent)*spa;  %50 steps per Angstrom
% fprintf(S,steps); 
arduino_message = sprintf('%i',steps1);
fprintf(S,arduino_message);  %ask the motor to move to the wavelength
timeD1=abs(steps1)/200/rpm*60+5; %wait for motor
pause(timeD1) %Wait for Motor 

%% begin to measure

wvls= wvlStart:wvlStep:wvlEnd;
scanlength = length(wvls);
spec = zeros(2,scanlength); 

figure;
 
axis([min(wvlStart, wvlEnd) max(wvlStart, wvlEnd)  0 1023]);
xlabel('Wavelength (Angstrom)');
ylabel('Intensity');
hold on;


for ind_wvl = 1:scanlength
     % move the motor to the next position
    if ind_wvl ~= 1
    steps2 = wvlStep*spa;
    arduino_message = sprintf('%i',steps2);
    fprintf(S,arduino_message);  %ask the motor to move to the wavelength
    
    timeD2=(abs(steps2)/200/rpm*60)*1.2; %wait for motor
    pause(timeD2) %Wait for Motor 
        
    Matlab_Position = wvlStart+(ind_wvl-1)*wvlStep;
    save('Wendy_Position','Matlab_Position');

    end
    
    % take data here!!!
    pause(shots/freq*3);      % wait for a specified number of shots to accumulate on the Boxcar
    fprintf(S,'PMT');       % ask the ADC to read the Boxcars
    data = fscanf(S);       % read the ADC reply
    
    [PMT1,b] = strtok(data,',');    % format for the reply is "<PMT1 value>,<PMT2 value>:"
    PMT1 = str2num(PMT1);
    PMT2 = str2num(strtok(b(2:length(b)),':'));    % blah
    
    % load the acquired values into the output matrix
    spec(1,ind_wvl) = PMT1;     % PMT 1 (LIF) --> column 1  
    spec(2,ind_wvl) = PMT2;     % PMT 2 (Laser intensity) --> column 2
       
    fprintf('Percentage Completed: %2.1f Ch1: %i Ch2: %i Iteration: %i \n', ...
             ind_wvl/scanlength*100,spec(1,ind_wvl),spec(2,ind_wvl),ind_wvl);

    plot(wvls(ind_wvl),PMT1,'.r');
    plot(wvls(ind_wvl),PMT2,'.g');
    
 
end
wvlCurrent3 = Matlab_Position;                                                              

%% save data
wvlsT=transpose(wvls);
specT=transpose(spec);

str1=datestr(clock);
str1( [15,18] ) = '_';
name=sprintf('YI3_%i_To_%iA_Time%s',wvlStart,wvlEnd,str1);
save(name,'wvlsT','specT');%save data

%% Move to scan begin wavelength, prepare for the second time scanning
% steps3=int16((wvlStart-wvlCurrent3)*spa);  %spa steps per Angstrom
% arduino_message = sprintf('%i',steps3);
% fprintf(S,arduino_message);  %ask the motor to move to the wavelength
% 
% timeD3=abs(steps3)/200/rpm*60+2; %wait for motor
% pause(timeD3) %Wait for Motor 
%Matlab_Position = wvlStart;

%% Save Spectrometer position
% Stopped_Safe = 1;

save('Wendy_Position','Matlab_Position');

%% Close connection to spectrometer
fclose(S);
delete(S);
clear S

%% Plot data
figure;
hSpectrum=plot(wvls,spec(1,:),'.-r',wvls,spec(2,:),'.-g');
xlabel('Wavelength (Angstrom)');
ylabel('Intensity');
title('YI  C')
%saveas(hSpectrum,name,'jpg');%save figure

