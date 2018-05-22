%Scantime calculation
freq = 20; % unit:Hz; external triggered by the DG535
shots = 3; % # of shots per data point
wvlStart = 4200; % Ang %don't scan in one dot, it will stops when trying to make the axis
wvlEnd = 4400; % Ang
wvlStep = 0.1; % Ang
intTime = 250; % ms
rpm=60; %revolution per minute,speed set for the motor
spa=50;%50 steps per Angstrom
Scantime=(abs(wvlEnd-wvlStart)/wvlStep*(shots/freq*3+(abs(wvlStep*spa)/200/rpm*60)*1.2))/60; % unit: miniute
fprintf('scan time: %f minutes= %f hours\n',Scantime,Scantime/60);