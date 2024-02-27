clc
clear all
close all

%% ASSIGNEMENT 4

% Creating a scenario:
startTime = datetime(2023,01,01,00,00,0);
stopTime = startTime + days(7);
sampleTime = 60;
sc = satelliteScenario(startTime,stopTime,sampleTime);

% Creating a satellite:
EOSat.semiMajorAxis = earthRadius + [500e3]; % [m]
EOSat.inclination = [97.4]; % [deg]
EOSat.eccentricity = [0];
EOSat.rightAscensionOfAscendingNode = [0];
EOSat.argumentOfPeriapsis = [0];
EOSat.trueAnomaly = [0];
EOSat.name = {'EO Sat'};
EOSat.satellite = satellite(sc, EOSat.semiMajorAxis, EOSat.eccentricity, EOSat.inclination, EOSat.rightAscensionOfAscendingNode, EOSat.argumentOfPeriapsis, EOSat.trueAnomaly ,"Name", "EOSat");
%show(EOSat.satellite)

%%

% SNR to EbNo:
Rc=9/10; %code rate
m=2;
Rb_target=375;  
Rb_gross=Rb_target/Rc;
Rs=Rb_gross/m; % in Msps
roll_off = 0.2;
B=(1+roll_off)*Rs; %equal to Rb_net
Rb_net=B;
Spectral_Efficiency=1.788612; 
Es_No_ideal = 6.42; 
link_margin = 3;

RequiredEbNo = Es_No_ideal - 10*log10(Spectral_Efficiency) + 10*log10(1+roll_off) + link_margin;

%Adding a gimbal (to steer the antenna):
EOSat.gimbal = gimbal(EOSat.satellite);

%Pointing a gimbal (or the satellite) towards another object:
%pointAt(Source, Target)

% Adding a transmitter
tx = transmitter(EOSat.gimbal,"Name", "Satellite Transmitter", 'Frequency', 8.2e9 ,"Power", 3,"SystemLoss", 1, "BitRate", Rb_net); % Bitrate in our case is the bandwidth, in MHz

% Adding a Gaussian antenna (approximation to speedup simulation, use parabolic reflector gain to determine parameters)
EOSat.antenna = gaussianAntenna(tx, "ApertureEfficiency", 0.55,"DishDiameter", 0.197); % parameters calculated through the parabolic reflector gain formula 


% Creating a Ground Station
gs1 = groundStation(sc, "Name","Inuvik", "Latitude",68.35, "Longitude",-133.72,"Altitude",15, "MinElevationAngle",10);
gs2 = groundStation(sc, "Name","Svalbard", "Latitude",78.22, "Longitude",15.38,"Altitude",440, "MinElevationAngle",10);
gs3 = groundStation(sc, "Name","Awarua", "Latitude",-46.52, "Longitude",168.48,"Altitude",0, "MinElevationAngle",10);
gs4 = groundStation(sc, "Name","Troll", "Latitude",-72.01, "Longitude",2.53,"Altitude",1200, "MinElevationAngle",10);

% Adding a Gimbal (on which the antennas will be mounted, mounting angles to point it towards sky)
gimbal(gs1, "MountingAngles", [0; 180; 0]);
pointAt(gs1.Gimbals,EOSat.satellite) % pointing the antenna towards the satellite
gimbal(gs2, "MountingAngles", [0; 180; 0]);
pointAt(gs2.Gimbals,EOSat.satellite)
gimbal(gs3, "MountingAngles", [0; 180; 0]);
pointAt(gs3.Gimbals,EOSat.satellite)
gimbal(gs4, "MountingAngles", [0; 180; 0]);
pointAt(gs4.Gimbals,EOSat.satellite)

G_rx_clear = 48.17; % dB calculated with the formula
T_sys_clear = 141.50 + 50; %system temperature

Gain_Temp = G_rx_clear - 10*log10(T_sys_clear);

% Adding a receiver
gs1_rx = receiver(gs1.Gimbals, "Name", "Ground Station 1 Receiver" ,"GainToNoiseTemperatureRatio", Gain_Temp,"SystemLoss", 1,"RequiredEbNo", RequiredEbNo,"PreReceiverLoss", 0);
gaussianAntenna(gs1_rx, "ApertureEfficiency", 0.65,"DishDiameter", 3.7)
gs2_rx = receiver(gs2.Gimbals, "Name", "Ground Station 2 Receiver" ,"GainToNoiseTemperatureRatio", Gain_Temp,"SystemLoss", 1,"RequiredEbNo", RequiredEbNo,"PreReceiverLoss", 0);
gaussianAntenna(gs2_rx, "ApertureEfficiency", 0.65,"DishDiameter", 3.7)
gs3_rx = receiver(gs3.Gimbals, "Name", "Ground Station 3 Receiver" ,"GainToNoiseTemperatureRatio", Gain_Temp,"SystemLoss", 1,"RequiredEbNo", RequiredEbNo,"PreReceiverLoss", 0);
gaussianAntenna(gs3_rx, "ApertureEfficiency", 0.65,"DishDiameter", 3.7)
gs4_rx = receiver(gs4.Gimbals, "Name", "Ground Station 4 Receiver" ,"GainToNoiseTemperatureRatio", Gain_Temp,"SystemLoss", 1,"RequiredEbNo", RequiredEbNo,"PreReceiverLoss", 0);
gaussianAntenna(gs4_rx, "ApertureEfficiency", 0.65,"DishDiameter", 3.7)

%% PART A

% Creating a link: 
pointAt(EOSat.satellite,gs1);
lnk_1 = link(tx,gs1_rx);
ac_1 = access(EOSat.satellite,gs1);
intervals_1 = linkIntervals(lnk_1);
[SNR_1, Time_1] = ebno(lnk_1);

RealEbNo_1 = SNR_1 - 10*log10(Spectral_Efficiency) + 10*log10(1+roll_off);

f=figure; plot(Time_1,RealEbNo_1,LineWidth=2);
title("SNR vs. Time, Inuvik");
ylim([0 35]);
xlabel("Time");
ylabel("SNR (dB)");
grid on;

PLcfgP618_1 = p618Config('Frequency', 8.2e9,'ElevationAngle', 10,'Latitude',68.35 ,'Longitude', -133.72,'TotalAnnualExceedance', 0.1,'AntennaDiameter',gs1_rx.Antenna.DishDiameter,'AntennaEfficiency',gs1_rx.Antenna.ApertureEfficiency);
[PL_1, ~, TSKY_1] = p618PropagationLosses(PLcfgP618_1);
SNR_1_at = RealEbNo_1 - PL_1.At;

hold on, plot(Time_1, SNR_1_at, 'r', LineWidth=2), yline(3,'o',LineWidth=2)
saveas(figure(1),"SNR_first_station",'epsc')

pointAt(EOSat.satellite,gs2);
lnk_2 = link(tx,gs2_rx);
ac_2 = access(EOSat.satellite,gs2);
intervals_2 = linkIntervals(lnk_2);
[SNR_2, Time_2] = ebno(lnk_2);

RealEbNo_2 = SNR_2 - 10*log10(Spectral_Efficiency) + 10*log10(1+roll_off);

f=figure; plot(Time_2,RealEbNo_2,LineWidth=2);
title("SNR vs. Time, Svalbard");
ylim([0 35]);
xlabel("Time");
ylabel("SNR (dB)");
grid on;

PLcfgP618_2 = p618Config('Frequency', 8.2e9,'ElevationAngle', 10,'Latitude',78.22 ,'Longitude', 15.38,'TotalAnnualExceedance', 0.1,'AntennaDiameter',gs2_rx.Antenna.DishDiameter,'AntennaEfficiency',gs2_rx.Antenna.ApertureEfficiency);
[PL_2, ~, TSKY_2] = p618PropagationLosses(PLcfgP618_2);
SNR_2_at = RealEbNo_2 - PL_2.At; % En_No considering rainy condition 

hold on, plot(Time_2, SNR_2_at, 'r', LineWidth=2), yline(3,'o',LineWidth=2)
saveas(figure(2),"SNR_second_station",'epsc')

pointAt(EOSat.satellite,gs3);
lnk_3 = link(tx,gs3_rx);
ac_3 = access(EOSat.satellite,gs3);
intervals_3 = linkIntervals(lnk_3);
[SNR_3, Time_3] = ebno(lnk_3);

RealEbNo_3 = SNR_3 - 10*log10(Spectral_Efficiency) + 10*log10(1+roll_off);

f=figure; plot(Time_3,RealEbNo_3,LineWidth=2)
title("SNR vs. Time, Awarua");
ylim([0 35]);
xlabel("Time");
ylabel("SNR (dB)");
grid on;

PLcfgP618_3 = p618Config('Frequency', 8.2e9,'ElevationAngle', 10,'Latitude',-46.52,'Longitude', 168.48,'TotalAnnualExceedance', 0.1,'AntennaDiameter',gs3_rx.Antenna.DishDiameter,'AntennaEfficiency',gs3_rx.Antenna.ApertureEfficiency);
[PL_3, ~, TSKY_3] = p618PropagationLosses(PLcfgP618_3);
SNR_3_at = RealEbNo_3 - PL_3.At;

hold on, plot(Time_3, SNR_3_at, 'r', LineWidth=2), yline(3,'o',LineWidth=2)
saveas(figure(3),"SNR_third_station",'epsc')

pointAt(EOSat.satellite,gs4);
lnk_4 = link(tx,gs4_rx);
ac_4 = access(EOSat.satellite,gs4);
intervals_4 = linkIntervals(lnk_4); 
[SNR_4, Time_4] = ebno(lnk_4);

RealEbNo_4 = SNR_4 - 10*log10(Spectral_Efficiency) + 10*log10(1+roll_off);

f=figure; plot(Time_4,RealEbNo_4,LineWidth=2);
title("SNR vs. Time, Troll");
ylim([0 35]);
xlabel("Time");
ylabel("SNR (dB)");
grid on;

PLcfgP618_4 = p618Config('Frequency', 8.2e9,'ElevationAngle', 10,'Latitude',-72.01 ,'Longitude', 2.53,'TotalAnnualExceedance', 0.1,'AntennaDiameter',gs4_rx.Antenna.DishDiameter,'AntennaEfficiency',gs4_rx.Antenna.ApertureEfficiency);
[PL_4, ~, TSKY_4] = p618PropagationLosses(PLcfgP618_4);
SNR_4_at = RealEbNo_4 - PL_4.At;

hold on, plot(Time_4, SNR_4_at, 'r', LineWidth=2), yline(3,'o',LineWidth=2)
saveas(figure(4),"SNR_fourth_station",'epsc')

%timetable with only numerical value 
t_1 = table2timetable(intervals_1(:,[3:4, 6:8]));
t_2 = table2timetable(intervals_2(:,[3:4, 6:8]));
t_3 = table2timetable(intervals_3(:,[3:4, 6:8]));
t_4 = table2timetable(intervals_4(:,[3:4, 6:8]));

%complete timetables
t_1_c = table2timetable(intervals_1); 
t_2_c = table2timetable(intervals_2);
t_3_c = table2timetable(intervals_3);
t_4_c = table2timetable(intervals_4);

table_1 = retime(t_1,'daily','mean'); %mean
table_1_2 = retime(t_1,'daily','count'); %count
table_1_3 = retime(t_1,'daily','sum'); %sum

% % Helper function to combine the tables of different stations
table_2 = retime(t_2,'daily','mean');
table_2_2 = retime(t_2,'daily','count');
table_2_3 = retime(t_2,'daily','sum');

table_3 = retime(t_3,'daily','mean');
table_3_2 = retime(t_3,'daily','count');
table_3_3 = retime(t_3,'daily','sum');

table_4 = retime(t_4,'daily','mean');
table_4_2 = retime(t_4,'daily','count');
table_4_3 = retime(t_4,'daily','sum');

tot_table_access= sortrows(vertcat(t_1_c, t_2_c, t_3_c, t_4_c));

tot_table_mean= sortrows(vertcat(table_1, table_2, table_3, table_4));
tot_table_count= sortrows(vertcat(table_1_2, table_2_2, table_3_2, table_4_2));
tot_table_sum= sortrows(vertcat(table_1_3, table_2_3, table_3_3, table_4_3));

MDC_1 = sum(table_1_3.Duration * Rb_net * 1e-6)/7; % mean daily capacity 
MDC_2 = sum(table_2_3.Duration * Rb_net * 1e-6)/7; % mean daily capacity 
MDC_3 = sum(table_3_3.Duration * Rb_net * 1e-6)/7; % mean daily capacity 
MDC_4 = sum(table_4_3.Duration * Rb_net * 1e-6)/7; % mean daily capacity 

MDC_EN = MDC_1 + MDC_2 + MDC_3 + MDC_4; %mean daily capacity for the entire network 

MDA_1 = sum(table_1_2.IntervalNumber)/7; %mean daily access
MDA_2 = sum(table_2_2.IntervalNumber)/7; %mean daily access
MDA_3 = sum(table_3_2.IntervalNumber)/7; %mean daily access
MDA_4 = sum(table_4_2.IntervalNumber)/7; %mean daily access

MDA_EN = MDA_1 + MDA_2 + MDA_3 + MDA_4; %mean daily access for the entire network 

AV_1 = sum(table_1_3.Duration)/(7*24*60*60)*100; % Availability
AV_2 = sum(table_2_3.Duration)/(7*24*60*60)*100; % Availability
AV_3 = sum(table_3_3.Duration)/(7*24*60*60)*100; % Availability
AV_4 = sum(table_4_3.Duration)/(7*24*60*60)*100; % Availability

AV_EN = AV_1 + AV_2 + AV_3 + AV_4; %mean daily access for the entire network 

% latency for station 1

lat_1 = [];

for i=1:(height(t_1_c)-1)
    latency = t_1_c.StartTime(i+1) - t_1_c.EndTime(i);
    lat_1 = [lat_1, latency];
end

lat_1_max = minutes(max(lat_1));
lat_1_mean = minutes(mean(lat_1));

% latency for station 2

lat_2 = [];

for i=1:(height(t_2_c)-1)
    latency = t_2_c.StartTime(i+1) - t_2_c.EndTime(i);
    lat_2 = [lat_2, latency];
end

lat_2_max = minutes(max(lat_2));
lat_2_mean = minutes(mean(lat_2));

% latency for station 3

lat_3 = [];

for i=1:(height(t_3_c)-1)
    latency = t_3_c.StartTime(i+1) - t_3_c.EndTime(i);
    lat_3 = [lat_3, latency];
end

lat_3_max = minutes(max(lat_3));
lat_3_mean = minutes(mean(lat_3));


% latency for station 4

lat_4 = [];

for i=1:(height(t_4_c)-1)
    latency = t_4_c.StartTime(i+1) - t_4_c.EndTime(i);
    lat_4 = [lat_4, latency];
end

lat_4_max = minutes(max(lat_4));
lat_4_mean = minutes(mean(lat_4));


% latency for the entire network

lat_EN = [];

for i=1:(height(tot_table_access)-1)
    latency = tot_table_access.StartTime(i+1) - tot_table_access.EndTime(i);
    lat_EN = [lat_EN, latency];
end

lat_EN_max = minutes(max(lat_EN));
lat_EN_mean = minutes(mean(lat_EN));

%% rain 

T_sys_r1 = 141.5 + TSKY_1; % new temperature with rain 
GTNTR_1 = G_rx_clear - 10*log10(T_sys_r1); %gain to noise temperature ratio 1 with rain

T_sys_r2 = 141.5 + TSKY_2; % new temperature with rain 
GTNTR_2 = G_rx_clear - 10*log10(T_sys_r2); %gain to noise temperature ratio 1 with rain

T_sys_r3 = 141.5 + TSKY_3; % new temperature with rain 
GTNTR_3 = G_rx_clear - 10*log10(T_sys_r3); %gain to noise temperature ratio 1 with rain

T_sys_r4 = 141.5 + TSKY_4; % new temperature with rain 
GTNTR_4 = G_rx_clear - 10*log10(T_sys_r4); %gain to noise temperature ratio 1 with rain

gs1_rx_r = receiver(gs1.Gimbals, "Name", "Ground Station 1 Receiver" ,"GainToNoiseTemperatureRatio", GTNTR_1,"SystemLoss", 1,"RequiredEbNo", RequiredEbNo,"PreReceiverLoss", 0);
gaussianAntenna(gs1_rx_r, "ApertureEfficiency", 0.65,"DishDiameter", 3.7)
gs2_rx_r = receiver(gs2.Gimbals, "Name", "Ground Station 2 Receiver" ,"GainToNoiseTemperatureRatio", GTNTR_2,"SystemLoss", 1,"RequiredEbNo", RequiredEbNo,"PreReceiverLoss", 0);
gaussianAntenna(gs2_rx_r, "ApertureEfficiency", 0.65,"DishDiameter", 3.7)
gs3_rx_r = receiver(gs3.Gimbals, "Name", "Ground Station 3 Receiver" ,"GainToNoiseTemperatureRatio", GTNTR_3,"SystemLoss", 1,"RequiredEbNo", RequiredEbNo,"PreReceiverLoss", 0);
gaussianAntenna(gs3_rx_r, "ApertureEfficiency", 0.65,"DishDiameter", 3.7)
gs4_rx_r = receiver(gs4.Gimbals, "Name", "Ground Station 4 Receiver" ,"GainToNoiseTemperatureRatio", GTNTR_4,"SystemLoss", 1,"RequiredEbNo", RequiredEbNo,"PreReceiverLoss", 0);
gaussianAntenna(gs4_rx_r, "ApertureEfficiency", 0.65,"DishDiameter", 3.7)

pointAt(EOSat.satellite,gs1);
lnk_1_r = link(tx,gs1_rx_r);
ac_1_r = access(EOSat.satellite,gs1);
intervals_1_r = linkIntervals(lnk_1_r);

pointAt(EOSat.satellite,gs2);
lnk_2_r = link(tx,gs2_rx_r);
ac_2_r = access(EOSat.satellite,gs2);
intervals_2_r = linkIntervals(lnk_2_r);

pointAt(EOSat.satellite,gs3);
lnk_3_r = link(tx,gs3_rx_r);
ac_3_r = access(EOSat.satellite,gs3);
intervals_3_r = linkIntervals(lnk_3_r);

pointAt(EOSat.satellite,gs4);
lnk_4_r = link(tx,gs4_rx_r);
ac_4_r = access(EOSat.satellite,gs4);
intervals_4_r = linkIntervals(lnk_4_r);

%timetable with only numerical value 
t_1_r = table2timetable(intervals_1_r(:,[3:4, 6:8]));
t_2_r = table2timetable(intervals_2_r(:,[3:4, 6:8]));
t_3_r = table2timetable(intervals_3_r(:,[3:4, 6:8]));
t_4_r = table2timetable(intervals_4_r(:,[3:4, 6:8]));

%complete timetables
t_1_c_r = table2timetable(intervals_1_r); 
t_2_c_r = table2timetable(intervals_2_r);
t_3_c_r = table2timetable(intervals_3_r);
t_4_c_r = table2timetable(intervals_4_r);

table_1_r = retime(t_1_r,'daily','mean'); %mean
table_1_2_r = retime(t_1_r,'daily','count'); %count
table_1_3_r = retime(t_1_r,'daily','sum'); %sum

% % Helper function to combine the tables of different stations
table_2_r = retime(t_2_r,'daily','mean');
table_2_2_r = retime(t_2_r,'daily','count');
table_2_3_r = retime(t_2_r,'daily','sum');

table_3_r = retime(t_3_r,'daily','mean');
table_3_2_r = retime(t_3_r,'daily','count');
table_3_3_r = retime(t_3_r,'daily','sum');

table_4_r = retime(t_4_r,'daily','mean');
table_4_2_r = retime(t_4_r,'daily','count');
table_4_3_r = retime(t_4_r,'daily','sum');

tot_table_access_r= sortrows(vertcat(t_1_c_r, t_2_c_r, t_3_c_r, t_4_c_r));

tot_table_mean_r= sortrows(vertcat(table_1_r, table_2_r, table_3_r, table_4_r));
tot_table_count_r= sortrows(vertcat(table_1_2_r, table_2_2_r, table_3_2_r, table_4_2_r));
tot_table_sum_r= sortrows(vertcat(table_1_3_r, table_2_3_r, table_3_3_r, table_4_3_r));

MDC_1_r = sum(table_1_3_r.Duration * Rb_net * 1e-6)/7; % mean daily capacity 
MDC_2_r = sum(table_2_3_r.Duration * Rb_net * 1e-6)/7; % mean daily capacity 
MDC_3_r = sum(table_3_3_r.Duration * Rb_net * 1e-6)/7; % mean daily capacity 
MDC_4_r = sum(table_4_3_r.Duration * Rb_net * 1e-6)/7; % mean daily capacity 

MDC_EN_r = MDC_1_r + MDC_2_r + MDC_3_r + MDC_4_r; %mean daily capacity for the entire network 

MDA_1_r = sum(table_1_2_r.IntervalNumber)/7; %mean daily accesses
MDA_2_r = sum(table_2_2_r.IntervalNumber)/7; %mean daily accesses
MDA_3_r = sum(table_3_2_r.IntervalNumber)/7; %mean daily accesses
MDA_4_r = sum(table_4_2_r.IntervalNumber)/7; %mean daily accesses

MDA_EN_r = MDA_1_r + MDA_2_r + MDA_3_r + MDA_4_r; %mean daily accesses for the entire network 

AV_1_r = sum(table_1_3_r.Duration)/(7*24*60*60)*100; % Availability
AV_2_r = sum(table_2_3_r.Duration)/(7*24*60*60)*100; % Availability
AV_3_r = sum(table_3_3_r.Duration)/(7*24*60*60)*100; % Availability
AV_4_r = sum(table_4_3_r.Duration)/(7*24*60*60)*100; % Availability

AV_EN_r = AV_1_r + AV_2_r + AV_3_r + AV_4_r; %mean daily access for the entire network 

% latency for station 1

lat_1_r = [];

for i=1:(height(t_1_c_r)-1)
    latency = t_1_c_r.StartTime(i+1) - t_1_c_r.EndTime(i);
    lat_1_r = [lat_1_r, latency];
end

lat_1_r_max = minutes(max(lat_1_r));
lat_1_r_mean = minutes(mean(lat_1_r));

% latency for station 2

lat_2_r = [];

for i=1:(height(t_2_c_r)-1)
    latency = t_2_c_r.StartTime(i+1) - t_2_c_r.EndTime(i);
    lat_2_r = [lat_2_r, latency];
end

lat_2_r_max = minutes(max(lat_2_r));
lat_2_r_mean = minutes(mean(lat_2_r));

% latency for station 3

lat_3_r = [];

for i=1:(height(t_3_c_r)-1)
    latency = t_3_c_r.StartTime(i+1) - t_3_c_r.EndTime(i);
    lat_3_r = [lat_3_r, latency];
end

lat_3_r_max = minutes(max(lat_3_r));
lat_3_r_mean = minutes(mean(lat_3_r));


% latency for station 4

lat_4_r = [];

for i=1:(height(t_4_c_r)-1)
    latency = t_4_c_r.StartTime(i+1) - t_4_c_r.EndTime(i);
    lat_4_r = [lat_4_r, latency];
end

lat_4_r_max = minutes(max(lat_4_r));
lat_4_r_mean = minutes(mean(lat_4_r));


% latency for the entire network

lat_EN_r = [];

for i=1:(height(tot_table_access_r)-1)
    latency = tot_table_access_r.StartTime(i+1) - tot_table_access_r.EndTime(i);
    lat_EN_r = [lat_EN_r, latency];
end

lat_EN_r_max = minutes(max(lat_EN_r));
lat_EN_r_mean = minutes(mean(lat_EN_r));

%% plot

%capacity

f=figure
subplot(2,2,1), bar(MDC_1,'b'), hold on, bar(2,MDC_1_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Avg. Daily Capacity [Tbit]'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 1')
subplot(2,2,2), bar(MDC_2,'b'), hold on, bar(2,MDC_2_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Avg. Daily Capacity [Tbit]'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 2')
subplot(2,2,3), bar(MDC_3,'b'), hold on, bar(2,MDC_3_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Avg. Daily Capacity [Tbit]'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 3')
subplot(2,2,4), bar(MDC_4,'b'), hold on, bar(2,MDC_4_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Avg. Daily Capacity [Tbit]'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 4')

sgtitle('Mean Daily Capacity')
saveas(figure(5),"Mean_Daily_Capacity_Stations",'epsc')


f=figure, bar(MDC_EN,'b'), hold on, bar(2,MDC_EN_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Avg. Daily Capacity [Tbit]'),  xticks([1 2]), xticklabels({'clear','rain'})
ylim([0 3])
title('Mean Daily Capacity, Entire Network')
saveas(figure(6),"Mean_Daily_Capacity_Network",'epsc')

% latencies

f=figure
subplot(2,2,1), bar(lat_1_mean,'b'), hold on, bar(2,lat_1_r_mean,'r'), bar(3,lat_1_max,'b'), bar(4,lat_1_r_max,'r'), grid on;
xlabel('Mean and Max'), ylabel('Minutes'),xticks([]), title('Station 1')
subplot(2,2,2), bar(lat_2_mean,'b'), hold on, bar(2,lat_2_r_mean,'r'), bar(3,lat_2_max,'b'), bar(4,lat_2_r_max,'r'), grid on;
xlabel('Mean and Max'), ylabel('Minutes'),xticks([]), title('Station 2')
subplot(2,2,3), bar(lat_3_mean,'b'), hold on, bar(2,lat_3_r_mean,'r'), bar(3,lat_3_max,'b'), bar(4,lat_3_r_max,'r'), grid on;
xlabel('Mean and Max'), ylabel('Minutes'),xticks([]), title('Station 3')
subplot(2,2,4), bar(lat_4_mean,'b'), hold on, bar(2,lat_4_r_mean,'r'), bar(3,lat_4_max,'b'), bar(4,lat_4_r_max,'r'), grid on;
xlabel('Mean and Max'), ylabel('Minutes'),xticks([]), title('Station 4')

sgtitle('Latency')
saveas(figure(7),"Latency_stations",'epsc')

f=figure, bar(lat_EN_mean,'b'), hold on, bar(2,lat_EN_r_mean,'r'), bar(3,lat_EN_max,'b'), bar(4,lat_EN_r_max,'r'), grid on;
xlabel('Mean and Max'), ylabel('Minutes'),xticks([]),title('Latency, Entire Network')
saveas(figure(8),"Latency_Network",'epsc')

% accesses

f=figure
subplot(2,2,1), bar(MDA_1,'b'), hold on, bar(2,MDA_1_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Number of times'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 1')
subplot(2,2,2), bar(MDA_2,'b'), hold on, bar(2,MDA_2_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Number of times'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 2'), ylim([0 12])
subplot(2,2,3), bar(MDA_3,'b'), hold on, bar(2,MDA_3_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Number of times'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 3')
subplot(2,2,4), bar(MDA_4,'b'), hold on, bar(2,MDA_4_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Number of times'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 4')

sgtitle('Mean Daily Accesses')
saveas(figure(9),"Mean_Daily_Access_Stations",'epsc')

f=figure, bar(MDA_EN,'b'), hold on, bar(2,MDA_EN_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Number of times'),  xticks([1 2]), xticklabels({'clear','rain'})
title('Mean Daily Accesses, Entire Network')
saveas(figure(10),"Mean_Daily_Access_Network",'epsc')

% Availability

f=figure
subplot(2,2,1), bar(AV_1,'b'), hold on, bar(2,AV_1_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Availability [%]'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 1')
subplot(2,2,2), bar(AV_2,'b'), hold on, bar(2,AV_2_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Availability [%]'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 2'), ylim([0 5])
subplot(2,2,3), bar(AV_3,'b'), hold on, bar(2,AV_3_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Availability [%]'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 3')
subplot(2,2,4), bar(AV_4,'b'), hold on, bar(2,AV_4_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Availability [%]'), xticks([1 2]), xticklabels({'clear','rain'}), title('Station 4')

sgtitle('Availability')
saveas(figure(11),"Availability_Stations",'epsc')

f=figure, bar(AV_EN,'b'), hold on, bar(2,AV_EN_r,'r'), grid on;
xlabel('Clear and Rainy Conditions'), ylabel('Availability [%]'),  xticks([1 2]), xticklabels({'clear','rain'})
title('Availability, Entire Network')
saveas(figure(12),"Availability_Network",'epsc')

%% SECOND PART

LDRS.planes = 3;
LDRS.satsPerPlane = 6;
LDRS.numSats = LDRS.planes*LDRS.satsPerPlane;
LDRS.eccentricity = 0; LDRS.inclination = 98.2;
LDRS.altitude = 1000e3;
LDRS.semiMajorAxis = earthRadius('meters') + LDRS.altitude;
LDRS.argOfPeriapsis = 0;
% Walker
LDRS.planePhase = 360/LDRS.planes;
LDRS.satsPhase = 360/LDRS.satsPerPlane;
LDRS.RAAN = zeros(1, LDRS.numSats);
LDRS.trueanomaly = zeros(1, LDRS.numSats);
idx = 0;
for plane = 1:LDRS.planes
    for sat=1:LDRS.satsPerPlane
        idx = idx + 1;
        LDRS.RAAN(idx)= plane * LDRS.planePhase;
        LDRS.trueanomaly(idx) = sat * LDRS.satsPhase;
    end
end

LDRS.Satellites = satellite(sc,LDRS.semiMajorAxis * ones(1, LDRS.numSats),LDRS.eccentricity * ones(1, LDRS.numSats),LDRS.inclination * ones(1, LDRS.numSats),LDRS.RAAN,LDRS.argOfPeriapsis * ones(1, LDRS.numSats),LDRS.trueanomaly);
%show(LDRS.Satellites)

% Add a conical sensor to gimbal
EOSat.LCTGimbal = gimbal(EOSat.satellite);
EOSat.LCT(1)= conicalSensor(EOSat.LCTGimbal, "MaxViewAngle", 100);

Max_range = 4e6;

pointAt(LDRS.Satellites,EOSat.satellite);
[az, el_ang, range] = aer(EOSat.satellite,LDRS.Satellites);
ac_new = access(EOSat.satellite,LDRS.Satellites);
ac_status = accessStatus(ac_new);

for i=1:size(range,1)
    for j=1:size(range,2)
        if range(i,j)>Max_range
            ac_status(i,j)=0;
        end
    end
end
connectivity=any(ac_status);
figure, bar(1:size(range,2),connectivity), xlabel('minutes'), ylabel('conn. status'), title('Connectivity, Entire System'), xlim([0 10081]),ylim([0 1])
saveas(figure(13),"Connectivity",'epsc')


AC_INT = accessIntervals(ac_new);

AC_INT_tt = sortrows(table2timetable(AC_INT)); % timetable

n=nnz(connectivity);
system_access_duration=n*sc.SampleTime;
scenarioDuration=seconds(sc.StopTime - sc.StartTime);
systemAccessPercentage=(system_access_duration/scenarioDuration)*100;

figure, bar(systemAccessPercentage), xticks(1), xticklabels({'system'}), ylabel('Availability [%]'), title('Availability, Entire LEO System'), ylim([0 100])
grid on 
saveas(figure(14),"SystemAccessPercentage",'epsc')

tot_latenza=[];

for i=1:length(connectivity)
    if connectivity(i)==0
        j=i;
        while connectivity(j)==0
            j=+j;
        end
    end
    latenza=(connectivity(j)-connectivity(i));
    tot_latenza=[tot_latenza,latenza];
end

latenza_mean=mean(tot_latenza);
latenza_max=max(tot_latenza);

figure, bar(1,latenza_mean), hold on, bar(2,latenza_max), xticks([1 2]), xticklabels({'mean','max'}), ylabel('Time [min]'), title('Latency, Entire LEO System'), ylim([0 100])
grid on 
saveas(figure(15),"Latency_LEO",'epsc')

%% Bonus: what happens if the satellite has to point in a certain direction (up (+90) or down (-90))?

time = startTime:seconds(sampleTime):stopTime;
steeringTable = timetable(time.', repmat([0 90], length(time), 1));
pointAt(EOSat.LCTGimbal,steeringTable)
[az, el_ang, range] = aer(EOSat.satellite,LDRS.Satellites);
ac_new = access(LDRS.Satellites,EOSat.LCT);
ac_status = accessStatus(ac_new);

for i=1:size(range,1)
    for j=1:size(range,2)
        if range(i,j)>Max_range
            ac_status(i,j)=0;
        end
    end
end
connectivity=any(ac_status);
figure, bar(1:size(range,2),connectivity), xlabel('minutes'), ylabel('conn. status'), title('Connectivity, Entire System'), xlim([0 10081]),ylim([0 1])
saveas(figure(16),"Connectivity2",'epsc')


AC_INT = accessIntervals(ac_new);

AC_INT_tt = sortrows(table2timetable(AC_INT)); % timetable

n=nnz(connectivity);
system_access_duration=n*sc.SampleTime;
scenarioDuration=seconds(sc.StopTime - sc.StartTime);
systemAccessPercentage=(system_access_duration/scenarioDuration)*100;

figure, bar(systemAccessPercentage), xticks(1), xticklabels({'system'}), ylabel('Availability [%]'), title('Availability, Entire LEO System'), ylim([0 100])
grid on 
saveas(figure(17),"SystemAccessPercentage2",'epsc')

tot_latenza=[];

for i=1:length(connectivity)
    if connectivity(i)==0
        j=i;
        while connectivity(j)==0
            j=j+1;
        end
    end
    latenza=(connectivity(j)-connectivity(i));
    tot_latenza=[tot_latenza,latenza];
end

latenza_mean=mean(tot_latenza);
latenza_max=max(tot_latenza);

figure, bar(1,latenza_mean), hold on, bar(2,latenza_max), xticks([1 2]), xticklabels({'mean','max'}), ylabel('Time [min]'), title('Latency, Entire LEO System'), ylim([0 10])
grid on 
saveas(figure(18),"Latency_LEO2",'epsc')