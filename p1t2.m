clc;
clear;
close all;
%% P1 Task 2. 

%% Constants
d_rod = 1*0.0254; % Diameter of rod, [m]
A_rod = 2*pi*d_rod/2; % Cross section of the rod, [m^2]
k = [130,130,115,115,16.2]; % Thermal Conductivity (k) [W/(m*K)]=[W/(m*C)];
Thermocouple_Err = 2; % Error in the thermocouples.

filename(1) = "Aluminum_21V_203mA.csv";
filename(2) = "Aluminum_30V_290mA.csv";
filename(3) = "Brass_21V_199mA.csv";
filename(4) = "Brass_30V_285mA.csv";
filename(5) = "Steel_21V_194mA.csv";

x_0 = (1+3/8).*0.0254;% Distance from x_0 to first thermocouple
spacing = 0.5.*0.0254; % Distance between thermocouples
pos_therm = linspace(x_0,x_0+(8*spacing),8); % [inches]

initialStates = zeros([5,9]);
M_exp = zeros(5,1);

figure('Position',[40 60 1100 700])
t = tiledlayout(2,3);
t.TileSpacing ="compact";
t.Padding = 'compact';
title(t,"Initial Temperature Gradient",'FontSize',16)

%% Go through each file
for i=1:length(filename)
    rawData(i) = importdata(filename(i));
    testData = rawData(i).data;
    initialStates(i,:) = testData(1,:);
    
    P = polyfit(pos_therm, initialStates(i,2:end),1);
    M_exp(i) = P(1);
    ys = polyval(P,pos_therm);

    titleFile = char (filename(i)); % Filename for the data
    
    % Individual PLots
    
    figure(1)
    nexttile;
    hold on;
    grid on;
    grid minor;
    
    plot(pos_therm,ys,'b',LineWidth=2)
    scatter(pos_therm,initialStates(i,2:end),25,'r','filled')
    plot(pos_therm,initialStates(i,2:end)+Thermocouple_Err,"--r", Color=[1,0,0])
    plot(pos_therm,initialStates(i,2:end)-Thermocouple_Err,"--r", Color=[1,0,0])

    ylim([0,25])
    % Titling Plots
    if (i==1 || i==2)
        volt = titleFile(1,10:11); % [V]
        curr = titleFile(1,14:16); % [mA]
        title(titleFile(1,1:8)+" "+volt+"V, "+curr+"mA",'FontSize',14);
    else
        volt = titleFile(1,7:8); % [V]
        curr = titleFile(1,11:13); % [mA]
        title(titleFile(1,1:5)+" "+volt+"V, "+curr+"mA",'FontSize',14);
    end
    xlabel("Distance [m]",'FontSize',14);
    ylabel("Temperature ["+char(176)+"C]",'FontSize',14)
end

ax = nexttile(1);
leg = legend("M_{exp}","Thermocouple Measurement","Thermocouple Error",'FontSize',13);
%leg.Position = [50,50,50,100];
leg.Layout.Tile = 6;



