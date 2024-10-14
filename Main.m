close all; clear; clc;

%% Constants
d_rod = 1*0.0254; % Diameter of rod, [m]
A_rod = pi*(d_rod^2)/4; % Cross section of the rod, [m^2]
k = [130,130,115,115,16.2]; % Thermal Conductivity (k) [W/(m*K)]=[W/(m*C)];

filename(1) = "Aluminum_21V_203mA.csv";
filename(2) = "Aluminum_30V_290mA.csv";
filename(3) = "Brass_21V_199mA.csv";
filename(4) = "Brass_30V_285mA.csv";
filename(5) = "Steel_21V_194mA.csv";
figure('Position', [40 60 1100 700]); hold on;
t2 = tiledlayout(2,3);
%t2.TileSpacing = 'compact';
t2.Padding = 'compact';
title(t2,"Steady State Temperatures along the Rod",'FontSize',16);


%% Go through each file
for i=1:length(filename)
    figure(1);
    nexttile; hold on; grid on; grid minor;
    %set(gca,"ColorOrder",parula(16))
    titleFile = char (filename(i)); % Filename for the data
    % Voltage and Current
    if (i==1 || i==2)
        volt = titleFile(1,10:11); % [V]
        curr = titleFile(1,14:16); % [mA]
    else
        volt = titleFile(1,7:8); % [V]
        curr = titleFile(1,11:13); % [mA]
    end

    rawData = importdata(filename(i));
    testData = rawData.data;

    xlabel("Location along Rod [m]",'FontSize',14);
    ylabel("Temperature ["+char(176)+"C]",'FontSize',14)
    
    if (i==1 || i==2)
        title(titleFile(1,1:8)+" "+volt+"V, "+curr+"mA",'FontSize',14);
    else
        title(titleFile(1,1:5)+" "+volt+"V, "+curr+"mA",'FontSize',14);
    end

    for j=2:9
        T_F(1,j-1)=testData(end,j);
    end
    x_0 = (1+3/8)*0.0254;% Distance from x_0 to first thermocouple
    spacing = 0.5*0.0254; % Distance between thermocouples
    pos_therm = linspace(x_0,x_0+(8*spacing),8); % [m]
    scatter(pos_therm,T_F,25,'r','filled');
    
    Coeff = polyfit(pos_therm,T_F,1);
    H_exp(i) = Coeff(1); % [C/m]
    T_0(i) = Coeff(2); % [C]
    H_an(i) = str2num(volt)*str2num(curr)*(10^-3)/k(i)/A_rod; % [C/m]
    v_exp = T_0(i) + H_exp(i)*pos_therm;
    v_an = T_0(i) + H_an(i)*pos_therm;
    [M,N] = size(testData(:,1));
    plot(pos_therm,H_exp(i)*pos_therm+Coeff(2),'r',LineWidth=2);
    plot(pos_therm,H_exp(i)*pos_therm+Coeff(2)+2,'--r',LineWidth=1);
    plot(pos_therm,H_exp(i)*pos_therm+Coeff(2)-2,'--r',LineWidth=1);
    plot(pos_therm,H_an(i)*pos_therm+Coeff(2),'b',LineWidth=2);
end
ax = nexttile(1);
lg  = legend('Experimental Data','H_e_x_p (Experimental Slope)','Thermocouple error','','H_a_n (Analytical Slope)','Orientation','Vertical','FontSize',13);
lg.Layout.Tile = 6;

