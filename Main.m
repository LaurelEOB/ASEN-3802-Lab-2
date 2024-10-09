close all; clear; clc;

%% Constants
d_rod = 1*0.0254; % Diameter of rod, [m]
A_rod = 2*pi*d_rod/2; % Cross section of the rod, [m^2]
k = [130,130,115,115,16.2]; % Thermal Conductivity (k) [W/(m*K)]=[W/(m*C)];

filename(1) = "Aluminum_21V_203mA.csv";
filename(2) = "Aluminum_30V_290mA.csv";
filename(3) = "Brass_21V_199mA.csv";
filename(4) = "Brass_30V_285mA.csv";
filename(5) = "Steel_21V_194mA.csv";

%% Go through each file
for i=1:length(filename)
    figure();
    hold on;
    set(gca,"ColorOrder",jet(25))

    titleFile = char (filename(i)); % Filename for the data
    % Voltage and Current
    if (i==1 || i==2)
        volt = titleFile(1,10:11); % [V]
        curr = titleFile(1,14:16); % [mA]
        title(titleFile(1,1:8)+" "+volt+"V, "+curr+"mA");
    else
        volt = titleFile(1,7:8); % [V]
        curr = titleFile(1,11:13); % [mA]
        title(titleFile(1,1:5)+" "+volt+"V, "+curr+"mA");
    end

    rawData = importdata(filename(i));
    testData = rawData.data;
    plot(testData(:,1),testData(:,2:end));
    xlabel("time [s]");
    ylabel("temperature ["+char(176)+"C]")
    legend("Ch1","Ch2","Ch3","Ch4","Ch5","Ch6","Ch7","Ch8");


    for j=2:9
        T_F(1,j-1)=testData(end,j);
    end
    x_0 = (1+3/8)*0.0254;% Distance from x_0 to first thermocouple
    spacing = 0.5*0.0254; % Distance between thermocouples
    pos_therm = linspace(x_0,x_0+(8*spacing),8); % [m]

    Coeff = polyfit(pos_therm,T_F,1);
    H_exp(i) = Coeff(1)/k(i)/A_rod;
    T_0(i) = Coeff(2);
    %%% Units!!
    H_an = str2num(volt)*str2num(curr)*(10^-3)/k(i)/A_rod;

    hold off;
end


