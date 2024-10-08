close all; clear; clc;

%% Constants
d_rod = 1*0.0254; % Diameter of rod, 1m
A_rod = 2*pi*d_rod/2; % Cross section of the rod, m^2
k = [130,130,115,115,16.2]; % Thermal Conductivity (k) [W/(m*K)];

filename(1) = "Aluminum_21V_203mA.csv";
filename(2) = "Aluminum_30V_290mA.csv";
filename(3) = "Brass_21V_199mA.csv";
filename(4) = "Brass_30V_285mA.csv";
filename(5) = "Steel_21V_194mA.csv";

%% Go through each file
for j=1:1%length(filename)
    figure();
    hold on;
    set(gca,"ColorOrder",jet(25))

    titleFile = char (filename(j)); % Filename for the data
    % Voltage and Current
    if (j==1 || j==2)
        volt = titleFile(1,10:11);
        curr = titleFile(1,14:16);
        title(titleFile(1,1:8)+" "+volt+"V, "+curr+"mA");
    else
        volt = titleFile(1,7:8);
        curr = titleFile(1,11:13);
        title(titleFile(1,1:5)+" "+volt+"V, "+curr+"mA");
    end

    rawData = importdata(filename(j));
    testData = rawData.data;
    plot(testData(:,1),testData(:,2:end));
    legend("Ch1","Ch2","Ch3","Ch4","Ch5","Ch6","Ch7","Ch8");


    T_F = zeros(1);
    for j=2:9
        T_F(1,j-1)=testData(end,j);
    end
    pos_therm = linspace(7.62/100,16.51/100,8);

    Coeff = polyfit(pos_therm,T_F,1);
    H_exp = Coeff(1);
    T_0 = Coeff(2);

    hold off;
end


