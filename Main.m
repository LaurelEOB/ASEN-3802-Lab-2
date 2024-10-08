close all; clear; clc;

d_rod = 1*0.0254; % Diameter of rod, 1m
A_rod = 2*pi*d_rod/2; % Cross section of the rod, m^2
k = [130,130,115,115,16.2]; % Thermal Conductivity (k) [W/(m*K)];

Alum21V = readtable("Aluminum_21V_203mA.csv");
Alum30V = readtable("Aluminum_30V_290mA.csv");
Brass21V = readtable("Brass_21V_199mA.csv");
Brass30V = readtable("Brass_30V_285mA.csv");
Steel21V = readtable("Steel_21V_194mA.csv");

figure();
plot(Alum21V,"Time_s_",["CH1_C_","CH2_C_","CH3_C_","CH4_C_","CH5_C_","CH6_C_","CH7_C_","CH8_C_"]);
title("Aluminum 21V,203mA");
legend("Time","CH1_C_")

T_F = zeros(1,5);
for i=1:8
    A = Alum21V(:,i+1);
    T_F(1,i)=A{end,1};
end
pos_therm = linspace(7.62/100,16.51/100,8);

CoeffAlum21V = polyfit(pos_therm,T_F,1);
H_exp = CoeffAlum21V(1);
T_0 = CoeffAlum21V(2);
figure();
plot(pos_therm,T_F);





% 
% figure();
% plot(Alum30V,"Time_s_",["CH1_C_","CH2_C_","CH3_C_","CH4_C_","CH5_C_","CH6_C_","CH7_C_","CH8_C_"]);
% title("Aluminum 30V,290mA");
% legend("Time","CH1_C_")
% 
% figure();
% plot(Brass21V,"Time_s_",["CH1_C_","CH2_C_","CH3_C_","CH4_C_","CH5_C_","CH6_C_","CH7_C_","CH8_C_"]);
% title("Brass 21V,199mA");
% legend("Time","CH1_C_")
% 
% figure();
% plot(Brass30V,"Time_s_",["CH1_C_","CH2_C_","CH3_C_","CH4_C_","CH5_C_","CH6_C_","CH7_C_","CH8_C_"]);
% title("Brass 30V,285mA");
% legend("Time","CH1_C_")
% 
% figure();
% plot(Steel21V,"Time_s_",["CH1_C_","CH2_C_","CH3_C_","CH4_C_","CH5_C_","CH6_C_","CH7_C_","CH8_C_"]);
% title("Steel 21V,194mA");
% legend("Time","CH1_C_")
% 
% 
% 
% 
% 
