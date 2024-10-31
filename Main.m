close all; clear; clc;

%% Constants
d_rod = 1*0.0254; % Diameter of rod, [m]
A_rod = pi*(d_rod^2)/4; % Cross section of the rod, [m^2]
k = [130,130,115,115,16.2]; % Thermal Conductivity (k) [W/(m*K)]=[W/(m*C)];
rho = [2810, 2810, 8500, 8500, 8000]; % Density [kg/m^3]
c_p = [960, 960, 380, 380, 500]; % Specific Heat Capacity [J/(kg*K)]
L = (6+(9/8))*0.0254; % Length of rod

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

%% Part 1
%% Go through each file
for i=1:length(filename)
    alpha(i) = k(i) / rho(i) / c_p(i);
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
    pos_therm = linspace(x_0,x_0+(7*spacing),8); % [m]
    scatter(pos_therm,T_F,25,'r','filled');
    
    T_init = testData(1,2:end);
    Coeff2 = polyfit(pos_therm,T_init, 1);
    H_m(i) = Coeff2(1);
    
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


%% Part 2
% Task 1
figure('Position', [40 60 600 300]); hold on; grid on; grid minor;

u1 = transientTemp_n(pos_therm(end), 1, T_0(1), H_an(1), L, k(1), rho(1), c_p(1), 10);
plot(linspace(0,10,11),u1,'.',LineWidth=6);
u2 = transientTemp_n(pos_therm(end), 1000, T_0(1), H_an(1), L, k(1), rho(1), c_p(1), 10);
plot(linspace(0,10,11),u2,".",LineWidth=6);
title({"Temperature at End of Aluminum 21V Rod"," General Solution Computed with Different N Values"})
xlabel("n value",'FontSize',14);
ylabel("Temperature ["+char(176)+"C]",'FontSize',14);
legend("t = 1s","t = 1000s",'Location','best');

hold off;

% Task 2
for i=1:length(filename)
    
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
    
    if (i==1 || i==2)
        titleArray = titleFile(1,1:8)+" "+volt+"V, "+curr+"mA";
    else
        titleArray = titleFile(1,1:5)+" "+volt+"V, "+curr+"mA";
    end

    figure('Position', [40 60 1000 400]); hold on; grid on; grid minor;
    xlabel("Time [min]",'FontSize',14);
    ylabel("Temperature ["+char(176)+"C]",'FontSize',14);
    title({titleArray,"Model 1A: H_a_n"},'FontSize',14);
    u3 = transientTemp(pos_therm,testData(:,1),T_0(i),H_an(i),L,k(i),rho(i),c_p(i),10);
    plot(testData(:,1)/60, u3,'-',"LineWidth",1.5);
    plot(testData(:,1)/60,testData(:,2:end),'--',"LineWidth",1.5)

   
    ax1 = gca;
    myColors1 = summer(13);
    myColors2 = winter(40);
    myColors = [myColors1(1:8,:);myColors2(15:23,:)];
    ax1.ColorOrder = myColors;
    legend("Exp Th_1","Exp Th_2","Exp Th_3","Exp Th_4","Exp Th_5","Exp Th_6","Exp Th_7","Exp Th_8","H_a_n Th_1","H_a_n Th_2","H_a_n Th_3","H_a_n Th_4","H_a_n Th_5","H_a_n Th_6","H_a_n Th_7","H_a_n Th_8","NumColumns",2,"location",'eastoutside');

    figure('Position', [40 60 1000 400]); hold on; grid on; grid minor;
    xlabel("Time [min]",'FontSize',14);
    ylabel("Temperature ["+char(176)+"C]",'FontSize',14);
    title({titleArray,"Model 1B: H_e_x_p"},'FontSize',14);
    u4 = transientTemp(pos_therm,testData(:,1),T_0(i),H_exp(i),L,k(i),rho(i),c_p(i),10);
    plot(testData(:,1)/60, u4,'-',"LineWidth",1.5);
    plot(testData(:,1)/60,testData(:,2:end),'--',"LineWidth",1.5)
    ax2 = gca;
    myColors1 = cool(30);
    myColors2 = winter(40);
    myColors = [myColors1(18:25,:);myColors2(15:23,:)];
    ax2.ColorOrder = myColors;
    legend("H_e_x_p Th_1","H_e_x_p Th_2","H_e_x_p Th_3","H_e_x_p Th_4","H_e_x_p Th_5","H_e_x_p Th_6","H_e_x_p Th_7","H_e_x_p Th_8","Exp Th_1","Exp Th_2","Exp Th_3","Exp Th_4","Exp Th_5","Exp Th_6","Exp Th_7","Exp Th_8","NumColumns",2,"location",'eastoutside');


    figure('Position', [40 60 1000 400]); hold on; grid on; grid minor;
    xlabel("Time [min]",'FontSize',14);
    ylabel("Temperature ["+char(176)+"C]",'FontSize',14);
    title({titleArray,"Model II: (M-H)x"},'FontSize',14);
    u4 = transientTemp_Hm(pos_therm,testData(:,1),T_0(i),H_exp(i),H_m(i),L,k(i),rho(i),c_p(i),10);
    plot(testData(:,1)/60, u4,'-',"LineWidth",1.5);
    plot(testData(:,1)/60,testData(:,2:end),'--',"LineWidth",1.5)
    ax2 = gca;
    myColors1 = hot(30);
    myColors2 = winter(40);
    myColors = [myColors1(7:14,:);myColors2(15:23,:)];
    ax2.ColorOrder = myColors;
    legend("H_e_x_p Th_1","H_e_x_p Th_2","H_e_x_p Th_3","H_e_x_p Th_4","H_e_x_p Th_5","H_e_x_p Th_6","H_e_x_p Th_7","H_e_x_p Th_8","Exp Th_1","Exp Th_2","Exp Th_3","Exp Th_4","Exp Th_5","Exp Th_6","Exp Th_7","Exp Th_8","NumColumns",2,"location",'eastoutside');
   
end

% This version of the function is for plotting u as n changes at a single x
% and t value. Used for p2t1.
function u = transientTemp_n(x,t,T_0,H,L,k,rho,c_p,n_max)
    total = 0; % Sets the initial transient state term (n=0) to be 0
    alpha = k./(rho.*c_p); % Gets alpha
    u = zeros(1,n_max+1); % Pre-allocates u for speed
    u(1) = T_0 + H.*x; % Establishes n=0 as the steady state case.
    for n = 1:n_max % Loops through values of n
        n;
        lambda = ((2.*n-1)*pi)./(2.*L); % Calculates lambda
        % The below calculates b
        if mod(n,2) == 0 % If n is even
            b = (8.*H.*L)./((2.*n-1).^2.*pi.^2);
        else
            b = (-8.*H.*L)./((2.*n-1).^2.*pi.^2);
        end
        % Updates the sum for each value of n
        total = total + b.*sin(lambda.*x).*exp(-1.*lambda.^2.*alpha.*t);
        % Gives a u value for each value of n
        u(n+1) = T_0 + H.*x + total;
    end

end


% This version of the function is for plotting u as x and t change, with a fixed value of n_max. Used for p2t2.
function u = transientTemp(x,t,T_0,H,L,k,rho,c_p,n_max)
    total = zeros(length(t),length(x)); % Sets the initial transient state term (n=0) to be 0
    alpha = k./(rho.*c_p); % Gets alpha
    u = zeros(length(t),length(x)); % Pre-allocates u for speed
    %u(1,:) = T_0 + H.*x; % Establishes n=0 as the steady state case.
    for n = 1:n_max % Loops through values of n
        n;
        lambda = ((2.*n-1)*pi)./(2.*L); % Calculates lambda
        % The below calculates b
        if mod(n,2) == 0 % If n is even
            b = (8.*H.*L)./((2.*n-1).^2.*pi.^2);
        else
            b = (-8.*H.*L)./((2.*n-1).^2.*pi.^2);
        end
        % Updates the sum for each value of n
        total = total + b.*sin(lambda.*x).*exp(-1.*lambda.^2.*alpha.*t);
        % Gives a u value for each value of n
        
    end
    % Adds the transient values for each x and t. 
    u(1:end,:) = T_0 + H.*x + total;
    %u(1,:) = T_0 + H.*x; % Establishes n=0 as the steady state case.
end

% This version of the function is for plotting u as x and t change, with a fixed value of n_max, but accounts for an initial T gradient. Used for p2t4.
function u = transientTemp_Hm(x,t,T_0,H,Hm,L,k,rho,c_p,n_max)
total = zeros(length(t),length(x)); % Sets the initial transient state term (n=0) to be 0
alpha = k./(rho.*c_p); % Gets alpha
u = zeros(length(t),length(x)); % Pre-allocates u for speed
%u(1,:) = T_0 + H.*x; % Establishes n=0 as the steady state case.
for n = 1:n_max % Loops through values of n
    n;
    lambda = ((2.*n-1)*pi)./(2.*L); % Calculates lambda
    % The below calculates b
    if mod(n,2) == 0 % If n is even
        b = (8.*(H-Hm).*L)./((2.*n-1).^2.*pi.^2);
    else
        b = (-8.*(H-Hm).*L)./((2.*n-1).^2.*pi.^2);
    end
    % Updates the sum for each value of n
    total = total + b.*sin(lambda.*x).*exp(-1.*lambda.^2.*alpha.*t);
    % Gives a u value for each value of n
    
end
% Adds the transient values for each x and t. 
u(1:end,:) = T_0 + H.*x + total;
%u(1,:) = T_0 + H.*x; % Establishes n=0 as the steady state case.
end
