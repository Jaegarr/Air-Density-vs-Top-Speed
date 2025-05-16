%% Track Data
Track = {'Bahrain'; 'Jeddah'; 'Melbourne'; 'Suzuka'; 'Shanghai';'Miami'; 'Imola'; 'Monaco'; 'Montreal'; 'Barcelona';'Spielberg'; 'Silverstone'; 'Budapest'; 'Spa'; 'Zandvoort';'Monza'; 'Baku'; 'Singapore'; 'Austin'; 'Mexico City';'São Paulo'; 'Las Vegas'; 'Lusail'; 'Yas Marina'};
Elevation = [7; 3; 31; 35; 3; 3; 35; 10; 10; 58; 660; 170; 100; 368; 5; 162; -28; 15; 225; 2240; 784; 630; 5; 5]; %In meters
AvgTemp  = [30; 32; 18; 25; 22; 29; 20; 28; 22; 26; 20; 17; 30; 18; 17; 24; 35; 30; 28; 16; 22; 35; 42; 32]; % In degree Celcius
tbl = table(Track, Elevation, AvgTemp);
%% Properties
%Universal constants
p0      = 101325; % Sea-level standard pressure [Pa]
H_scale = 8434;  % Scale height [m] for barometric formula
g       = 9.81;      % Gravity [m/s^2]
R       = 287;       % Specific gas constant for air [J/(kg·K)]
gamma = 1.4;     % Heat capacity ratio
% Engine parameters
Vd      = 0.0016;        % Displacement volume [m^3], Article xxx
RPM_max = 12500;    % Maximum engine speed [RPM], used the rev limiter value from Article xxx
VE      = 0.95;            % Volumetric efficiency (assumed constant)
AFR     = 14.7;           % Stoichiometric AFR 
LHV     = 44e6;          % Lower heating value of fuel [J/kg]
n_combustion= 0.50; % Combustion efficiency (assumed), this number is thrown around a lot so why not?
n_MGUK = 0.88;          % This efficiency value is taken from the efficiency map I have used in my SAE publication for 12500 RPM
inlet_d_2025 = 345e-4 % Inlet valve diameter, Article 5.6.13
inlet_d_2026 = 345e-4 % Inlet valve diameter, Article 5.6.13
% MGU-K
P_MGUK_2025 = 120e3;    % MGU-K power limit 2025 [W]
P_MGUK_2026 = 350e3;    % MGU-K power limit 2026 [W]

% Vehicle Characteristics 
Cd     = 0.78; % Drag coefficient
Af = 1.5         % frontal area [m^2]
Cr = 0.01;    % Rolling resistance coefficient, this value is taken from my SAE publication
mass = 795;  % 2026 Vehicle mass including driver [kg], Article 4.1 and 4.6.2 

%% 3) Preallocate results
n = height(tbl);
rho              = zeros(n,1);
mdot_air         = zeros(n,1);
mdot_fuel        = zeros(n,1);
P_ICE            = zeros(n,1);
P_tot_2025       = zeros(n,1);
P_tot_2026       = zeros(n,1);
v_max_2025       = zeros(n,1);
v_max_2026       = zeros(n,1);

%% 4) Loop through each track and compute
for i = 1:n
    % Ambient conditions
    T_K = tbl.AvgTemp(i) + 273.15;
    p_track = p0 * exp(-tbl.Elevation(i) / H_scale);
    rho(i) = p_track / (R * T_K);
    % Air mass flow [kg/s]
    mdot_air(i) = rho(i) * Vd * RPM_max/60 * VE;  % RPM→rev/s division
    % Fuel mass flow [kg/s]
    mdot_fuel(i) = mdot_air(i) / AFR;
    % ICE power [W]
    P_ICE(i) = mdot_fuel(i) * LHV * n_combustion;
    % Total peak power [W]
    P_tot_2025(i) = P_ICE(i) + P_MGUK_2025*n_MGUK;
    P_tot_2026(i) = P_ICE(i) + P_MGUK_2026*n_MGUK;
    % Solve for top speed [m/s]: P = F_drag + F_rr * v
    f_speed = @(v,P) 0.5 * rho(i) * (Cd*Af) * v.^3 + Cr * mass * g * v - P;
    v_max_2025(i) = fzero(@(v) f_speed(v, P_tot_2025(i)), 100);
    v_max_2026(i) = fzero(@(v) f_speed(v, P_tot_2026(i)), 100);
end
%% Results
tbl.Rho_kg_m3 = rho;
tbl.MdotAir_kg_s = mdot_air;
tbl.MdotFuel_kg_s = mdot_fuel;
tbl.Power_ICE_kW = P_ICE/1e3;
tbl.Power_total2025_kW = P_tot_2025/1e3;
tbl.Power_total2026_kW = P_tot_2026/1e3;
tbl.Vmax2025_kph = v_max_2025 * 3.6;
tbl.Vmax2026_kph = v_max_2026 * 3.6;
disp(tbl);
subplot(2,1,1);
plot(tbl.Power_total2025_kW)
hold on
plot(tbl.Power_total2026_kW)
subplot(2,1,2);
plot(tbl.Vmax2025_kph)
hold on
plot(tbl.Vmax2026_kph)