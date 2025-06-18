%% Data Input
tracks = { ...
    'Bahrain','Jeddah','Melbourne','Suzuka','Shanghai','Miami', ...
    'Imola','Monaco','Montreal','Barcelona','Spielberg','Silverstone', ...
    'Budapest','Spa','Zandvoort','Monza','Baku','Singapore','Austin', ...
    'Mexico City','São Paulo','Las Vegas','Lusail','Yas Marina'}';
elevation_m = [7;3;31;35;3;3;35;10;10;58;660;170;100;368;5;162;-28;15;225;2240;784;630;5;5];
avgTemp_C = [30;32;18;25;22;29;20;28;22;26;20;17;30;18;17;24;35;30;28;16;22;35;42;32];
data = table(tracks, elevation_m, avgTemp_C, 'VariableNames', ...
    {'Track','Elevation_m','AvgTemp_C'});
%% Properties
% Universal constants (SI)
p0_Pa       = 101325; % Sea-level pressure
Hscale_m    = 8434;   % Scale height
g_m_s2      = 9.81;    % Gravity
R_J_kgK     = 287;     % Gas constant
gamma       = 1.4;     % Heat capacity ratio
% Power-unit parameters
Vd_m3       = 0.0016;
RPM_max     = 12500;
AFR_g         = 14.7;
AFR_e         = 9;
LHV_g_J_kg  = 44e6;
LHV_e_J_kg  = 26.8e6;
eta_comb    = 0.50;
eta_MGUK    = 0.88;
inlet_d_m   = 345e-4;
P_MGUK25_W  = 120e3;
P_MGUK26_W  = 350e3;
% Vehicle characteristics
Cd          = 0.78;
Af_m2       = 1.5;
Cr          = 0.01;
mass25_kg   = 800;
mass26_kg   = 768;
%% Drag + Rolling-Resistance function
function F = dragRR(v,rho,Cd,Af,Cr,mass,g)
    F = 0.5*rho*(Cd*Af)*v.^3 + Cr*mass*g.*v;
end
%% Rho and Air Mass Flow
T_K = data.AvgTemp_C + 273.15;
p_Pa = p0_Pa .* exp(-data.Elevation_m ./ Hscale_m);
rho = p_Pa ./ (R_J_kgK .* T_K);
mdot_air = rho .* (pi/4*inlet_d_m^2) .* sqrt(gamma * R_J_kgK .* T_K);
mdot_fuel = mdot_air ./ AFR_g;
mdot_fuel_e = mdot_air./AFR_e;
%% Power calculations
P_ICE25 = mdot_fuel .* LHV_g_J_kg .* eta_comb;
P_ICE26 = mdot_fuel .* LHV_g_J_kg .* eta_comb;
P_ICE26e= mdot_fuel_e .* LHV_e_J_kg .* eta_comb;
P_tot25 = P_ICE25 + P_MGUK25_W .* eta_MGUK;
P_tot26 = P_ICE26 + P_MGUK26_W .* eta_MGUK;
P_tot26e= P_ICE26e + P_MGUK26_W .* eta_MGUK;
%% Top-speed solver
solveVmax = @(Ptot,mass) arrayfun(@(rho_i,P_i) ...
    fzero(@(v) dragRR(v,rho_i,Cd,Af_m2,Cr,mass,g_m_s2) - P_i, 100), ...
    rho, Ptot);
v25 = solveVmax(P_tot25, mass25_kg);
v26 = solveVmax(P_tot26, mass26_kg);
v26e= solveVmax(P_tot26e,mass26_kg);
%% Results
data.Rho = rho;
data.MdotAir = mdot_air;
data.MdotFuel= mdot_fuel;
data.Power25_kW = P_tot25/1e3;
data.Power26_kW = P_tot26/1e3;
data.Power26e_kW= P_tot26e/1e3;
data.Vmax25_kph = v25*3.6;
data.Vmax26_kph = v26*3.6;
data.Vmax26e_kph= v26e*3.6;
%Validation
% Top speed data from 2024 season (fill in what you have)
%Top speed data: Las Vegas: 356.4, Monza: 353.5, Mexico City: 370, Miami:
%343.5, , Singapore: 314.8
realTrackNames = {'Spa'; 'Monza'; 'Las Vegas'; 'Jeddah'; 'Austin','Miami'};
realVmax = [357.2; 353.5; 356.4; 341; 329.4;343.5]; % [km/h]
% Step 1: Initialize result arrays
n_real = numel(realTrackNames);
matchedVcalc = zeros(n_real, 1);
accuracy = zeros(n_real, 1);

% Step 2: Loop and find matches in your dataset
for i = 1:n_real
    idx = find(strcmp(data.Track, realTrackNames{i}));
    if ~isempty(idx)
        matchedVcalc(i) = data.Vmax2025_kph(idx); % or use Vmax2026_kph, etc.
        accuracy(i) = 1 - abs((matchedVcalc(i) - realVmax(i)) / realVmax(i));
    else
        matchedVcalc(i) = NaN;
        accuracy(i) = NaN;
    end
end
% Step 3: Build results table
comparison = table(realTrackNames, realVmax, matchedVcalc, accuracy, ...
    'VariableNames', {'Track', 'ObservedVmax_kph', 'CalculatedVmax_kph', 'Accuracy'})
% Step 4 (optional): Plot with ±10% tolerance bands
figure;
hold on
x = categorical(realTrackNames);
x = reordercats(x, realTrackNames);
bar(x, realVmax, 'FaceAlpha', 0.3); % Real
plot(x, matchedVcalc, 'ro-', 'LineWidth', 2); % Calculated
% Plot tolerance bands
upper = realVmax * 1.1;
lower = realVmax * 0.9;
plot(x, upper, 'k--');
plot(x, lower, 'k--');
ylabel('Top Speed [km/h]')
legend('2024 Observed', 'Calculated', '+10% Band', '-10% Band', 'Location', 'best');
title('Top Speed Comparison with ±10% Tolerance')
grid on
% Plotting

trackCat = categorical(data.Track);
trackCat = reordercats(trackCat, data.Track);  % preserve original order
figure;
tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile;
plot(trackCat,data.Rho);ylabel('Air Density[kg/m3]');
nexttile;
plot(trackCat, data.Power25_kW, 'DisplayName','2025'); 
hold on
plot(trackCat, data.Power26_kW, 'DisplayName','2026');
plot(trackCat, data.Power26e_kW,'DisplayName','2026 Ethanol');
ylabel('Total Power [kW]'); legend('Location','southwest');
nexttile;
plot(trackCat, data.Vmax25_kph, 'DisplayName','2025'); 
hold on
plot(trackCat, data.Vmax26_kph, 'DisplayName','2026');
plot(trackCat, data.Vmax26e_kph,'DisplayName','2026 Ethanol');
ylabel('Top Speed [kph]'); legend('Location','southwest');

