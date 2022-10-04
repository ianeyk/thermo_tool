% create thermo_process object to model the engine
engine = thermo_processes;

% set up constants, initial conditions
engine.m = 0.01; % kg
engine.R = 287; % J/kg/K
engine.cv = 716; % J/kg/K
engine.P0 = 1 .* 10 .^ 6; % Pa
engine.T0 = 450; % K
engine = engine.startCycle();

% define the cycle processes
engine = engine.isentropic(T = 300);
engine = engine.isothermal(V = engine.stateProperties.V(1), Qout = true);
engine = engine.isochoric(T = 450, Qin = true, stateName = "State '0");

% print out state variables and changes in states
disp(engine.stateProperties)
disp(engine.pathProperties)

Qin = engine.Qin % Heat transfer into the engine through Path "State 2 -> State '0"
Wnet = engine.Wnet % Net work done by the engine
efficiency = engine.Wnet ./ engine.Qin % Efficiency of the engine


% for plot readability, zoom out on the axes slightly
figure(1);
subplot(1, 2, 1);
expandAxes(1.1);

figure(1);
subplot(1, 2, 2);
expandAxes(1.1);