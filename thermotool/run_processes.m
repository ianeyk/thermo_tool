carnot = thermo_processes;

carnot.m = 0.01; % kg
carnot.P0 = 1 .* 10 .^ 5; % Pa
carnot.T0 = carnot.TK - 5; % K
carnot.get_V0

% carnot = carnot.isothermal(P = 3 .* 10 .^ 5);
% q(1) = carnot.Qnet
% w(1) = carnot.Wnet
% carnot = carnot.isochoric(P = 1 .* 10 .^ 5);
% q(2) = carnot.Qnet
% w(2) = carnot.Wnet
% % carnot = carnot.isobaric(T = carnot.TK - 5);
% carnot = carnot.isobaric(S = 0);
% q(3) = carnot.Qnet
% w(3) = carnot.Wnet

carnot = carnot.isentropic(P = 3 .* 10 .^ 5);
q(1) = carnot.Qnet
w(1) = carnot.Wnet
carnot = carnot.isobaric(T = carnot.TK + 25);
q(2) = carnot.Qnet
w(2) = carnot.Wnet
carnot = carnot.isentropic(P = 1 .* 10 .^ 5);
q(3) = carnot.Qnet
w(3) = carnot.Wnet
carnot = carnot.isobaric(T = carnot.TK - 5);
q(4) = carnot.Qnet
w(4) = carnot.Wnet