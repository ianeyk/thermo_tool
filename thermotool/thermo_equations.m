syms Esys Q W Strans Sgen Esys0 m R cv cp T0 P0 V0 S0 T P V S x;

% Esys = 8;  % J
% Q = 4;     % J
% T = 100;   % K
% Sgen = 0;  % J/K
% R = 287;   % J/kg/K
% cv = 716;  % J/kg/K
% cp = cv + R; % J/kg/K
% m = 0.01;     % kg


TK = 273.15;

first_law = Esys - Esys0 == sum(Q) - sum(W);
Strans = sum(Q ./ T);
second_law = S - S0 == Strans + Sgen;
internal_energy = Esys - Esys0 == m .* cv .* (T - T0);
entropy_eqn = S - S0 == m .* (cv .* log(P ./ P0) + cp .* log(V ./ V0));
gas_law = P .* V == m .* R .* T;
gas_law_0 = P0 .* V0 == m .* R .* T0;

% T0 = TK - 5; % K
% P0 = 10 .^ 5; % N/m^2
% Esys0 = 0;
% S0 = 0; % unitless

% var_vals = [T0 == TK - 5;
%             P0 == 10 .^ 5;
% %             Esys0 == 0;
% %             S0 == 0;
% %             Sgen == 0;
%             R == 287;
%             cv == 716;
%             cp == cv + R;
%             m == 0.01;
%             S == S0;
%             P == 3 .* 10 .^ 5
%             ].';

var_vals = [
    {T0,    TK - 5      };
    {P0,    10 .^ 5     };
    {Esys0, 0           };
    {S0,    0           };
    {Sgen,  0           };
    {cp,    cv + R      };
    {R,     287         };
    {cv,    716         };
    {m,     0.01        };
    {S,     S0          };
    {P,     3 .* 10 .^5}
];

%% make plots
all_eqns = [[first_law, second_law, internal_energy, entropy_eqn, gas_law, gas_law_0], []];%var_vals];
new_eqn = Esys0 == 0;
for eqn = all_eqns
    eqn
    new_eqn = subs(new_eqn, eqn)
end
for var_i = 1:size(var_vals, 1)
    var = var_vals(var_i, 1)
    new_eqn = subs(new_eqn, var_vals(var_i, 1), var_vals(var_i, 2))
end

