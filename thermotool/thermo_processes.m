classdef thermo_processes
    properties
        m {mustBeNumeric}
        R {mustBeNumeric}
        cv {mustBeNumeric}
        TK {mustBeNumeric} % Celcius to Kelvin conversion
        nPoints {mustBeNumeric}

        T0 {mustBeNumeric}
        P0 {mustBeNumeric}
        V0 {mustBeNumeric}
        Snet {mustBeNumeric}

        Esys_net {mustBeNumeric}
        Esys {mustBeNumeric}
        Qnet {mustBeNumeric}
        Wnet {mustBeNumeric}
        Strans {mustBeNumeric}
        Sgen {mustBeNumeric}
        cp {mustBeNumeric}

        Qin {mustBeNumeric}
        Qout {mustBeNumeric}
        Win {mustBeNumeric}
        Wout {mustBeNumeric}

        stateProperties
        pathProperties
        currentStateNum {mustBeNumeric}
        currentStateName
        prevStateName
    end

    % properties (Dependent)
    % end

    methods
        function obj = thermo_processes()
            % any of these may be overwritten; these are just good default values
            obj.m = 1; % kg
            obj.R = 287; % J/kg/K
            obj.cv = 716; % J/kg/K
            obj.TK = 273.15; % K
            obj.nPoints = 100;

            obj.Esys_net = 0; % J
            obj.Esys = NaN;
            obj.Qnet = 0;
            obj.Wnet = 0;
            obj.Strans = NaN;
            obj.Sgen = NaN;
            obj.cp = NaN; % will be re-assigned upon being called
            obj.T0 = NaN;
            obj.P0 = NaN;
            obj.V0 = NaN;
            obj.Snet = 0; % J/kg/K

            obj.Qin = 0;
            obj.Qout = 0;
            obj.Win = 0;
            obj.Wout = 0;

            obj.stateProperties = array2table(zeros(0, 8));
            obj.pathProperties = array2table(zeros(0, 8));

            obj.currentStateNum = 0;
            obj.currentStateName = NaN;
            obj.prevStateName = "State 0";

            figure(1);
            clf;
        end

        function obj = startCycle(obj, NameValueArgs) % run this function after defining all the initial states in the script
            arguments
                obj
                NameValueArgs.stateName = ""
            end
            if NameValueArgs.stateName == ""
                stateName = obj.prevStateName;
            else
                stateName = NameValueArgs.stateName;
            end
            obj.stateProperties.Properties.VariableNames = {'State',  'T',        'ΔSnet',  'P',        'V',        'ΔEnet',      'ΔQnet',  'ΔWnet'};
            obj.stateProperties = [obj.stateProperties;    {stateName, obj.get_T0, obj.Snet, obj.get_P0, obj.get_V0, obj.Esys_net, obj.Qnet, obj.Wnet}];
        end

        function value = get.cp(obj)
            value = obj.cv + obj.R;
        end

        function value = get_T_from_PV(obj, P, V)
            value = P .* V ./ (obj.m .* obj.R);
        end

        function value = get_P_from_TV(obj, T, V)
            value = obj.m .* obj.R .* T ./ V;
        end

        function value = get_V_from_PT(obj, P, T)
            value = obj.m .* obj.R .* T ./ P;
        end

        function value = get_T0(obj)
            if isnan(obj.T0)
                value = obj.P0 .* obj.V0 ./ (obj.m ./ obj.R);
            else
                value = obj.T0;
            end
        end

        function value = get_P0(obj)
            if isnan(obj.P0)
                value = obj.m .* obj.R .* obj.T0 ./ obj.V0;
            else
                value = obj.P0;
            end
        end

        function value = get_V0(obj)
            if isnan(obj.V0)
                value = obj.m .* obj.R .* obj.T0 ./ obj.P0;
            else
                value = obj.V0;
            end
        end

        function obj = get_initial_state(obj) % solves for missing values using PV = mRT, if any
            obj.T0 = obj.get_T0;
            obj.P0 = obj.get_P0;
            obj.V0 = obj.get_V0;
        end


        function obj = isentropic(obj, NameValueArgs)
            arguments
                obj
                NameValueArgs.P = NaN
                NameValueArgs.V = NaN
                NameValueArgs.T = NaN
                NameValueArgs.Qin = false
                NameValueArgs.Qout = false
                NameValueArgs.Win = false
                NameValueArgs.Wout = false
                NameValueArgs.stateName = ""
            end
            obj = obj.get_initial_state();
            S = linspace(obj.Snet, obj.Snet, obj.nPoints); % constant value

            if ~isnan(NameValueArgs.P)
                P = linspace(obj.P0, NameValueArgs.P, obj.nPoints);
                T = obj.T0 .* (P ./ obj.P0) .^ (obj.R ./ (obj.cv + obj.R));
                V = obj.V0 .* (P ./ obj.P0) .^ (-obj.cv ./ (obj.cv + obj.R));
            elseif ~isnan(NameValueArgs.V)
                V = linspace(obj.V0, NameValueArgs.V, obj.nPoints);
                T = obj.T0 .* (V ./ obj.V0) .^ (-obj.R ./ obj.cv);
                P = obj.P0 .* (V ./ obj.V0) .^ (-(obj.cv + obj.R) ./ obj.cv);
            elseif ~isnan(NameValueArgs.T)
                T = linspace(obj.T0, NameValueArgs.T, obj.nPoints);
                P = obj.P0 .* (T ./ obj.T0) .^ ((obj.cv + obj.R) ./ obj.R);
                V = obj.V0 .* (T ./ obj.T0) .^ (-obj.cv ./ obj.R);
            end

            Esys = obj.Esys_net + obj.m .* obj.cv .* (T(end) - obj.T0);
            % Esys = Q - W
            Q = obj.Qnet;
            W = obj.Wnet + (Q - obj.Qnet) - (Esys - obj.Esys_net);

            obj = obj.set_initial_state(P, V, T, S, Q, W, Esys, NameValueArgs);
        end

        function obj = isobaric(obj, NameValueArgs)
            arguments
                obj
                NameValueArgs.V = NaN
                NameValueArgs.T = NaN
                NameValueArgs.S = NaN
                NameValueArgs.Qin = false
                NameValueArgs.Qout = false
                NameValueArgs.Win = false
                NameValueArgs.Wout = false
                NameValueArgs.stateName = ""
            end
            obj = obj.get_initial_state();
            P = linspace(obj.P0, obj.P0, obj.nPoints); % constant value

            if ~isnan(NameValueArgs.S)
                S = linspace(obj.Snet, NameValueArgs.S, obj.nPoints);
                T = obj.T0 .* exp((S - obj.Snet) ./ (obj.m .* (obj.cv + obj.R)));
                V = obj.V0 .* exp((S - obj.Snet) ./ (obj.m .* (obj.cv + obj.R)));
            elseif ~isnan(NameValueArgs.V)
                V = linspace(obj.V0, NameValueArgs.V, obj.nPoints);
                T = obj.get_T_from_PV(P, V);
                S = obj.Snet + obj.m .* (obj.cv + obj.R) .* log(V ./ obj.V0);
            elseif ~isnan(NameValueArgs.T)
                T = linspace(obj.T0, NameValueArgs.T, obj.nPoints);
                V = obj.get_V_from_PT(P, T);
                S = obj.Snet + obj.m .* (obj.cv + obj.R) .* log(T ./ obj.T0);
            end

            Esys = obj.Esys_net + obj.m .* obj.cv .* (T(end) - obj.T0);
            W = obj.Wnet + P(end) .* (V(end) - obj.V0); % constant volume; int(PdV) = 0J
            Q = obj.Qnet + (Esys - obj.Esys_net) + (W - obj.Wnet);

            obj = obj.set_initial_state(P, V, T, S, Q, W, Esys, NameValueArgs);
        end

        function obj = isothermal(obj, NameValueArgs)
            arguments
                obj
                NameValueArgs.V = NaN
                NameValueArgs.P = NaN
                NameValueArgs.S = NaN
                NameValueArgs.Qin = false
                NameValueArgs.Qout = false
                NameValueArgs.Win = false
                NameValueArgs.Wout = false
                NameValueArgs.stateName = ""
            end
            obj = obj.get_initial_state();
            T = linspace(obj.T0, obj.T0, obj.nPoints); % constant value

            if ~isnan(NameValueArgs.S)
                S = linspace(obj.Snet, NameValueArgs.S, obj.nPoints);
                P = obj.P0 .* exp(-(S - obj.Snet) ./ (obj.m .* obj.R));
                V = obj.V0 .* exp((S - obj.Snet) ./ (obj.m .* obj.R));
            elseif ~isnan(NameValueArgs.V)
                V = linspace(obj.V0, NameValueArgs.V, obj.nPoints);
                P = obj.get_P_from_TV(T, V);
                S = obj.Snet + obj.m .* obj.R .* log(V ./ obj.V0);
            elseif ~isnan(NameValueArgs.P)
                P = linspace(obj.P0, NameValueArgs.P, obj.nPoints);
                V = obj.get_V_from_PT(P, T);
                S = obj.Snet + obj.m .* obj.R .* -1 .* log(P ./ obj.P0);
            end

            Esys = obj.Esys_net + obj.m .* obj.cv .* (T(end) - obj.T0); % 0J
            % Esys = Q - W
            Q = obj.Qnet + obj.T0 .* (S(end) - obj.Snet);
            W = obj.Wnet + (Q - obj.Qnet) - (Esys - obj.Esys_net);

            obj = obj.set_initial_state(P, V, T, S, Q, W, Esys, NameValueArgs);
        end

        function obj = isochoric(obj, NameValueArgs)
            arguments
                obj
                NameValueArgs.T = NaN
                NameValueArgs.P = NaN
                NameValueArgs.S = NaN
                NameValueArgs.Qin = false
                NameValueArgs.Qout = false
                NameValueArgs.Win = false
                NameValueArgs.Wout = false
                NameValueArgs.stateName = ""
            end
            obj = obj.get_initial_state();
            V = linspace(obj.V0, obj.V0, obj.nPoints); % constant value

            if ~isnan(NameValueArgs.S)
                S = linspace(obj.Snet, NameValueArgs.S, obj.nPoints);
                T = obj.T0 .* exp((S - obj.Snet) ./ (obj.m .* obj.cv));
                P = obj.P0 .* exp((S - obj.Snet) ./ (obj.m .* obj.cv));
            elseif ~isnan(NameValueArgs.T)
                T = linspace(obj.T0, NameValueArgs.T, obj.nPoints);
                P = obj.get_P_from_TV(T, V);
                S = obj.Snet + obj.m .* obj.cv .* log(T ./ obj.T0);
            elseif ~isnan(NameValueArgs.P)
                P = linspace(obj.P0, NameValueArgs.P, obj.nPoints);
                T = obj.get_T_from_PV(P, V);
                S = obj.Snet + obj.m .* obj.cv .* log(P ./ obj.P0);
            end

            Esys = obj.Esys_net + obj.m .* obj.cv .* (T(end) - obj.T0);
            W = obj.Wnet; % constant volume; int(PdV) = 0J
            Q = obj.Qnet + (Esys - obj.Esys_net) + (W - obj.Wnet);

            obj = obj.set_initial_state(P, V, T, S, Q, W, Esys, NameValueArgs);
        end

        function PVTS(obj, P, V, T, S, currentStateName, prevStateName)
            figure(1);
            subplot(1, 2, 1);
            hold on;
            plot(S, T, 'g-');
            % plot(S(1), T(1), 'ks');
            text(S(1), T(1), prevStateName, HorizontalAlignment = 'center');
            % plot(S(end), T(end), 'ks');
            % text(S(end), T(end), currentStateName, HorizontalAlignment = 'center');
            xlabel("Entropy relative to S_0 (J/K)");
            ylabel("Absolute Temperature (K)");
            subplot(1, 2, 2);
            hold on;
            plot(V, P, 'b-');
            % plot(V(1), P(1), 'ro');
            text(V(1), P(1), prevStateName, HorizontalAlignment = 'center');
            % plot(V(end), P(end), 'ro');
            % text(V(end), P(end), currentStateName, HorizontalAlignment = 'center');
            xlabel("Volume (m^3)");
            ylabel("Pressure (Pa)");
        end

        function obj = set_initial_state(obj, P, V, T, S, Q, W, Esys, NameValueArgs)


            obj.currentStateNum = obj.currentStateNum + 1;
            if NameValueArgs.stateName == ""
                obj.currentStateName = "State " + obj.currentStateNum;
            else
                obj.currentStateName = NameValueArgs.stateName;
            end
            pathName = obj.prevStateName + " -> " + obj.currentStateName;
            obj.PVTS(P, V, T, S, obj.currentStateName, obj.prevStateName);
            obj.prevStateName = obj.currentStateName;

            Pnew = P(end);
            Vnew = V(end);
            Tnew = T(end);
            Snew = S(end);

            if NameValueArgs.Qin
                obj.Qin = obj.Qin + (Q - obj.Qnet);
            end
            if NameValueArgs.Qout
                obj.Qout = obj.Qout + (Q - obj.Qnet);
            end
            if NameValueArgs.Win
                obj.Win = obj.Win + (W - obj.Wnet);
            end
            if NameValueArgs.Wout
                obj.Wout = obj.Wout + (W - obj.Wnet);
            end

            obj.stateProperties.Properties.VariableNames = {'State',             'T',    'Snet', 'P',    'V',  'Esys_net', 'Qnet', 'Wnet'};
            obj.stateProperties = [obj.stateProperties;    {obj.currentStateName, Tnew,   Snew,   Pnew,   Vnew, Esys,       Q,      W     }];
            obj.pathProperties.Properties.VariableNames = {'Path',  'ΔT',          'ΔS',            'ΔP',          'ΔV',          'ΔE_sys',            'ΔQ',         'ΔW'         };
            obj.pathProperties = [obj.pathProperties;     {pathName, Tnew - obj.T0, Snew - obj.Snet, Pnew - obj.P0, Vnew - obj.V0, Esys - obj.Esys_net, Q - obj.Qnet, W - obj.Wnet}];

            obj.Esys_net = Esys;
            obj.P0 = Pnew;
            obj.T0 = Tnew;
            obj.V0 = Vnew;
            obj.Snet = Snew;
            obj.Qnet = Q;
            obj.Wnet = W;
        end
    end
end