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
        S0 {mustBeNumeric}

        T {mustBeNumeric}
        P {mustBeNumeric}
        V {mustBeNumeric}
        S {mustBeNumeric}

        Esys0 {mustBeNumeric}
        Esys {mustBeNumeric}
        Q0 {mustBeNumeric}
        Qnet {mustBeNumeric}
        W0 {mustBeNumeric}
        Wnet {mustBeNumeric}
        Strans {mustBeNumeric}
        Sgen {mustBeNumeric}
        cp {mustBeNumeric}
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

            obj.Esys0 = 0; % J
            obj.Esys = NaN;
            obj.Q0 = 0;
            obj.Qnet = NaN;
            obj.W0 = 0;
            obj.Wnet = NaN;
            obj.Strans = NaN;
            obj.Sgen = NaN;
            obj.cp = NaN; % will be re-assigned upon being called
            obj.T0 = NaN;
            obj.P0 = NaN;
            obj.V0 = NaN;
            obj.S0 = 0; % J/kg/K
            obj.T = NaN;
            obj.P = NaN;
            obj.V = NaN;
            obj.S = NaN;

            figure(1);
            clf;
        end

        function value = get.cp(obj)
            value = obj.cv + obj.R;
        end

        function value = get_T_from_PV(obj, P, V)
            value = P .* V ./ (obj.m .* obj.R);
        end

        function value = get_P_vrom_TV(obj, T, V)
            value = obj.m .* obj.R .* T ./ V;
        end

        function value = get_V_from_PT(obj, P, T)
            value = obj.m .* obj.R .* T ./ P;
        end

        function value = get_T0(obj)
            if isnan(obj.T0)
                value = obj.P0 .* obj.V0 ./ (obj.m ./ obj.R);
            else
                value = obj.T0
            end
        end

        function value = get_P0(obj)
            if isnan(obj.P0)
                value = obj.m .* obj.R .* obj.T0 ./ obj.V0;
            else
                value = obj.P0
            end
        end

        function value = get_V0(obj)
            if isnan(obj.V0)
                value = obj.m .* obj.R .* obj.T0 ./ obj.P0;
            else
                value = obj.V0
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
            end
            obj = obj.get_initial_state()
            S = linspace(obj.S0, obj.S0, obj.nPoints); % constant value

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

            Esys = obj.Esys0 + obj.m .* obj.cv .* (T(end) - obj.T0);
            % Esys = Q - W
            Q = obj.Q0;
            W = obj.W0 + (Q - obj.Q0) - (Esys - obj.Esys0);

            obj = obj.set_initial_state(P, V, T, S, Q, W, Esys);
        end

        function obj = isobaric(obj, NameValueArgs)
            arguments
                obj
                NameValueArgs.V = NaN
                NameValueArgs.T = NaN
                NameValueArgs.S = NaN
            end
            obj = obj.get_initial_state()
            P = linspace(obj.P0, obj.P0, obj.nPoints); % constant value

            if ~isnan(NameValueArgs.S)
                S = linspace(obj.S0, NameValueArgs.S, obj.nPoints);
                T = obj.T0 .* exp((S - obj.S0) ./ (obj.m .* (obj.cv + obj.R)));
                V = obj.V0 .* exp((S - obj.S0) ./ (obj.m .* (obj.cv + obj.R)));
            elseif ~isnan(NameValueArgs.V)
                V = linspace(obj.V0, NameValueArgs.V, obj.nPoints);
                T = obj.get_T_from_PV(P, V);
                S = obj.S0 + obj.m .* (obj.cv + obj.R) .* log(V ./ obj.V0);
            elseif ~isnan(NameValueArgs.T)
                T = linspace(obj.T0, NameValueArgs.T, obj.nPoints);
                V = obj.get_V_from_PT(P, T);
                S = obj.S0 + obj.m .* (obj.cv + obj.R) .* log(T ./ obj.T0);
            end

            Esys = obj.Esys0 + obj.m .* obj.cv .* (T(end) - obj.T0);
            W = obj.W0 + P(end) .* (V(end) - obj.V0); % constant volume; int(PdV) = 0J
            Q = obj.Q0 + (Esys - obj.Esys0) + (W - obj.W0);

            obj = obj.set_initial_state(P, V, T, S, Q, W, Esys);
        end

        function obj = isothermal(obj, NameValueArgs)
            arguments
                obj
                NameValueArgs.V = NaN
                NameValueArgs.P = NaN
                NameValueArgs.S = NaN
            end
            obj = obj.get_initial_state()
            T = linspace(obj.T0, obj.T0, obj.nPoints); % constant value

            if ~isnan(NameValueArgs.S)
                S = linspace(obj.S0, NameValueArgs.S, obj.nPoints);
                P = obj.P0 .* exp(-(S - obj.S0) ./ (obj.m .* obj.R));
                V = obj.V0 .* exp((S - obj.S0) ./ (obj.m .* obj.R));
            elseif ~isnan(NameValueArgs.V)
                V = linspace(obj.V0, NameValueArgs.V, obj.nPoints);
                P = obj.get_P_from_TV(T, V);
                S = obj.S0 + obj.m .* obj.R .* log(V ./ obj.V0);
            elseif ~isnan(NameValueArgs.P)
                P = linspace(obj.P0, NameValueArgs.P, obj.nPoints);
                V = obj.get_V_from_PT(P, T);
                S = obj.S0 + obj.m .* obj.R .* -1 .* log(P ./ obj.P0);
            end

            Esys = obj.Esys0 + obj.m .* obj.cv .* (T(end) - obj.T0); % 0J
            % Esys = Q - W
            Q = obj.Q0 + obj.T0 .* (S(end) - obj.S0);
            W = obj.W0 + (Q - obj.Q0) - (Esys - obj.Esys0);

            obj = obj.set_initial_state(P, V, T, S, Q, W, Esys);
        end

        function obj = isochoric(obj, NameValueArgs)
            arguments
                obj
                NameValueArgs.T = NaN
                NameValueArgs.P = NaN
                NameValueArgs.S = NaN
                NameValueArgs.Tinterface = NaN
            end
            obj = obj.get_initial_state()
            V = linspace(obj.V0, obj.V0, obj.nPoints); % constant value

            if ~isnan(NameValueArgs.S)
                S = linspace(obj.S0, NameValueArgs.S, obj.nPoints);
                T = obj.T0 .* exp((S - obj.S0) ./ (obj.m .* obj.cv));
                P = obj.P0 .* exp((S - obj.S0) ./ (obj.m .* obj.cv));
            elseif ~isnan(NameValueArgs.T)
                T = linspace(obj.T0, NameValueArgs.T, obj.nPoints);
                P = obj.get_P_from_TV(T, V);
                S = obj.S0 + obj.m .* obj.cv .* log(T ./ obj.T0);
            elseif ~isnan(NameValueArgs.P)
                P = linspace(obj.P0, NameValueArgs.P, obj.nPoints);
                T = obj.get_T_from_PV(P, V);
                S = obj.S0 + obj.m .* obj.cv .* log(P ./ obj.P0);
            end

            Esys = obj.Esys0 + obj.m .* obj.cv .* (T(end) - obj.T0);
            W = obj.W0; % constant volume; int(PdV) = 0J
            Q = obj.Q0 + (Esys - obj.Esys0) + (W - obj.W0);

            obj = obj.set_initial_state(P, V, T, S, Q, W, Esys);
        end

        function PVTS(obj, P, V, T, S)
            [P; T; V; S].' % print
            figure(1);
            subplot(1, 2, 1);
            hold on;
            plot(S, T, 'g-');
            plot(S(1), T(1), 'ks');
            plot(S(end), T(end), 'ks');
            xlabel("Entropy relative to S_0 (J/K)");
            ylabel("Absolute Temperature (K)");
            subplot(1, 2, 2);
            hold on;
            plot(V, P, 'b-');
            plot(V(1), P(1), 'ro');
            plot(V(end), P(end), 'ro');
            xlabel("Volume (m^3)");
            ylabel("Pressure (Pa)");
        end

        function obj = set_initial_state(obj, P, V, T, S, Q, W, Esys)
            obj.PVTS(P, V, T, S);
            obj.P = P(end);
            obj.S = S(end);
            obj.T = T(end);
            obj.V = V(end);
            obj.Esys0 = Esys;
            obj.Qnet = Q; % for external display
            obj.Wnet = W; % for external display

            obj.P0 = obj.P;
            obj.T0 = obj.T;
            obj.V0 = obj.V;
            obj.S0 = obj.S;
            obj.Q0 = Q;
            obj.W0 = W;
        end
    end
end