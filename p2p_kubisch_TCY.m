function [ S, dot_S, ddot_S, T ] = p2p_kubisch_TCY( W_stuetz, T_ges, delta_T )
% Erzeugt aus N_I Stuetzpunkten in W_stuetz je mit Anfangs- und Endpunkt N_I-1 Trajektorien je in Form eines kubischen Polynoms
% S         := Trajektorie auf Positionsebene
% dot_S     := Trajektorie auf Geschwindigkeitsebene
% ddot_S    := Trajektorie auf Beschleunigungsebene
% T         := Zeitvektor der Trajektorie

% W_stuetz  := Stuetzpunkte
% T_ges     := Dauer der Bewegung/Interpolation
% delta_T   := Taktzeit

% Anzahl der Freiheitsgrade
N_Q       = size( W_stuetz,1 );

% Anzahl der Stuetzpunkte
N_I       = size( W_stuetz,2 );

% Zeitintervall fuer Interpolation
T_I       = 0:delta_T:(T_ges/(N_I-1));  % Zeitintervall fuer ein Teilstueck
N_T_I     = length(T_I);          % Anzahl der Zeitpunkte eines Teilstuecks

%% Berechnung der Trajektorie

% Initialisierung fuer Teilstuecke
S_I       = zeros( N_Q, N_T_I );
dot_S_I   = zeros(size(S_I));
ddot_S_I  = zeros(size(S_I));

% Initialisierung fuer Gesamttrajektorie
S         = zeros( N_Q, (N_T_I-1)*(N_I-1)+1 );
dot_S     = zeros(size(S));
ddot_S    = zeros(size(S));
T         = 0:delta_T:((N_I-1)*T_I(N_T_I));

%% --- ARBEITSBEREICH: ------------------------------------------------
% Erzeuge Trajektorie 

for i=1:1:N_I-1   % Schleife ueber Stuetzpunktepaare
    
    for m=1:3;   %Schleife ueber Koordinaten x,y,z
        
        % Erzeuge 3.Ordnung also 4x4 Matrix V
        T_ab=[T_I(1) T_I(N_T_I)]; % Start- und Endzeit fuer Teilstueck
        n=length(T_ab);
        T_ab=T_ab(:);
        v=ones(n,1);
        V(1:2*n,2*n)=[v;zeros(n,1)];
    
        for j = 1:2*n-1
            V(n+1:end,2*n-j)=j*v;
            v=T_ab.*v;
            V(1:n,2*n-j)=v;
        end
        
        % Koeffizientvektor fuer stueckweise Interpolation mit Randbedingung
        H = (V\[W_stuetz(m,i);W_stuetz(m,i+1);0;0]).';  
        
        % Berechnung der Trajektorie fuer Teilstuecke
        for p=1:4
            delta_S_I=H(p)*T_I.^(4-p);
            S_I(m,:)=S_I(m,:)+delta_S_I;
        end
        % Berechnung der Geschwingdigkeit fuer Teilstuecke
        for p=1:3
            delta_dot_S_I=(4-p)*H(p)*T_I.^(3-p);
            dot_S_I(m,:)=dot_S_I(m,:)+delta_dot_S_I;
        end
        % Berechnung der Beschleunigung fuer Teilstuecke
        for p=1:2
            delta_ddot_S_I=(4-p)*(3-p)*H(p)*T_I.^(2-p);
            ddot_S_I(m,:)=ddot_S_I(m,:)+delta_ddot_S_I;
        end
    
    end
    
    
    % Berechnung der Trajektorie fuer Gesamttrajektorie
    S(1:3,(i-1)*54+1:i*54+1)=S_I;
    % Berechnung der Geschwindigkeit fuer Gesamttrajektorie
    dot_S(1:3,(i-1)*54+1:i*54+1)=dot_S_I;
    % Berechnung der Beschleunigung fuer Gesamttrajektorie
    ddot_S(1:3,(i-1)*54+1:i*54+1)=ddot_S_I;

    % Rueckstellung fuer Teilstuecke
    S_I       = zeros( N_Q, N_T_I );
    dot_S_I   = zeros(size(S_I));
    ddot_S_I  = zeros(size(S_I));
end

%% --- ENDE ARBEITSBEREICH --------------------------------------------
end % function
