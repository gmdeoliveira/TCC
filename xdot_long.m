%% Modelo longitudinal do dirigível
function xd = xdot_long(x, u)

%% Estados
U = x(1); % Velocidade linear no eixo x
W = x(2); % Velocidade linear no eixo y
theta = x(3); % Ângulo de arfagem
q = x(4); % Velocidade angular de arfagem
xe = x(5);
h = x(6);

%% Entradas
T = u(1); % Força de tração de cada motor (Newtons)
mi = u(2); % Ângulo da vetorização (graus)
de = u(3); % Deflexão do profundor (graus)
d2r = pi/180; % Conversão de graus para radianos
mi_r = mi*d2r;
de_r = de*d2r;

%% Parâmetros gerais
g = 9.81; % Aceleração da gravidade
rho = 1.225; % Densidade do ar

%% Parâmetros do dirigível
a1 = 3.13; %Semi-eixo menor
a2 = 3.37; %Semi-eixo maior
L = a1 + a2; %Comprimento do envelope
D = 1.65; %Maior diâmetro do envelope
Iy = 50.29; %Momento de inércia em torno do eixo Y
b = D/2; % Variável auxiliar
Vol = (2/3)*pi*(a1+a2)*(b^2); %Volume do envelope
m = rho*Vol; %Massa do dirigível
az = D/4; %Posição do CG no eixo Z
dz = D/2; %Posição dos motores no eixo Z

%% Matriz de massa
% Variáveis auxiliares
a = (a1 + a2)/2;
e = sqrt(1 - (b/a)^2);
a0 = 2*((1-e^2)/e^3)*(.5*log((1+e)/(1-e))-e);
b0 = 1/e^2 - ((1-e^2)/(2*e^3))*log((1+e)/(1-e));
% Fatores de massa e inércia virtual
k1 = a0/(2-a0);
k2 = b0/(2-b0);
k3 = ((e^4)*(b0-a0))/((2-e^2)*(2*e^2-(2-e^2)*(b0-a0)));
Vol = (2/3)*pi*(a1+a2)*(b^2); %Volume do envelope
m_ar = rho*Vol; % Massa de ar equivalente ao volume do envelope
Iy_ar = m_ar*(L^2 + D^2)/20; % Momento de inércia da massa de ar entorno do eixo Y
X_udot = -k1*m_ar; mx = m - X_udot;
Z_wdot = -k2*m_ar; mz = m - Z_wdot;
M_qdot = -k3*Iy_ar; Jy = Iy - M_qdot;
M = [mx 0 m*az; 0 mz 0; m*az 0 Jy]; %Matriz de massa

%% Forças e momento dinâmicos
f1 = -mz*W*q;
f2 = mx*q*U + m*az*q^2;
f3 = -az*W*q;
Fd = [f1; f2; f3];

%% Forças e momento aerodinâmicos
VT = sqrt(U^2 + W^2); %Velocidade total
alpha = atan(W/U); %Ângulo de ataque
%Coeficientes aerodinâmicos
CX1 = -0.4054;
CX2 = 1.0624;
CZ1 = 1.0624;
CZ2 = -1.3380;51
CZ3 = -5.5606;
CZ4 = -0.2896;
CM1 = -4.6638;
CM2 = -7.1851;
CM3 = -21.4411;
CM4 = -1.5549;
Cmq = -0.259;
Xa = .5*rho*VT^2*(CX1*cos(alpha)^2 + CX2*sin(2*alpha)*sin(alpha/2));
Za = .5*rho*VT^2*(CZ1*(cos(alpha/2))*sin(2*alpha) + CZ2*sin(2*alpha) + CZ3*sin(alpha
Ma = .5*rho*VT^2*(CM1*(cos(alpha/2))*sin(2*alpha) + CM2*sin(2*alpha) + CM3*sin(alpha
A = [Xa; Za; Ma];

%% Vetor de peso e empuxo
Wgt = m*g; %Peso do dirigível
B = rho*Vol*g; %Empuxo do dirigível
G1 = -(Wgt-B)*sin(theta);
G2 = (Wgt-B)*cos(theta);
G3 = -az*Wgt*sin(theta);
G = [G1; G2; G3];

%% Forças e momento propulsivos
Xp = 2*T*cos(mi_r);
Zp = -2*T*sin(mi_r);
Mp = 2*T*dz*cos(mi_r);
P = [Xp; Zp; Mp];

%% Equações de trajetória
gamma = theta - alpha;
xe_dot = VT*cos(gamma);
hdot = VT*sin(gamma);
thetadot = q; % Derivada do ângulo de arfagem
uwqdot = M\(A + G + P + Fd);
udot = uwqdot(1);
wdot = uwqdot(2);
qdot = uwqdot(3);
xd = [udot; wdot; thetadot; qdot; xe_dot; hdot];

end
