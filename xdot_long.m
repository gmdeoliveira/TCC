function xd = xdot_long(x, u)
%% Modelo longitudinal do dirig�vel

%% Estados
U = x(1); % Velocidade linear no eixo x
W = x(2); % Velocidade linear no eixo y
theta = x(3); % �ngulo de arfagem
q = x(4); % Velocidade angular de arfagem
xe = x(5);
h = x(6);

%% Entradas
T = u(1); % Tra��o fornecida por cada motor
% Limite da tra��o de cada motor
% if T > 20
%     T = 20;
% end

mi = u(2); % �ngulo da vetoriza��o (graus)
de = u(3); % Deflex�o do profundor (graus)

% Batente superior do profundor
% if de >= 28
%     de = 28;
% end
% % Batente inferior do profundor
% if de <= -28
%     de = -28;
% end

d2r = pi/180; % Convers�o de graus para radianos
mi_r = mi*d2r;
de_r = de*d2r;

%% Par�metros gerais
g = 9.81; % Acelera��o da gravidade
rho = 1.225; % Densidade do ar

%% Par�metros do dirig�vel
a1 = 3.13; %Semi-eixo menor
a2 = 3.37; %Semi-eixo maior 
L = a1 + a2; %Comprimento do envelope
D = 1.65; %Maior di�metro do envelope
Iy = 50.29; %Momento de in�rcia em torno do eixo Y 
b = D/2; % Vari�vel auxiliar
Vol = (2/3)*pi*(a1+a2)*(b^2); %Volume do envelope
m = rho*Vol; %Massa do dirig�vel  
az = D/4; %Posi��o do CG no eixo Z
dz = D/2; %Posi��o dos motores no eixo Z

%% Matriz de massa
% Vari�veis auxiliares 
a = (a1 + a2)/2;
e = sqrt(1 - (b/a)^2);
a0 = 2*((1-e^2)/e^3)*(.5*log((1+e)/(1-e))-e);
b0 = 1/e^2 - ((1-e^2)/(2*e^3))*log((1+e)/(1-e));

% Fatores de massa e in�rcia virtual
k1 = a0/(2-a0);
k2 = b0/(2-b0);
k3 = ((e^4)*(b0-a0))/((2-e^2)*(2*e^2-(2-e^2)*(b0-a0)));

Vol = (2/3)*pi*(a1+a2)*(b^2); %Volume do envelope
m_ar = rho*Vol; % Massa de ar equivalente ao volume do envelope
Iy_ar = m_ar*(L^2 + D^2)/20; % Momento de in�rcia da massa de ar entorno do eixo Y

X_udot = -k1*m_ar;  mx = m - X_udot;
Z_wdot = -k2*m_ar;  mz = m - Z_wdot;
M_qdot = -k3*Iy_ar; Jy = Iy - M_qdot; 

M = [mx 0 m*az; 0 mz 0; m*az 0 Jy]; %Matriz de massa

%% For�as e momento din�micos
f1 = -mz*W*q;
f2 = mx*q*U + m*az*q^2;
f3 = -az*W*q;
Fd = [f1; f2; f3];

%% For�as e momento aerodin�micos
VT = sqrt(U^2 + W^2); %Velocidade total
alpha = atan(W/U); %�ngulo de ataque

%Coeficientes aerodin�micos
CX1 = -0.4054;
CX2 = 1.0624;

CZ1 = 1.0624;
CZ2 = -1.3380;
CZ3 = -5.5606;
CZ4 = -0.2896;

CM1 = -4.6638;
CM2 = -7.1851;
CM3 = -21.4411;
CM4 = -1.5549;
Cmq = -0.259;

Xa = .5*rho*VT^2*(CX1*cos(alpha)^2 + CX2*sin(2*alpha)*sin(alpha/2));
Za = .5*rho*VT^2*(CZ1*(cos(alpha/2))*sin(2*alpha) + CZ2*sin(2*alpha) + CZ3*sin(alpha)*sin(abs(alpha)) + 2*CZ4*de_r);
Ma = .5*rho*VT^2*(CM1*(cos(alpha/2))*sin(2*alpha) + CM2*sin(2*alpha) + CM3*sin(alpha)*sin(abs(alpha)) + 2*CM4*de_r + L*Cmq*q);
A = [Xa; Za; Ma];

%% Vetor de peso e empuxo
Wgt = m*g; %Peso do dirig�vel
B = rho*Vol*g; %Empuxo do dirig�vel

G1 = -(Wgt-B)*sin(theta);
G2 = (Wgt-B)*cos(theta);
G3 = -az*Wgt*sin(theta);
G = [G1; G2; G3];

%% For�as e momento propulsivos
Xp = 2*T*cos(mi_r);
Zp = -2*T*sin(mi_r);
Mp = 2*T*dz*cos(mi_r);
P = [Xp; Zp; Mp];

%% Equa��es de trajet�ria 
gamma = theta - alpha;
xe_dot = VT*cos(gamma);
hdot = VT*sin(gamma);

thetadot = q; % Derivada do �ngulo de arfagem

uwqdot = M\(A + G + P + Fd);
udot = uwqdot(1);
wdot = uwqdot(2);
qdot = uwqdot(3);

xd = [udot; wdot; thetadot; qdot; xe_dot; hdot];
end