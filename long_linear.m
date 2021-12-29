clear all; close all
%%Linearização do modelo longitudinal
% Executar o código long_trim.m antes deste código
r2d = 180/pi;
global utrim xtrim

% Condições iniciais dos estados
xref = xtrim;

% Condições iniciais das entradas
uref = utrim;

% Função utilizada para a linearização em torno do ponto de equilíbrio
fx = @(x) xdot_long(x,uref);
fu = @(u) xdot_long(xref,u);

% Matrizes jacobianas
A = jac(fx,xref);
B = jac(fu,uref);

% Matriz reduzida p/ o controlador LQR, excluindo a vetorização
B_lqr = [B(:,1), B(:,3)];

% Calculo da matriz de ganhos do controlador do tipo LQR
Q = eye(6); % Matriz de pesos dos estados
Q(6,6) = 1000;
Q(3,3) = 10;
Q(4,4) = 10;
R = eye(2); % Matriz de pesos dos controles
R(1,1) = 0.01;
R(2,2) = 10;
[K_lqr] = lqr(A,B_lqr,Q,R)
