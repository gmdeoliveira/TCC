clear all; close all
%% Linearização do modelo longitudinal
run('C:\Users\gmdeo\Documents\Gabriel\UFMG\TCC\Modeling an RC airship\long_trim.m')
r2d = 180/pi;
global utrim xtrim %A_long B_long
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

% Matrizes reduzidas
A_long = A(1:4,1:4); 
% for i = 1:5
%     A_long(5,i) = A(6,i);
% end
% B_long = B(1:4,:);

% % M - Autovetores | D - Autovalores
[M,D] = eig(A_long); 
damp(A_long)

M_inv = inv(M);

% Variáveis auxiliares
r1 = M(1,:);
r2 = M(2,:);
r3 = M(3,:);
r4 = M(4,:);

c1 = M_inv(:,1);
c2 = M_inv(:,2);
c3 = M_inv(:,3);
c4 = M_inv(:,4);

C1 = zeros(4,4);
for i = 1:4
    C1(i,i) = c1(i);
end

C2 = zeros(4,4);
for i = 1:4
    C2(i,i) = c1(i);
end

C3 = zeros(4,4);
for i = 1:4
    C3(i,i) = c1(i);
end

C4 = zeros(4,4);
for i = 1:4
    C4(i,i) = c1(i);
end

S = zeros(4,4); % Matriz de sensibilidade
S(1,:) = r1*C1;
S(2,:) = r2*C2;
S(3,:) = r3*C3;
S(4,:) = r4*C4;

% Módulo dos números reais e complexos
for i = 1:4
    for j = 1:4
        S(i,j) = abs(S(i,j));
    end
end

% Soma dos valores por linha para normalização
soma = zeros(4,1);
for i = 1:4
    for j = 1:4
        soma(i) = S(i,j) + soma(i);
    end
end

% Normalização
for i = 1:4
    for j = 1:4
        S(i,j) = S(i,j)/soma(i);
    end
end
% 
% %%% Funcoes de transferencia
% % Em relacao a velocidade linear U
% c_u = [1 0 0 0];
% d = [0 0 0];
% 
% [num_T,den_T] = ss2tf(A_long,B_long,c_u,d,1);
% u_T = tf(num_T,den_T);
% [Z1,P1,K1] = tf2zpk(num_T,den_T);
% %pidTuner(u_T)
% 
% % Em relacao ao angulo de arfagem
% c_theta = [0 0 1 0];
% [num_theta,den_theta] = ss2tf(A_long,B_long,c_theta,d,3);
% theta_de = tf(num_theta,den_theta);
% [Z2,P2,K2] = tf2zpk(num_theta,den_theta);
% 
% % Em relacao a velocidade de arfagem
% c_q = [0 0 0 1];
% [num_q,den_q] = ss2tf(A_long,B_long,c_q,d,1);
% q_de = tf(num_q,den_q);
% [Z3,P3,K3] = tf2zpk(num_q,den_q);
% 
% % % Em relacao a altitude
% c_h = [0 0 0 0 0 1];
% [num_h,den_h] = ss2tf(A,B,c_h,d,3);
% h_de = tf(num_h,den_h);
% [Z4,P4,K4] = tf2zpk(num_h,den_h);

% Calculo da matriz de ganhos do controlador do tipo LQR
%Q = eye(6); % Matriz de pesos dos estados (ver nomenclatura correta)
%Q(6,6) = 1000; 
%Q(3,3) = 10;
%Q(4,4) = 10;
%R = eye(2); % Matriz de pesos dos controles (ver nomenclatura correta)
%R(1,1) = 0.01;
%R(2,2) = 10;
%[K_lqr] = lqr(A,B_lqr,Q,R)