clear all
%% Programa de trimagem do modelo longitudinal
global xtrim utrim 

dec0 = zeros(3, 1);
dectrim = fsolve(@trim_fun, dec0);
xtrim = dec2x(dectrim);
utrim = dec2u(dectrim);

function x = dec2x(dec)
alpha = pi/180; % VARIAVEL DE DECISÃO
gamma = 0; % VARIAVEL DE DECISÃO (VALOR PADRÃO: ZERO!)
theta =  alpha + gamma;
VT = 6; % VARIAVEL DE DECISÃO
U = VT*cos(alpha);
W = VT*sin(alpha);
q = 0;
xe = 0; % Nao importa
h = 0;  % Nao importa

x = [U; W; theta; q; xe; h];
end

function u = dec2u(dec)
T = dec(1);
mi = dec(2);
de = dec(3);

u = [T; mi; de];
end

function err = xdot2err(xdot)
err = xdot([1,2,4]);
end

function err = trim_fun(dec)
x = dec2x(dec);
u = dec2u(dec);
xdot = xdot_long(x, u);
err = xdot2err(xdot);
end