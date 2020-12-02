clear ; close all; clc

load('pattern_letters.mat');
load('answers_data.mat');
Q = size(P, 1);
P = P';

n1 = 40;  %Numero de neuronas en la capa oculta
ep = 1;   % Ventana de valores iniciales
% Valores iniciales
W1 = ep*(2*rand(n1,63)-1);
b1 = ep*(2*rand(n1,1)-1);
W2 = ep*(2*rand(4,n1)-1);
b2 = ep*(2*rand(4,1)-1);
alfa = 0.01;

a2 = zeros(12,4);
ec = zeros(1,12);
for Epocas = 1:600
    for q = 1:Q
        % Propagación de la entrada hacia la salida
        a1 = tg_sig(W1*P(:,q) + b1);
        a2 = tg_sig(W2*a1 + b2);
        % Retropropagación de la sensibilidades
        e = T(:,q) - a2;
        s2 = -2*diag(1-a2.^2)*e;
        s1 = diag(1-a1.^2)*W2'*s2;
         % Actualización de pesos sinapticos y polarizaciones
        W2 = W2 - alfa*s2*a1';
        b2 = b2 - alfa*s2;
        W1 = W1 - alfa*s1*P(:,q)';
        b1 = b1 - alfa*s1; 
        % Error Cuadrático
        ec(q) = e'*e;
    end
    ecm(Epocas)= sum(ec)/Q;
    if mod(Epocas,50) == 0
        Epocas
        plot(ecm)
        pause(0.1)
    end
end

a = zeros(12,4);
% Salidas de la RNA de verificacion.
for q=1:Q
    a1 = tg_sig(W1*P(:,q) + b1);
    a2 = hardlims (W2*a1 + b2);
    b = size(a2, 1);
    for i = 1:b
        a(q, i) = a2(i);
    end
end
a