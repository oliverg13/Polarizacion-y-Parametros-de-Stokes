clear, clc, close all, format compact;

%% Definir campo óptico gaussiano
n=7; % Resolución 2^n puntos
w0=2; % Anchura
m =0;  % Carga topologica
lado=linspace(-3*w0,3*w0,2^n); % se considera un sensor cuadrado de lado 2*(3*w0)
[x,y]= meshgrid(lado);
[phi,r]=cart2pol(x,y);
E=(r/w0).^abs(m).*exp(-(r/w0).^2).*exp(1i*m*phi);

figure(1)
hold on
surf(x/w0,y/w0,real(E)), shading interp, axis square, colormap(hot), view(25,30)
xlabel('X'), ylabel('Y'), zlabel('Z'), title("Perfil de amplitud para m=" + num2str(m))
hold off

%% Definición de estados de polarización en vectores de Jones
H=[1;0]; % horizontal
V=[0;1]; % vertical
R=[1;-1i]/sqrt(2); % circular derecho
L=[1;1i]/sqrt(2); % circular izquierdo
D=[1;1]/sqrt(2); % diagonal
AD=[1;-1]/sqrt(2); % antidiagonal

% cualquier otra combinación 
E0=1;% Amplitud
angulo_de_inclinacion=pi/5;
E0x=E0*cos(angulo_de_inclinacion); %
E0y=E0*sin(angulo_de_inclinacion); % 
delta=pi/3; % FaseEnY - FaseEnX
ep= [E0x/E0;E0y*exp(1i*delta)/E0]; % vector de polarización del campo

%% Polarización utilizada para este trabajo
ep=H; 
%ep=H;
Ex=E*ep(1);
Ey=E*ep(2);

%% Grafica vectorial de estados de polarización
figure(2)
hold on
quiver(x,y,Ex,Ey)
xlabel('x'), ylabel('y'), title("Grafica vectorial de estados de polarización")
hold off

%% Definición de un polarizador lineal (matriz de Jones)
na=16;%Numero de angulos diferentes de inclinación del polarizador
a_pl=linspace(0,pi,na); %Vector de angulos de inclinación dle polarizador

I=zeros(1,na); %Preasignación de memoria
for  c=1:na
    %a_pl(c) c-esimo ángulo de inclinación del polarizador lineal
    ep_pl= [cos(a_pl(c));sin(a_pl(c))]; % vector de polarización del polarizador
    M_pl=(ep_pl)*(ep_pl'); % matriz de polarizador lineal angulo a_pl
    Ex_out=Ex*M_pl(1,1)+Ey*M_pl(1,2);
    Ey_out=Ex*M_pl(2,1)+Ey*M_pl(2,2);
    E_out=sqrt(abs(Ex_out).^2 + abs(Ey_out).^2);

%% Validación de ley de Malus
h=(max(lado)-min(lado))/(2^n-1);
I(c)=sum(sum(E_out.^2))*h^2;
end



I0=sum(sum(abs(E).^2))*h^2;
figure(3)
hold on
theta=linspace(0,pi,1000);
plot(a_pl,I/I0,theta,(cos(theta)).^2)
%plot(theta,(cos(theta)).^2,a_pl,I/I0)
xlabel('Angulo del polarizador \theta [rad]')
ylabel('I_{resultante}/I_0')
title('Verificación de Ley de Malus')
legend(["Valor númerico" "Resultado teorico esperado"],'Location','northeast')
grid on
hold off

%% De aquí en adelante, se utilizará un estado de polarización vertical
ep=V;
Ex=E*ep(1);
Ey=E*ep(2);

%% Definición de un retardador de 1/4 de lambda
delta=pi/2; %desfase introducido de y - x
%M=[1 0; 0 exp(-1i*delta)];
M=[exp(+1i*delta/2) 0; 0 exp(-1i*delta/2)];

%% Ajustar retardador de 1/4 de onda con polarizador lineal para minimizar variación de la amplitud de la señal de potencia
q=40; %Numero de angulos diferentes de inclinación
phi=linspace(0,pi,q); %angulo de inclinación del retardador 
[k1, k2]=meshgrid(phi);
for cc=1:q %contador para angulo de rotacion del polarizador lineal
for c=1:q %contador para angulo de rotacion de retardador
    rotar_alfa=[cos(k1(c,cc)) sin(k1(c,cc)); -sin(k1(c,cc)) cos(k1(c,cc))]; %rotacion horaria si alfa>0
    ep_out=(rotar_alfa^-1)*M*rotar_alfa*ep;
    Ex_out=E*ep_out(1);
    Ey_out=E*ep_out(2);
    
    %% Utilizando una proyección para el polariador lineal
    Up=[cos(k2(c,cc)) sin(k2(c,cc))];
    Proy = Ex_out*Up(1) + Ey_out*Up(2);
    Ex_out = Proy*Up(1);
    Ey_out = Proy*Up(2);
    
    E_out=sqrt(abs(Ex_out).^2 + abs(Ey_out).^2);
    I(c,cc)=sum(sum(E_out.^2))*h^2;
end
end

%% Grafica 3D cuando el polarizador y el retardador se rotan simultaneamente
figure (14)
hold on
surf(k1,k2,I/I0), shading interp, axis square, colormap(hot),view(2)
xlabel('Angulo del retardador'), ylabel('Angulo de polarizador'), zlabel('Intensidad')
hold off

%% Definición de retardador de 1/2 de lambda
delta=pi; %desfase introducido de y - x
%M=[1 0; 0 exp(-1i*delta)];
M=[exp(+1i*delta/2) 0; 0 exp(-1i*delta/2)];

%% Ajustar retardador de 1/2 de onda con polarizador lineal para minimizar variación de la amplitud de la señal de potencia
q=40; %Numero de angulos diferentes de inclinación
phi=linspace(0,pi,q); %angulo de inclinación del retardador 
[k1, k2]=meshgrid(phi);
for cc=1:q %contador para angulo de rotacion del polarizador lineal
for c=1:q %contador para angulo de rotacion de retardador
    rotar_alfa=[cos(k1(c,cc)) sin(k1(c,cc)); -sin(k1(c,cc)) cos(k1(c,cc))]; %rotacion horaria si alfa>0
    ep_out=(rotar_alfa^-1)*M*rotar_alfa*ep;
    Ex_out=E*ep_out(1);
    Ey_out=E*ep_out(2);
    
    %% Utilizando una proyección para el polariador lineal
    Up=[cos(k2(c,cc)) sin(k2(c,cc))]; 
    Proy = Ex_out*Up(1) + Ey_out*Up(2);
    Ex_out = Proy*Up(1);
    Ey_out = Proy*Up(2);
    
    E_out=sqrt(abs(Ex_out).^2 + abs(Ey_out).^2);
    I(c,cc)=sum(sum(E_out.^2))*h^2;
end
end

%% Grafica 3D cuando el polarizador y el retardador se rotan simultaneamente
figure (15)
hold on
surf(k1,k2,I/I0), shading interp, axis square, colormap(hot),view(10,45)
xlabel('Angulo del retardador'), ylabel('Angulo de polarizador'), zlabel('Intensidad')
hold off

%% Retardador de media onda para polarizacion circular
%% Definición de retardador de 1/2 de lambda
delta=pi; %desfase introducido de y - x
%M=[1 0; 0 exp(-1i*delta)];
M=[exp(+1i*delta/2) 0; 0 exp(-1i*delta/2)];

%% Ajustar retardador de 1/2 de onda con polarizador lineal para minimizar variación de la amplitud de la señal de potencia
q=40; %Numero de angulos diferentes de inclinación
phi=linspace(0,pi,q); %angulo de inclinación del retardador 
[k1, k2]=meshgrid(phi);
I=zeros(q,q); S0=I; S1=I; S2=I; S3=I;
ep=R;
for cc=1:q %contador para angulo de rotacion del polarizador lineal
for c=1:q %contador para angulo de rotacion de retardador
    rotar_alfa=[cos(k1(c,cc)) sin(k1(c,cc)); -sin(k1(c,cc)) cos(k1(c,cc))]; %rotacion horaria si alfa>0
    ep_out=(rotar_alfa^-1)*M*rotar_alfa*ep;
    Ex_out=E*ep_out(1);
    Ey_out=E*ep_out(2);
    
    %% Utilizando una proyección para el polariador lineal
    Up=[cos(k2(c,cc)) sin(k2(c,cc))];    
    Proy = Ex_out*Up(1) + Ey_out*Up(2);
    Ex_out = Proy*Up(1);
    Ey_out = Proy*Up(2);
    
    %% Proyeccion con matrices
    %rotar_beta=[cos(k2(c,cc)) sin(k2(c,cc)); -sin(k2(c,cc)) cos(k2(c,cc))];
    %ep_out=(rotar_beta^-1)*[1;0].*(rotar_beta*(rotar_alfa^-1)*M*rotar_alfa*ep);
    %Ex_out=E*ep_out(1);
    %Ey_out=E*ep_out(2);
    
    %M_pl=(ep_pl)*(ep_pl'); %matriz de polarizador lineal angulo a_pl
    %ep_out=M_pl*ep_out; 
    %Ex_out=E*ep_out(1);
    %Ey_out=E*ep_out(2);
    
    E_out=sqrt(abs(Ex_out).^2 + abs(Ey_out).^2);
    I(c,cc)=sum(sum(E_out.^2))*h^2;
    
    %%Parametros de stokes respecto k1,k2 (dependiente de 2 dimensiones)
    S0(c,cc)=sum(sum(abs(Ex_out).^2))*h^2+sum(sum(abs(Ey_out).^2))*h^2;
    S1(c,cc)=sum(sum(abs(Ex_out).^2))*h^2-sum(sum(abs(Ey_out).^2))*h^2;    %calculo de epsilon = fase_y-fase_x
    %efx=angle(Ex_out);
    %efy=angle(Ey_out);
    %ef=efy-efx;
    S21(c,cc)=2*sum(sum(real(Ex_out.*conj(Ey_out))))*h^2;
    S31(c,cc)=-2*sum(sum(imag(Ex_out.*conj(Ey_out))))*h^2;
     S2(c,cc)=2*sum(sum(real(Ey_out.*conj(Ex_out))))*h^2;
    S3(c,cc)=2*sum(sum(imag(Ey_out.*conj(Ex_out))))*h^2;
    
end
end
figure(50)
surf(k1,k2,S31-S3)
figure(51)
surf(k1,k2,S21-S2)
%% Grafica 3D cuando el polarizador y el retardador se rotan simultaneamente
figure (16)
hold on
surf(k1,k2,I/I0), shading interp, axis square, colormap(hot),view(10,45)
xlabel('Angulo del retardador'), ylabel('Angulo de polarizador'), zlabel('Intensidad')
hold off

%% Grafica 3D de grado P
P=sqrt(S1.^2+S2.^2+S3.^2)./S0;
figure (17)
hold on
surf(k1,k2,P), shading interp, axis square, colormap(hot),view(10,45)
xlabel('Angulo del retardador'), ylabel('Angulo de polarizador'), zlabel('Grado P')
hold off

%% Obtencion de grado P usando medicion perpendicular y paralela del polarizador
I_min=max(I);
I_max=min(I);
P=(I_max-I_min)/(I_max+I_min)
