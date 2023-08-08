clc
clear all;

%Inicialización
Nx = 100; %Tamaño de la cuerda
dx = 1;
x(:,1) = (0:Nx-1)*dx;

T = 1000; %Tiempo
dt = 0.0005;
t(:,1)= (0:T-1)*dt;

c = 800; %Velocidad de las ondas 
C = c*(dt/dx); 

y = zeros(T,Nx);  %y(tiempo, espacio)



%Condiciones iniciales
f = 50; 
liminf = 50;
limsup = 90;
rang = limsup - liminf;

y(1,(liminf:limsup)) = sin(2*pi*f.*t(1:rang+1));
y(2,(liminf:limsup)) = sin(2*pi*f.*t(1:rang+1));

y(:,1) = 0;
y(:,Nx) = 0;


% y(1,(1:Nx/2)) = linspace(0,1, Nx/2);
% y(1,(Nx/2+1:Nx)) = linspace(1,0, Nx/2);
% y(2,(1:Nx/2)) = linspace(0,1, Nx/2);
% y(2,(Nx/2+1:Nx)) = linspace(1,0, Nx/2);

for n = 3:T
    for k = 2:Nx-1
        y(n,k) = 2*y(n-1,k) - y(n-2,k) + (C^2)*(y(n-1,k-1)-2*y(n-1,k)+y(n-1,k+1));    
    end                   
end

for n=1:T
    plot(x,y(n,:));
    title('Propagación de una onda en una cuerda')
    xlabel('x')
    ylabel('y')
    ylim([-1 1])
    pause(0.000001)
end

Sonido = zeros(Nx*T,1);

for n=1:T
    for i=1:Nx
        Sonido((n-1)*(Nx)+i) = 0.5*y(n,i);
    end
    
end

Fs = 44100;
audiowrite('Audio2.wav',Sonido,Fs);

