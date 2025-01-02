% Driven Cavity flow in slender cavities with unsteady lid velocity
%
% Codice per la simulazione delle equazioni di Navier Stokes 2D 
% incompressibili nella formulazione vorticità-funzione di corrente
% (modello 'psi-zita').
% Il problema affrontato in questo codice è quello del flusso dentro una 
% cavità RETTANGOLARE con un lato che scorre parallelamente a se stesso con 
% velocità NON uniforme u_w (problema della 'lid-driven cavity').
%
% Il modello utilizzato è quello 'psi-zita' 
%
%  @zita/@t + 0.5[V*grad(zita)+div(V*zita)] = 1/Re*nabla^2(zita)
%  nabla^2 psi = zita
%
% Le condizioni al contorno sono date alla Dirichlet per la ellittica per
% psi, e dalle condizioni alla 'Thom' per la zita, che stima la produzione
% di vorticita' alla parete tramite:
%           zita_wall = 2*(psi_int - psi_wall + uw*h)/h^2
%
% La collocazione delle variabili sul mesh è di tipo 'collocated', ovvero
% tutte le variabili sono collocate nei nodi di griglia, che cadono anche
% direttamente sul bordo.

clc; clear all; close all;
global Nx Ny V U h hq Re 
% DOMINIO
Ly  = 1;        % Estenzione dominio lungo y 
AR  = 3;        % AR = Aspect Ratio (per lo studio parametrico)
Lx = AR*Ly;     % Estenzione dominio lungo x 
Ny  = 80;   Nx  = AR*Ny;         % Numero di nodi lungo i lati
x  = linspace(0,Lx,Nx);          % Mesh nella direzione x
y  = linspace(0,Ly,Ny);          % Mesh nella direzione y
hx = x(2)-x(1);                  % Passo spaziale lungo x
hy = y(2)-y(1);                  % Passo spaziale lungo y
h  = hx;       hq  = h*h;        % Impongo un mesh uniforme

% Matrice 'G' topologica per la definizione del LAPLACIANO
G = zeros(Nx,Ny); % Inizializzazione     
k = 0;  % Contatore --> indice progressivo di numerazione dei nodi
for j = 2:Ny-1      % Riempimento matrice G punti interni 
    for i = 2:Nx-1
            k = k+1;     % Indice progressivo per i punti interni (fluido) 
            G(i,j) = k;  % Numerazione punti interni
    end
end
Lap  = -delsq(G)/hq;     % Generazione del Laplaciano 

% Inizializzazione matrici e Reynolds
PSI  = zeros(Nx,Ny);     ZITA = zeros(Nx,Ny);
U = zeros(Nx,Ny);        V = zeros(Nx,Ny);
Re = 1200;

% TEMPO
T = 40;                        % Tempo finale di simulazione
beta = 0.1;  Dtd = hq*beta*Re; % Dt basato su beta   
Dt = Dtd;     
Nt = round(T/Dt);              % Numero di timestep temporali 

% Inizializzazione array estrofia e velocita della parete
enstrophy = NaN(Nt,1);
u_w = NaN(Nt,1);

% Ciclo instazionario
it = 0;
for it = 1:Nt
    t = it*Dt;
    velocita = 'lin'; % Parametro per scegliere il tipo di velocita
    switch velocita
        case 'sin'    % Velocita parete sinusuidale
         A = 1;       % Ampiezza sinusoide
         om = 0.6;    % Pulsazione
         phi0 = 0;    % Sfasamento iniziale
         u_w(it) = A*sin(om*t + phi0);
     
         case 'lin'   % Velocita parete instazionaria (Variazione lineare)
         u_min = -1;
         u_max =  1;
         u_w(it) = ((u_max - u_min)/T)*t + u_min;
        
        case 'cost'  % Velocita costante (serve solo per confronto ed estrofia)
        u_w(it) = 1;
    end    

    % Velocita parete NORD 
    U(:,Ny) = u_w(it); 
    
    % Calcolo della velocità dalla psi (punti interni)
    i = 2:Nx-1;   j = 2:Ny-1;
    U(i,j) =  (PSI(i,j+1)-PSI(i,j-1))/(2*h);
    V(i,j) = -(PSI(i+1,j)-PSI(i-1,j))/(2*h);

    % Formula di Thom - CC per l'equazione della vorticita 
    %CC Parete SUD
    i = 2:Nx-1;        j = 1;               
    ZITA(i,j) = 2*(PSI(i,j+1) - PSI(i,j))/hq; 
    %CC parete EST
    j = 1:Ny-1;        i = 1;                
    ZITA(i,j) = 2*(PSI(i+1,j) - PSI(i,j))/hq;
    %CC parete WEST
    j = 1:Ny-1;        i = Nx;              
    ZITA(i,j) = 2*(PSI(i-1,j) - PSI(i,j))/hq;
    %CC parete NORD - unica con u_w non nulla
    i = 2:Nx-1;        j = Ny;       
    if u_w < 0
       ZITA(i,j) = 2*(PSI(i,j-1) - PSI(i,j) - u_w(it)*h)/hq;
    else
       ZITA(i,j) = 2*(PSI(i,j-1) - PSI(i,j) + u_w(it)*h)/hq;
    end

    % Integrazione temporale mediante il metodo di Runge-Kutta 
    ZETA1 = ZITA;              F1 = RHS_Zita_skwsim(t,ZETA1);
    ZETA2 = ZITA + 0.5*Dt*F1;  F2 = RHS_Zita_skwsim(t,ZETA2);
    ZETA3 = ZITA + 0.5*Dt*F2;  F3 = RHS_Zita_skwsim(t,ZETA3);
    ZETA4 = ZITA + Dt*F3;      F4 = RHS_Zita_skwsim(t,ZETA4);

    ZITA  = ZITA + Dt*((1/6)*(F1+F4)+(1/3)*(F2+F3));
   
    zita = ZITA(G>0); % Zita nei punti interni messo in un vettore colonna 
  
    % Calcolo Psi con metodo diretto - Formulazione con Backslash
    Psi = Lap\zita;
    
    % Riversiamo la Psi nell'array 2D 
    PSI(G>0)=Psi;

    % Calcolo dell'enstrofia globale 
    enstrophy(it) = sum(ZITA(:).^2);
    
    % Grafica
   if mod(it,5) == 0
        
        % Andamento della velocita della parete
        figure(1);
        plot(Dt:Dt:Nt*Dt,u_w) 
        title('Velocity'); 
        axis([0 T -3 3]); grid on;
        xlabel('t'); ylabel('u_w(t)'); hold on
        
        % Evoluzione del campo di moto
        figure(2);clf;
        pcolor(x,y,ZITA'); shading interp; colormap jet; hold on
        vneg = linspace(min(min(PSI)),0,10);
        vpos = linspace(0,max(max(PSI)),10);
        contour(x,y,PSI',vneg,'k'); axis square; hold on;
        contour(x,y,PSI',vpos,'r'); 
        title(['Psi-zita. Lid Driven Cavity. t = ',num2str(t)]);
        axis image; drawnow; hold on
        
        % Evoluzione temporale dell'enstrofia adimensionalizzata
        figure(3)
        plot(Dt:Dt:Nt*Dt,enstrophy./enstrophy(1) -1)
        title('Enstrophy')
        xlabel('t','Interpreter','latex');
        ylabel('$\Omega$','Interpreter','latex');
        hold off 
    end
end


        
