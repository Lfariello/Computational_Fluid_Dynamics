function  F  = RHS_Zita_skwsim(~,ZITA) 
% Funzione per il calcolo dei termini convettivi secondo la formulazione
% Skew-Symmetric e dei termini diffusivi per l'equazione Psi-Zita.
% Il termine diffusivo e' ottenuto discretizzando 1/Re*nabla^2(zita)
% Il termine convettivo e' combinazione lineare delle discretizzazioni 
% della forma advective V*grad(zita) e della forma divergence div(V*zita).

global Nx Ny V U h hq Re

F       = zeros(Nx,Ny);
for i = 2:Nx-1 
    for j = 2:Ny-1
% TERMINE CONVETTIVO 
        adv = ( U(i,j)*(ZITA(i+1,j)-ZITA(i-1,j))/(2*h)+...
               V(i,j)*(ZITA(i,j+1)-ZITA(i,j-1))/(2*h) );
        
        div =  (U(i+1,j)*ZITA(i+1,j)-U(i-1,j)*ZITA(i-1,j))/(2*h)+...
               (V(i,j+1)*ZITA(i,j+1)-V(i,j-1)*ZITA(i,j-1))/(2*h) ;
        
        F(i,j) = - 0.5*(div + adv); 
        
% TERMINE DIFFUSIVO    
        F(i,j) = F(i,j) + (1/Re)*(ZITA(i+1,j)+ZITA(i-1,j)+ZITA(i,j+1)+...
                                  ZITA(i,j-1) - 4*ZITA(i,j))/hq;                       
    end
end

end