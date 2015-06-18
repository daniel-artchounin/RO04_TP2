N=2; // 10
A=toeplitz([2,-1,zeros(1,N-2)]);
// b = rand(N,1);
b = [0.1531217;0.6970851]; 

function y=J(x,A,b)
    // calcul d'une forme quadratique
    y = (1/2)* x'*A*x - b'*x;
endfunction

function dy=gradJ(x,A,b)
        // calcul d'un gradient
    dy = A*x-b;
endfunction

function dy=dk(x,A,b)
        // calcul d'un gradient
    dy = -gradJ(x,A,b);
endfunction

function dy=phi(x,A,b,rho)
    // calcul de phi
    dx = -gradJ(x,A,b);
    dy = J(x+rho*dx,A,b);
endfunction

function dy=phiPrime(x,A,b,rho)
    // calcul de phi
    dx = -gradJ(x,A,b);
    dy = gradJ(x+rho*dx,A,b)*dx;
endfunction

function rhok=Wolfe(J,dJ,Phi,PhiPrime,xk,Dk,m1,m2)
    // calcule le pas rhok d’apr\‘es la r\‘egle de Goldstein
    rho = 0,1;
    rit = 0;
    alpha = 0;
    _beta = %inf; 
    while ~( (Phi(xk,A,b,rho)<= Phi(xk,A,b,0) + m1*PhiPrime(x,A,b,0)*rho )  & ( PhiPrime(x,A,b,rho)>=m2*PhiPrime(x,A,b,0) ) )
        if  Phi(xk,A,b,rho) > (Phi(xk,A,b,0) + m1*PhiPrime(x,A,b,0)*rho) then
            _beta = rho;
            rho = (alpha + _beta)/2;
        elseif  (Phi(xk,A,b,rho)<= Phi(xk,A,b,0) + m1*PhiPrime(x,A,b,0)*rho ) & ( PhiPrime(x,A,b,rho)<m2*PhiPrime(x,A,b,0) ) then
            alpha = rho;
            if (_beta  <  %inf) then
                rho = (alpha + _beta)/2
            else
                rho = 2*rho;
            end           

        end
        rit = rit + 1;
    end
    rhok = rho;
endfunction

function [c,xn,Jxn,errors]=gradient(x0,J,dJ,stop,nmax,x_opt,m1,m2)
    // a la fin de l’algorithme c est un vecteur colonne contenant toutes les valeurs J(xk)
    // a la fin de l’algorithme errors est un vecteur colonne contenant toutes les erreurs
    // ||xk-x_opt|| pour la norme euclidienne
    // xn et Jxn sont les valeurs finales obtenues
    // stop : valeur num\’erique du crit\‘ere d’arret
    // nmax : nombre d’it\’erations maximales
    // rho : le pas est calcul\’e \‘a chaque it\’eration \‘a l’aide de la fonction pas_var_opt
    i = 1;    
    xk = x0;
    c=[J(x0,A,b)];
    xkplus1 = xk + Goldstein(J,dJ,phi,xk,dk,m1,m2)*dk(xk,A,b);
    errors=[norm((xkplus1-x_opt),2)];
    c = [c,J(xkplus1,A,b)];
    while i<nmax & norm((xkplus1-xk),2)> stop    
        xk = xkplus1;
        xkplus1 = xk + Goldstein(J,dJ,phi,xk,dk,m1,m2)*dk(xk,A,b);  
        errors=[errors, norm((xkplus1-x_opt),2)];        
        c = [c, J(xkplus1,A,b)];
        i = i + 1;
    end
    Jxn = J(xkplus1,A,b);
    xn = xkplus1;    
endfunction

x0=zeros(N,1);
nmax = 20000;  // 500
// stop = 1e-3; // 0.1 1e-2; 1e-4; 1e-6; 1e-9;disp('nb iterations :');
// stop = 1e-2; // 
// stop = 0.1 ; 
//stop = 1e-4;
// stop = 1e-6;
// stop = 1e-9;
// stop = 1e-3;
disp('A=');
disp(A);
disp('x_opt_1');
x_opt_1 = inv(A)*b;
disp(x_opt_1);
// m1 = 1e-4;
// m2 = 0.9;
m1 = 0.1;
m2 = 0.7;

[c,xn,Jxn,errors]=gradient(x0,J,gradJ,stop,nmax,x_opt_1,m1,m2);
// disp('c');
// disp(c);
Nit = length(c)-1;
disp(size([0:Nit]));
disp('nb iterations ');
disp(Nit);
disp(size(c));
clf;
figure(2);
disp('nb iterations :');
disp(Nit);
plot([0:Nit],c,'k.');
legend('Wolfe'); // affichage de la légende
xlabel('$i$'); // titre de l'abscisse
ylabel('$f(x_i)$'); // titre de l'ordonnée
title('Méthode du gradient avec calcul du pas via la méthode de wolfe'); // titre du graphique

