function y=rastrigin(x)  // the function to optimize
    n=max(size(x));
    y=n+sum(x.^2-cos(2*%pi*x));
endfunction

//-----------------------------------------------------
function y=rastrigingrad(x)  // the gradient of the function to optimize
    y=2*x+2*%pi*sin(2*%pi*x);
endfunction

function dy=dk(x)
    // direction de descente
    dy = -rastrigingrad(x);
endfunction

function dy=phi(x,rho)
    // calcul de phi
    dx = dk(x);
    dy = rastrigin(x+rho*dx);
endfunction

function dy=phiPrime(x,rho)
    // calcul de phi
    dx = dk(x);
    dy = rastrigingrad(x+rho*dx)'*dx;
endfunction

function rhok=Wolfe(rastrigin,rastrigingrad,Phi,PhiPrime,xk,Dk,m1,m2)
    // calcule le pas rhok d’apr\‘es la r\‘egle de Goldstein
    rho = rand()/1000+0,000001; // calcul de rho_0 -> de départ joue un rôle prépondérant -> un rho_0 trop grand entraine un arrêt de l'algorithme
    disp(rho);
    rit = 0;
    alpha = 0;
    _beta = %inf; 
    while ~( (Phi(xk,rho)<= Phi(xk,0) + m1*PhiPrime(xk,0)*rho )  & ( PhiPrime(xk,rho)>=m2*PhiPrime(xk,0) ) )
        if  Phi(xk,rho) > (Phi(xk,0) + m1*PhiPrime(xk,0)*rho) then
            _beta = rho;
            rho = (alpha + _beta)/2;
        elseif  (Phi(xk,rho)<= Phi(xk,0) + m1*PhiPrime(xk,0)*rho ) & ( PhiPrime(xk,rho) < m2*PhiPrime(xk,0) ) then
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

function [c,xn,rastriginxn,errors]=gradient(x0,rastrigin,rastrigingrad,stop,nmax,x_opt,m1,m2)
    // a la fin de l’algorithme c est un vecteur colonne contenant toutes les valeurs rastrigin(xk)
    // a la fin de l’algorithme errors est un vecteur colonne contenant toutes les erreurs
    // ||xk-x_opt|| pour la norme euclidienne
    // xn et rastriginxn sont les valeurs finales obtenues
    // stop : valeur num\’erique du crit\‘ere d’arret
    // nmax : nombre d’it\’erations maximales
    // rho : le pas est calcul\’e \‘a chaque it\’eration \‘a l’aide de la fonction pas_var_opt
    i = 1;    
    xk = x0;
    c=[rastrigin(x0)];
    xkplus1 = xk + Wolfe(rastrigin,rastrigingrad,phi,phiPrime,xk,dk,m1,m2)*dk(xk);
    errors=[norm((xkplus1-x_opt),2)];
    c = [c,rastrigin(xkplus1)];
    while i<nmax & norm((xkplus1-xk),2)> stop    
        xk = xkplus1;
        xkplus1 = xk + Wolfe(rastrigin,rastrigingrad,phi,phiPrime,xk,dk,m1,m2)*dk(xk);  
        errors=[errors, norm((xkplus1-x_opt),2)];        
        c = [c, rastrigin(xkplus1)];
        i = i + 1;
    end
    rastriginxn = rastrigin(xkplus1);
    xn = xkplus1;    
endfunction

// programme principal
N=10;
x0= ones(N,1); // le x0 joue un rôle important ici (si on prend x0=0 cela ne converge pas)
nmax = 500;
stop = 1e-7;

m1 = 1e-4;
m2 = 0.9;

[c,xn,rastriginxn,errors]=gradient(x0,rastrigin,rastrigingrad,stop,nmax,x0,m1,m2); // /!\ on n'a pas l'optimal -> calcul d'erreur pas possible
Nit = length(c)-1; // -1 pour que la taille du vecteur des abscisses coincide avec celui des ordonnées 
clf; // on efface le contenu de la figure
figure(0);
plot([0:Nit],c,'k.'); // affichage du graphique dans la figure 0
legend('Wolfe'); // affichage de la légende
xlabel('$i$'); // titre de l'abscisse
ylabel('$f(x_i)$'); // titre de l'ordonnée
title('Minimum de la fonction de rastrigin avec la méthode de Wolfe'); // titre du graphique

