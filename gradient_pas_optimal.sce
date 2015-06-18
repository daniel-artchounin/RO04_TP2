N=2; // 10
A=toeplitz([2,-1,zeros(1,N-2)]);
// b = rand(N,1);
// b = ones(N,1).*0.5;
b = [0.1531217;0.6970851]; 

function y=J(x,A,b)
    // calcul d'une forme quadratique
    y = (1/2)* x'*A*x - b'*x;
endfunction

function dy=gradJ(x,A,b)
        // calcul d'un gradient
    dy = A*x-b;
endfunction

function rho_opt=pas_var_opt(dk, A)
    // pas variable optimal pour la matrice Aendfunction
    rho_opt = (dk'*dk)/(dk'*A*dk);
endfunction

function [c,xn,Jxn,errors]=gradient_pasoptimal(x0,J,dJ,stop,nmax,x_opt)
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
    dk = -dJ(xk,A,b);
    xkplus1 = xk + pas_var_opt(dk, A)*dk;
    errors=[norm((xkplus1-x_opt),2)];
    c = [c,J(xkplus1,A,b)];
    while i<nmax & norm((xkplus1-xk),2)> stop    
        xk = xkplus1;
        dk = -dJ(xk,A,b)
        disp('****');
        disp(dk);
        xkplus1 = xk + pas_var_opt(dk, A)*dk;        
        errors=[errors, norm((xkplus1-x_opt),2)];        
        c = [c, J(xkplus1,A,b)];
        i = i + 1;
    end
    Jxn = J(xkplus1,A,b);
    xn = xkplus1;    
endfunction


x0=zeros(N,1);
nmax =  20000; // 500;

// stop = 1e-2;
// stop = 0.1 ; 
// stop = 1e-4;
// stop = 1e-6;
stop = 1e-9;
// stop = 1e-3;
disp('A=');
disp(A);
disp('x_opt_1');
x_opt_1 = inv(A)*b;
disp(x_opt_1);

[c,xn,Jxn,errors]=gradient_pasoptimal(x0,J,gradJ,stop,nmax,x_opt_1);
disp('c');
disp(c);
Nit = length(c)-1;
disp(size([0:Nit]));
disp('nb iterations (pas optimal)');
disp(Nit);
disp(size(c));
figure(2);
plot([0:Nit],c,'k.');
legend('Pas optimal'); // affichage de la légende
xlabel('$i$'); // titre de l'abscisse
ylabel('$f(x_i)$'); // titre de l'ordonnée
title('Méthode du gradient à pas optimal'); // titre du graphique
