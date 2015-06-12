N=10;
A=toeplitz([2,-1,zeros(1,N-2)]);
b = rand(N,1);

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

function [c,xn,Jxn,errors,nbIterations]=gradient_pasoptimal(x0,J,dJ,stop,nmax,x_opt)
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
    dk = -dJ(xk,A,b)
    xkplus1 = xk + pas_var_opt(dk, A)*dk;
    errors=[norm((xkplus1-x_opt),2)];
    c = [c,J(xkplus1,A,b)];
    while i<nmax & norm((xkplus1-xk),2)> stop    
        xk = xkplus1;
        dk = -dJ(xk,A,b)
        xkplus1 = xk + pas_var_opt(dk, A)*dk;        
        errors=[errors, norm((xkplus1-x_opt),2)];        
        c = [c, J(xkplus1,A,b)];
        i = i + 1;
    end
    Jxn = J(xkplus1,A,b);
    xn = xkplus1;   
    nbIterations = i; 
endfunction

function rho_opt=pas_cst_opt(A)
    // pas optimal pour la matrice A
    vp_max = max(spec(A));
    vp_min = min(spec(A));
    rho_opt= 2/(vp_max+vp_min);
endfunction

function [c,xn,Jxn,errors,nbIterations]=gradient_pasconstant_new(x0,J,dJ,rho,stop,nmax,x_opt)
    // a la fin de l’algorithme c est un vecteur colonne contenant toutes les valeurs J(xk)
    // a la fin de l’algorithme errors est un vecteur colonne contenant toutes les erreurs
    // ||xk-x_opt|| pour la norme euclidienne
    // xn et Jxn sont les valeurs finales obtenues
    // stop : valeur num\’erique du crit\‘ere d’arret
    // nmax : nombre d’it\’erations maximales
    // rho : choisir une valeur strictement positive mais inf\’erieur \‘a rho_max
    i = 1;    
    xk = x0;
    c=[J(x0,A,b)];
    xkplus1 = xk - rho*dJ(xk,A,b);
    errors=[norm((xkplus1-x_opt),2)];
    c = [c,J(xkplus1,A,b)];
    while i<nmax & norm((xkplus1-xk),2)> stop    
        xk = xkplus1;
        xkplus1 = xk -rho*dJ(xk,A,b);        
        errors=[errors, norm((xkplus1-x_opt),2)];        
        c = [c, J(xkplus1,A,b)];
        i = i + 1;
    end
    Jxn = J(xkplus1,A,b);
    xn = xkplus1;   
    nbIterations =  i;
endfunction

x0=zeros(N,1);
nmax = 500;
stop = 1e-3;


disp('A=');
disp(A);
disp('x_opt_1');
x_opt_1 = inv(A)*b;
disp(x_opt_1);

[c,xn,Jxn,errors,iterations]=gradient_pasoptimal(x0,J,gradJ,stop,nmax,x_opt_1);
// disp('c');
// disp(c);
Nit = length(errors)-1;
disp(size([0:Nit]));
disp(size(errors));
figure(0);
plot([0:Nit],errors,'k.');


[c_bis,xn_bis,Jxn_bis,errors_bis,iterations_bis]=gradient_pasconstant_new(x0,J,gradJ,pas_cst_opt(A),stop,nmax,x_opt_1);
Nit = length(errors_bis)-1;
disp(size([0:Nit]));
disp(size(errors_bis));
plot([0:Nit],errors_bis,'b.');

disp('nb iterations pas optimal');
disp(iterations);
disp('nb iterations pas constant');
disp(iterations_bis);


     
legend('Pas optimal','Pas constant'); // legend
xlabel('$i$'); // titre de l'abscisse
ylabel('$\epsilon(x_i)$'); // titre de l'ordonnée
title('Comparaison entre les erreurs de la méthode du gradient à pas constant et les erreurs du gradient à pas optimal'); // titre du graphique

