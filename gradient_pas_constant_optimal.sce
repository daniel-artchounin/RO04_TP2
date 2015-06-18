N= 2;// 10;
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

function rho_max=rmax(A)
    // valeur propre maximale de A
    alpha = max(spec(A));
    // contrainte sur le pas de gradient
    rho_max= 2/alpha;
endfunction

function rho_opt=pas_cst_opt(A)
    // pas optimal pour la matrice A
    vp_max = max(spec(A));
    vp_min = min(spec(A));
    rho_opt= 2/(vp_max+vp_min);
endfunction

function [c,xn,Jxn,errors]=gradient_pasconstant_new(x0,J,dJ,rho,stop,nmax,x_opt)
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
endfunction


x0=zeros(N,1);
nmax = 20000; // 500;
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


figure(1);
rho_tab = linspace(0,rmax(A)-0.1,5);
disp(rho_tab);

[c_bis,xn_bis,Jxn_bis,errors_bis]=gradient_pasconstant_new(x0,J,gradJ,rho_tab(1),stop,nmax,x_opt_1);
disp('c_bis');
disp(c_bis);
Nit = length(c_bis)-1;
disp(size([0:Nit]));
disp(size(c_bis));
plot([0:Nit],c_bis,'y.');

[c_bis,xn_bis,Jxn_bis,errors_bis]=gradient_pasconstant_new(x0,J,gradJ,rho_tab(2),stop,nmax,x_opt_1);
Nit = length(c_bis)-1;
disp(size([0:Nit]));
disp(size(c_bis));
plot([0:Nit],c_bis,'g.');

[c_bis,xn_bis,Jxn_bis,errors_bis]=gradient_pasconstant_new(x0,J,gradJ,rho_tab(3),stop,nmax,x_opt_1);
Nit = length(c_bis)-1;
disp(size([0:Nit]));
disp(size(c_bis));
plot([0:Nit],c_bis,'m.');

[c_bis,xn_bis,Jxn_bis,errors_bis]=gradient_pasconstant_new(x0,J,gradJ,rho_tab(4),stop,nmax,x_opt_1);
Nit = length(c_bis)-1;
disp(size([0:Nit]));
disp(size(c_bis));
plot([0:Nit],c_bis,'b.');

[c_bis,xn_bis,Jxn_bis,errors_bis]=gradient_pasconstant_new(x0,J,gradJ,rho_tab(5),stop,nmax,x_opt_1);
Nit = length(c_bis)-1;
disp(size([0:Nit]));
disp(size(c_bis));
figure(1);
plot([0:Nit],c_bis,'r.');

[c_bis,xn_bis,Jxn_bis,errors_bis]=gradient_pasconstant_new(x0,J,gradJ,pas_cst_opt(A),stop,nmax,x_opt_1);
Nit = length(c_bis)-1;
disp('nb iterations (pas constant optimal)');
disp(Nit);
disp(size([0:Nit]));
disp(size(c_bis));
plot([0:Nit],c_bis,'k.');

     
legend('$\rho = 0$','$\rho = 0.1275840$','$\rho = 0.2551680$','$\rho = 0.3827521$','$\rho = 0.5103361$', 'pas constant optimal'); // affichage de la légende
xlabel('$i$'); // titre de l'abscisse
ylabel('$f(x_i)$'); // titre de l'ordonnée
title('Méthode du gradient à pas constant'); // titre du graphique


