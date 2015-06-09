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


// Programmer un algorithme de gradient avec un pas constant.
// Les donn ́ees sont la fonction coˆ
// ut J, son gradient gradJ, le premier terme x0 de la suite (x k ) k et les
// crit`eres d’arrˆet : un nombre d’it ́eration maximal nmax ou une borne maximale pour l’erreur entre deux
// it ́er ́es stop. Le pas rho devra ˆetre donn ́e dans le script de la fonction.
// Remarque: Il faut  ́evaluer num ́eriquement la borne maximale du pas ρ (cf exercice 1 du TD3).


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

function [c,xn,Jxn]=gradient_pasconstant(x0,J,dJ,rho,stop,nmax)
    // a la fin de l’algorithme c est un vecteur colonne contenant toutes les valeurs J(xk)
    // xn et Jxn sont les valeurs finales obtenues
    // stop : valeur num\’erique du crit\‘ere d’arret
    // nmax : nombre d’it\’erations maximales
    // rho : choisir une valeur strictement positive mais inf\’erieur \‘a rho_max
    i = 2;    
    xk = x0;
    c=[J(x0,A,b)];
    xkplus1 = xk - rho*dJ(xk,A,b);
    c = [c,J(xkplus1,A,b)];
    while i<nmax & norm((xkplus1 -xk),2)> stop      
        xk = xkplus1;  
        xkplus1 = xk - rho*dJ(xk,A,b);
        c = [c,J(xkplus1,A,b)];
        i = i + 1;
    end
    Jxn = J(xkplus1,A,b);
    xn = xkplus1;     
endfunction

function [c,xn,Jxn,errors]=gradient_pasconstant_new(x0,J,dJ,rho,stop,nmax,x_opt)
    // a la fin de l’algorithme c est un vecteur colonne contenant toutes les valeurs J(xk)
    // a la fin de l’algorithme errors est un vecteur colonne contenant toutes les erreurs
    // ||xk-x_opt|| pour la norme euclidienne
    // xn et Jxn sont les valeurs finales obtenues
    // stop : valeur num\’erique du crit\‘ere d’arret
    // nmax : nombre d’it\’erations maximales
    // rho : choisir une valeur strictement positive mais inf\’erieur \‘a rho_max
    i = 2;    
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
    i = 2;    
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
endfunction

x0=zeros(N,1);
nmax = 500;
stop = 1e-3;
rho = rmax(A);
disp('A=');
disp(A);
disp('x_opt_1');
x_opt_1 = inv(A)*b;
disp(x_opt_1);
x_opt_2 = A \ b;
disp('x_opt_2');
disp(x_opt_2);
// x_opt_3 = linsolve(A , b, x0);
// disp('x_opt_3');
// disp(x_opt_3);
[c,xn,Jxn]=gradient_pasconstant(x0,J,gradJ,rho,stop,nmax);
disp('xn');
disp(xn);
disp('erreur : ');
disp(norm(xn-x_opt_1,2));
Nit = length(c)-1;
disp(size([0:Nit]));
disp(size(c));
clf; // on efface le contenu des figures
figure(0);
plot([0:Nit],c,'k.');
legend('pas maximal'); // affichage de la légende
xlabel('$i$'); // titre de l'abscisse
ylabel('$f(x_i)$'); // titre de l'ordonnée
title('Méthode du gradient à pas constant'); // titre du graphique

figure(1);
rho_tab = linspace(0,rmax(A),5);
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
disp(size([0:Nit]));
disp(size(c_bis));
plot([0:Nit],c_bis,'k.');

     
legend('$\rho = 0$','$\rho = 0.1275840$','$\rho = 0.2551680$','$\rho = 0.3827521$','$\rho = 0.5103361$', 'pas constant optimal'); // affichage de la légende
xlabel('$i$'); // titre de l'abscisse
ylabel('$f(x_i)$'); // titre de l'ordonnée
title('Méthode du gradient à pas constant'); // titre du graphique



[c,xn,Jxn,errors]=gradient_pasoptimal(x0,J,gradJ,stop,nmax,x_opt_1);
disp('c');
disp(c);
Nit = length(c)-1;
disp(size([0:Nit]));
disp(size(c));
figure(2);
plot([0:Nit],c,'k.');
legend('Pas optimal'); // affichage de la légende
xlabel('$i$'); // titre de l'abscisse
ylabel('$f(x_i)$'); // titre de l'ordonnée
title('Méthode du gradient à pas optimal'); // titre du graphique
