//
//------------------------------------------------------
// BFGS for Rastrigin
disp('BFGS for Rastrigin');
//
//-----------------------------------------------------
//
function y=rastrigin(x)  // the function to optimize
n=max(size(x));
y=n+sum(x.^2-cos(2*%pi*x));
endfunction
//-----------------------------------------------------
function y=rastrigingrad(x)  // the gradient of the function to optimize
y=2*x+2*%pi*sin(2*%pi*x);
endfunction
//--------------------------------------------------------
function z=dotprod(x,y); // computes the dot product of x and y
z=sum(x.*y);
endfunction
//-----------------------------------------------------
function [d,gx]=descentdirection(Bk,gx);   // descent direction: gradient
d=linsolve(Bk,gx);
endfunction
//----------------------------------------------------
function [xnew,fnew]=backtracking(f,x,fx,gx,d);// line search by backtracking until Armijo condition
tau=0.3;
bet=0.0001;
alphainit=1;
alpha=alphainit;
xnew=x+alpha*d;
fnew=f(xnew);
while(fnew>fx+bet*alpha*dotprod(gx,d))
  alpha=tau*alpha;
  xnew=x+alpha*d;
  fnew=f(xnew);
end
endfunction


//-------------------------------------------------
// main program
//
timer();
n=2;   
epsilon=1E-5;
//
xmin=-5*ones(1,n);
xmax=5*ones(1,n);
u=rand(1,n);
x0=xmin+(xmax-xmin).*u;
x=x0;fx=rastrigin(x);gx=rastrigingrad(x);
itgrad=1;
Xbest=x;Fbest=fx;
Bk=eye(n,n);
//
while (norm(gx)>epsilon)
  x0=x;gx0=gx;
  d=descentdirection(Bk,gx');  
  [x,fx]=backtracking(rastrigin,x,fx,gx,d');
  Xbest=[Xbest;x];
  Fbest=[Fbest;fx];
  gx=rastrigingrad(x);
  itgrad=itgrad+1;
   yk=gx'-gx0';sk=x'-x0';
  Bk=Bk+(yk*yk')/dotprod(sk,yk)-(Bk*sk*sk'*Bk)/dotprod(sk,Bk*sk);
 // disp(mini(spec(Bk)));
end
//------------------------------------------------------
// results display
//
disp('function evaluation number:');disp(itfct);
disp('gradient evaluation number:');disp(itgrad);
//
disp('minimum obtained:');disp(x);
disp('corresponding value by f:');disp(fx);
//disp('computational time:');disp(timer());
//
// case of the rastrigin function with 2 parameters (n=2)------
if (n==2)
xmin=-5.12;xmax=5.12;N=300;
xplot=xmin:((xmax-xmin)/(N-1)):xmax;
yplot=xplot;
zplot=zeros(N,N);
for i=1:N
for j=1:N
zplot(i,j)=rastrigin([xplot(i),yplot(j)]);
end
end
xset('window',0)
clf()
plot2d(Xbest(:,1),Xbest(:,2),rect=[-5.12,-5.12,5.12,5.12]); 
contour2d(xplot,yplot,zplot,[0:0.01:0.1,0.2:1,1:10]);
xset('window',1)
clf()
plot2d(Xbest(:,1),Xbest(:,2),rect=[x(1)-0.1,x(2)-0.1,x(1)+0.1,x(2)+0.1]); 
contour2d(xplot,yplot,zplot,[fx:0.1:(fx+1)]);
end








  










  
