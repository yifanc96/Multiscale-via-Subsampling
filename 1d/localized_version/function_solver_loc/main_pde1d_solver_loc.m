clear
N=2^12; hg=1/(N+1);
t=(0:1:N)+0.5; v=kappa(t*hg)'; 
A=spdiags([-v(2:N+1),v(1:N)+v(2:N+1),-v(1:N)],-1:1,N,N)/hg^2; clear v
f=force_f((1:N)*hg)';
u=A\f; u=[0;u;0];

H=2^(-7); h=H;
sizeH=N*H; sizeh=N*h;
tic;
uh=pde1d_subsample_solver_loc(A,f,N,sizeH,sizeh,1);
toc
err_H=norm(gradient(u-uh,hg),2)*sqrt(hg);
err_L2=norm(u-uh,2)*sqrt(hg);

fprintf('H^1 error is %g, L^2 error is %g\n',err_H,err_L2);

function [y]=kappa(x)
    k=100;
    W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
    tmp_cos=cos((1:k)'*x);
    tmp_sin=sin((1:k)'*x);
    y=1+0.5*sin(W1'*tmp_cos+W2'*tmp_sin);%row vector
end
function [y]=force_f(x)
    k=100;
    W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
    tmp_cos=cos((1:k)'*x);
    tmp_sin=sin((1:k)'*x);
    y=1+0.5*sin(W1'*tmp_cos+W2'*tmp_sin);%row vector
end