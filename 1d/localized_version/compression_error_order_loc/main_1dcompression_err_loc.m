clear
N=2^11; hg=1/(N+1);

t=(0:1:N)+0.5; v=kappa(t*hg)'; 
A=spdiags([-v(2:N+1),v(1:N)+v(2:N+1),-v(1:N)],-1:1,N,N)/hg^2; clear v
f=force_f((1:N)*hg)';
u=A\f; u=[0;u;0];
B=inv(A);


arrH=2.^(-1:-1:-8); nH=length(arrH);
ratio=[1,1/2,1/4]; n_ratio=length(ratio);
[arr_errH,arr_errL2,arr_err_comp]=deal(zeros(n_ratio,nH));
over_samp=1;

tic
for ii=1:n_ratio
    for jj=1:nH
        H=arrH(jj);h=H*ratio(ii);
        t=floor(over_samp*log(1/H)/log(2));
        sizeH=N*H; sizeh=N*h;
        pde1d_compression_err_loc;
        arr_errH(ii,jj)=err_H;
        arr_errL2(ii,jj)=err_L2;
        arr_err_comp(ii,jj)=err_comp;
        fprintf('H=%g,h=%g,completed, running time %g s \n',H,h, toc);
    end
end



%% figure
figure
plot(0:h:1,kappa(0:h:1))

figure
loglog(arrH',arr_errH','-o');
% legend('H/h=1','H/h=2','H/h=4')
xlabel('H')
ylabel('H^1 error')

figure
loglog(arrH',arr_errL2','-o');
% legend('H/h=1','H/h=2','H/h=4')
xlabel('H')
ylabel('L^2 error')

figure
loglog(arrH',arr_err_comp','-o');
% legend('H/h=1','H/h=2','H/h=4')
xlabel('H')
ylabel('compression error')


function [y]=kappa(x)
    k=100;
    W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
    tmp_cos=cos((1:k)'*x);
    tmp_sin=sin((1:k)'*x);
    y=1+0.5*sin(W1'*tmp_cos+W2'*tmp_sin);%row vector
end
function [y]=force_f(x)
    k=50;
    W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
    tmp_cos=cos((1:k)'*x);
    tmp_sin=sin((1:k)'*x);
    y=sin(W1'*tmp_cos+W2'*tmp_sin);%row vector
end