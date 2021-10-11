clear
N=2^12; hg=1/(N+1);

t=(0:1:N)+0.5; v=kappa(t*hg)'; 
A=spdiags([-v(2:N+1),v(1:N)+v(2:N+1),-v(1:N)],-1:1,N,N)/hg^2; clear v
% A=sparse(full(A)^(0.1));
% A=A*A;
f=force_f((1:N)*hg)';
u=A\f; u=[0;u;0];


arrH=2.^(-2:-1:-10); nH=length(arrH);
ratio=1; n_ratio=length(ratio);
[arr_errH,arr_errL2,arr_err_H2,arr_err_H3,arr_err_W11,arr_err_L1,arr_err_W21,arr_err_W31]=deal(zeros(n_ratio,nH));

tic
for ii=1:n_ratio
    for jj=1:nH
        H=arrH(jj);
%         h=1/N*2;
        h=H*ratio(ii);
        sizeH=N*H; sizeh=N*h;
        pde1d_compression_err;
        factor1=1; factor2=1;
        arr_errH(ii,jj)=err_H/factor1;
        arr_errL2(ii,jj)=err_L2/factor1;
        err_H2=sqrt(sum(abs(gradient(gradient(u-uh,hg),hg)).^2))*sqrt(hg);
        err_H3=sqrt(sum(abs(gradient(gradient(gradient(u-uh,hg),hg),hg)).^2))*sqrt(hg);
        arr_err_H2(ii,jj)=err_H2/factor1;
        arr_err_H3(ii,jj)=err_H3/factor1;
        err_W11=sum(abs(gradient(u-uh,hg)))*hg;
        err_L1=sum(abs(u-uh))*hg;
        err_W21=sum(abs(gradient(gradient(u-uh,hg),hg)))*hg;
        err_W31=sum(abs(gradient(gradient(gradient(u-uh,hg),hg),hg)))*hg;
        arr_err_L1(ii,jj)=err_L1/factor2;
        arr_err_W11(ii,jj)=err_W11/factor2;
        arr_err_W21(ii,jj)=err_W21/factor2;
        arr_err_W31(ii,jj)=err_W31/factor2;
%         arr_err_comp(ii,jj)=err_comp;
        fprintf('H=%g,h=%g,completed, running time %g s \n',H,h, toc);
    end
end



%% figure
% figure
% plot(0:h:1,kappa(0:h:1))

% figure
% loglog(1./arrH',arr_errH','-o');
% % legend('H/h=1','H/h=2','H/h=4')
% xlabel('1/H')
% ylabel('H^1 error')
% 
% figure
% loglog(1./arrH',arr_errL2','-o');
% % legend('H/h=1','H/h=2','H/h=4')
% xlabel('1/H')
% ylabel('L^2 error')
% 
% figure
% loglog(1./arrH',arr_err_comp','-o');
% % legend('H/h=1','H/h=2','H/h=4')
% xlabel('1/H')
% ylabel('compression error')


function [y]=kappa(x)
    y=ones(1,length(x));
%     k=100;
%     W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
%     tmp_cos=cos((1:k)'*x);
%     tmp_sin=sin((1:k)'*x);
%     y=1+0.5*sin(W1'*tmp_cos+W2'*tmp_sin);%row vector
end
function [y]=force_f(x)
%     y=sin(x);
    k=1000;
    tmp_sin=sin(2*pi*(1:k)'*x);
    s=0; %f in H^s
    coef=randn(1,k)./(1:k).^(s+1/2);
    y=coef*tmp_sin;

%     k=100;
%     W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
%     tmp_cos=cos((1:k)'*x);
%     tmp_sin=sin((1:k)'*x);
%     y=1+0.5*sin(W1'*tmp_cos+W2'*tmp_sin);%row vector
end