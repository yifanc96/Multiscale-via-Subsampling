clear

d=2; N=2^6; hg=1/(N+1); %size N, total N^2 points
% stiffness matrix
t_mid=(0:1:N)'+0.5; t_grid=1:N; [tt,ss]=meshgrid(t_mid*hg,t_grid*hg);
a1=kappa([tt(:),ss(:)]); a2=kappa([ss(:),tt(:)]); 
a1=reshape(a1,N,N+1);a2=reshape(a2,N+1,N);
a_diag=reshape(a1(:,1:N)+a1(:,2:N+1)+a2(1:N,:)+a2(2:N+1,:),[],1);
temp1=reshape([a2(2:N,:);zeros(1,N)],[],1); temp1=temp1(1:N^2-1);
temp2=reshape(a1(:,2:N),[],1);
a_sub1=[temp1;0]; a_super1=[0;temp1];
a_sub2=[temp2;zeros(N,1)]; a_super2=[zeros(N,1);temp2];
A=spdiags([-a_sub2,-a_sub1,a_diag,-a_super1,-a_super2], [-N,-1,0,1,N],N^2,N^2)/hg^2; %stiffness matrix
clear temp1 temp2 a1 a2 a_sub2 a_sub1 a_diag a_super1 a_super2
B=inv(A);

f=reshape(force_f([tt(:),ss(:)]),N,N+1); f=f(:,1:end-1); f= f(:);
u=A\f;

arrH=2.^(-1:-1:-4);nH=length(arrH);
ratio=[1,0]; n_ratio=length(ratio);
[arr_errH,arr_errL2,arr_err_comp]=deal(zeros(n_ratio,nH));

over_samp=1;
tic
for ii=1:n_ratio
    for jj=1:nH
        H=arrH(jj);h=H*ratio(ii);
        t=floor(over_samp*log(1/H)/log(2));
        sizeH=N*H; sizeh=N*h;
        if sizeh<1
            sizeh=1;
        end
        pde2d_compression_err_loc;
        arr_errH(ii,jj)=err_H;
        arr_errL2(ii,jj)=err_L2;
        arr_err_comp(ii,jj)=err_comp;
        fprintf('H=%g,h=%g,completed, running time %g s \n',H,h, toc);
    end
end

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
    eps=[1/5,1/13,1/17,1/31];
    y=1/6*((1.1+sin(2*pi*x(:,1)/eps(1)))./(1.1+cos(2*pi*x(:,2)/eps(1)))+...
        (1.1+sin(2*pi*x(:,2)/(eps(2))))./(1.1+cos(2*pi*x(:,1)/eps(2)))+...
        (1.1+cos(2*pi*x(:,1)/(eps(3))))./(1.1+sin(2*pi*x(:,2)/eps(3)))+...
        (1.1+sin(2*pi*x(:,2)/(eps(4))))./(1.1+cos(2*pi*x(:,1)/eps(4)))+sin(4*(x(:,1).^2).*(x(:,2)).^2)+1);
end

function [y]=force_f(x)
    k=100;
    W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
    tmp_cos=cos((1:k)'*x(:,1)');
    tmp_sin=sin((1:k)'*x(:,2)');
    y=1*sin(W1'*tmp_cos+W2'*tmp_sin); %row vector
end