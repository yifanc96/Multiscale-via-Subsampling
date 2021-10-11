% 1d compression error

N=2^11; h=1/(N+1);
num=8;
arr_H=zeros(3,num);
arr_L2=zeros(3,num);
arr_comp=zeros(3,num);

t=(0:1:N)+0.5; v=kappa(t*h)'/h;
A=spdiags([-v(2:N+1),v(1:N)+v(2:N+1),-v(1:N)],-1:1,N,N)/h;
B=inv(full(A));
% A=A'*A;

% f=force_f((1:N)*h)'*h;
% f=randn(N,1)*h;
f=force_f((1:N)*h)';

for jj=0:2 % subsampling rate 2^(-jj)
    jj
    for ii=1:8  %H=2^(-ii)
        dimH=N*2^(-ii); H=dimH*h; num=N/dimH;
        rate=2^(-jj); %subsampling rate
        Hs=H*rate; dimHs=dimH*rate;
        
        v=ones(dimHs,1); u=house(v); U=eye(dimHs)-2*(u*u');
        Us=blkdiag(eye((dimH-dimHs)/2),U(:,2:dimHs),eye((dimH-dimHs)/2));
        Us=kron(eye(num),Us);
        Ah=Us'*A*Us;
        
        i=(dimH-dimHs)/2+1:(dimH+dimHs)/2; i=reshape(bsxfun(@plus,i',0:dimH:N-dimH),[],1);
        j=reshape(repmat(1:num,dimHs,1),[],1); phi=sparse(i,j,1,N,num);
        psi=Us*(Ah\(Us'*(A*phi))); psi=phi-psi;
        
        %solution error
        u=A\f; u=[0;u;0];
        uh=psi*((psi'*A*psi)\(psi'*f)); uh=[0;uh;0];
        arr_H(jj+1,ii)=norm(gradient(u-uh,h),2)*sqrt(h);
        arr_L2(jj+1,ii)=norm(u-uh,2)*sqrt(h);
        
        %compression error
        arr_comp(jj+1,ii)=norm(B-full(psi*((psi'*A*psi)\psi')),2);
    end
end

figure
plot(0:h:1,kappa(0:h:1))


figure
Harray=2.^(-8:1:-1);
loglog(Harray',flip(arr_H'),'-o');
legend('H/h=1','H/h=2','H/h=4')
xlabel('H')
ylabel('H^1 error')

figure
Harray=2.^(-8:1:-1);
loglog(Harray',flip(arr_L2'),'-o');
legend('H/h=1','H/h=2','H/h=4')
xlabel('H')
ylabel('L^2 error')

figure
Harray=2.^(-8:1:-1);
loglog(Harray',flip(arr_comp'),'-o');
legend('H/h=1','H/h=2','H/h=4')
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
    k=500;
    W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
    tmp_cos=cos((1:k)'*x);
    tmp_sin=sin((1:k)'*x);
    y=1+0.5*sin(W1'*tmp_cos+W2'*tmp_sin);%row vector
end


function [u]=house(v)
	n=size(v,1);
	sgn=sign(v(1));
	u=(v+sgn*norm(v)*eye(n,1));
	u=u./norm(u);
end

