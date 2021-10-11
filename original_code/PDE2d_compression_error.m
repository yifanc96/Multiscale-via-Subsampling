%compression error
clear
N=2^6; h=1/(N+1); %total N^2*N^2
d=2;

t_mid=(0:1:N)'+0.5; t_grid=1:N; [tt,ss]=meshgrid(t_mid*h,t_grid*h);
a1=kappa([tt(:),ss(:)]); a2=kappa([ss(:),tt(:)]); 
a1=reshape(a1,N,N+1);a2=reshape(a2,N+1,N);
a_diag=reshape(a1(:,1:N)+a1(:,2:N+1)+a2(1:N,:)+a2(2:N+1,:),[],1);
temp1=reshape([a2(2:N,:);zeros(1,N)],[],1); temp1=temp1(1:N^2-1);
temp2=reshape(a1(:,2:N),[],1);
a_sub1=[temp1;0]; a_super1=[0;temp1];
a_sub2=[temp2;zeros(N,1)]; a_super2=[zeros(N,1);temp2];
clear temp1 temp2
A=spdiags([-a_sub2,-a_sub1,a_diag,-a_super1,-a_super2], [-N,-1,0,1,N],N^2,N^2)/h^2; %stiffness matrix
B=inv(full(A));

f=reshape(force_f([tt(:),ss(:)]),N,N+1); f=f(:,1:end-1); f= f(:);
num=4;
arr_H=zeros(3,num);
arr_L2=zeros(3,num);
arr_comp=zeros(3,num);

for jj=0:2
    for ii=1:4
        ii
        sizeH=N*2^(-ii); H=sizeH*h; dimH=sizeH^d; num=N^d/dimH;
        rate=2^(-jj); %subsampling rate
        Hs=H*rate; dimHs=dimH*rate^d; sizeHs=sizeH*rate;
        v=ones(dimHs,1); u=house(v); U=eye(dimHs)-2*(u*u'); U=U(:,2:dimHs);
        
        %subsampled location
        temp=zeros(sizeH,sizeH); 
        if dimHs>1
            idx_loc=(sizeH-sizeHs)/2+1:(sizeH+sizeHs)/2; 
        elseif dimHs==1
            idx_loc=sizeH/2+1;
        end
        temp(idx_loc,idx_loc)=1;
        temp=temp(:);[idx_loc,~,~]=find(temp); [idx_loc2,~,~]=find(~temp);
        idx_global=reshape(bsxfun(@plus,(1:sizeH)',0:N:N*(sizeH-1)),[],1);
        i=[repmat(idx_global(idx_loc),dimHs-1,1);idx_global(idx_loc2)];
        j=[reshape(repmat(1:dimHs-1,dimHs,1),[],1);(dimHs:dimH-1)'];
        k=[U(:);ones(dimH-dimHs,1)];
        shift=reshape(bsxfun(@plus,(0:sizeH:(N-sizeH))',0:sizeH*N:(N-sizeH)*N),[],1);
        i=reshape(bsxfun(@plus,i,shift'),[],1);
        j=reshape(bsxfun(@plus,j,0:max(j):max(j)*(num-1)),[],1);
        k=repmat(k,num,1);
        Us=sparse(i,j,k,N^2,max(j));
        Ah=Us'*A*Us;
        
        %phi
        i=reshape(bsxfun(@plus,idx_global(idx_loc),shift'),[],1);
        j=reshape(repmat(1:num,dimHs,1),[],1); phi=sparse(i,j,1,N^2,num);
        psi=Us*(Ah\(Us'*(A*phi))); psi=phi-psi;
            
        %solution error
        u=A\f; 
        Ap=psi'*A*psi;
        uh=psi*(Ap\(psi'*f)); 
        arr_H(jj+1,ii)=sqrt(sum(sum(gradient(reshape(u-uh,N,N),h,h).^2)))*h;
        arr_L2(jj+1,ii)=sqrt(sum((u-uh).^2))*h;
        
        %compression error
        arr_comp(jj+1,ii)=norm(B-full(psi*(Ap\psi')),2);
    end
end

% xx=t_grid*h; yy=t_grid*h;
% figure;
% s=surf(xx,yy,reshape(psi(:,num/2),N,N),'EdgeColor','none');

% figure
% Harray=2.^(-5:1:-1);
% loglog(Harray',flip(arr_H'),'-o');
% legend('H/h=1','H/h=2')
% xlabel('H')
% ylabel('H^1 error')
% 
% figure
% Harray=2.^(-5:1:-1);
% loglog(Harray',flip(arr_L2'),'-o');
% legend('H/h=1','H/h=2')
% xlabel('H')
% ylabel('L^2 error')
% 
% figure
% Harray=2.^(-5:1:-1);
% loglog(Harray',flip(arr_comp'),'-o');
% legend('H/h=1','H/h=2')
% xlabel('H')
% ylabel('compression error')

function [y]=kappa(x)
    eps=[1/5,1/13,1/17,1/31];
    y=1/6*((1.1+sin(2*pi*x(:,1)/eps(1)))./(1.1+cos(2*pi*x(:,2)/eps(1)))+...
        (1.1+sin(2*pi*x(:,2)/(eps(2))))./(1.1+cos(2*pi*x(:,1)/eps(2)))+...
        (1.1+cos(2*pi*x(:,1)/(eps(3))))./(1.1+sin(2*pi*x(:,2)/eps(3)))+...
        (1.1+sin(2*pi*x(:,2)/(eps(4))))./(1.1+cos(2*pi*x(:,1)/eps(4)))+sin(4*(x(:,1).^2).*(x(:,2)).^2)+1);
end

function [y]=force_f(x)
    k=50;
    W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
    tmp_cos=cos((1:k)'*x(:,1)');
    tmp_sin=sin((1:k)'*x(:,2)');
    y=1+0.5*sin(W1'*tmp_cos+W2'*tmp_sin); %row vector
end

function [u]=house(v)
	n=size(v,1);
	sgn=sign(v(1));
	u=(v+sgn*norm(v)*eye(n,1));
	u=u./norm(u);
end