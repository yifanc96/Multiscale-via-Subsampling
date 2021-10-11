clear

%--------------------------------------------------------------------------
% fine grid size, coarse grid size, dim of coarse grid, num of coarse grid
d=2; N=2^6; hg=1/(N+1); %size N, total N^2 points
sizeH=2^3; H=sizeH*hg; dimH=sizeH^d; num=N^d/dimH; 
fprintf('d=%.0f, grid points %.0f^2, loc-dom size %.0f\n',d,N,sizeH);


%--------------------------------------------------------------------------
% stiffness matrix in hg (finite difference)
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


%--------------------------------------------------------------------------
% subsampling
rate=2^(-3); %subsampling rate
h=H*rate; dimh=dimH*rate^d; sizeh=sizeH*rate;
fprintf('subsampling ratio is %.2f, subsampled dom size %.3f\n',rate,sizeh);
coarse_patch=zeros(sizeH,sizeH);
if dimh>1
    idx_loc=(sizeH-sizeh)/2+1:(sizeH+sizeh)/2; 
elseif dimh==1
    idx_loc=sizeH/2+1;
elseif dimh<1
    fprintf('subsampling rate is too small!\n');
    return
end
coarse_patch(idx_loc,idx_loc)=1;


%--------------------------------------------------------------------------
% construct phi psi, from left-bottom patch to the whole domain
coarse_patch=coarse_patch(:);
[idx_loc,~,~]=find(coarse_patch); [idx_loc2,~,~]=find(~coarse_patch);
idx_global=reshape(bsxfun(@plus,(1:sizeH)',0:N:N*(sizeH-1)),[],1);
idx_global1=idx_global(idx_loc); index_global2=idx_global(idx_loc2);
i=[repmat(idx_global1,dimh-1,1);index_global2]; 
j=[reshape(repmat(1:dimh-1,dimh,1),[],1);(dimh:dimH-1)'];
v=ones(dimh,1); u=house(v); U=eye(dimh)-2*(u*u'); U=U(:,2:dimh); %householder
k=[U(:);ones(dimH-dimh,1)]; 
shift=reshape(bsxfun(@plus,(0:sizeH:(N-sizeH))',0:sizeH*N:(N-sizeH)*N),[],1);
i=reshape(bsxfun(@plus,i,shift'),[],1);
j=reshape(bsxfun(@plus,j,0:max(j):max(j)*(num-1)),[],1);
k=repmat(k,num,1);
Us=sparse(i,j,k,N^2,max(j));
Ah=Us'*A*Us;

i=reshape(bsxfun(@plus,idx_global(idx_loc),shift'),[],1);
j=reshape(repmat(1:num,dimh,1),[],1); phi=sparse(i,j,1,N^2,num); %whole phi
phi_f=phi(:,num/2-4);
psi=Us*(Ah\(Us'*(A*phi_f))); psi=phi_f-psi;
phi_f=reshape(phi_f,N,N);
psi=reshape(psi,N,N);


%--------------------------------------------------------------------------
% figure 
xx=t_grid*hg; yy=t_grid*hg;
figure
subplot(2,2,1)
surf(xx,yy,phi_f,'EdgeColor','none');
subplot(2,2,2)
contourf(xx,yy,phi_f)
subplot(2,2,3)
s=surf(xx,yy,psi,'EdgeColor','none');
subplot(2,2,4);
surf(xx,yy,log(abs(psi)),'EdgeColor','none')
shading interp



function [y]=kappa(x)
    eps=[1/5,1/13,1/17,1/31];
    y=1/6*((1.1+sin(2*pi*x(:,1)/eps(1)))./(1.1+cos(2*pi*x(:,2)/eps(1)))+...
        (1.1+sin(2*pi*x(:,2)/(eps(2))))./(1.1+cos(2*pi*x(:,1)/eps(2)))+...
        (1.1+cos(2*pi*x(:,1)/(eps(3))))./(1.1+sin(2*pi*x(:,2)/eps(3)))+...
        (1.1+sin(2*pi*x(:,2)/(eps(4))))./(1.1+cos(2*pi*x(:,1)/eps(4)))+sin(4*(x(:,1).^2).*(x(:,2)).^2)+1);
end
function [u]=house(v)
	n=size(v,1);
	sgn=sign(v(1));
	u=(v+sgn*norm(v)*eye(n,1));
	u=u./norm(u);
end