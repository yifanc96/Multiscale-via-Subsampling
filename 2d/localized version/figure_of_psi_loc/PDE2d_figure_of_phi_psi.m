clear

%--------------------------------------------------------------------------
% fine grid size, coarse grid size, dim of coarse grid, num of coarse grid
d=2; N=2^6; hg=1/(N+1); %size N, total N^2 points
sizeH=2^2; H=sizeH*hg; dimH=sizeH^d; num=N^d/dimH; 
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
rate=2^(-2); %subsampling rate
h=H*rate; dimh=dimH*rate^d; sizeh=sizeH*rate;
fprintf('subsampling ratio is %.2f, subsampled dom size %.3f\n',rate,sizeh);
coarse_patch=zeros(sizeH,sizeH);
idx_loc=floor((sizeH-sizeh)/2)+1:floor((sizeH+sizeh)/2);
coarse_patch(idx_loc,idx_loc)=1;


%--------------------------------------------------------------------------
% construct phi psi, from left-bottom patch to the whole domain
coarse_patch=coarse_patch(:);
[idx_loc,~,~]=find(coarse_patch); [idx_loc2,~,~]=find(~coarse_patch);
over_samp=1; t=floor(over_samp*log(1/H)/log(2));
v=ones(dimh,1); v=house(v); U=eye(dimh)-2*(v*v'); U=U(:,2:dimh); %householder

psi_I=zeros((2*t+1)^2*sizeH^2,1);psi_J=psi_I;psi_K=psi_I; %sparse matrix
index=1;
for patch=1:num
   [Iy,Ix]=ind2sub([N/sizeH,N/sizeH],patch);
   Ix_l=max([1,Ix-t]); Ix_r=min([N/sizeH,Ix+t]);
   Iy_u=min([N/sizeH,Iy+t]); Iy_d=max([1,Iy-t]);
   iy=(Iy_d-1)*sizeH+1:1:(Iy_u)*sizeH; ly=length(iy);
   ix=(Ix_l-1)*sizeH+1:1:(Ix_r)*sizeH; lx=length(ix); dim_loc=lx*ly;
   iy_tot=repmat(iy',lx,1); ix_tot=reshape(repmat(ix,ly,1),[],1);
   ind_tot=sub2ind([N,N],iy_tot,ix_tot);
   
   %construct Us
   Lx=lx/sizeH; Ly=ly/sizeH; num_loc=Lx*Ly;
   idx_global=reshape(bsxfun(@plus,(1:sizeH)',0:ly:ly*(sizeH-1)),[],1);
   idx_global1=idx_global(idx_loc); index_global2=idx_global(idx_loc2);
   i=[repmat(idx_global1,dimh-1,1);index_global2];
   j=[reshape(repmat(1:dimh-1,dimh,1),[],1);(dimh:dimH-1)'];
   
   k=[U(:);ones(dimH-dimh,1)];
   shift=reshape(bsxfun(@plus,(0:sizeH:(ly-sizeH))',0:sizeH*ly:(lx-sizeH)*ly),[],1);
   i=reshape(bsxfun(@plus,i,shift'),[],1);
   j=reshape(bsxfun(@plus,j,0:max(j):max(j)*(num_loc-1)),[],1);
   k=repmat(k,num_loc,1);
   Us=sparse(i,j,k,dim_loc,max(j));
   A_loc=A(ind_tot,ind_tot);
   Ah=Us'*A_loc*Us;
   
   % construct phi
   i=reshape(bsxfun(@plus,idx_global(idx_loc),shift'),[],1);
   j=reshape(repmat(1:num_loc,dimh,1),[],1); phi=sparse(i,j,1,dim_loc,num_loc);
   temp=(Ix-Ix_l)*Ly+Iy-Iy_d+1; phi=phi(:,temp);
   psi_patch=Us*(Ah\(Us'*(A_loc*phi))); psi_patch=phi-psi_patch;
   
   len=length(ind_tot);
   psi_I(index:index+len-1)=ind_tot; psi_J(index:index+len-1)=patch*ones(length(ind_tot),1);
   psi_K(index:index+len-1)=psi_patch;
   index=index+len;
%    psi(:,patch)=sparse(ind_tot,ones(length(ind_tot),1),psi_patch,N^2,1);
end
psi=sparse(psi_I,psi_J,psi_K,N^2,num);


idx_global=reshape(bsxfun(@plus,(1:sizeH)',0:N:N*(sizeH-1)),[],1);
shift=reshape(bsxfun(@plus,(0:sizeH:(N-sizeH))',0:sizeH*N:(N-sizeH)*N),[],1);
i=reshape(bsxfun(@plus,idx_global(idx_loc),shift'),[],1);
j=reshape(repmat(1:num,dimh,1),[],1); phi=sparse(i,j,1,N^2,num); %whole phi
phi_f=phi(:,num/2-4);
psi_f=psi(:,num/2-4);
phi_f=reshape(phi_f,N,N);
psi_f=reshape(psi_f,N,N);

%--------------------------------------------------------------------------
% figure 
xx=t_grid*hg; yy=t_grid*hg;
figure(1);
surf(xx,yy,phi_f,'EdgeColor','none');
figure(2);
pcolor(xx,yy,phi_f)
figure(3);
surf(xx,yy,psi_f,'EdgeColor','none');
figure(4);
surf(xx,yy,log(abs(psi_f)),'EdgeColor','none')
%shading interp

figure;
pcolor(xx,yy,reshape(sum(phi,2),N,N))

% figure
% load force_data
% u=A\f;
% ind=find(sum(phi,2)==1);
% [mx,my]=meshgrid(xx,yy);
% mesh(mx,my,reshape(u,N,N)); hold on
% plot3(mx(ind),my(ind),u(ind),'.');
% save subsample_plot mx my ind

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