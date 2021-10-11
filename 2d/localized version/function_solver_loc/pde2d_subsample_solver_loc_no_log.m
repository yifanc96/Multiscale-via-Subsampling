function [ uh ] = pde2d_subsample_solver_loc_no_log( A,f,N,sizeH,sizeh, over_samp)
d=2; 
dimH=sizeH^d; num=N^d/dimH; dimh=sizeh^d; H=sizeH/N;

% construct phi psi, from left-bottom patch to the whole domain
coarse_patch=zeros(sizeH,sizeH);
idx_loc=floor((sizeH-sizeh)/2)+1:floor((sizeH+sizeh)/2);
coarse_patch(idx_loc,idx_loc)=1;
coarse_patch=coarse_patch(:);
[idx_loc,~,~]=find(coarse_patch); [idx_loc2,~,~]=find(~coarse_patch);
t=over_samp;

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
%solution error
Ap=psi'*A*psi;
uh=psi*(Ap\(psi'*f));
end

function [u]=house(v)
n=size(v,1);
sgn=sign(v(1));
u=(v+sgn*norm(v)*eye(n,1));
u=u./norm(u);
end

