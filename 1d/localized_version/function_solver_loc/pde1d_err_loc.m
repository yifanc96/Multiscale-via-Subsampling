%--------------------------------------------------------------------------
% fine grid size, coarse grid size, dim of coarse grid, num of coarse grid
d=1; dimH=sizeH^d; num=(N/sizeH)^d;

%--------------------------------------------------------------------------
% subsampling
dimh=(sizeh)^d;
idx_loc=floor((sizeH-sizeh)/2)+1:floor((sizeH+sizeh)/2);

%--------------------------------------------------------------------------
% construct phi psi
v=ones(dimh,1); v=house(v); U=eye(dimh)-2*(v*v');
psi=zeros(N,num);
for iii=1:num
    % for the ii-th patch, the coarse idx will be max(ii-t,1) to
    % min(ii+t,num), and the fine index will be (..-1)*sizeH+1
    % num_c1=l=num of left patch; num_cr=num of right patch
    num_cl=iii-max([iii-t,1]); num_cr=min([iii+t,num])-iii; 
    num_c=num_cl+num_cr+1; dim_ii=num_c*sizeH;
    index_ii=(iii-num_cl-1)*sizeH+1:1:(iii+num_cr)*sizeH;
    Aii=A(index_ii,index_ii); %localized A
    Us=blkdiag(eye(ceil((dimH-dimh)/2)),U(:,2:dimh),eye(floor((dimH-dimh)/2)));
    Us=kron(eye(num_c),Us);
    i=idx_loc'+num_cl*sizeH; j=ones(dimh,1);
    phi_ii=sparse(i,j,1,dim_ii,1); %localized phi
    Ah=Us'*Aii*Us; psi_ii=phi_ii-Us*(Ah\(Us'*(Aii*phi_ii)));
    psi(:,iii)=sparse(index_ii,ones(dim_ii,1),psi_ii,N,1); %lift to the global domain
end

%--------------------------------------------------------------------------
% solution
Ap=psi'*A*psi; uh=psi*(Ap\(psi'*f)); uh=[0;uh;0];

err_H=norm(gradient(u-uh,hg),2)*sqrt(hg);
err_L2=norm(u-uh,2)*sqrt(hg);

function [u]=house(v)
n=size(v,1);
sgn=sign(v(1));
u=(v+sgn*norm(v)*eye(n,1));
u=u./norm(u);
end
