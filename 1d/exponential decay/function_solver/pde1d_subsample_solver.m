function [uh] = pde1d_subsample_solver(A,f,N,sizeH,sizeh)

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
Us=blkdiag(eye(ceil((dimH-dimh)/2)),U(:,2:dimh),eye(floor((dimH-dimh)/2)));
Us=kron(eye(num),Us);
i=reshape(bsxfun(@plus,idx_loc',0:dimH:N-dimH),[],1);
j=reshape(repmat(1:num,dimh,1),[],1);
phi=sparse(i,j,1,N,num);
Ah=Us'*A*Us; psi=phi-Us*(Ah\(Us'*(A*phi)));

%--------------------------------------------------------------------------
% solution
Ap=psi'*A*psi; uh=psi*(Ap\(psi'*f)); uh=[0;uh;0];
end

function [u]=house(v)
n=size(v,1);
sgn=sign(v(1));
u=(v+sgn*norm(v)*eye(n,1));
u=u./norm(u);
end
