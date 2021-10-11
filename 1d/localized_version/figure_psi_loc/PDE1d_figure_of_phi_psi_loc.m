clear

%--------------------------------------------------------------------------
% fine grid size, coarse grid size, dim of coarse grid, num of coarse grid
d=1; N=2^10; hg=1/(N+1);
sizeH=2^5; dimH=sizeH^d; H=sizeH*hg; num=(N/sizeH)^d;
fprintf('d=%.0f, grid points %g, loc-dom size %g\n',d,N,sizeH);

%--------------------------------------------------------------------------
% stiffness matrix in hg (finite difference)
t=(0:1:N)+0.5; v=kappa(t*hg)'; 
A=spdiags([-v(2:N+1),v(1:N)+v(2:N+1),-v(1:N)],-1:1,N,N)/hg^2;
% A=A'*A;  
clear v

%--------------------------------------------------------------------------
% subsampling
rate=2^(0); %subsampling rate
h=H*rate; sizeh=sizeH*rate; dimh=dimH*rate^d;
fprintf('subsampling ratio is %g, subsampled dom size %g\n',rate,sizeh);
idx_loc=floor((sizeH-sizeh)/2)+1:floor((sizeH+sizeh)/2);

%--------------------------------------------------------------------------
% construct phi psi
H=sizeH/N; t=floor(0.5*log(1/H)/log(2));
fprintf('oversampling coarse region size %g\n',t);
v=ones(dimh,1); v=house(v); U=eye(dimh)-2*(v*v');
psi=zeros(N,num);
for ii=1:num
    % for the ii-th patch, the coarse idx will be max(ii-t,1) to
    % min(ii+t,num), and the fine index will be (..-1)*sizeH+1
    % num_c1=l=num of left patch; num_cr=num of right patch
    num_cl=ii-max([ii-t,1]); num_cr=min([ii+t,num])-ii; 
    num_c=num_cl+num_cr+1; dim_ii=num_c*sizeH;
    index_ii=(ii-num_cl-1)*sizeH+1:1:(ii+num_cr)*sizeH;
    Aii=A(index_ii,index_ii); %localized A
    Us=blkdiag(eye(ceil((dimH-dimh)/2)),U(:,2:dimh),eye(floor((dimH-dimh)/2)));
    Us=kron(eye(num_c),Us);
    i=idx_loc'+num_cl*sizeH; j=ones(dimh,1);
    phi_ii=sparse(i,j,1,dim_ii,1); %localized phi
    Ah=Us'*Aii*Us; psi_ii=phi_ii-Us*(Ah\(Us'*(Aii*phi_ii)));
    psi(:,ii)=sparse(index_ii,ones(dim_ii,1),psi_ii,N,1); %lift to the global domain
end

i=reshape(bsxfun(@plus,idx_loc',0:dimH:N-dimH),[],1);
j=reshape(repmat(1:num,dimh,1),[],1); 
phi=sparse(i,j,1,N,num); 
phi_f=phi(:,num/2);
psi_f=psi(:,num/2);
%--------------------------------------------------------------------------
% figure 
figure
plot(0:hg:1,kappa(0:hg:1))
figure;
plot(0:hg:1,[0;phi_f;0]);
figure
plot(0:hg:1,[0;psi_f;0]);


function [y]=kappa(x)
    k=100;
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

