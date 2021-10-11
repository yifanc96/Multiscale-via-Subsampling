clear

%--------------------------------------------------------------------------
% fine grid size, coarse grid size, dim of coarse grid, num of coarse grid
d=1; N=2^11; hg=1/(N+1);
sizeH=2^4; dimH=sizeH^d; H=sizeH*hg; num=(N/sizeH)^d;
fprintf('d=%g, grid points %g, loc-dom size %g\n',d,N,sizeH);

%--------------------------------------------------------------------------
% stiffness matrix in hg (finite difference), rhs
t=(0:1:N)+0.5; v=kappa(t*hg)'; 
A=spdiags([-v(2:N+1),v(1:N)+v(2:N+1),-v(1:N)],-1:1,N,N)/hg^2; clear v
% A=A'*A;  
f=force_f((1:N)*hg)';

%--------------------------------------------------------------------------
% subsampling
rate=2^(-1); %subsampling rate
h=H*rate; sizeh=sizeH*rate; dimh=dimH*rate^d;
fprintf('subsampling ratio is %g, subsampled dom size %g\n',rate,sizeh);
coarse_patch=zeros(sizeH); 
if dimh>1
    idx_loc=(sizeH-sizeh)/2+1:(sizeH+sizeh)/2; 
elseif dimh==1
    idx_loc=sizeH/2+1;
elseif dimh<1
    fprintf('subsampling rate is too small!\n');
    return
end
coarse_patch(idx_loc)=1; 

%--------------------------------------------------------------------------
% construct phi psi
v=ones(dimh,1); u=house(v); U=eye(dimh)-2*(u*u');
Us=blkdiag(eye(ceil((dimH-dimh)/2)),U(:,2:dimh),eye(floor((dimH-dimh)/2)));
Us=kron(eye(num),Us);
i=reshape(bsxfun(@plus,idx_loc',0:dimH:N-dimH),[],1);
j=reshape(repmat(1:num,dimh,1),[],1); 
phi=sparse(i,j,1,N,num); 
Ah=Us'*A*Us; psi=phi-Us*(Ah\(Us'*(A*phi)));

%--------------------------------------------------------------------------
% solution
u=A\f; u=[0;u;0];
Ap=psi'*A*psi;
uh=psi*(Ap\(psi'*f)); uh=[0;uh;0];
errH=norm(gradient(u-uh,hg),2)*sqrt(hg);
errL2=norm(u-uh,2)*sqrt(hg);
B=inv(full(A));
err_comp=norm(B-full(psi*(Ap\psi')),2);
fprintf('H^1 err %g, L^2 err %g, compression err %g\n',errH,errL2,err_comp);

%figure
figure
plot(0:hg:1,kappa(0:hg:1))
figure
plot(0:hg:1,u);
hold on
plot(0:hg:1,uh);
legend('true solution','multiscale method solution');


function [y]=kappa(x)
    k=100;
    W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
    tmp_cos=cos((1:k)'*x);
    tmp_sin=sin((1:k)'*x);
    y=1+0.5*sin(W1'*tmp_cos+W2'*tmp_sin);%row vector
end

function [y]=force_f(x)
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

