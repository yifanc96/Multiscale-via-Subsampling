clear

%--------------------------------------------------------------------------
% fine grid size, coarse grid size, dim of coarse grid, num of coarse grid
d=1; N=2^12; hg=1/(N+1);
sizeH=2^7; dimH=sizeH^d; H=sizeH*hg; num=(N/sizeH)^d;
fprintf('d=%.0f, grid points %.0f, loc-dom size %.0f\n',d,N,sizeH);

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
fprintf('subsampling ratio is %.2f, subsampled dom size %.3f\n',rate,sizeh);
coarse_patch=zeros(sizeH); 
% if dimh>1
%     idx_loc=(sizeH-sizeh)/2+1:(sizeH+sizeh)/2; 
% elseif dimh==1
%     idx_loc=sizeH/2+1;
% elseif dimh<1
%     fprintf('subsampling rate is too small!\n');
%     return
% end
idx_loc=floor((sizeH-sizeh)/2)+1:floor((sizeH+sizeh)/2);
coarse_patch(idx_loc)=1; 

%--------------------------------------------------------------------------
% construct phi psi
v=ones(dimh,1); v=house(v); U=eye(dimh)-2*(v*v');
Us=blkdiag(eye((dimH-dimh)/2),U(:,2:dimh),eye((dimH-dimh)/2));
Us=kron(eye(num),Us);

i=reshape(bsxfun(@plus,idx_loc',0:dimH:N-dimH),[],1);
j=reshape(repmat(1:num,dimh,1),[],1); 
phi=sparse(i,j,1,N,num); phi_f=phi(:,num/2);
Ah=Us'*A*Us; psi=phi_f-Us*(Ah\(Us'*(A*phi_f)));

%--------------------------------------------------------------------------
% figure 
figure
plot(0:hg:1,kappa(0:hg:1))
figure;
plot(0:hg:1,[0;phi_f;0]);
figure
plot(0:hg:1,[0;psi;0]);


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

