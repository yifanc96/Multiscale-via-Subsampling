clear

%--------------------------------------------------------------------------
% fine grid size, coarse grid size, dim of coarse grid, num of coarse grid
d=2; N=2^7; hg=1/(N+1); %size N, total N^2 points
sizeH=2^3; H=sizeH*hg; dimH=sizeH^d; num=N^d/dimH; 
fprintf('d=%.0f, grid points %.0f^2, loc-dom size %.0f\n',d,N,sizeH);


%--------------------------------------------------------------------------
A=stiff_fem_q(N+1);

%--------------------------------------------------------------------------
% subsampling
rate=2^(-3); %subsampling rate
h=H*rate; dimh=dimH*rate^d; sizeh=sizeH*rate;
fprintf('subsampling ratio is %g, subsampled dom size %g\n',rate,sizeh);
coarse_patch=zeros(sizeH,sizeH);
idx_loc=floor((sizeH-sizeh)/2)+1:floor((sizeH+sizeh)/2);
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

ind_phi=num/2+N/(2*sizeH)+1;  % index of selected phi
phi_f=phi(:,ind_phi);
position=idx_global+shift(ind_phi);
position=reshape(position,sizeH,sizeH);
%index of covered elements

psi=Us*(Ah\(Us'*(A*phi_f))); psi=phi_f-psi;
phi_f=reshape(phi_f,N,N);
psi=reshape(psi,N,N);


%--------------------------------------------------------------------------
% figure 
% t_grid=1:N;
% xx=t_grid*hg; yy=t_grid*hg;
% figure
% subplot(2,2,1)
% surf(xx,yy,phi_f,'EdgeColor','none');
% subplot(2,2,2)
% contourf(xx,yy,phi_f)
% subplot(2,2,3)
% surf(xx,yy,psi,'EdgeColor','none');
% subplot(2,2,4);
% surf(xx,yy,log(abs(psi)),'EdgeColor','none')
% shading interp

%quantify exponential decay rate
norm_sum=psi(:)'*A*(psi(:));
element_index_loc=elem_ind(position,N);
A_loc=stiff_fem_q_interior(N+1,element_index_loc);
A_locc=A_loc(position,position);

% norm_loc=sum(sum(norm_all(position)));
norm_loc=psi(position(:))'*A_locc'*psi(position(:));
decay_rate=1-norm_loc/norm_sum;
fprintf('Energy norm decay rate is %g\n', decay_rate); 

patch=ind_phi;
[Iy,Ix]=ind2sub([N/sizeH,N/sizeH],patch);
% t=floor(log(1/H)/log(2));
t=1;
Ix_l=max([1,Ix-t]); Ix_r=min([N/sizeH,Ix+t]);
Iy_u=min([N/sizeH,Iy+t]); Iy_d=max([1,Iy-t]);
iy=(Iy_d-1)*sizeH+1:1:(Iy_u)*sizeH; ly=length(iy);
ix=(Ix_l-1)*sizeH+1:1:(Ix_r)*sizeH; lx=length(ix); dim_loc=lx*ly;
iy_tot=repmat(iy',lx,1); ix_tot=reshape(repmat(ix,ly,1),[],1);
ind_tot=sub2ind([N,N],iy_tot,ix_tot);
position2=ind_tot;
position2=reshape(position2,3*sizeH,3*sizeH);
element_index_loc2=elem_ind(position2,N);
A_loc2=stiff_fem_q_interior(N+1,element_index_loc2);
A_locc2=A_loc2(position2,position2);
norm_loc2=psi(position2(:))'*A_locc2'*psi(position2(:));
decay_rate2=1-norm_loc2/norm_sum;
fprintf('One oversampling: Energy norm decay rate is %g\n', decay_rate2); 

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

function [elem_index]=elem_ind(position_nodes,N)
% return the element index of the patch_i-th patch domain
% input: position_nodes will be the node index of patch_i-th
 position=Dirichlet_in2all(position_nodes,N); 
 %position is the index in the whole nodes
 [m,n]=size(position);
 position=position(:);
 elem_index=[];
 for iter=1:m*(n-1)
     if mod(iter,m)~=0
         tmp=floor(position(iter)/(N+2));
         res=position(iter)-tmp*(N+2);
         elem_index=[elem_index;tmp*(N+1)+res];
     end
 end
end


function [ind_all]=Dirichlet_in2all(ind_in,N)
% input the index of interior nodes, output its index in all nodes 
% including the boundary nodes
% N is the number of interior nodes in each side
t=floor(ind_in/N);
ll=(t<ind_in/N);
ind_all=ind_in+((N+2)+2*t+1).*ll+((N+2)+2*t-1).*(1-ll);
end


function A=stiff_fem_q(n)
%reference rectangle: (0,0),(1,0),(1,1),(0,1)
num=4; %total num of element in the rectangle

%analytic element
phi1=@(x) (1-x(:,1)).*(1-x(:,2));
phi2=@(x) x(:,1).*(1-x(:,2));
phi3=@(x) x(:,1).*x(:,2);
phi4=@(x) (1-x(:,1)).*x(:,2);

%analytic gradient
grad1=@(x) [x(:,2)-1,x(:,1)-1];
grad2=@(x) [1-x(:,2),-x(:,1)];
grad3=@(x) [x(:,2),x(:,1)];
grad4=@(x) [-x(:,2),1-x(:,1)];

%numerical element and gradient inner-product
ref_x=0:0.5:1; ref_y=0:0.5:1; ref_length=length(ref_x); ref_dim=ref_length^2;
ref_x=repmat(ref_x',ref_length,1); ref_y=repmat(ref_y,ref_length,1); ref_y=ref_y(:);
ref_xy=[ref_x,ref_y]; %form the reference coordinate
ref_phi=zeros(ref_dim,num);
ref_phi(:,1)=phi1(ref_xy);
ref_phi(:,2)=phi2(ref_xy);
ref_phi(:,3)=phi3(ref_xy);
ref_phi(:,4)=phi4(ref_xy);
ref_grad=zeros(ref_dim,2,num);
ref_grad(:,:,1)=grad1(ref_xy);
ref_grad(:,:,2)=grad2(ref_xy);
ref_grad(:,:,3)=grad3(ref_xy);
ref_grad(:,:,4)=grad4(ref_xy);
ref_inprod=zeros(ref_dim,num^2);
for i = 1:num
    for j = 1:num
        ref_inprod(:,num*(i-1)+j)=sum(ref_grad(:,:,i).*ref_grad(:,:,j),2);
    end
end

% finite element method
h=1/n; N=n+1; %total nodes each coordinate
[x,y]=meshgrid((0:n)/n,(0:n)/(n));
p=[x(:),y(:)];
element=[1,N+1,N+2,2]; %the first element, counter-clockwise orientation
element=kron(element,ones(N-1,1))+kron(ones(size(element)),(0:N-2)');
element=kron(element,ones(N-1,1))+kron(ones(size(element)),(0:N-2)'*N); %form the whole elements
b=[1:N,N+1:N:N*N,2*N:N:N*N,N*N-N+2:N*N-1]; %boundary point
interior=[];
for ii=1:N^2
    if all(ii-b)
        interior=[interior,ii];
    end
end

Ne=size(element,1); %number of elements
I=zeros(16*Ne,1);J=I;K=I; %index of sparse stiffness matrix
                        %each element leads to 16 entries
w=[1/6,2/3,1/6]; w=w'*w; w=w(:); w=w';
for e=1:Ne
    nodes=element(e,:);
    points=p(nodes,:);
    new_points=bsxfun(@plus,points(1,:),h*ref_xy); %include the midpoint
    kappa_e=kappa(new_points);
    for i=1:4
        for j=1:4
            index=16*(e-1)+4*(i-1)+j;
            I(index)=nodes(i);
            J(index)=nodes(j);
            K(index)=w*(kappa_e.*ref_inprod(:,(i-1)*4+j));
        end
    end
end
A=sparse(I,J,K,N^2,N^2);
A(b,:)=0; A(:,b)=0; A(b,b)=speye(length(b),length(b));
A=A(interior,interior);

end


function A=stiff_fem_q_interior(n,elem_index)
%reference rectangle: (0,0),(1,0),(1,1),(0,1)
num=4; %total num of element in the rectangle

%analytic element
phi1=@(x) (1-x(:,1)).*(1-x(:,2));
phi2=@(x) x(:,1).*(1-x(:,2));
phi3=@(x) x(:,1).*x(:,2);
phi4=@(x) (1-x(:,1)).*x(:,2);

%analytic gradient
grad1=@(x) [x(:,2)-1,x(:,1)-1];
grad2=@(x) [1-x(:,2),-x(:,1)];
grad3=@(x) [x(:,2),x(:,1)];
grad4=@(x) [-x(:,2),1-x(:,1)];

%numerical element and gradient inner-product
ref_x=0:0.5:1; ref_y=0:0.5:1; ref_length=length(ref_x); ref_dim=ref_length^2;
ref_x=repmat(ref_x',ref_length,1); ref_y=repmat(ref_y,ref_length,1); ref_y=ref_y(:);
ref_xy=[ref_x,ref_y]; %form the reference coordinate
ref_phi=zeros(ref_dim,num);
ref_phi(:,1)=phi1(ref_xy);
ref_phi(:,2)=phi2(ref_xy);
ref_phi(:,3)=phi3(ref_xy);
ref_phi(:,4)=phi4(ref_xy);
ref_grad=zeros(ref_dim,2,num);
ref_grad(:,:,1)=grad1(ref_xy);
ref_grad(:,:,2)=grad2(ref_xy);
ref_grad(:,:,3)=grad3(ref_xy);
ref_grad(:,:,4)=grad4(ref_xy);
ref_inprod=zeros(ref_dim,num^2);
for i = 1:num
    for j = 1:num
        ref_inprod(:,num*(i-1)+j)=sum(ref_grad(:,:,i).*ref_grad(:,:,j),2);
    end
end

% finite element method
h=1/n; N=n+1; %total nodes each coordinate
[x,y]=meshgrid((0:n)/n,(0:n)/(n));
p=[x(:),y(:)];
element=[1,N+1,N+2,2]; %the first element, counter-clockwise orientation
element=kron(element,ones(N-1,1))+kron(ones(size(element)),(0:N-2)');
element=kron(element,ones(N-1,1))+kron(ones(size(element)),(0:N-2)'*N); %form the whole elements
b=[1:N,N+1:N:N*N,2*N:N:N*N,N*N-N+2:N*N-1]; %boundary point
interior=[];
for ii=1:N^2
    if all(ii-b)
        interior=[interior,ii];
    end
end

w=[1/6,2/3,1/6]; w=w'*w; w=w(:); w=w';
Ne=length(elem_index);
I=zeros(16*Ne,1);J=I;K=I; %index of sparse stiffness matrix
                        %each element leads to 16 entries
for e=1:Ne
    nodes=element(elem_index(e),:);
    points=p(nodes,:);
    new_points=bsxfun(@plus,points(1,:),h*ref_xy); %include the midpoint
    kappa_e=kappa(new_points);
    for i=1:4
        for j=1:4
            index=16*(e-1)+4*(i-1)+j;
            I(index)=nodes(i);
            J(index)=nodes(j);
            K(index)=w*(kappa_e.*ref_inprod(:,(i-1)*4+j));
        end
    end
end
A=sparse(I,J,K,N^2,N^2);
A(b,:)=0; A(:,b)=0; A(b,b)=speye(length(b),length(b));
A=A(interior,interior);
end