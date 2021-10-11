clear

d=2; N=2^6; hg=1/(N+1); %size N, total N^2 points
% stiffness matrix
t_mid=(0:1:N)'+0.5; t_grid=1:N; [tt,ss]=meshgrid(t_mid*hg,t_grid*hg);
% a1=kappa([tt(:),ss(:)]); a2=kappa([ss(:),tt(:)]); 
% a1=reshape(a1,N,N+1);a2=reshape(a2,N+1,N);
% a_diag=reshape(a1(:,1:N)+a1(:,2:N+1)+a2(1:N,:)+a2(2:N+1,:),[],1);
% temp1=reshape([a2(2:N,:);zeros(1,N)],[],1); temp1=temp1(1:N^2-1);
% temp2=reshape(a1(:,2:N),[],1);
% a_sub1=[temp1;0]; a_super1=[0;temp1];
% a_sub2=[temp2;zeros(N,1)]; a_super2=[zeros(N,1);temp2];
% A=spdiags([-a_sub2,-a_sub1,a_diag,-a_super1,-a_super2], [-N,-1,0,1,N],N^2,N^2)/hg^2; %stiffness matrix
% clear temp1 temp2 a1 a2 a_sub2 a_sub1 a_diag a_super1 a_super2

A=stiff_fem_q(N+1)/hg^2;

B=inv(A);

f=reshape(force_f([tt(:),ss(:)]),N,N+1); f=f(:,1:end-1); f= f(:);
u=A\f;

arrH=2.^(-2:-1:-4);nH=length(arrH);
ratio=[1,1/2,1/4]; n_ratio=length(ratio);
[arr_errH,arr_errL2,arr_err_comp]=deal(zeros(n_ratio,nH));

over_samp=2;
tic
for ii=1:n_ratio
    for jj=1:nH
        H=arrH(jj);h=H*ratio(ii);
%         t=over_samp;
        t=floor(over_samp*log(1/H)/log(2));
        sizeH=N*H; sizeh=N*h;
        pde2d_compression_err_loc;
        arr_errH(ii,jj)=err_H;
        arr_errL2(ii,jj)=err_L2;
        arr_err_comp(ii,jj)=err_comp;
        fprintf('H=%g,h=%g,completed, running time %g s \n',H,h, toc);
    end
end

figure
loglog(arrH',arr_errH','-o');
legend('H/h=1','H/h=2','H/h=4')
xlabel('H')
ylabel('H^1 error')

figure
loglog(arrH',arr_errL2','-o');
legend('H/h=1','H/h=2','H/h=4')
xlabel('H')
ylabel('L^2 error')

figure
loglog(arrH',arr_err_comp','-o');
legend('H/h=1','H/h=2','H/h=4')
xlabel('H')
ylabel('compression error')

function [y]=kappa(x)
    eps=[1/5,1/13,1/15,1/25];
    y=1/6*((1.1+sin(2*pi*x(:,1)/eps(1)))./(1.1+cos(2*pi*x(:,2)/eps(1)))+...
        (1.1+sin(2*pi*x(:,2)/(eps(2))))./(1.1+cos(2*pi*x(:,1)/eps(2)))+...
        (1.1+cos(2*pi*x(:,1)/(eps(3))))./(1.1+sin(2*pi*x(:,2)/eps(3)))+...
        (1.1+sin(2*pi*x(:,2)/(eps(4))))./(1.1+cos(2*pi*x(:,1)/eps(4)))+sin(4*(x(:,1).^2).*(x(:,2)).^2)+1);
end

function [y]=force_f(x)
    k=30;
    W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
    tmp_cos=cos((1:k)'*x(:,1)');
    tmp_sin=sin((1:k)'*x(:,2)');
    y=1*sin(W1'*tmp_cos+W2'*tmp_sin); %row vector
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