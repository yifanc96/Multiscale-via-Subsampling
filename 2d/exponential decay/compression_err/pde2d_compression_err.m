d=2;
dimH=sizeH^d; num=N^d/dimH; dimh=sizeh^d; 
v=ones(dimh,1); v=house(v); U=eye(dimh)-2*(v*v'); U=U(:,2:dimh);

%subsampled location
temp=zeros(sizeH,sizeH);
idx_loc=floor((sizeH-sizeh)/2)+1:floor((sizeH+sizeh)/2);
temp(idx_loc,idx_loc)=1;
temp=temp(:);[idx_loc,~,~]=find(temp); [idx_loc2,~,~]=find(~temp);

idx_global=reshape(bsxfun(@plus,(1:sizeH)',0:N:N*(sizeH-1)),[],1);
i=[repmat(idx_global(idx_loc),dimh-1,1);idx_global(idx_loc2)];
j=[reshape(repmat(1:dimh-1,dimh,1),[],1);(dimh:dimH-1)'];
k=[U(:);ones(dimH-dimh,1)];
shift=reshape(bsxfun(@plus,(0:sizeH:(N-sizeH))',0:sizeH*N:(N-sizeH)*N),[],1);
i=reshape(bsxfun(@plus,i,shift'),[],1);
j=reshape(bsxfun(@plus,j,0:max(j):max(j)*(num-1)),[],1);
k=repmat(k,num,1);
Us=sparse(i,j,k,N^2,max(j));
Ah=Us'*A*Us;

%phi
i=reshape(bsxfun(@plus,idx_global(idx_loc),shift'),[],1);
j=reshape(repmat(1:num,dimh,1),[],1); phi=sparse(i,j,1,N^2,num);
psi=Us*(Ah\(Us'*(A*phi))); psi=phi-psi;

%solution error
Ap=psi'*A*psi;
uh=psi*(Ap\(psi'*f));
err_H=sqrt(sum(sum(gradient(reshape(u-uh,N,N),hg,hg).^2)))*hg;
err_L2=norm(u-uh,2)*hg;
%compression error
err_comp=norm(B-full(psi*(Ap\psi')),2);

function [u]=house(v)
n=size(v,1);
sgn=sign(v(1));
u=(v+sgn*norm(v)*eye(n,1));
u=u./norm(u);
end