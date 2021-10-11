clear

d=2; 
N=2^5;
hg=1/(N+1); %size N, total N^2 points
sizeH=2^3; 
sizeh=sizeH/2^3;


% stiffness matrix
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

a1=kappa_weighted([tt(:),ss(:)]); a2=kappa([ss(:),tt(:)]); 
a1=reshape(a1,N,N+1);a2=reshape(a2,N+1,N);
a_diag=reshape(a1(:,1:N)+a1(:,2:N+1)+a2(1:N,:)+a2(2:N+1,:),[],1);
temp1=reshape([a2(2:N,:);zeros(1,N)],[],1); temp1=temp1(1:N^2-1);
temp2=reshape(a1(:,2:N),[],1);
a_sub1=[temp1;0]; a_super1=[0;temp1];
a_sub2=[temp2;zeros(N,1)]; a_super2=[zeros(N,1);temp2];
A_weighted=spdiags([-a_sub2,-a_sub1,a_diag,-a_super1,-a_super2], [-N,-1,0,1,N],N^2,N^2)/hg^2; %stiffness matrix
clear temp1 temp2 a1 a2 a_sub2 a_sub1 a_diag a_super1 a_super2


f=reshape(force_f([tt(:),ss(:)]),N,N+1); f=f(:,1:end-1); f= f(:);
% load force_data
u=A\f;
tic
uh=interpolation_2d_subsample_loc(A,u,N,sizeH,sizeh,1);
toc
disp(toc)
u_weighted=interpolation_2d_subsample_loc(A_weighted,u,N,sizeH,sizeh,1);
toc


err_H=sqrt(sum(sum(gradient(reshape(u-uh,N,N),hg,hg).^2)))*hg;
err_L2=norm(u-uh,2)*hg;
fprintf('H^1 error is %g, L^2 error is %g, time %g \n',err_H,err_L2,toc);
xx=t_grid*hg; yy=t_grid*hg;

% figure;
% contourf(xx,yy,reshape(u,N,N))
figure;
surf(xx,yy,reshape(u,N,N))
shading interp
% figure;
% contourf(xx,yy,reshape(uh,N,N))
figure;
surf(xx,yy,reshape(uh,N,N))
shading interp
% figure;
% contourf(xx,yy,reshape(u_weighted,N,N))
figure;
surf(xx,yy,reshape(u_weighted,N,N))
shading interp

function [y]=kappa(x)
%     eps=[1/5,1/13,1/17,1/31];
%     y=1/6*((1.1+sin(2*pi*x(:,1)/eps(1)))./(1.1+cos(2*pi*x(:,2)/eps(1)))+...
%         (1.1+sin(2*pi*x(:,2)/(eps(2))))./(1.1+cos(2*pi*x(:,1)/eps(2)))+...
%         (1.1+cos(2*pi*x(:,1)/(eps(3))))./(1.1+sin(2*pi*x(:,2)/eps(3)))+...
%         (1.1+sin(2*pi*x(:,2)/(eps(4))))./(1.1+cos(2*pi*x(:,1)/eps(4)))+sin(4*(x(:,1).^2).*(x(:,2)).^2)+1);
    y=ones(size(x,1),1);
end

function [y]=kappa_weighted(x)
N=2^5;
hg=1/(N+1); %size N, total N^2 points
sizeH=2^3;
sizeh=sizeH/2^3;
    cellnum=N/sizeH; cell_H=sizeH*hg; cell_h=sizeh*hg*0.5;
%     eps=[1/5,1/13,1/17,1/31];
%     y=1/6*((1.1+sin(2*pi*x(:,1)/eps(1)))./(1.1+cos(2*pi*x(:,2)/eps(1)))+...
%         (1.1+sin(2*pi*x(:,2)/(eps(2))))./(1.1+cos(2*pi*x(:,1)/eps(2)))+...
%         (1.1+cos(2*pi*x(:,1)/(eps(3))))./(1.1+sin(2*pi*x(:,2)/eps(3)))+...
%         (1.1+sin(2*pi*x(:,2)/(eps(4))))./(1.1+cos(2*pi*x(:,1)/eps(4)))+sin(4*(x(:,1).^2).*(x(:,2)).^2)+1);
    y=ones(size(x,1),1);
    for ii=1:length(y)
        ix=floor(x(ii,1)/cell_H); iy=floor(x(ii,2)/cell_H);
        central_point=[(ix+1/2)*cell_H, (iy+1/2)*cell_H];
        ratio=cell_H/max([norm(central_point-x(ii,:)),cell_h]);
        y(ii)=ratio*log(ratio)*y(ii);
    end
end


function [y]=force_f(x)
    k=50;
    W1=rand(k,1)-0.5; W2=rand(k,1)-0.5;
    tmp_cos=cos((1:k)'*x(:,1)');
    tmp_sin=sin((1:k)'*x(:,2)');
    y=0.4*sin(W1'*tmp_cos+W2'*tmp_sin); %row vector
end