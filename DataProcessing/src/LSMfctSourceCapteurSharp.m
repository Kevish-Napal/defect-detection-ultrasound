%function [Feta, g, HEST,ignormLSM,ignormGLSM]=LSMfct(index,sortie, params,h)
%function [Feta, g, HEST,ignormLSM,ignormGLSM,rhsSave]=LSMfct(index, params)
%[Feta, g1,g2,g,T1,T2,Morozov1, Morozov2,t, HEST]
function  [Fsharp, g,HEST,asharp]=LSMfctSourceCapteurSharp(params) %, F_sharp,aSharp,aLsm

global GUP  SQ HEST Feta  SIGMAN rhs aTK
nsource=params.nsource;
ncapteur=params.ncapteur;

source=params.source;
capteur=params.capteur;
%load(index)
%F=index;
k=params.k;
Nx=params.Nx;
xmin=params.xmin;
xmax=params.xmax;
Ny=params.Ny;
ymin=params.ymin;
ymax=params.ymax;

Dx=(-xmin+xmax)/(Nx-1);
Dy=(ymax-ymin)/(Ny-1);
[X,Y]=meshgrid(xmin:Dx:xmax,ymin:Dy:ymax);

% eta=params.eta; %0.001;
% noise=1+eta*(rand(size(F))*2-1)+sqrt(-1)*eta*(rand(size(F))*2-1);
% Feta=F.*noise;

[Vr,Dr]=eig((Feta+Feta')/2);
[Vi,Di]=eig((Feta-Feta')/(2*sqrt(-1)));
Fsharp=Vr*abs(Dr)*Vr^-1+Vi*abs(Di)*Vi^-1;


% HEST=norm(Feta-F); % le delta de Morozov
% 
% if HEST < 1.e-02
%     HEST= 1.e-08;
% end
% 
% nF=norm(Feta);
% disp(['Bruit= ',num2str(HEST/nF*100)])
options=optimset('Display','off');

nFsharp = norm(Fsharp);
asharp = zeros(Nx,Ny);
for ix=1:Nx%52 %Nx

    disp(['Doing column ix= ',num2str(ix)])

    for iy=1:Ny%52%Ny

        xz=X(ix,iy);
        yz=Y(ix,iy);

        dx=cos(capteur);
        dy=sin(capteur);
        rhs=zeros(length(dx),1);

        rhs=(exp(-sqrt(-1).*k.*(dx.*xz+dy.*yz)))/norm(exp(-sqrt(-1).*k.*(dx*xz+dy*yz)));
        a=aTK(ix,iy)./(nFsharp+0.05);
        asharp(ix,iy)=a;
        g(:,ix,iy)=(a*Fsharp+a*HEST*eye(size(Fsharp))+Feta'*Feta)^-1*(Feta'*rhs);
    end
end
