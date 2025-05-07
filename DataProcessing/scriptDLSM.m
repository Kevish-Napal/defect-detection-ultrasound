% This script computes the DLSM indicator function.
% Note that it uses Fsharp_1 instead of Fsharp_0, where Fsharp_1 is the penalization used for the damaged material:
% asharp_0 *<Fsharp_1 g,g> +  asharp_0 * HEST_0 * |g|^2 


clear all

global GUP  SQ HEST Feta  SIGMAN rhs aTK


Dossier =  '/home/napal/Documents/numerique/3_Close_Disk_lambda0.5_Noise0.05_Avril2018/'


								
%% Choose Ref materal (default 0) and damaged material (Modif) %%

Ref ='0'
Modif ='1'



%% Load data folders
dossierRef = [Dossier Ref];
dossierModif = [Dossier Modif];
														 						 
 %% Dossier de sortie %%
 DLSM = ['/DLSM_Fsharp_geo' Ref '&' Modif '.mat'];
 

load([dossierRef '/GLSM_Fsharp.mat'],'asharp','HEST');
load([dossierModif '/LSM.mat'],'Feta');
load([dossierModif '/GLSM_Fsharp.mat'],'Fsharp');
load([dossierModif '/dataMat.mat'],'capteur');


lambda = 0.5;
nbcapteur = length(capteur);
k=2.0*pi/lambda;
Nx=200;
Ny=200;
xmin=-3;
xmax=3;
ymin=-3;
ymax=3;


Dx=(-xmin+xmax)/(Nx-1);
Dy=(ymax-ymin)/(Ny-1);
[X,Y]=meshgrid(xmin:Dx:xmax,ymin:Dy:ymax); 


options=optimset('Display','off');


for ix=1:Nx%52 %Nx

    disp(['Doing column ix= ',num2str(ix)])

    for iy=1:Ny%52%Ny

        xz=X(ix,iy);
        yz=Y(ix,iy);

        dx=cos(capteur);
        dy=sin(capteur);
        rhs=zeros(length(dx),1);

        rhs=(exp(-sqrt(-1).*k.*(dx.*xz+dy.*yz)))/norm(exp(-sqrt(-1).*k.*(dx*xz+dy*yz)));
        a=asharp(ix,iy);
        gsharp(:,ix,iy)=(a*Fsharp+a*HEST*eye(size(Fsharp))+Feta'*Feta)^-1*(Feta'*rhs);
    end
end


save([dossierModif DLSM],'gsharp','Fsharp');



