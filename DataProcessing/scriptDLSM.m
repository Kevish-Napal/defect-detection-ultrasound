% la pénalisation utilisée initialement diffère legèrement de celle effectuée avec la donnée de F0;
% on utilise comme pénalisation
% asharp_0 *<Fsharp_1 g,g> +  asharp_0 * HEST_0 * |g|^2 
% on a donc utilisée Fsharp_1 au lieu de Fsharp_0 dans le premier terme de la pénalisation.
% le fichier scriptDLSMbis utilise la même pénalisation que celle utilisée pour la première GLSM.

clear all

global GUP  SQ HEST Feta  SIGMAN rhs aTK


Dossier =  '/home/napal/Documents/numerique/3_Close_Disk_lambda0.5_Noise0.05_Avril2018/'


								
								%% A modifier %%
								
								Ref ='0'
								Modif ='1'



%% Dossier de données 
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


%%%
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
%%%
save([dossierModif DLSM],'gsharp','Fsharp');



