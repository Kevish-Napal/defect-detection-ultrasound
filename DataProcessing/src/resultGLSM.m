clear all

Dossier =  '/home/napal/Documents/numerique/Dtntest';%3_Close_Disk_lambda0.5_Noise0.05_Avril2018/'
								
								%% A modifier %%
								
								Ref ='0'
				
%% Dossier de données 
DossierRef = [Dossier Ref];

	%%%%%%%%%%% INITIAL BACKGROUND DATA %%%%%%%%%%%%%%%%%%%%%%
	
    InitialBackground = load([Dossier '/GLSM_Fsharp.mat']);
    g = InitialBackground.('gsharp');
    F = InitialBackground.('Fsharp');
    a = InitialBackground.('asharp');
    HEST = InitialBackground.('HEST');
    
Nx=200;
Ny=200;
xmin=-3;
xmax=3;
ymin=-3;
ymax=3;
Dx=(-xmin+xmax)/(Nx-1);
Dy=(ymax-ymin)/(Ny-1);
 

parfor ii = 1:Nx
    for jj =1:Ny
    	A(ii,jj)=abs(g(:,ii,jj)'*(F*g(:,ii,jj))) + HEST * norm(g(:,ii,jj))^2;
    end
end
        
	%/////////////////////////////////////////////////////////////%
	%/////////////////////////////////////////////////////////////%

close all
%figure 1 - ImageT
figure, imagesc(xmin:Dx:xmax,ymin:Dy:ymax,1./A) % ImageT
set(gca,'Ydir','Normal')
colorbar
colormap jet 
