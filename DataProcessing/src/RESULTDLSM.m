

clear all

	%3DisksConfig_lambda0.5_Noise0.05_Mars2018
    %3_Close_Disk_lambda0.5_Noise0.05_Avril2018
    
Dossier =  '/home/napal/Documents/numerique/3DisksConfig_lambda0.5_Noise0.05_Mars2018/'


								
								%% A modifier %%
								
								Ref ='0'
								Modif ='1'



%% Dossier de données 
dossierRef = [Dossier Ref];
dossierModif = [Dossier Modif];
											

	%%%%%%%%%%% INITIAL BACKGROUND DATA %%%%%%%%%%%%%%%%%%%%%%
	
    InitialBackground = load([dossierRef '/GLSM_Fsharp.mat']);
    g = InitialBackground.('gsharp');
    F = InitialBackground.('Fsharp');
    a = InitialBackground.('asharp');
    HEST = InitialBackground.('HEST');
    
    
    %%%%%%%%%%%%%%%%%%%% MODIFIED BACKGROUND %%%%%%%%%%%%%%%%%%%%%%%%%
    ModifiedBackground = load([dossierModif '/DLSM_Fsharp_geo' Ref '&' Modif '.mat']);
    g1 = ModifiedBackground.('gsharp');
    F1 = ModifiedBackground.('Fsharp');
   

    
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
    	A1(ii,jj)=abs(g1(:,ii,jj)'*(F1*g1(:,ii,jj)))+HEST * norm(g1(:,ii,jj))^2;
        A(ii,jj)=abs(g(:,ii,jj)'*(F*g(:,ii,jj))) + HEST * norm(g(:,ii,jj))^2;
        D(ii,jj) = abs((g(:,ii,jj)-g1(:,ii,jj))'*(F*(g(:,ii,jj)-g1(:,ii,jj)))) + HEST * norm(g(:,ii,jj)-g1(:,ii,jj))^2; % Ou B à la place de Feta
      end
end
        
	%/////////////////////////////////////////////////////////////%
	%/////////////////////////////////////////////////////////////%


	%warningmessage = warndlg(['saving images on ' dossier1]);
	%drawnow     
	%waitfor(warningmessage);

     ImageT = 1./sqrt(A1.*(1+A1*1./D));
     ImageMB = 1./sqrt(A+A1.*(1+A./D));
     ImageT_2 = 1./sqrt(A.*(1+A./D));
     
close all
%figure 1 - ImageT
for c = 2:0.1:3
figure, imagesc(xmin:Dx:xmax,ymin:Dy:ymax,1./sqrt(A1.*(1+c*A1*1./D))) % ImageT
set(gca,'Ydir','Normal')
title('$$\frac{1}{\sqrt{A*(1+\frac{A}{D})}}, avec  ~ g_{GLSM}$$','interpreter','latex','FontSize',20)
colorbar
end
%saveas(gcf,[dossier1 label 'ImageT.jpg']);

close all
%figure 2 - 
for c = 0:0.1:2
figure, imagesc(xmin:Dx:xmax,ymin:Dy:ymax,1./(A+c*A1.*(1+A./D))) % ImageMB
set(gca,'Ydir','Normal')
title('$$\frac{1}{\sqrt{A0+A*(1+\frac{A0}{D})}}, avec ~  g_{GLSM}$$','interpreter','latex','FontSize',20)
colorbar
end
%saveas(gcf,[dossier1 label 'ImageT_2N.jpg']);


close all
figure, imagesc(xmin:Dx:xmax,ymin:Dy:ymax,1./(A1.*(1+A1./D))) % ImageT2
set(gca,'Ydir','Normal')
title('$$\frac{1}{\sqrt{A_1*(1+\frac{A_1}{D_0})}}$$','interpreter','latex','FontSize',20)
colormap jet 
colorbar
figure, imagesc(xmin:Dx:xmax,ymin:Dy:ymax,1./(A.*(1+A./D))) % ImageT2
set(gca,'Ydir','Normal')
title('$$\frac{1}{\sqrt{A_0*(1+\frac{A_0}{D_0})}}$$','interpreter','latex','FontSize',20)
colormap jet 
colorbar



close all
X = 0.75*cos([0:2*pi/100:2*pi]);
Y = 0.75*sin([0:2*pi/100:2*pi]);
c=0.5;
figure, imagesc(xmin:Dx:xmax,ymin:Dy:ymax,1./(A.*(1+c*A./D))) % ImageT2
set(gca,'Ydir','Normal')
hold on,plot(X-2,Y+2,'black--')
hold on,plot(X+2,Y+2,'black--')
hold on,plot(X,Y-2,'black--')
colormap jet 
colorbar

hold on, plot([1.5,2.5],[2,2],'*black','markers',12) %long crack droit

hold on, plot([-1.5,-2.5],[2,2],'*black','markers',12) %long crack gauche

hold on, plot([-.5,.5],[-2,-2],'*black','markers',12) %long crack bas

hold on, plot([-.5,.5],[-2.3,-2.3],'*black','markers',12) %long crack plus bas

hold on, plot([-0.5,-0.17],[-2,-2],'*black','markers',12) % mini crack bas

hold on, plot([-0.5,0.17],[-2,-2],'*black','markers',12) % medium crack bas


figure, imagesc(xmin:Dx:xmax,ymin:Dy:ymax,1./A1) % ImageT2





