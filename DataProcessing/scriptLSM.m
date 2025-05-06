clear all

global GUP  SQ HEST Feta  SIGMAN rhs aTK


dossier = '/home/napal/Documents/numerique/3_Close_Disk_lambda0.5_Noise0.05_Avril2018/1/';
load([dossier 'dataMat.mat'],'F','source','capteur');
lambda = 0.5;
nbcapteur = length(capteur);
params.k=2.0*pi/lambda;
params.Nx=200;
params.Ny=200;
params.xmin=-3;
params.xmax=3;
params.ymin=-3;
params.ymax=3;

Nx=200;
Ny=200;
xmin=-3;
xmax=3;
ymin=-3;
ymax=3;


params.selectsource=1:nbcapteur;%1:49;%1:50
params.selectcapteur=1:nbcapteur;%201:249;%51:100;
params.nsource=length(params.selectsource);
params.ncapteur=length(params.selectcapteur);

params.source=source(params.selectsource);
params.capteur=capteur(params.selectcapteur);


params.eta=0.06; %0.001;

Feta=zeros(params.ncapteur,params.nsource);
g=zeros(params.nsource,params.Nx,params.Ny);

Fdata=F(params.selectcapteur,params.selectsource);

Dx=(-params.xmin+params.xmax)/(params.Nx-1);
Dy=(params.ymax-params.ymin)/(params.Ny-1);
[X,Y]=meshgrid(params.xmin:Dx:params.xmax,params.ymin:Dy:params.ymax); 

[Feta,gTK, HEST,aTK]=LSMfctSourceCapteurTK(Fdata, params);

save([dossier 'LSM.mat'],'HEST','Feta','gTK', 'aTK');

params.TKdata='LSM.mat';
[Fsharp, gsharp,HEST,asharp]=LSMfctSourceCapteurSharp(params);


save([dossier 'GLSM_Fsharp.mat'],'gsharp','Fsharp','HEST','asharp');



