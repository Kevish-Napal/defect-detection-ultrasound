clear all;

chemin = '/home/napal/Documents/numerique/3_Close_Disk_lambda0.5_Noise0.05_Avril2018/1/';

Source = importdata([chemin 'thetaSource.dat'],'\n');
source = Source(2:101);

Capteur = importdata([chemin 'thetaCapteur.dat'],'\n');
capteur = Capteur(2:101);

aRe=dlmread([chemin 'FF_Real.dat']);
ARe = aRe(2:101,:);
aIm=dlmread([chemin 'FF_Image.dat']);
AIm = aIm(2:101,:);

F=ARe + sqrt(-1)*AIm;

save([chemin 'dataMat.mat'],'F','source','capteur');
