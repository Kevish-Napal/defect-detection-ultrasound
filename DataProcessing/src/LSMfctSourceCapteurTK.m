%function [Feta, g, HEST,ignormLSM,ignormGLSM]=LSMfct(index,sortie, params,h)
%function [Feta, g, HEST,ignormLSM,ignormGLSM,rhsSave]=LSMfct(index, params)
%[Feta, g1,g2,g,T1,T2,Morozov1, Morozov2,t, HEST]
function  [Feta, g,HEST,aa]=LSMfctSourceCapteurTK(index, params) %, F_sharp,aSharp,aLsm

global GUP  SQ HEST Feta  SIGMAN rhs 
nsource=params.nsource;
ncapteur=params.ncapteur;               % nombre de source et de capteurs

source=params.source;
capteur=params.capteur;                 % vecteurs des positions des sources et capteurs en radiant
%load(index)
F=index;
k=params.k;
Nx=params.Nx;                           % Nx et Ny sont la résolution pour le maillage
xmin=params.xmin;                       % Xmin, Xmax, Ymin, Ymax sont les sommet du carré contenant l'objet et qui sera découpé selon un maillage de petits carré.
xmax=params.xmax;

Ny=params.Ny;
ymin=params.ymin;
ymax=params.ymax;                                    

Dx=(-xmin+xmax)/(Nx-1);       
Dy=(ymax-ymin)/(Ny-1);                  % pas du maillage en abscisse et en ordonnée.
[X,Y]=meshgrid(xmin:Dx:xmax,ymin:Dy:ymax);

eta=params.eta; %0.001;
noise=1+eta*(rand(size(F))*2-1)+sqrt(-1)*eta*(rand(size(F))*2-1);
Feta=F.*noise;

HEST=norm(Feta-F)%; % le delta de Morozov

if HEST < 1.e-02
    HEST= 1.e-02;                       % Pour que la perturbation soit significative.
end

nF=norm(Feta);
 disp(['Bruit en pourcentage = ',num2str(HEST/nF*100)]) % simplement une commande d'affichage, sans num2str(le nombre), le nombre ne s'affiche pas.
[Ueta,Seta,Veta]=svd(Feta);             % Ueta et Veta sont des matrices de rotations, Seta une matrice diagonale contenant les valeurs singulières de Feta.
SIGMAN=diag(Seta);                      % SIGMAN est un vecteur contenant la diagonale de Seta.
SQ=SIGMAN.^2;                           % Vecteur contenant les valeurs propres de (Feta)*Feta.
Us=Ueta';

options=optimset('Display','off');                                           % Qu'est-ce que c'est?


for ix=1:Nx%52 %Nx

    disp(['Doing column ix= ',num2str(ix)])

    for iy=1:Ny%52%Ny

        xz=X(ix,iy);
        yz=Y(ix,iy);                   % coordonnées des différents points du maillage.

        dx=cos(capteur);
        dy=sin(capteur);               % coordonnées cartésiennes des capteurs.
        rhs=zeros(length(dx),1);

        rhs=(exp(-sqrt(-1).*k.*(dx.*xz+dy.*yz)))/norm(exp(-sqrt(-1).*k.*(dx*xz+dy*yz)));
        %rhs=transpose(rhs);
        %             if isempty(h)
        sc=Us*rhs;
        %             else
        %            sc=Us*rhs-Veta*h(:,ix,iy);
        %             end
        GUP=conj(sc).*sc;
        SIGMAN=diag(Seta);
        SQ=SIGMAN.^2;
        gammamin=HEST*min(SIGMAN);    % borne inf pour la recherche de alpha* de morozov.
        gammamax=HEST*max(SIGMAN);    % le alpha(delta) de morozov se trouve dans l'intervalle [gammamin,gammamax].

        [a,fval,exitflag]=fzero('dfun',[gammamin,gammamax],options); %Faire par dichotomie pour voir (nombre de dichotomie?)
% a =morozov,     fval = dfun(a),     exitflag < 0 => aucuns zero 
      
        if a<0||exitflag<0
            disp(['Problem: a= ',num2str(a),' Resetting a=0'])
            disp(['Problem: exitflag= ' num2str(exitflag)])
            a=0;
        end

        aa(ix,iy)=a;

        Valxi=(SIGMAN.*sc)./(SQ+a);
        
        g(:,ix,iy)=Veta*Valxi;
    end
end
