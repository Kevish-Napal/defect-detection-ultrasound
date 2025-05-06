function y=dfun(x)
global GUP SQ HEST %GUP=conj(sc).*sc       SQ = SIGMAN.^2     HEST = delta 
h=HEST;
xs=x^2;
hs=h^2;
nn=size(SQ);
SQA=(SQ+x).^2;
TMP=GUP.*(xs-(hs*SQ));
y =  sum(TMP./SQA);


% fonction utilisé dans le morozov: le zero de cette fonction donne le
% delta optimal (le x de la fonction)
