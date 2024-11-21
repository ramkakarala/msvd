% batch file to compute statistics for experiments in IEEE Transactions
% paper, "signal approximation using a multiresolution form of the SVD"
% by R. Kakarala and P. Ogunbona
%
%                                 kakarala@labs.agilent.com
% 
% Agilent GIVES NO EXPRESS OR IMPLIED WARRANTY OF ANY KIND AND 
% ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR PURPOSE ARE DISCLAIMED.
% Agilent SHALL NOT BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
% OR CONSEQUENTIAL DAMAGES ARISING OUT OF ANY USE OF THIS SOFTWARE.
clear;
L=4;
p=4;
et=0.5*ones(4,1);
for k=1:L Nl(L-k+1)=512*512/(2^(2*k)); end;
Nl=Nl;
load barb;
disp('barb');
[X,Cbarb,mu]=mpc2wcov(barb,L,1);
for k=1:L
  [Vbarb(:,(k*4-3):4*k),Sbarb(:,(k*4-3):4*k),U]=svd(Cbarb(:,(k*4-3):4*k));
  Ebarb(:,k)=diag(Sbarb(:,(k*4-3):4*k));
end;
for k=1:L
  cbarb(k)=p*log(cg(Ebarb(:,k))); 
  Vlbarb(:,k)=Vbarb(:,4*k-3);
  C=Cbarb(:,(k*4-3):4*k);
  hbarb(k)=Nl(k)*(Ebarb(1,k)*et'*inv(C)*et+1/Ebarb(1,k)*et'*C*et-2);
  V=Vbarb(:,(L*4-3):4*L);
  F=V'*C*V;
  cpcbarb(k)=Nl(k)*log(det(diag(diag(F)))/det(F));
end;
cbarb=cbarb.*(Nl-17/6);;
clear barb;
load boats;
disp('boats');
[X,Cboat,mu]=mpc2wcov(boats,L,1);
for k=1:L
  [Vboat(:,(k*4-3):4*k),Sboat(:,(k*4-3):4*k),U]=svd(Cboat(:,(k*4-3):4*k));
  Eboat(:,k)=diag(Sboat(:,(k*4-3):4*k));
end;
for k=1:L
  cboat(k)=p*log(cg(Eboat(:,k))); 
  Vlboat(:,k)=Vboat(:,4*k-3);
    C=Cboat(:,(k*4-3):4*k);
  hboat(k)=Nl(k)*(Eboat(1,k)*et'*inv(C)*et+1/Eboat(1,k)*et'*C*et-2);
    V=Vboat(:,(L*4-3):4*L);
  F=V'*C*V;
  cpcboat(k)=Nl(k)*log(det(diag(diag(F)))/det(F));
end;
cboat=cboat.*(Nl-17/6);;
clear boats;
load grass;
disp('grass');
[X,Cgrass,mu]=mpc2wcov(grass,L,1);
for k=1:L
  [Vgrass(:,(k*4-3):4*k),Sgrass(:,(k*4-3):4*k),U]=svd(Cgrass(:,(k*4-3):4*k));
  Egrass(:,k)=diag(Sgrass(:,(k*4-3):4*k));
end;
for k=1:L
  cgrass(k)=p*log(cg(Egrass(:,k))); 
  Vlgrass(:,k)=Vgrass(:,4*k-3);
    C=Cgrass(:,(k*4-3):4*k);
  hgrass(k)=Nl(k)*(Egrass(1,k)*et'*inv(C)*et+1/Egrass(1,k)*et'*C*et-2);
    V=Vgrass(:,(L*4-3):4*L);
  F=V'*C*V;
  cpcgrass(k)=Nl(k)*log(det(diag(diag(F)))/det(F));
end;
cgrass=cgrass.*(Nl-17/6);;
clear grass;
load noisy_checker;
disp('checker');
[Xc,Cchecker,mu]=mpc2wcov(checker,L,1);
for k=1:L
  [Vchecker(:,(k*4-3):4*k),Schecker(:,(k*4-3):4*k),U]=svd(Cchecker(:,(k*4-3):4*k));
  Echecker(:,k)=diag(Schecker(:,(k*4-3):4*k));
end;
for k=1:L
  cchecker(k)=p*log(cg(Echecker(:,k))); 
  Vlchecker(:,k)=Vchecker(:,4*k-3);
  C=Cchecker(:,(k*4-3):4*k);
  hchecker(k)=Nl(k)*(Echecker(1,k)*et'*inv(C)*et+1/Echecker(1,k)*et'*C*et-2);
  Vo=Vchecker(:,(L*4-3):4*L);
  F=V'*C*V;
  cpcchecker(k)=Nl(k)*log(det(diag(diag(F)))/det(F));
end;
cchecker=cchecker.*(Nl-17/6);;
clear checker,Xc;
for k=1:L-1
  Vl1=Vboat(:,4*k+1:4*k+4);
  C=Cboat(:,(k*4-3):4*k);
  F=Vl1'*C*Vl1;
  rpcboat(k)=Nl(k)*log(det(diag(diag(F)))/det(F));
  Vl1=Vbarb(:,4*k+1:4*k+4);
  C=Cbarb(:,(k*4-3):4*k);
  F=Vl1'*C*Vl1;
  rpcbarb(k)=Nl(k)*log(det(diag(diag(F)))/det(F));
  Vl1=Vgrass(:,4*k+1:4*k+4);
  C=Cgrass(:,(k*4-3):4*k);
  F=Vl1'*C*Vl1;
  rpcgrass(k)=Nl(k)*log(det(diag(diag(F)))/det(F));
  Vl1=Vchecker(:,4*k+1:4*k+4);
  C=Cchecker(:,(k*4-3):4*k);
  F=Vl1'*C*Vl1;
  rpcchecker(k)=Nl(k)*log(det(diag(diag(F)))/det(F));
end;
  
save fourres ;




