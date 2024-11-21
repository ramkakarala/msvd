function x=imsvd2d(X,V,mu)
%
% usage
%        x=imsvd2d(X,V,mu)
% Performs the inverse L level multiresoltion analysis using principal 
% components on two-dimensional data x; inverse of twodmaupc.
% SEE ALSO MSVD2D
%
%
%                      Written by Ramakrishna Kakarala 17/5/99
%                                 kakarala@labs.agilent.com
% 
% Agilent GIVES NO EXPRESS OR IMPLIED WARRANTY OF ANY KIND AND 
% ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR PURPOSE ARE DISCLAIMED.
% Agilent SHALL NOT BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
% OR CONSEQUENTIAL DAMAGES ARISING OUT OF ANY USE OF THIS SOFTWARE.

x=X;
[M,N]=size(x);
[p,L]=size(mu);

Lindex=1;
totalsize=M*N;
while L>=1
   LsizeM=M/2^L; LsizeN=N/2^L;
   Nb=LsizeM*LsizeN;
   Xh=zeros(4,Nb);
   Xh(1,:)=reshape(x(1:LsizeM,1:LsizeN),1,Nb);
   Xh(2,:)=reshape(x(1:LsizeM,(1:LsizeN)+LsizeN),1,Nb);
   Xh(3,:)=reshape(x((1:LsizeM)+LsizeM,(1:LsizeN)+LsizeN),1,Nb);
   Xh(4,:)=reshape(x((1:LsizeM)+LsizeM,1:LsizeN),1,Nb);
   Xh=V(:,4*Lindex-3:4*Lindex)*Xh; % inverse transform
   Xh=Xh+mu(:,Lindex)*ones(1,Nb); 
   x(1:2:2*LsizeM,1:2:2*LsizeN)=reshape(Xh(1,:),LsizeM,LsizeN);
   x(1:2:2*LsizeM,2:2:2*LsizeN)=reshape(Xh(2,:),LsizeM,LsizeN);
   x(2:2:2*LsizeM,1:2:2*LsizeN)=reshape(Xh(3,:),LsizeM,LsizeN);
   x(2:2:2*LsizeM,2:2:2*LsizeN)=reshape(Xh(4,:),LsizeM,LsizeN);
   % Xt is now 4 x Nb
   Lindex=Lindex+1;
   L=L-1;
end;

