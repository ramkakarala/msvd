function x=imsvd1d(X,U,mu)
%
% usage
%        x=imsvd2d(X,U,mu);
% Performs the inverse L level multiresolution svd on transform X,
% using eigenvectors U and means mu.  If
% [X,S,U,mu]=msvd1d(x,L,meancorrect)
% then the inverse is 
%          x=imsvd1d(X,U,mu);
%
%                      Written by Ramakrishna Kakarala 17/5/99
%                                 kakarala@labs.agilent.com
% 
% Agilent GIVES NO EXPRESS OR IMPLIED WARRANTY OF ANY KIND AND 
% ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR PURPOSE ARE DISCLAIMED.
% Agilent SHALL NOT BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
% OR CONSEQUENTIAL DAMAGES ARISING OUT OF ANY USE OF THIS SOFTWARE.
%
% see also msvd1d.m

%force transform into row vector form if needed
N=length(X);
X=reshape(X,1,N);

x=X;
[p,L]=size(mu);  % L is the number of levels

Lindex=1;

% compute inverse transform level by level
while L>=1
   Nl=N/2^L;
   Xh=zeros(2,Nl);
   
   Xh(1,:)=reshape(x(1:Nl),1,Nl);
   Xh(2,:)=reshape(x((1:Nl)+Nl),1,Nl);
   
   Xh=U(:,2*Lindex-1:2*Lindex)*Xh;  % inverse transform
   Xh=Xh+mu(:,Lindex)*ones(1,Nl);   % add the mean back 
 
   x(1:2:2*Nl)=reshape(Xh(1,:),1,Nl);
   x(2:2:2*Nl)=reshape(Xh(2,:),1,Nl);

   Lindex=Lindex+1;
   L=L-1;
end;

