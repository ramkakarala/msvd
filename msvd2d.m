function [X,S,V,mu]=msvd2d(x,L,meancorrect)
%
% usage
%        [X,S,V,mu]=msvd2d(x,L,meancorrect);
% Performs the L level multiresoltion analysis using principal 
% components on two-dimensional data x
% meancorrect=1 means substract mean, 
% meancorrect=0 means do no subtract mean
% outputs
%        X = transform, same size as x, organised as
%            with smoothest comp in upper left, with decreaseing
%            magnitude principal components in clockwise order
%        S = singular values matrix, dimensions 4 x L, organised 
%            into s.v.'s of increasing size images going to the right.
%        V = eigenvector matrix, dimensions 4 x 4L
%        mu = mean matrix dimensions 4xL
%
%                      Written by Ramakrishna Kakarala 17/5/99
%                                 kakarala@labs.agilent.com
% 
% Agilent GIVES NO EXPRESS OR IMPLIED WARRANTY OF ANY KIND AND 
% ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR PURPOSE ARE DISCLAIMED.
% Agilent SHALL NOT BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
% OR CONSEQUENTIAL DAMAGES ARISING OUT OF ANY USE OF THIS SOFTWARE.

[M,N]=size(x);
X=zeros(M,N);
S=zeros(4,L);
V=zeros(4,4*L);
mu=zeros(4,L);
if (mod(N,2^L)>0) | (mod(M,2^L)>0)
   disp('Dimensions of x must be divisible by 2^L!');
   return;
end;
xt=x;
Lindex=L;
while Lindex>=1
   M=M/2; N=N/2;
   Nb=M*N;
   Xt=zeros(4,Nb);
   for q=1:N
      A=xt(:,2*q-1:2*q)'; % two columns, transposed into two rows
                          % dimensions are 2 x 2*M = M blocks
      Xt(:,M*(q-1)+1:M*q) = reshape(A,4,M);
   end;                    
   % Xt is now 4 x Nb
  mu(:,Lindex)=meancorrect*mean(Xt,2); % mean
  Xt=Xt-mu(:,Lindex)*ones(1,Nb); 
  T=Xt*Xt';
  [Vl,D,Vlt]=svd(T);
  D=sqrt(D);  % get singular values 
              % they are already sorted in decreasing order
  S(:,Lindex)=diag(D);
  V(:,4*(Lindex-1)+1:4*Lindex)=Vl;
  Xh=Vl'*Xt;  % compute transform
  X(1:M,1:N)=reshape(Xh(1,:),M,N); % strongest svalue
  X(1:M,(1:N)+N)=reshape(Xh(2,:),M,N); % next strongest sval
  X((1:M)+M,(1:N)+N)=reshape(Xh(3,:),M,N); % third strongest sval
  X((1:M)+M,1:N)=reshape(Xh(4,:),M,N); % weakest sval
  xt=X(1:M,1:N);
  Lindex=Lindex-1;
end;

