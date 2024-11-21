function [X,S,U,mu]=msvd1d(x,L,meancorrect)
%
% usage
%        [X,S,U,mu]=msvd1d(x,L,meancorrect);
% Performs the L level multiresoltion singular value decomposition on
% one dimensional data x
% meancorrect=1 means substract mean, 
% meancorrect=0 means do no subtract mean
% outputs
%        X = transform, same size as x, organised as
%            with smoothest comp at left, with decreasing
%            magnitude detail components as you go to the right
%        S = singular values matrix, dimensions 2 x L, organised 
%            into s.v.'s of increasing size vectors going to the right.
%        U = eigenvector matrix, dimensions 2 x 2L
%        mu = mean matrix dimensions 2xL
%
%                      Written by Ramakrishna Kakarala 17/5/99
%                                 kakarala@labs.agilent.com
% 
% Agilent GIVES NO EXPRESS OR IMPLIED WARRANTY OF ANY KIND AND 
% ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR PURPOSE ARE DISCLAIMED.
% Agilent SHALL NOT BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
% OR CONSEQUENTIAL DAMAGES ARISING OUT OF ANY USE OF THIS SOFTWARE.
% 
% see also imsvd1d.m--this is the inverse transform

% make sure input is a row vector
N=length(x);
x=reshape(x,1,N);

X=zeros(1,N);
S=zeros(2,L);
U=zeros(2,2*L);
mu=zeros(2,L);
if mod(N,2^L)>0
   disp('Length of x must be divisible by 2^L!');
   return;
end;
xt=x;
Lindex=L;
while Lindex>=1
   N=N/2;
   Nb=N;
   
   % arrange odd indices on top row, even indices on bottom
   Xt=zeros(2,N);
   Xt(1,:)=xt(1:2:end);
   Xt(2,:)=xt(2:2:end);
  
   %compute mean and remove it from each row
   mu(:,Lindex)=meancorrect*mean(Xt,2); % mean
   Xt=Xt-mu(:,Lindex)*ones(1,N); 
   T=Xt*Xt';
   [Ut,D,Ut_trans]=svd(T);
   
   % make sure 1st element of each eigenvector is positive
   % this is stated as "without loss of generality" assumption 
   % to force unique result in paper.
   % msvd works without this, just that results will have negative
   % sign.
   if Ut(1,1)<0
    Ut(:,1)=-Ut(:,1);
   end;
   if Ut(1,2)<0
    Ut(:,2)=-Ut(:,2);
   end;
   
   D=sqrt(D);  % get singular values 
              % they are already sorted in decreasing order
   S(:,Lindex)=diag(D);
   U(:,2*(Lindex-1)+1:2*Lindex)=Ut;  % save transform matrix
   Xh=Ut'*Xt;  % compute transform
   X(1:N)=Xh(1,:); %stronger singular value
   X((1:N)+N)=Xh(2,:); % weater singular val
   xt=X(1:N);
   Lindex=Lindex-1;
end;