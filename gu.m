%example 1
X=[1 2 3 4 5 6 7 8];
len1=length(X);
%1st level
%X1=[1 3 5 7; 2 4 6 8]
X1=[X(1:2:len1); X(2:2:len1)];
[M1 N1]=size(X1);
% remove the mean in each row
X1a=X1*(eye(N1)-ones(N1,1)*ones(N1,1)'/N1)
% it also calculated by matlab
X1a=X1- mean(X1')'*ones(1,N1);
T1=X1a*X1a'
[U1,S1,V1]=svd(T1)
X1hat=U1'*X1a  %*sqrt(2)

%next level
F1=X1hat(1,:); %*sqrt(2)
len2=length(F1);
X2=[F1(1:2:len2); F1(2:2:len2)]
[M2 N2]=size(X2);
X2a=X2*(eye(N2)-ones(N2,1)*ones(N2,1)'/N2)
T2=X2a*X2a'
[U2,S2,V2]=svd(T2)
X2hat=U2'*X2a