function f=dymconv_rec(v,loc,y)
%%
%%  In this function,
%%  given an n-by-n matrix as the operator needed in dynamical sampling problem. Sampling set is round(n/2). 
%%  This function will show the figure the relationship between the error of the signal from the dynamical sampling and the true signal and # of sample levels
%%% check whether length of v and f are same as n???
% generate n-by n convolution matrix
% A=rand(n,n);

% set b(i)=0 if abs(b(i))<nv*sigma
nv=2;

n=length(v);
if size(v,1)==1
    v=v';
end
A=amat(v);
%A=convf(v,n);
% I is n-by-n identity matrix
I=eye(n);
% sampling locations
%loc=[1  3 5 7 9  11];

nloc=length(loc);

loc1=int8(linspace(1,n,n));
loc1(loc)=0;
loc1=find(loc1~=0)

num=length(y)/nloc;
f=zeros(n,1);
%f(loc)=y(1:nloc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To generate an oprator with all the singular values not larger than 1
%[U,D,V]=svd(A);
%%%%%%%%%%%%

%
%% starting from here
tempA=I(loc,:);
B=tempA;
%num1=ceil(num/2)+3;
for i=1:num-1
    tempA=tempA*A;
    B=[B;tempA];
end
%f=B\y;
[U,D,V]=svd(B,'econ');
size(U)
size(D)
size(V)
d=length(diag(D));
figure
d=diag(D);
plot(d)
legend('sigular values')

temp=0.0;
temp1=length(d);
for i=1:length(d)
    if d(i)<temp
        %temp=d(i-1);
        for j=i:length(d)
            d(j)=0;
        end
        temp1=i;
        break;
    end
end

% for i=1:length(d)
%     D(i,i)=d(i);
% end

f0=U'*y;
for i=1:temp1
    f0(i)=f0(i)/D(i,i);
end
for i=temp1+1:length(d)
    f0(i)=0;
end

f=V*f0;
%f(loc)=y(1:nloc);
% B=U*D*V';
% d=diag(D)
% 
% plot(d)
%f=B\y;
%[D(1,1) D(d,d)]
norm(y-B*f)/norm(f);

