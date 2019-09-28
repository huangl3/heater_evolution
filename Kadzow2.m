function avv = Kadzow2(A,K,k)
%Cadzow denoising of the Hankel matrix A
%K is the presumed rank of the denoised matrix
%eps is the tolerance parameter
% [U,S,V] = svd(A);
% s1=flip(sort(diag(S)));
% eps=s1(k);
% for i=1:size(S,1)
%     if S(i,i)<eps
%         S(i,i)=0;
%     end
% end
% A=U*S*V';
% A=Hankel(A);
flag = 0;

% for one hotspot falg<10 with 1:m:d
% for two hotspots falg<7 with 2:m:d
while flag<12
 [U,S,V] = svd(A);
 for i=k+1:size(S,1)
     S(i,i)=0;
 end
%  s1=flip(sort(diag(S)));
%  eps=s1(k);
%  for i=1:size(S,1)
%     if S(i,i)<eps
%         S(i,i)=0;
%     end
% end
A=U*S*V';
A=Hankel(A);
flag=flag+1;
end

% err = S(K+1,K+1)/S(K,K);
% flag = 0;
% while (err > eps) && (flag < 100) 
%     dd=size(S);
%     d=min(dd);
%     for j = (K+1):d
%         S(j,j)=0;
%     end
%     A = U*S*V';
%     A = Hankel(A); %projection onto Hankel matrices
%     [U,S,V] = svd(A);
%     err = S(K+1,K+1)/S(K,K);
%     flag = flag+1;
% end
avv = A;