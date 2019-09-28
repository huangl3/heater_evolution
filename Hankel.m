function avv = Hankel(A)
%projection onto Hankel matrices
m = size(A);
c = A(:,1);
r = A(m(1),:);
for k = 2:m(1)
    kdiag=[]; %k-th anti-diagonal
    for i = 1:k
        kdiag=[kdiag;A(i, k-i+1)];
    end
    c(k) = mean(kdiag);
end
for k = 2:(m(2)-1)
    kdiag=[];
    for i = 1:(m(1)-k+1)
      kdiag=[kdiag;A(m(1)-i+1,k-1+i)];
    end
%     for i=m(1):-1:max(k+m(1)-m(2),1)
%         kdiag=[kdiag;A(i,k+m(1)-i)];
%     end
        
    r(k) = mean(kdiag);
end
r(1) = c(m(1));
B = hankel(c,r); %built in routine for constructing a hankel matrix from a column and row
avv = B;