function B = measmat(u)
%converts a vector of measurements at a single location into a Hankel
%matrix
%[~,H] = size(u);
H=length(u);
H = ceil(H/2);
H1=length(u)-H+1;
A = zeros(H,H1);   
for k=1:H
  A(k,:) = u(k:k+H1-1);    
end;
B=A;

