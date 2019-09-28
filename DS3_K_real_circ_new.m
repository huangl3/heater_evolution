function  DS3_K_real_circ_new(m,h)
%spectral recovery with Kadzow denoising from the circular rim
%m - subsampling factor
%L - oversampling parameter (signal is sampled at 2(m+L) time levels)
%h - time step

 %[v, ~] = stempread('wheel_cooldown_one_hotspot_15_sensors.csv');
 %[baseline,~] = stempread('wheel_cooldown_one_hotspot_15_sensors.csv');
[v, ~] = stempread('wheel_cooldown_two_hotspots_15_sensors.csv');
[baseline,~] = stempread('wheel_cooldown_two_hotspots_15_sensors.csv');

%% delete the first 19 rows data
v=v(20:end,:);

rv=size(v,1);% the # of rows in v
avv=10;% Take every avv rows of v and average the rows
r_avv=rem(rv,avv);% the left rows
n_avv=(rv-r_avv)/avv; % the # of rows after we do the average
w=[];
% The following are the average process
for i=1:n_avv
    w=[w;mean(v((i-1)*avv+1:i*avv,:))];
end
w=[w;mean(v(n_avv*avv+1:size(v,1),:),1)];
v=w;

size(v)
%%%%%%
figure('DefaultAxesFontSize',16)
%%%%%%%% plot the evolved data at different locations
for i=1:length(v(1,:))
    plot((v(:,i))')
    hold on
end
locations=(1:15)';
legend(strcat('location=',num2str(locations)))
hold off
%%%%%%%%%%%%%%%%%
figure('DefaultAxesFontSize',16)
%%%%%% plot the evolved data
for i=1:length(v(:,1))
    plot(v(i,:))
    hold on
end
hold off

% To get the data in fourier domain
y = fft(v');
%
[~,ry]=size(y);
y1y2=[];

for j=1:size(y,1)
     A = measmat(y(j,:));
     %A=hankel(y(j,1:n-m),y(j,n-m:n));
     %size(A)
     C = Kadzow2(A,1,5); %for one-hotspot
     b = -C(:,2);
     r = [C(:,1)\b;1]; 
%         r= fliplr(r');
    r= fliplr(transpose(r));
    sr=roots(r);
    y1y2=[y1y2 sr];
end


% for i=1:ry-1
%     y1y2=[y1y2;abs(y(i+1,:)./y(i,:))];
% end
% y1y2(1,:)=[];
% y1y2(end,:)=[];
cvf=y1y2
% draw the graph of the estimated eigenvalues of the filter

figure('DefaultAxesFontSize',16)
% plot(sort((abs(cvf))'),'r--') 
% hold on

d = size(v,2);
d = floor(d/m)*m;
%if d/2 == floor(d/2)
%    d = d-m;
%end
J = d/m;
% L = size(v,1);
% L = floor(L/2)-m;
z = v';
%u = meastaken5(xv, fil, m, L);
%u = z(3:m:d,1:2*(m+L));
%u = z(1:m:d,:);
u = z(2:m:d,:);
y = (fft(u));

%r=zeros(d,1);
L=-1;
D=[];

n=length(y(1,:));

% Cadzow Denoising
for j=1:J
     A = measmat(y(j,:));
     %A=hankel(y(j,1:n-m),y(j,n-m:n));
     %size(A)
     if j==1
        if L < 0
         %C = Kadzow2(A,m,(1+m)/2); %for one-hotspot locations 1:m:d with
         %10 cadzow iteration
         %C = Kadzow2(A,m,(1+m)/2+3); %for one-hotspot locations 2:m:d with 1 cadzow iteration
         %C = Kadzow2(A,m,(1+m)/2+3); %for one-hotspot locations 3:m:d with 3 cadzow iteration 
         C = Kadzow2(A,m,(1+m)/2+3); %for two-hotspots locations 2:m:d
         %C = Kadzow2(A,m,(1+m)/2+2); %for two-hotspots locations 3:m:d
           
        else
         C = A;
        end;
%    C=A;
    %C=C(1:m,:);
    
        b = -C(:,(m+3)/2);
        r = [C(:,1:(m+1)/2)\b;1]; 
%         r= fliplr(r');
        r= fliplr(transpose(r));
        sr=roots(r);
        sr=sort(sr);
        D=[D sr' (sr(1:end-1))'];
     else
        if L < 0
         %C = Kadzow2(A,m,m+3); %for one-hotspot locations 1:m:d
         %C = Kadzow2(A,m,m+6); %for one-hotspot locations 2:m:d
         %C = Kadzow2(A,m,m+6); %for one-hotspot locations 3:m:d
         %C = Kadzow2(A,m,m+6); %for two-hotspots when the uniform sample start from 1
         C = Kadzow2(A,m,m+4); %for two-hotspots when the uniform sample start from 2  m+4
         %C = Kadzow2(A,m,m+3); %for two-hotspots when the uniform sample start from 3
        else
         C = A;
        end;
        b = -C(:,m+1);
        r = [C(:,1:m)\b;1]; 
        %r= fliplr(r');
        r= fliplr(transpose(r));
        D=[D roots(r)'];
    end
end

D=sort(abs(D));
index=find(D>1);
if ~isempty(index)
    D(index)=1-(1e-6);
end

%plot(D,'b-.')
scatter((1:length(D))',D,'b')
title('Spectrum','FontSize',16)
%legend('eigenvalues by averaging','eigenvalues by down sampled data')
%hold off

% 
% Enforce the eigenvalues to be pairs.
len=length(D);
D1=[];

for i=1:(len-1)/2
    D1=[D1 (D(2*i-1)+D(2*i))/2];
end
D1=[D1 D(len)];
D=[flip(D1) D1(1:end-1)];

y=D;
F=dftmtx(d);
A=1/d*F*diag(D)*F';
A=real(A);
%A(:,1)
omega=2:m:d;
%omega=1:m:d;
if m~=1
    omega=[omega 3 10];
    %omega=[omega 3 15];
end
omega
inif=v(1,:);
v=(v(:,omega))';
v=reshape(v,[size(v,1)*size(v,2) 1]);
f1=dymconv_recN(A(:,1),omega,v);
f1=real(f1');
x=1:1:length(f1);
%f1=f1'+baseline(1,1:2:end);
figure('DefaultAxesFontSize',16)
scatter(x,f1,'r*');
%plot(x,f1,'r-');
hold on

f1(1,omega);
(v(1:length(omega)))';
 scatter(x,inif,'b');
% plot(x,inif,'b--')
 legend('recovered signal','original signal')
 hold off
norm(inif(omega)-f1(omega))/norm(inif)
norm(inif-f1)/norm(inif)


