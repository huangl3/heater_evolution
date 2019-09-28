function mat = amat(column)
%Builds a circular matrix from the column
d = size(column,1);
mat = zeros(d);
for j = 1:d
    for k = 1:d
        ind = mod(-k+j+1,d);
        if ind == 0
            ind = d;
        end
        %if abs(column(ind))>0.001
        mat(j,k) = column(ind);
        %end
    end
end
        