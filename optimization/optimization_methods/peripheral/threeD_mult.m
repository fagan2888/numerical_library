% vector / 3D matrix multiplication
function product = threeD_mult(v,M)

    n = size(M,1);
    m = size(M,3);

    sum = zeros(n);
    for i = 1:m
        sum = sum + v(i)*M(:,:,i);
    end
    
    product = sum;

end