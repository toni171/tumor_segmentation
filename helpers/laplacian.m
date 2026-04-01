function L = laplacian(phi)
    kernel = [1 4 1;4 -20 4;1 4 1];
    L = conv2(phi, kernel, 'same');
end