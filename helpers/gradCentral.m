function [Ix,Iy] = gradCentral(I)
    [rows,cols] = size(I);
    Ix = zeros(size(I));
    Iy = zeros(size(Ix));
    Ix(2:rows-1,:) = (I(3:rows,:) - I(1:rows-2,:)) / 2;
    Iy(:,2:cols-1) = (I(:,3:cols) - I(:,1:cols-2)) / 2;
    Ix(1,:) = I(2,:) - I(1,:);
    Ix(rows,:) = I(rows,:) - I(rows-1,:);
    Iy(:,1) = I(:,2) - I(:,1);
    Iy(:,cols) = I(:,cols) - I(:,cols-1);
end