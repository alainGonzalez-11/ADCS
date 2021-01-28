function [x] = CrossMatrix(vector)
%CROSSMATRIX Returns the Cross Product Matrix
%   Creates the matrix [Ax] defined by:
%   a x b = [Ax]b
%   [Ax] = |  0   -a3   a2  |
%          |  a3   0   -a1  |
%          | -a2   a1   0   |
%   This matrix is also referred as the tilde matrix.
x = [0 -vector(3) vector(2);
    vector(3) 0 -vector(1);
    -vector(2) vector(1) 0];
end

