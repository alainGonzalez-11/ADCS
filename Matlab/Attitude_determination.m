function [q] = Attitude_determination(a_i, v_i, s_i)
%Attitude_determination Determination of CRPs applying QUEST method.
%   Arguments:
%   - a_i: Relative weights vector of size (n)
%   - v_i: Reference unit vectors matrix of size (3, n)
%   - s_i: Observed unit vectors matrix of size (3, n)
%   Returns:
%   - q: CRPs vector of size (3)
%   Notes:
%   - "n" is the number of observed vectors
    B = (a_i .* s_i) * v_i';
    S = B + B';
    sigma = trace(B);
    Z = [B(2,3) - B(3,2); B(3,1) - B(1,3); B(1,2) - B(2,1)];
    K = [sigma Z';
        Z, (S-sigma*eye(3))];
    A = sum(a_i);
    error = 1;
    while error > 1e-10
        A_old = A;
        r = det(K - A * eye(4))
        r1 = det(K - A * eye(4)) * trace(inv(K - A * eye(4)) * (-eye(4)));
        A = A - r/r1;
        error = abs(A_old-A);
    end
    q = (inv((A + sigma) * eye(3) - S))* Z;
end
