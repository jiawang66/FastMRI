function a = matnorm(A, v)
%
% function a = matnorm(A, v)
%
% Compute the norm of matrix A
%
% Input - 
% A: MxN matrix
% v: 0 -- L_infinity norm
%    1 -- L1 norm
%    2 -- L2 norm
%    3 -- F norm
%
% Output -
% a: the norm of matrix A

%% Parameter check
if nargin ~=2 
    error('Error! Two parameter needed.')
end

if ~ismatrix(A)
    error('Error! The first parameter should be a numerical matrix.')
end

%% Compute norm according to parameter v
switch v
    case 0
        a = max(sum(abs(A),2));
    case 1
        a = max(sum(abs(A),1));
    case 2
        [~,d] = eig(A'*A);
        a = sqrt(max(diag(d)));
    case 3
        A = A.*A;
        a = sqrt(sum(A(:)));
    otherwise
        error('The second inputed parameter should be number 0, 1, 2, or 3.')
end

a = abs(a);

end