function [K, R, t]= proj_decomp(P)
%
% function [K, R, t]= proj_decomp(P)
%
% Decompose a projection matrix P (3x4)
% in an intrinsic parameters matrix K
% and a rigid transformation [R t], such that
%   P = K * [R | t]
%
% P : 3x4 : projection matrix
% K : 3x3 : intrinsic parameters
% R : 3x3 : rotation matrix from world to camera coords
% t : 3x1 : translation vector from world to camera coords
%

% Textbook reference: R. Hartley, A. Zisserman, "Multiple View Geometry",
%   Cambridge University Press 2000, page 150.

% 24/3/01, J. Gaspar

P = P ./ norm(P(3,1:3));

E= [0 0 1; 0 1 0; 1 0 0];
[Q,U]= qr(P(:,1:3)'*E);

K= -E*U'*E;

%R= inv(K) * P(:,1:3);
R= -E*Q';

t= inv(K) * P(:,4);

% "The ambiguity in the decomposition is removed by requiring that
%  K have positive diagonal entries", in [Hartley & Zisserman]
%
D= diag( sign(diag(K)) );
%
% note that inv(D)==D
%
K= K*D;
R= D*R;
t= D*t;

if det(R) <0
    [K, R, t]= camera.models.pinhole.utils.proj_decomp(-P);
end

