function u = myTikho_inv (A, b, lambda, option)
% function u = myTikhonov (A, b, lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tikhonov regularization to minimize:
% arg min(u) = ||Au-b||_2^2 +  ||u||_2^2
    % (AtA + lambda*I)u = At.b
    % u = (AtA + lambda*I)^-1 At.b
    % P u = r;
    % u = cgs(P, r)
% INPUTS:
%       - A     : Linear System Matrix A
%       - b     : Observed data vector b
%       - lambda: Regularization pararameter
%       - option: (1)cgs, (2)inv, (3)QR solver
% OUTPUS:
%       - u     : u solved by Tikhonov
% Developed by: Edmundo Arom Miranda
% Version 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    At = A';
    AtA = At*A;
    [m, n] = size(AtA);
    switch option
        case 1 % option 1          
            u = cgs(AtA + lambda*speye(m,n), At*b);
        case 2 % option 2 
            u = inv(AtA + lambda*speye(m,n))*At*b;
        case 3 % option 3 
            u = (AtA + lambda*speye(m,n))\(At*b);

        case 11 
            u = tikho_personalized(A, b, lambda );

    end
    
return
function x = tikho_personalized (A, b, lambda)
Param.k = 1;
Param.N_window = 20; 
Param.beta = 1/100000;
Param.tolerance = 1/1000;
Param.operator = 'G';
Param.alpha=2;

switch (Param.operator)
    case 'I'
        L=speye(size(A,2)); 
    case 'G'
        L=speye(size(A,2))-circshift(speye(size(A,2)),[0 1]); 
        L(end,:)=0; %Tikhonov matrix L
end
At = A';
AtA = A'*A;
B =b';
x0 = zeros(size(A,2),1);
err = 1;
while err > Param.tolerance
    W = sparse(size(A,2));
    Lx = L*x0; 
    
%     for i = 1:numel(x0)
%         W(i,i) = Param.k/2*(abs(Lx(i))^2+Param.beta)^((Param.k/2)-1);
%     end

    W_diagonal = Param.k/2 * (abs(Lx).^2 + Param.beta).^((Param.k/2) - 1);
    W = spdiags(W_diagonal, 0, size(A,2), size(A,2));

    x = (AtA+Param.alpha^2*L'*W*L)\(At*b); %del Paper A Regularized Inverse Approach to Ultrasonic
    %Pulse-Echo Imaging página 3
    err = norm(x-x0)^2/norm(x)^2; %error relativo
    x0 = x;
end
return
