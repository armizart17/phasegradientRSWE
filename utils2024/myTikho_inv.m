function u = myTikho_inv (A, b, lambda, option, pars)
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
% Developed by: Edmundo Arom Miranda, based on R-WAVE code
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

        case 6 % option Generalized
            u = tikho_personalized(A, b, pars);
        case 11 % new CGPT
            u = tikho_gpt (A, b, pars);

    end
    
end
function x = tikho_personalized (A, b, Param)

    Param.k = 1;
    Param.beta = 1/100000; % 1e5;
    Param.beta = 1; % 1e5;
    Param.tolerance = 1e-5;
    Param.operator = 'G';
    Param.alpha = 2;
    
    switch (Param.operator)
        case 'I'
            L = speye(size(A,2)); 
        case 'G'
            L = speye(size(A,2))-circshift(speye(size(A,2)),[0 1]); 
            L(end,:) = 0; %Tikhonov matrix L
    end
    At = A';
    AtA = A'*A;
    B = b';
    x0 = zeros(size(A,2),1);
    err = 1;
    while err > Param.tolerance
    
        Lx = L*x0; 
        
        %%%%%%%%%%% ORIGINAL SLOW %%%%%%%%%%%
    %     W = sparse(size(A,2));
    %     for i = 1:numel(x0)
    %         W(i,i) = Param.k/2*(abs(Lx(i))^2+Param.beta)^((Param.k/2)-1);
    %     end
    
        W_diagonal = Param.k/2 * (abs(Lx).^2 + Param.beta).^((Param.k/2) - 1);
        W = spdiags(W_diagonal, 0, size(A,2), size(A,2));
    
        x = (AtA+Param.alpha^2*L'*W*L)\(At*b); %del Paper A Regularized Inverse Approach to Ultrasonic
        %Pulse-Echo Imaging página 3
        err = norm(x-x0)^2 / norm(x)^2; % relative error 
        x0 = x;
    end
    disp(x)
end

function x_reg = tikho_gpt(A, b, pars)
    % Tikhonov regularization to solve Ax = b
    % A: coefficient matrix
    % b: result vector
    % lambda: regularization parameter

    % Construct the Tikhonov regularization matrix (L)
    % For simplicity, using identity matrix as L, assuming a zero-order Tikhonov regularization
    lambda = pars.lambda;
    lambda = 0.11;
    L = speye(size(A,2)); 
    
    % Solve the regularized linear system
    % The regularization term is lambda^2 * L' * L
    % The modified system to solve is: (A'*A + lambda^2 * L'*L) x_reg = A'*b
    x_reg = (A'*A + lambda^2 * (L'*L)) \ (A'*b);
    
    % x_reg contains the regularized solution
end
