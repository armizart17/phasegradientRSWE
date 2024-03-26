function [x, cost, error, fide, regul] = pdo_TVV_single(Y, A, lambda, tau, maxIter, tol, numberEstimators, stableIter)
% function [x, cost, error, fide, regul] = pdo_TVV_single(Y, A, lambda, tau, maxIter, tol, numberEstimators, stableIter)
% SINGLE CHANNEL TOTAL VARIATION for Ax = Y, based on Primal Dual Algorithm
% Condat, 
% Edmundo A. Miranda
	
    rho = 1.99;		% relaxation parameter, in [1,2]
	sigma = 1/tau/8; % proximal parameter

    global H W numEstim;
    
    [H,W]=size(Y);
    numEstim = numberEstimators;

    y = Y(:);
    At = A';
    AtA = At*A;
    
%%%% INITIALIZATION %%%%
%     x_ini = A\y;                  % INVERSE PROBLEM  
%     x_ini = cgs(A'*A, A'*y);      & CGS SOLUTION

%     x_ini = zeros([H, W]);  % ZEROS INITIALIZATION
    x_vect_ini = cgs(A'*A, A'*y); 
    x_ini = reshape(x_vect_ini, H, W);
 %%%% INITIALIZATION %%%%  
 
    %%% JUST ONCE %%%%
    global izq 
    [~, n] = size(A);   
    izq = ( sparse(eye(n)) + tau*(AtA) )^-1;
    %%% JUST ONCE %%%%

    x2 = reshape(x_ini, [H,W]);
   
    x2(:,:,1) = Y(:,:,1 ); % 1st estimator medium value of SLD_TERM

    u2 = zeros([size(x2) 2]);

	ee = 1;	 
    iter = 1;
    
    fide = zeros(1,maxIter);
    cost = fide; regul = fide; error = fide;
    
    error(1) = 1;
    cost(1) = 100;
    fide(1) = 1;
    regul(1) = 1;
   
    disp('Itn   |    Cost   |   Delta(Cost)');
    disp('----------------------------------');
    
    while (iter < maxIter) && (ee > tol)    

        x = prox_tau_f_generic(x2-tau*opDadj(u2),Y,tau, A, AtA, At);
		u = prox_sigma_g_conj(u2+sigma*opD(2*x-x2), lambda);
        
		x2 = x2 + rho*(x-x2);
		u2 = u2 + rho*(u-u2);
        
  %         cost(iter+1) = 0.5*sum(sum(sum((x(:) - y).^2)))+lambda*TVnuclear_EMZ(x);
%       cost(iter+1) = 0.5* norm(A*x(:) - y,2).^2 + lambda*TVnuclear_EMZ(x);
        
        fide(iter+1) = 0.5*sum(sum(sum((reshape(A*x(:), size(Y)) - Y).^2)));
%         regul(iter+1) = lambda*sum(sum(sqrt(sum(opD(x).^2,3))));
        regul(iter+1) = lambda*TVV(opD(x));
        
        cost(iter+1) = fide(iter+1) + regul(iter+1);
%         cost(iter+1) = 0.5*sum(sum(sum((reshape(A*x(:), size(Y)) - Y).^2)))+lambda*TVnuclear_EMZ(x);
        
        ee = abs(cost(iter+1) - cost(iter));     
        % NO NORMALIZADO
%         error(iter+1) = ee / cost(iter+1);
        
        % NORMALIZADO
        ee = ee / cost(iter+1);
        error(iter+1) = ee;
        
        if (iter < stableIter) % en las primeras 100 iteraciones inestable 
            ee = 1; % poner error grande (cualquier valor) para que no pare iteracion            
        end
    
		if mod(iter,5)==0
            
%             primalcost = 0.5* norm(A*x(:) - y,2).^2 + lambda*TVnuclear_EMZ(x);
                       
            fprintf('%4d  | %f | %e\n',iter, cost(iter+1), error(iter+1));

        end
        iter = iter+1;
    end
return

function [u_3d] = prox_tau_f_generic(v, B, tau, A, AtA, At)

    global H W ;
    
    global izq
%     [~, n] = size(A);   
%     izq = sparse(eye(n)) + tau*(AtA);

    der = v(:) + tau*At*B(:);
    u = ( ( izq ) ) * ( der );
     
    u_3d = reshape(u, [H W]);
    
return
    
function [u] = prox_sigma_g_conj(v, lambda)
    
%     bsxfun(@rdivide,u,max(sqrt(sum(u.^2,3))/lambda,1));
    u = v ./ max( sqrt( sum (v.^2, 3) ) / lambda, 1);
    
    
return

function [u] = opD(v)
    % Description: Dy(:,:,:,1) Dx (:,:,:,2), 3rd dimension is channel
    [H,W] = size(v);

%     u = cat(4,[diff(v,1,1);zeros(1,W,C)],[diff(v,1,2) zeros(H,1,C)]);
      u = cat(3,[diff(v,1,1);zeros(1,W)],[diff(v,1,2) zeros(H,1)]);
return

function [u] = opDadj(v)
    % u = Dy Dx
%     u = -[v(1,:,:,1);diff(v(:,:,:,1),1,1)]-[v(:,1,:,2) diff(v(:,:,:,2),1,2)];	
      u = -[v(1,:,1);diff(v(:,:,1),1,1)]-[v(:,1,2) diff(v(:,:,2),1,2)];	 
return
    
function [sumF] = TVV(opD)
    % Description: Dy(:,:,:,1) Dx (:,:,:,2), 3rd dimension is channel
%     Dx = opD(:,:,2);
%     Dy = opD(:,:,1);
%     sum2D = 0;
%     for  ii = 1:size(opD,3) %% number of channels
%         sum2D = sum2D + Dx(:,:,ii).^2 + Dy(:,:,ii).^2;
%     end
    
%     sum2D = sqrt(sum2D);
    
%     sumF = sum(sum(sum2D));

    sumF = sum(sum(sqrt(sum(opD.^2,3))));
return

