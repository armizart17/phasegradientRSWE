function [grad_z,grad_x,k,sws_matrix] = phase_estimator_TikhoVersion(u, w_kernel,f_v,dinf,og_size,constant,lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that yields the shear wave speed of a region with the 
% phase gradient method with QR solver MATLAB first version. 
% 
% Inputs:  
%          u           : 2D region of interest to evaluate (previously mirror padding)
%          w_kernel    : vector with the size of the window kernel
%          f_v         : vibration frequency
%          dinf        : structure that contains the spatial resolutions
%          og_size     : vector containing the original size of the data
%          lambda      : regularization coefficient

% Outputs: 
%          grad_z       : Gradient matrix for the axial direction
%          grad_x       : Gradient matrix for the lateral direction
%          k            : Total Wavenumber matrix 
%          sws_matrix   : Shear wave speed matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    res_z = dinf.dz; % Axial resolution
    res_x = dinf.dx; % Lateral resolution
    
     % Axis of the kernels
    z_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(1))*res_z; 
    x_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(2))*res_x; 
    
    %% Initializing vectors and matrixes
    grad_z = zeros(og_size); 
    grad_x = grad_z; % matrixes containing the estimated wavenumber along each direction
    
    angle_u = angle(u);
    %angle_u_unwrp = unwrap(angle_u, [], direction);
    %angle_u_unwrp = unwrap(angle_u, [], 2);   
    

    [X,Z] = meshgrid(x_axis,z_axis);
    A = [X(:) Z(:) ones(length(x_axis)*length(z_axis),1)];

	 [U,S,V] = svd(A); % SVD for Tikho
     s = diag(S);
	 
    angle_z = unwrap(angle_u,[],1);
    for ii = 1:og_size(1)
        for jj=1:og_size(2) %% for faster computing pararell toolbox
            area = angle_z(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            %% QR solver
            b = area(:);
            
			[results,rho,eta] = tikhonov(U,s,V,b,lambda);
            grad_z(ii,jj) = results(2);
        end
    end
    
    angle_x = unwrap(angle_u,[],2);
    for ii = 1:og_size(1)
        for jj=1:og_size(2) %% for faster computing pararell toolbox
            area = angle_x(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            %% QR solver
            b = area(:);
            [results,rho,eta] = tikhonov(U,s,V,b,lambda);
            grad_x(ii,jj) = results(1);
        end
    end
    
    k = ((grad_x.^2 + grad_z.^2)/constant).^0.5;

    sws_matrix = (2*pi*f_v)./k;
end