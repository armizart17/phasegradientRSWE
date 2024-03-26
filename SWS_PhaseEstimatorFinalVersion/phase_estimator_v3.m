function [k_z,k_x,c0,k,sws_matrix] = phase_estimator_v3(u, w_kernel,f_v,dinf,og_size, constant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that yields the shear wave speed of a region with the 
% phase gradient method. 
% 
% Inputs:  
%          u           : 2D region of interest to evaluate (previously mirror padding)
%          w_kernel     : vector with the size of the window kernel
%          f_v         : vibration freqeuncy
%          dinf        : structure that contains the spatial resolutions
%          og_size     : vector containing the original size of the data
%          constant    : constant from equations (chech paper 
%                        - typically 15*sqrt(pi/2)

% Outputs: 
%          k_z         : Wavenumber matrix for the axial direction
%          k_x         : Wavenumber matrix for the lateral direction
%          k         : Total Wavenumber matrix 
%          sws_matrix  : Shear wave speed matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    res_z = dinf.dz; % Axial resolution
    res_x = dinf.dx; % Lateral resolution
    
     % Axis of the kernels
    z_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(1))*res_z; 
    x_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(2))*res_x; 
    
    step = 2;
    
    %% Initializing vectors and matrixes
    k_z = zeros(length(1:step:og_size(1)),length(1:step:og_size(2))); 
    k_x = k_z; % matrixes containing the estimated wavenumber along each direction
    c0 = k_z;
    angle_u = angle(u);
    %angle_u_unwrp = unwrap(angle_u, [], direction);
    %angle_u_unwrp = unwrap(angle_u, [], 2);   
    
    % ----- Sebastian ----
    [X,Z] = meshgrid(x_axis,z_axis);
    A = [X(:) Z(:) ones(length(x_axis)*length(z_axis),1)];
    
    [U,S,V] = svd(A); % SVD for Tikho
    s = diag(S);
    lambda = 0.0001e-4; 
    lambda = 0;
    search_z = 1:step:og_size(1); % To save time, the windowed matrix is evaluated in points separated by a step of 2
    
    % ------------------
    for ii = 1:length(search_z)
        search_x = 1:step:og_size(2);
        
        for jj = 1:length(search_x) %% for faster computing pararell toolbox
            area = angle_u(search_z(ii):search_z(ii)+w_kernel(1)-1,search_x(jj):search_x(jj)+w_kernel(2)-1); % Window kernel
%             area = unwrap(unwrap(area), [],2); % doble unwraping
%             area = unwrap((area), [],1);
            area = unwrap((area), [],2);
%             area = unwrap((area), [],2);
            %% Least squares
            b = area(:);
%             results = LS_mat * b;
            
%             [x_lambda,rho,eta] = tikhonov(U,s,V,b,lambda);
  
            x_lambda=A\b;
% chvr ahi_vF
           % x_lambda = cgs(A,b);
            k_z(ii,jj) = x_lambda(2);
            k_x(ii,jj) = x_lambda(1);
            c0(ii,jj) = x_lambda(3);
        end
    end
%     k = (constant*(k_x.^2+k_z.^2)).^0.5;
    k = (1.*(k_x.^2+k_z.^2)).^0.5;
    sws_matrix = (2*pi*f_v)./k;
end
    