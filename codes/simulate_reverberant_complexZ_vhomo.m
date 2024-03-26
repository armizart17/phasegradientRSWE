function [pv,x,z] = simulate_reverberant_complexZ_vhomo(numWaves, freq, c_back)
%------------------------------------------------------------------------------------
%   Inputs
%       numWaves        Numero de ondas, puede ser 1000 - 10 0000
%       freq            Vibration frequency
%       c_back          SWS background
%
%   Outputs
%       pv              Particle velocity en el componente axial
%       x,z,t           ejes coordenados
%------------------------------------------------------------------------------------

% Field dimensions
Xmin = -2.5e-2; Xmax = 2.5e-2; %5cm or 50mm
resX = 0.1E-3;
x = Xmin:resX:Xmax;
Zmin = -2.5e-2; Zmax = 2.5e-2;
resZ= resX;
z = Zmin:resZ:Zmax;
[X,Z]=meshgrid(x,z);

% Time parameters
w = 2*pi*freq;

% Elasticity parameters
k = 2*pi*freq/c_back; % background
K = k*ones(size(X));

% Random number generation
Alpha = rand(numWaves,1)*2*pi;
vqa = randn(numWaves,1);
shift = randn(3,numWaves);
shift = shift./sqrt(sum(shift.^2)) * 50E-3;
nq = randn(3,numWaves);
nq = nq./sqrt(sum(nq.^2));

pv = zeros(size(K));

    for q = 1:numWaves    
        npa_z = -sin(Alpha(q))*sqrt(1 - nq(3,q)^2); % La raiz es sin(arcos(Theta))
        phase =K.*( nq(1,q)*(X-shift(1,q)) + nq(2,q)*(0-shift(2,q)) + nq(3,q)*(Z-shift(3,q)) );
        pv = pv + npa_z * vqa(q) * exp(1j*( phase ));   
    end

end