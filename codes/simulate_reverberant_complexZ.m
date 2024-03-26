function [pv,x,z,t] = simulate_reverberant_complexZ(numOndas, freq, c_back, c_inc, vid_duration)
%------------------------------------------------------------------------------------
%   Inputs
%       numOndas        Numero de ondas, puede ser 1000 - 10 0000
%       freq            Vibration frequency
%       c_back          SWS background
%       c_inc           SWS inclusion
%       vid_duration    Duracion del video, recomiendo 10 ms aunque si se
%                       puede  50ms mejor
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
resT = 1/(3600); % 3600 Hz: PRF
t(1,1,:) = 0:resT:vid_duration;

% Elasticity parameters
k = 2*pi*freq/c_back; % background
K = k*ones(size(X));
rr = sqrt(X.^2+Z.^2);
K(rr<10E-3) = 2*pi*freq/c_inc;

% Random number generation
Alpha = rand(numOndas,1)*2*pi;
vqa = randn(numOndas,1);
shift = randn(3,numOndas);
shift = shift./sqrt(sum(shift.^2)) * 50E-3;
nq = randn(3,numOndas);
nq = nq./sqrt(sum(nq.^2));

pv = zeros([size(K),length(t)]);
for q = 1:numOndas    
    npa_z = -sin(Alpha(q))*sqrt(1 - nq(3,q)^2); % La raiz es sin(arcos(Theta))
    phase =K.*( nq(1,q)*(X-shift(1,q)) + nq(2,q)*(0-shift(2,q)) + nq(3,q)*(Z-shift(3,q)) );
    pv = pv + npa_z * vqa(q) * exp(1j*( phase +w*t));   
end

end