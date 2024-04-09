% A.1 Angle Compounding

% This function processes angle-compounded IQ data.
% A moving average across all frames is returned.
% =========================================================================
function IQ = A1_compound_angle_IQ( IQ_raw, num_angles )

IQ = nan( size(IQ_raw) - [0 0 num_angles-1] );
% IQ = nan([size(IQ_raw,1) size(IQ_raw,2) round(size(IQ_raw,3)/num_angles)] );
for frm = 1 : size(IQ,3)
    % averaging every num_angles consecutive frames
    % the index 3 indicates an averaging over the 3rd dimension
    IQ(:,:,frm) = mean( IQ_raw(:,:,frm:frm+num_angles-1), 3);
%     IQ(:,:,frm) = mean( IQ_raw(:,:,num_angles*(frm-1)+1:num_angles*frm), 3);
end
