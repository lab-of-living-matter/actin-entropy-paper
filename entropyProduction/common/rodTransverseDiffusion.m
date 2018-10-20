% function dRot = rodRotationalDiffusion(L, r, kT, eta)
%
% Calculates the transverse diffusion coefficient of a rigid rod
%   D_trans = kBT / gamma_trans
% with
%
%                  4*pi*eta*L
% gamma_trans = -----------------
%               ln(L/(2r)) + 0.84
%
%
% From De la Torre, J., & Bloomfield, V. (1981). Quarterly Reviews of Biophysics, Eq 73
%
% Parameters
% ----------
% L : float
%     length of filament in um
% r : float
%     radius of filament in um
% kT : float (optional)
%     thermal energy in untis of pN um. Defaults T = 25 C
%     0.0041 pN*um
% eta : float (optional)
%     shear modulus of medium in pN s / um^2. Defaults to that of water
%     0.91 * 10^-3 pN s / um^2
%
% Returns
% -------
% dRot : float
%     rotational diffusion constant for a rod, calculated as given above.
function dTrans = rodTransverseDiffusion(L, r, kT, eta)

    switch nargin
        case 2
            eta = 0.00091;
            kT = 0.0041;
        case 3
            eta = 0.00091;
        otherwise
            error('Need at least 2 inputs (L & r) or at most 4 inputs (L, r, kT, eta). \n%s',...
                  ['Current number of inputs: ', num2str(nargin)])
    end

    drag = 4 * pi * eta * L / ((log(L / (2*r)) - 0.84));
    dTrans = kT / drag;
end
