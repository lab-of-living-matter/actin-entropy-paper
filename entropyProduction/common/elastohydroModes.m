% Finds mode decomposition of a vector into elastohydrodynamic modes of a free
% slender rod
%
% INPUT     theta : 1xN vector of tangent angles from filament positions
%           ds    : 1xN vector of displcements between each filament data point
%           L     : length of filament
%           nmax  : max number of modes to use
%           bcs   : boundary conditions for the mode (right now, only use free ends)
%
% OUTPUT    a     : (nmax)x1 vector with the coefficients

function a = elastohydroModes(theta, ds, L, nmax, bcs)

    if ~strcmp(bcs,'free')
        error('I don''t know how to do anything but free boundary conditions yet...')
    end
    % Get positions to evaluate the function
    s = linspace(-L/2,L/2,numel(theta));
    a = zeros(nmax,1);
    for n = 1:2:nmax
        phi = oddModes(n,L,s);
%         keyboard
        a(n) = -L^(-1/2) .* trapz(ds, phi.*theta);
    end
    for n = 2:2:nmax
        phi = evenModes(n,L,s);
        a(n) = -L^(-1/2) .* trapz(ds, phi.*theta);
    end
end

function k = waveVec(n)
    k = (n+1/2)*pi;
end

function phi = oddModes(n,L,s)
    k = waveVec(n);
    phi =  L^(1/2)/k .* (sin((k/L) .* s) ./ cos(k/2) + sinh((k/L) .* s) ./ cosh(k/2));
end

function phi = evenModes(n,L,s)
    k = waveVec(n);
    phi =  L^(1/2)/k .* (cos((k/L) .* s) ./ sin(k/2) + cosh((k/L) .* s) ./ sinh(k/2));
end
