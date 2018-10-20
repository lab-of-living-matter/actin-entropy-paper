% function aa = getAxonemeModes(tangents, nmax, px2um, dt)
%
% This function rotates the filaments so that they lie horizontally,
% with the first frame end points of the filament determining what is
% 'horizontal'. Get tangent angle data along the arc length of the filament.
% Plots both of these things together. Then also decomposes the tangent angle
% into elastohydrodynamic modes of thin rods using elastohydroModes
%
% Parameters
% ----------
% tangents : array
%     Array of filament tangent angles. Each row is a position
%     along the arc length, each column is a time point
% nmax : int
%     Maximum number of modes to calculate for the filament
% px2um : float
%     Conversion factor of pixels to microns. Units: [um/pixel]
% dt : float
%     time between frames. Used to subtract linear increase of tangent angle
%     measured as axoneme rotates continuously
%
%
% Returns
% -------
% aa : array
%     (N)x(nmax) array of the mode coefficients for every
%     each time point contained in a single row
%
% This version is made for axoneme data from the following Dryad database:
% http://datadryad.org/resource/doi:10.5061/dryad.0529j
%
% For the following paper:
%
% Dynamic curvature regulation accounts for the symmetric and asymmetric beats of
% Chlamydomonas flagella
% Sartori P, Geyer VF, Scholich A, Jülicher F, Howard J
% Date Published: May 12, 2016
% DOI: http://dx.doi.org/10.5061/dryad.0529j
%
% Created by Daniel Seara, 07/27/2017
function aa = getAxonemeModes(tangents, nmax, px2um, dt)
    aa = [];

    % First subtract off linear portion of tangents
    t = cumsum(ones(size(tangents, 2), 1) .* dt);
	p = polyfit(t, tangents(1,:)',1);
	tangents = tangents - repmat((p(1) * t' + p(2)), [size(tangents, 1), 1]);

    nframes = size(tangents, 2);

    for jj = 1:nframes
        theta = tangents(:,jj);
        ds = [0; ones(numel(theta) - 1, 1)] .* px2um;
        L = (numel(theta)-1) * px2um;
        % keyboard
        elastoCoeffs = elastohydroModes(theta', cumsum(ds)' - L/2, L, nmax, bcs);
        aa = [aa; elastoCoeffs'];
    end % End loop over frames
end
