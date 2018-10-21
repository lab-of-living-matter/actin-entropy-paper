% function aa = getJFilamentModes(filament,nmax,bcs,savestuff)
% This function rotates the filaments so that they lie horizontally,
% with the first frame end points of the filament determining what is
% 'horizontal'. Get tangent angle data along the arc length of the filament.
% Decomposes the tangent angle into elastohydrodynamic modes of thin rods
%
% Parameters
% ----------
% filament : array
%     Mx5 array filament positions over time, output of loadJFilamentData. Array structure:
%     frame | segment # | x | y | all zeros (?)
% nmax : int
%     Maximum number of modes to calculate for the filaments
% px2um : float
%     Conversion factor of pixels to microns. Units: [um/pixel]
%
% Returns
% -------
% aa : array
%     Mx(nmax) array of the mode coefficients for every filament at every time point,
%     each time point contained in a single row
% aaEach Filament : cell array
%     cell array with mode amplitudes for each individual filament in each element
%     of the cell array
% avgL : array
%     array where nth element is the average length measured over time of the nth
%     filament.
%
% Created by Daniel Seara, 05/10/2017
function [aaAll, aaEachFilament, avgL] = getJFilamentModes(filament, nmax, px2um)
    aaAll = [];
    avgL = [];
	for ii = 1:numel(filament)
        lArr = [];
        aa = [];
		filaData = filament{ii};
		frames = unique(filaData(:,1));
		for jj = 1:numel(frames)
			f = frames(jj);
			xy = [filaData(filaData(:,1)==f,3)'.*px2um;...
				  filaData(filaData(:,1)==f,4)'.*px2um];
			npoints = size(xy,2);
			if jj==1
				theta = atan2(xy(2,end)-xy(2,1), xy(1,end)-xy(1,1));
				% rotate clockwise by angle theta
		    	R = [cos(theta), sin(theta);-sin(theta), cos(theta)];
		    end % End if statement to get rotation angle

		    % Rotate and subtract starting position
	        xy = R*(xy - [xy(1,1).*ones(1,size(xy,2));xy(2,1).*ones(1,size(xy,2))]);
	        displacements = diff(xy, 1,2); % take difference along columns
	        ds = sqrt(displacements(1,:).^2 + displacements(2,:).^2);
	        L = sum(ds);
	        lArr = [lArr, L];
	        % minus sign because of flip in y-axis between image and plot
	        tangents = -atan2(displacements(2,:), displacements(1,:));
	        elastoCoeffs = elastohydroModes(tangents, cumsum(ds) - L/2, L, nmax, bcs);
	        aa = [aa; elastoCoeffs'];
		end % End loop over frames
	    aaAll = [aaAll;aa];
        aaEachFilament{ii} = aa;
        avgL = [avgL; mean(lArr)];
	end % End loop over filaments
end
