% function aa = getJFilamentModes(filament,nmax,bcs,savestuff)
% This function rotates the filaments so that they lie horizontally,
% with the first frame end points of the filament determining what is
% 'horizontal'. Get tangent angle data along the arc length of the filament.
% Decomposes the tangent angle into elastohydrodynamic modes of thin rods
%
% Parameters
% ----------
% simData : array
%     The data format in files is as follows:
%     x-position, y-position, FilamentID, FrameNumber
% nmax : int
%     Maximum number of modes to calculate for the filaments
% bcs : string
%     String of boundary conditions of the rods (use 'free')
%
% Returns
% -------
% aa : array
%     Mx(nmax) array of the mode coefficients for every filament at every time point,
%     each time point contained in a single row
%
% Created to read cytosim data found in the following location:
% llmStorage203/Vikrant/DannyData
%
% Created by Daniel Seara, 06/07/2018
function [aaAll, avgL] = getModes(simData, nmax, bcs)
    fids = unique(simData(:,3));
    frames = unique(simData(:,4));
    aaAll = zeros(numel(fids) .* numel(frames), nmax);
    avgL = zeros(numel(fids), 1);
    for ii = 1:numel(fids)
        lArr = zeros(numel(frames),1);
        aa = zeros(numel(frames), nmax);
        for jj = 1:numel(frames)
            disp(['Filament ', num2str(ii), ' of ', num2str(numel(fids)), '. Frame ', num2str(jj), ' of ', num2str(numel(frames))])
            filaData = simData(simData(:,3) == fids(ii) & simData(:,4) == frames(jj), [1,2]);
            xy = [filaData(:, 1)'; filaData(:, 2)'];
            npoints = size(xy,2);

            if npoints < 3
                % skip filaments with only 2 points in them
                continue
            end

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
            lArr(jj) = L;
            % minus sign because of flip in y-axis between image and plot
            tangents = -atan2(displacements(2,:), displacements(1,:));
            % keyboard
            elastoCoeffs = elastohydroModes(tangents, cumsum(ds), L, nmax, bcs);
            aa(jj, :) = elastoCoeffs';
        end % End loop over frames

        aaAll((ii-1)*jj + 1 : ii*jj, :) = aa;
        % aaEachFilament{ii} = aa;
        avgL(ii) = mean(lArr);
    end % End loop over filaments
end
