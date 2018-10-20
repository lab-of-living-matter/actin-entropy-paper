% function filament = loadJFilamentData(snakes)
%
% This function reads in the data given by JFilament
% Load the txt file that JFilament spits out as a table with NaNs
% replacing the useless elements.
% 
% INPUTS  
%			snakes: matlab table of imported JFilament .txt output file
%
% OUTPUTS 
%		  filament: a cell array, where each cell element is a Nx5 array
%					where N=(number of segments of filament)x(number of frames)
%					with columns:
%					frame | segment # | x | y | all zeros (?)
%
% Created by Daniel Seara, 05/10/2017
function filament = loadJFilamentData(snakes)

	snakes.Properties.VariableNames = {'frame','segmentNumber','x','y','var5'};

	% First, get all the data out of the table made from the output txt file.
	% All the empty spaces are replaced with NaNs.
	% Every new filament is tracked starting with a row of all NaNs,
	% followed by a row with a 0 and all NaNs, then all the numeric stuff
	inds = find(isnan(snakes.frame));
	inds = inds(10:end); % First 9 rows just configuration info

	% Split up all the filament data
	nFila = numel(inds);
	for ii = 1:nFila-1
		% This won't work for the last filament, can't + 1 if ii goes to end
		filament{ii} = table2array(snakes(inds(ii)+2:inds(ii+1)-1,:));
	end
	% Last filament
	filament{nFila} = table2array(snakes(inds(end)+2:end,:));