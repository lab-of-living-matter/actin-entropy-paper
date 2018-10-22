% function cytosimDataExtractor(input_path, output_path, size_xy, nframes, first_frame, n_corrections)
%
% Function to take cytosim outputs and put them into a flat csv file of positions.
%
% Parameters
% ----------
% input_path : str
%     full path to where cytosim output lives
% output_path : str
%     full path to where parsed data will be saved
% size_xy : [float, float]
%     2 element array that indicates the size of the simulation box in x and y direction
% nframes : int
%     the total number of frames to analyze
% first_frame : int
%     the first frame to analyze, can be used to skip early, equilibration frames simulation
% n_corrections : int
%     some filaments get pushed outside size_xy. Perform a correction for this n_corrections
%     times to fix the problem. 4 is usually a good choice
%
% Written by Vikrant Yadav, 10/22/2018
function cytosimDataExtractor(input_path, output_path, size_xy, nframes, first_frame, n_corrections)

    % read cytosim output file
    posfile = posreader3(input_path);

    % parameters of interest
    sx = size_xy(1); % size of simulation dimension along x
    sy = size_xy(2); % size of simulation dimension along y
    % num = 1;

    xdim = sx ./ 2; % parameter to take care of spacing in cytosim. Half of sx
    ydim = sy ./ 2; % half of sy

    fnum = nframes; % number of frames to analyze
    fini = first_frame; % first frame to analyze. Depends on number of frames over which filaments in simulation are growing
    % lensc = 1; % length scale of interest (to coarse grain)
    finres = [0 0 0 0];

    % Analysis
    tic;
    for fnum = fini:fini + fnum - 1;
        disp(fnum);
        fr = find(posfile(:, 2) == 'frame'); % find locations where the text 'frame' occurs
        data1 = posfile(fr(fnum, 1):fr(fnum + 1, 1), :); % extract all data between consecutive appearence of frames
        listf = find(data1(:, 2) == 'fiber'); % find locations where the text 'fiber' occurs in that frame
        nf = size(listf);
        ctr = 1;
        for i = 1:nf(1) - 1
            p = data1(listf(i, 1) + 1:listf(i + 1, 1) - 1, :); % extract data between conscutive occurence of term fiber
            p = str2double(p(:, :)); % convert text to number for doing mathsy stuff later.
            l = size(p);
            for j = 1:l(1) - 1 % do successive substraction here to get local orientation information
                res1(ctr, 1) = p(j, 2);
                res1(ctr, 2) = p(j, 3);
                res1(ctr, 3) = p(j, 1);
                res1(ctr, 4) = fnum; %p(j,2);
                ctr = ctr+1;
            end
        end

        % Some filaments are not parsed correctly by cytosim, this block will take care of it.
        counter = 0;
        while counter < n_corrections
            fixlist = find(res1(:, 1) > xdim);
            res1(fixlist, 1) = res1(fixlist, 1) - sx;
            fixlist = find(res1(:, 1) < -xdim);
            res1(fixlist, 1) = res1(fixlist, 1) + sx;
            fixlist = find(res1(:, 2) > ydim);
            res1(fixlist, 2) = res1(fixlist, 2) - sy;
            fixlist = find(res1(:, 2) < -ydim);
            res1(fixlist, 2) = res1(fixlist, 2) + sy;

            counter = counter + 1;
        end

        finres = cat(1, finres, res1);

    end
    finres(1, :) = [];
    toc;
    disp('Saving output...');
    dlmwrite(output_path, finres); %Change content inside quotes to output file name
    disp('... saved');
end
