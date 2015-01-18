classdef NeighborhoodVideoReader < handle
    %Reads a neighborhood of frames into a buffer
    
    properties
        inputVideoReader
        neighborhoodSize
        neighborhood % frame buffer
        frame
        index
    end
    
    methods
        function NVR = NeighborhoodVideoReader(inputVideoName, neighborhoodSize)
            NVR.inputVideoReader = VideoReader(inputVideoName);
            NVR.neighborhood = zeros(NVR.inputVideoReader.Height, NVR.inputVideoReader.Width, 3,neighborhoodSize, 'uint8');
            NVR.neighborhoodSize = neighborhoodSize;
            for frame = 1 : NVR.neighborhoodSize;
                NVR.neighborhood(:,:,:,frame) = read(NVR.inputVideoReader, frame);
            end
            NVR.index = NVR.neighborhoodSize / 2;
            NVR.frame = frame;
        end
        
        function read(NVR)
            if(NVR.frame > NVR.inputVideoReader.NumberOfFrames - 1)
                err = MException('ResultChk:OutOfRange', 'Attempted reading beyond file.');
                throw(err);
            end
            NVR.neighborhood = cat(4, NVR.neighborhood, read(NVR.inputVideoReader, NVR.frame));
            NVR.frame = NVR.frame + 1;
            NVR.index = NVR.index + 1;

            if(mod(NVR.frame, NVR.neighborhoodSize) == 0)
               NVR.neighborhood(:,:,:,1:end - NVR.neighborhoodSize) = [];
               NVR.index = NVR.neighborhoodSize / 2;
            end 
        end
    end    
end

function videocompletion(inputVideoName, outputVideoName, maskVideoName)
%% videocompletion: loads input video file and copies to output AVI
%  videocompletion('../data/crossing_ladies_input.avi', '../data/crossing_ladies_myoutput.avi', '../data/crossing_ladies_mask_sequence.avi')
    % Open input & output video files
    starttime = cputime;
    NEIGHBORHOOD_SIZE = 30;
    inputNVR = NeighborhoodVideoReader(inputVideoName, NEIGHBORHOOD_SIZE);
    neighborhood = inputNVR.neighborhood;
    index = NEIGHBORHOOD_SIZE / 2;%inputNVR.index;
    maskVideo = VideoReader(maskVideoName);
	outputVideo = VideoWriter(outputVideoName);
	open(outputVideo);
    
    % Define the hole H from a mask input video
    H = [];
    for t = 49 : maskVideo.NumberOfFrames
        mask = read(maskVideo, t);
        imshow(mask);
        [Hx, Hy] = find(rgb2gray(mask));
        H = cat(1, H, [Hx, Hy, ones(size(Hx, 1), 1) .* t, zeros(size(Hx,1),1)]);
    end

    % Compute the dist from boundary beforehand
    % Hunt over mask and read for each hole pixel, what is the closest
    % non-hole pixel.
    % NOTE: Very parallel
    
    % Distance transform?
    [nHx, nHy] = find(mask == 0);
    for x = 1 : size(Hx)
        mindist = inf;
        h1 = H(x,:,:,:);
        for nx = 1 : 5%size(nHx)
            nx1 = nHx(nx);
            ny1 = nHy(nx);
            % Todo: include time in dist!
            dist = sqrt((h1(1) - nx1) .* (h1(1) - nx1) + (h1(2) - ny1) .* (h1(2) - ny1));
            if(dist < mindist)
                mindist = dist;
            end
        end
        % save mindist
        h1(4) = mindist;
        H(x, :, :, :) = h1;
    end

	% Read neighborhood range of input frames into buffer, outputting to
	% video exactly
    
    %TOOD: impyramid
    cnt = 0;
    horizontal = fspecial('sobel');
    Hrow = 1;
    Ht = H(Hrow,:,:,:); % Start filling on first frame in time valid
    lasttime = H(end,:,:,:);
    lasttime = lasttime(3);

	for frame = inputNVR.frame : (lasttime + NEIGHBORHOOD_SIZE / 2)
        % Read neighborhood
        neighborhood = cat(4, neighborhood, read(inputNVR.inputVideoReader, frame));
        index = index + 1;

        % Advance neighborhood
        if(mod(frame, NEIGHBORHOOD_SIZE) == 0)
           neighborhood(:,:,:,1:end - NEIGHBORHOOD_SIZE) = [];
            index = NEIGHBORHOOD_SIZE / 2;
        end

        % Get neighborhood centered on hole time index (confirmed working)
        if(frame == Ht(3) + (NEIGHBORHOOD_SIZE / 2))
            S = neighborhood(:,:,:,index:index + NEIGHBORHOOD_SIZE / 2);
            nextHrow = Hrow;
            % zero it out to show what pixels get filled
            for row = Hrow : size(H,1)
                p = H(row, :);
                if(p(3) == Ht(3))
                    neighborhood(p(1), p(2), :, index) = [0,0,0];
                else
                    nextHrow = row;
                    break;
                end
            end

            % 2.2 LOCAL SIMILARITY MEASURE
            % Local implementation of wexler for single hole pixel at t for
            % implementation ease
            SGray = S(:,:,1,:) + S(:,:,2,:) + S(:,:,3,:);
            HGradS = imfilter(SGray(:,:,:,:), horizontal);
            VGradS = imfilter(SGray(:,:,:,:), horizontal');
            for row = 1 : (nextHrow - Hrow)
                p = H(Hrow + row, :);
                p(3) = size(S, 4) / 2;
                % Compute distance from closest boundary of hole to weight
                alpha = power(1.3, -p(4));
                c = equation3(p, alpha, S, HGradS, VGradS);
                neighborhood(p(1), p(2), :, index) = c;
            end
            % Advance H to next frame
            Hrow = nextHrow;
            Ht = H(Hrow,:,:,:);
        end
        
        % Write to output
        writeVideo(outputVideo, neighborhood(:, :, :, index));
        cnt = cnt + 1;
	end
	close(outputVideo);
    fprintf(1,'100%% complete.\n%d seconds spent in processing.\n%d frames written out of %d total.\n', round(cputime - starttime), cnt, inputNVR.inputVideoReader.NumberOfFrames);
    end

    % equation 3 - Implement improvements such as mean shift algorithm later.
% TODO: Check video borders to verify patch exists, for now, limit p
% minimize the variance of the color
function c = equation3(p, alpha, frames, HGradFrames, VGradFrames)
    % TODO: Optimize patch iteration!
	%neighbors = reshape(meshgrid(-2:2, -2:2, -2:2), 1, []); %? Maybe...
	%for all neighbors i of p
    max_sim = 0;
    
    c = [1; 1; 1];
	for ix = -2  : 2
		for iy = -2 : 2
			for it = -2 : 2
				ip = [p(1) + ix, p(2) + iy, p(3) + it];
				Wip = frames(...
                    ip(1) - 2 : ip(1) + 2,...
                    ip(2) - 2 : ip(2) + 2,...
                    :,...
                    ip(3) - 2 : ip(3) + 2);
                HGradWip = HGradFrames(...
                    ip(1) - 2 : ip(1) + 2,...
                    ip(2) - 2 : ip(2) + 2,...
                    :,...
                    ip(3) - 2 : ip(3) + 2);
                VGradWip = VGradFrames(...
                    ip(1) - 2 : ip(1) + 2,...
                    ip(2) - 2 : ip(2) + 2,...
                    :,...
                    ip(3) - 2 : ip(3) + 2);
                for jx = -1  : 1
                    for jy = -1 : 1
                        jp = [ip(1) + jx, ip(2) + jy, ip(3)];
                        Vi = frames(...
                        jp(1) - 2 : jp(1) + 2,...
                        jp(2) - 2 : jp(2) + 2,...
                        :,...
                        jp(3) - 2 : jp(3) + 2);
                        HGradVi = HGradFrames(...
                        jp(1) - 2 : jp(1) + 2,...
                        jp(2) - 2 : jp(2) + 2,...
                        :,...
                        jp(3) - 2 : jp(3) + 2);
                        VGradVi = VGradFrames(...
                        jp(1) - 2 : jp(1) + 2,...
                        jp(2) - 2 : jp(2) + 2,...
                        :,...
                        jp(3) - 2 : jp(3) + 2);
                        current_sim = sum(sim(Wip, Vi, alpha, HGradWip, VGradWip, HGradVi, VGradVi));
                        if(current_sim > max_sim)
                            max_sim = current_sim;
                            omega = alpha .* max_sim;
                            c = Vi(2, 2, :, 2);
                        end
                    end
                end
			end
		end
	end
    % TODO: Max likelihood solution
end


% equation 2: at least 3x3x3 Wp and Vq
function res = sim(Wp, Vq, alpha, HGradWp, VGradWp, HGradVq, VGradVq)
    % sigma = 75th percentile of all distances in current search in all locations

    sigma = 5.0 / 255.0; % From short paper
    Wpm = measurement(Wp, alpha, HGradWp, VGradWp);
    Vqm = measurement(Vq, alpha, HGradVq, VGradVq);
    d = Wpm - Vqm;
    s = sum(d .* d);
    res = sum(s);
    res = exp(-res ./ (2.0 * sigma * sigma));
end

% similarity measure: at least 3x3x3 Wp and Vq
function sum = d(Wp, Vq, alpha, HGradWp, VGradWp, HGradVp, VGradVp)
    sum = 0.0;
    for x = 2 : size(Wp, 1) - 1
        for y = 2 : size(Wp, 2) - 1
            for t = 2 : size(Wp, 4) - 1
                Wpm = measurement(Wp(x-1:x+1, y-1:y+1, :, t), alpha, HGradWp(x-1:x+1, y-1:y+1, :, t), VGradWp(x-1:x+1, y-1:y+1, :, t));
                Vqm = measurement(Vq(x-1:x+1, y-1:y+1, :, t), alpha, HGradVp(x-1:x+1, y-1:y+1, :, t), VGradVp(x-1:x+1, y-1:y+1, :, t));
                
                dif = double(Wpm - Vqm);
                sum = sum + dif .* dif;
            end
        end
    end
end

% Measurement: takes in an at least 3x3x1 RGB patch
% HGradY is horizontal gradient of gray version of S
function m = measurement(S, alpha, HGradY, VGradY)
    % Yx = dY/dx : Can use sobel
    % Yy = dY/dy
    % Yt = frame difference between points being compared.
    % TODO: u and v are NORM flow vectors:
   % alpha = 5;
    %u = Yt/Yx;
    %v = Yt/Yy;
    % Get Yt = frame difference between points being compared.,
    % should be distance from t of hole.
    %measurement = [Yr, Yg, Yb, alpha * u, alpha *v];
    m = [S(:,:,1,:), S(:,:,2,:), S(:,:,3,:), alpha .* HGradY, alpha .* VGradY];
end
