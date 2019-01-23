classdef get_Input
    %GET_INPUT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        % #################################################################
        % make some artificial patches
        function [times, indices] = artificial_rand_patches(num, method)
            if(nargin < 2)
                method = "deterministic";
            end
            % base patches for repeat
            L = 5;
            patches(:,:,1) = [0 0 0 0 0;...
                              0 0 0 0 0;...
                              1 1 1 1 1;...
                              0 0 0 0 0;...
                              0 0 0 0 0];
            patches(:,:,2) = [0 0 0 0 0;...
                              1 1 0 0 0;...
                              0 0 1 0 0;...
                              0 0 0 1 1;...
                              0 0 0 0 0];
            patches(:,:,3) = [1 0 0 0 0;...
                              0 1 0 0 0;...
                              0 0 1 0 0;...
                              0 0 0 1 0;...
                              0 0 0 0 1];
            patches(:,:,4) = [0 1 0 0 0;...
                              0 1 0 0 0;...
                              0 0 1 0 0;...
                              0 0 0 1 0;...
                              0 0 0 1 0];
            patches(:,:,5) = [0 0 1 0 0;...
                              0 0 1 0 0;...
                              0 0 1 0 0;...
                              0 0 1 0 0;...
                              0 0 1 0 0];
            patches(:,:,6) = [0 0 0 1 0;...
                              0 0 0 1 0;...
                              0 0 1 0 0;...
                              0 1 0 0 0;...
                              0 1 0 0 0];
            patches(:,:,7) = [0 0 0 0 1;...
                              0 0 0 1 0;...
                              0 0 1 0 0;...
                              0 1 0 0 0;...
                              1 0 0 0 0];
            patches(:,:,8) = [0 0 0 0 0;...
                              0 0 0 1 1;...
                              0 0 1 0 0;...
                              1 1 0 0 0;...
                              0 0 0 0 0];                          
            [~,~,f] = size(patches);
            c = 0;
            times   = [];
            indices = [];
            param.binNumber     = 1000;
            param.maxPercSpikes = 0.25; % min time for first spike
            param.tau           = 0.2; 
            for r = 1:num
                for i=randperm(f,1)
                    switch(method)
                        case "deterministic"
                            [imT_times, imT_index] = mat2inSpike.deterministic(patches(:,:,i),param);
                        case "poisson"
                            [imT_times, imT_index] = mat2inSpike.poisson(patches(:,:,i),param);
                    end
                    times   = [times   imT_times+c + rand(size(imT_times))/40];
                    indices = [indices imT_index];
                    c = c + 1;
                end
            end
        end
        % #################################################################
        % make some artificial patches
        function [times, indices] = artificial_patches(L, repeat)
            % base patches for repeat
            patches = zeros(L,L,2);
            [~,~,f] = size(patches);
            patches(1:L^2)       = [ones(1,L) zeros(1,L^2-L)];
            patches(L^2+1:2*L^2) = [zeros(1,L^2-L) ones(1,L)];
            
            c = 0;
            times   = [];
            indices = [];
            param.binNumber     = 1000;
            param.maxPercSpikes = 0.5; % max number of spiking bins for max value
            param.tau           = 0.2; 
            for r = 1:repeat
                for i=1:f
                    [imT_times, imT_index] = mat2inSpike.deterministic(patches(:,:,i),param);
                    times   = [times   imT_times+c];
                    indices = [indices imT_index];
                    c = c + 1;
                end
            end
        end
        % #################################################################
        % get spikes from picture patches. The indices are indices of the
        % patch pixels
        function [times, indices] = pic_patches(patchSize)
            im2D  = imread('cameraman.tif');
            [imT] = Dat2Dconvolver.get_dat2Dconv(im2D, 3,'centerON');
            % reshape to get patches
            [r,~] = size(imT);

            pShift= patchSize;
            reduceSize = mod(r(1)-patchSize,pShift);
            imT(end-reduceSize+1:end,:) = [];
            imT(:,end-reduceSize+1:end) = [];
            [r,~] = size(imT);
            % generate patches from picture
            patches = zeros(patchSize,patchSize,(r/patchSize)^2);
            patchNum = r/patchSize;
            for patchRow = 1:patchNum
                patches(:,:,1+patchNum*(patchRow-1):patchNum*patchRow) = ...
                    reshape(imT(1+(patchRow-1)*patchSize:patchSize*patchRow,:),...
                            patchSize,patchSize,[]);
            end

            % make spike signals from patches
            pixelIndex = 1:numel(patches(:,:,1));
            times = zeros(1,pixelIndex(end)*length(patches));
            indices = zeros(1,pixelIndex(end)*length(patches));
            pStart = 1;
            for i=1:length(patches)
                [imT_times, imT_index] = mat2inSpike.deterministic(patches(:,:,i));
                L = length(imT_times);
                times(pStart:pStart+L-1) = imT_times + 1*(i-1);
                indices(pStart:pStart+L-1) = imT_index;
                pStart = pStart+L;
            end
        end       
    end
    
end

