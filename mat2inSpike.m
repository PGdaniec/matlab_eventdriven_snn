classdef mat2inSpike
    %T2SPIKE Convert value matrix input to spike time
    % For matrix time data:
    % data2D(Rvalue, Cvalue, timeIdx)
    %                   [1 4 7]                              
    % data2D(:,:,1)   = [2 5 8]    <- time index 1
    %                   [3 6 9]   
    %                   [10 13 16]                              
    % data2D(:,:,2)   = [11 14 17] <- time index 2
    %                   [12 15 18] 
    %
    % or single data matrices, where size(timeIdx) = 1 and therfore
    % size(data2D) = m x m
    %
    % resort by element as time-series
    %
    % linMat = [1 10]
    %          [2 11]
    %          [....]
    %          [9 18]
    %
    % time_s   = {t1_elem1 t2_elem1 ... tend_elem1}
    %            {t1_elem2 ............ tend_elem2} <- may have different L
    properties
    end
    
    methods (Static)
        % Transform continuous value to spike time poisson train ##########
        % Get ON / OFF spikes for increase / decrease of value over time
        function [timesONspike_s, timesOFFspike_s, indexONspikes, indexOFFspikes] = ...
                timeContrast(data2D_frames, Ts, treshold)
            % Will generate ON spike for relative increase of value over time
            % Will generate OFF spike for relative decrease of value over time
            % make linear index of 2DData   
            [r,c,f] = size(data2D_frames);
            if(f == 1)
               data2D_frames = reshape(data2D_frames,1,1,numel(data2D_frames)); 
               [r,c,f] = size(data2D_frames);
            end
            % linMat(1,:) contains diff to last, which is per def. 0.
            % Diff current to last
            if verLessThan('matlab','9.3')
                linMat = zeros(r*c,f);
                for i=1:f
                    frame = data2D_frames(:,:,i);
                    linMat(:,i) = frame(1:end)';
                end
            else
                linMat = reshape(data2D_frames(1:end)',r*c,f);
            end
                
            vMin = min(min(min(linMat)));
            vMax = max(max(max(linMat)));
            offsetForLog = 1; % avoids log(0)
            linMat = (linMat-vMin)./(vMax-vMin)+offsetForLog;
            
            logSignal       = log(linMat+offsetForLog);
            % signals with pos. spikes
            times_s         = (0:f-1)*Ts;
            
            timesONspike_s  = nan(1,length(times_s)*r*c);
            timesOFFspike_s = nan(1,length(times_s)*r*c); 
            indexONspikes   = nan(1,length(times_s)*r*c);
            indexOFFspikes  = nan(1,length(times_s)*r*c); 
            ONidx  = 1;
            OFFidx = 1;
            for i=1:r*c
                % pos spikes
                logSignal_row_sig = logSignal(i,:);
                
                lastSP_idx = 1;
                for j = 2:f
                    if((logSignal_row_sig(j)-logSignal_row_sig(lastSP_idx))/...
                       logSignal_row_sig(lastSP_idx) > treshold)
                        timesONspike_s(ONidx) = times_s(j);
                        indexONspikes(ONidx)  = i;
                        ONidx                 = ONidx + 1;
                        lastSP_idx            = j;
                    elseif((logSignal_row_sig(j)-logSignal_row_sig(lastSP_idx))/...
                           logSignal_row_sig(lastSP_idx) < -treshold)
                        timesOFFspike_s(OFFidx) = times_s(j);
                        indexOFFspikes(OFFidx)  = i;
                        OFFidx                  = OFFidx + 1;
                        lastSP_idx              = j;
                    end
                end
            end
            timesONspike_s (isnan(timesONspike_s )) = [];
            timesOFFspike_s(isnan(timesOFFspike_s)) = []; 
            indexONspikes  (isnan(indexONspikes  )) = []; 
            indexOFFspikes (isnan(indexOFFspikes )) = []; 
            [timesONspike_s,  tON_sort_idx]  = sort(timesONspike_s );
            [timesOFFspike_s, tOFF_sort_idx] = sort(timesOFFspike_s);
            indexONspikes                    = indexONspikes (tON_sort_idx);
            indexOFFspikes                   = indexOFFspikes(tOFF_sort_idx);
        end
        
        % Transform continuous value to spike time poisson train ##########
        % Get poisson spike train with high probability for high value
        function [times_s, index] = ...
                poisson(data2D, param)
            % Will map values to time range 0 ... 1 with binNumber of dT
            % bins. Max value will produce on average max number of  spikes 
            % for the given bins. Min value will produce no spikes.
            if(nargin < 2)
                param.binNumber     = 1000;
                param.maxPercSpikes = 0.5; % max number of spiking bins for max value
            end
            
            spikeProb_scaled_with_time = ones(1,param.binNumber);
            [times_s, index] = mat2inSpike.transform_val2RandSpikeTrain(...
                                 data2D, spikeProb_scaled_with_time, param); 
        end
        % Transform continuous value to exp decay spike  train ############
        % Get poisson spike train with high probability for high value
        function [times_s, index] = ...
                expDec(data2D, param)
            % Will map values to time range 0 ... 1 with binNumber of dT
            % bins. Max value will produce on average max number of  spikes 
            % for the given bins. Min value will produce no spikes.
            if(nargin < 2)
                param.binNumber     = 1000;
                param.maxPercSpikes = 0.5; % max number of spiking bins for max value
                param.tau           = 0.2; 
            end
            
            minT   = 0;
            maxT   = 1;
            
            spikeProb        = exp(-linspace(minT,maxT,param.binNumber)/param.tau);
            [times_s, index] = mat2inSpike.transform_val2RandSpikeTrain(...
                                 data2D, spikeProb, param); 
        end
        
        % Transform continuous value to spike time ########################
        % get first poisson spike per data value
        function [times_s, index] = ...
                poisson_first(data2D, param)
            if(nargin < 2)
                param.binNumber     = 1000;
                param.maxPercSpikes = 0.5; % max number of spiking bins for max value
            end
            param.take_first = true;
            [times_s, index] = mat2inSpike.poisson(data2D, param); 
        end
        % Transform continuous determ. spike time #########################
        % get deterministic value hither comes first, 0 does not come
        function [times_s, index] = ...
                deterministic(data2D, param)
            if(nargin < 2)
                param.binNumber     = 1000;
                param.maxPercSpikes = 0.5; % max number of spiking bins for max value
            end
            param.take_first    = true;
            param.deterministic = true;
            
            spikeProb_scaled_with_time = linspace(1,0,param.binNumber);
            [times_s, index] = mat2inSpike.transform_val2RandSpikeTrain(...
                                 data2D, spikeProb_scaled_with_time, param); 
        end
    end
    methods(Static, Access = private)
        % helper funcs ####################################################
        function [times_s, index] = transform_val2RandSpikeTrain(data2D, spikeProb_scaled_with_time, param)
            % Will map values to time range 0 ... 1 with binNumber of dT
            % bins. Max value will produce on average max number of  spikes 
            % for the given bins. Min value will produce no spikes.
            binNumber     = param.binNumber;
            maxPercSpikes = param.maxPercSpikes; 
            if(isfield(param,'take_first') && param.take_first == true)
                take_first = true;
            else
                take_first = false;
            end
            if(isfield(param,'deterministic') && param.deterministic == true)
                deterministic = true;
            else
                deterministic = false;
            end
             
            [data2D]     = mat2inSpike.rescale01(data2D, param);
            [r,c]        = size(data2D);
            matValues    = data2D(1:end);
            % high value leads to high spike prop per bin
            allTSamples_s = linspace(0, 1, binNumber);
            times_s = nan(1,r*c*binNumber);
            index   = nan(1,r*c*binNumber);
            nextStartidx = 1;
            
            if(deterministic)
                rand_div_spikeProb = kron(ones(r*c,1), spikeProb_scaled_with_time);
            else
                if verLessThan('matlab','9.3')
                    % If no cude, use:
                    rand_div_spikeProb = ...
                        rand(r*c, param.binNumber) ./ ...
                        spikeProb_scaled_with_time;  
                else
                    parallel.gpu.rng('shuffle', 'Philox4x32-10');
                    rand_div_spikeProb = ...
                        gpuArray.rand([r*c param.binNumber],'single') ./ ...
                        spikeProb_scaled_with_time;   
                end
            end
            
            for i=1:r*c  
                t = allTSamples_s(matValues(i)*maxPercSpikes >= rand_div_spikeProb(i,:) & matValues(i) ~= zeros(size(rand_div_spikeProb(i,:))));
                if(take_first)
                   t = min(t);
                end
                L = length(t);
                if(L ~= 0)
                    times_s(nextStartidx:nextStartidx+L-1) = t;
                    index(nextStartidx:nextStartidx+L-1)   = i;
                end
                nextStartidx            = nextStartidx + L; 
            end
            index(isnan(index))        = [];
            times_s(isnan(times_s))    = [];
            [times_s, sortIdx]         = sort(times_s);
            index                      = index(sortIdx);            
        end
        % #################################################################
        function [minV, maxV] = getMinMaxVal(data2D, param)
            if(isfield(param,'minmax'))
                minV = param.minV;
                maxV = param.maxV;
            else
                minV = min(min(data2D));
                maxV = max(max(data2D));
            end
        end
        function [resc_data2D] = rescale01(data2D, param)
            [minV, maxV] = mat2inSpike.getMinMaxVal(data2D, param);
            resc_data2D  = (data2D-minV)./(maxV-minV);
        end
    end
    
end

