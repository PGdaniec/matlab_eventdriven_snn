classdef mapMatSegments
    %mapMatSegments Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        % Get mapping of 2D-data to neuron index ##########################
        % Structure of
        % data2DIdx_2_neuronIdx is a cell vector
        %   { [target indices of neurons for first lin. index of data2D matrix] } 
        %   { ...                                                               }
        %   { ...                                                               }
        %   { [target indices of neurons for last lin.  index of data2D matrix] }
        % Each index from the same segment in data2D projects to every index in 
        % the target neuron segment. For example data2D has a segment
        % consisting of index 1,2,3,4 and the target segment of neurons
        % consists of indices 5,6,7,8 then there will be the projection
        %                         { [5 6 7 8 ...] }
        % data2DIdx_2_neuronIdx = { [5 6 7 8 ...] }
        %                         { [5 6 7 8 ...] }
        %                         { [5 6 7 8 ...] }
        %                         { ...           }
        function [data2DIdx_2_neuronIdx, neuronIdxSize] = ...
                one2one(data2D_size, ...
                segment_size_data2D, segment_shift_data2D, ...
                segment_size_neuron, segment_shift_neuron)
            % segment size and shift shall be symetric
            sSize    = segment_size_data2D;
            sShift   = segment_shift_data2D; 
            nsSize   = segment_size_neuron;  % segment size in 2D neuron formation 
            nsShift  = segment_shift_neuron;
            % number of neurons per segment
            nNeuronsPerSeg = nsSize^2;
            
            % size of whole data
            dR                    = data2D_size(1);
            dC                    = data2D_size(2);
            data2DIdx_2_neuronIdx = zeros(dR*dC,1);
            
            %          [s1 s4 s7]
            % data2D = [s2 s5 s8] segments can also overlap
            %          [s3 s6 s9]
            
            numSegments_row = (dR-sSize)/sShift + 1;
            if(rem(numSegments_row,1) ~= 0)
               error('(dR-sSize)/sShift is not whole number!'); 
            end
            numSegments_col = (dC-sSize)/sShift + 1;
            if(rem(numSegments_col,1) ~= 0)
               error('(dC-sSize)/sShift is not whole number!'); 
            end
            % Each segment from data2D projects to segment in neuron index
            numSegments_row_N = numSegments_row;
            numSegments_col_N = numSegments_col;
            
            numSegments_data2D = numSegments_row   * numSegments_col;
            numSegments_neuron = numSegments_row_N * numSegments_col_N;
            % calculate number of neurons
            nR              = (numSegments_row_N-1)*nsShift + nsSize;
            nC              = (numSegments_col_N-1)*nsShift + nsSize;
            neuronIdx       = (1:nR*nC)';
            neuronIdxSize   = [nR nC];
            
            for colShiftNum=0:numSegments_col-1
                
                % data column-indices of segment
                colIndices   = 1+colShiftNum*sShift:sSize+colShiftNum*sShift;
                N_colIndices = 1+colShiftNum*nsShift:nsSize+colShiftNum*nsShift;
                for rowShiftNum=0:numSegments_row-1
                    
                    % data row indices of segment
                    rowIndices   = 1+rowShiftNum*sShift:sSize+rowShiftNum*sShift;  
                    N_rowIndices = 1+rowShiftNum*nsShift:nsSize+rowShiftNum*nsShift;
                    % calculate neuron projection and store indices
                    neuronIndices = sub2ind([nR nC],...
                        kron(ones(1,length(N_colIndices)), N_rowIndices),...
                        kron(N_colIndices                , ones(1,length(N_rowIndices))) );
                    % calculate linear index

                    Idx2D = sub2ind(data2D_size,...
                        kron(ones(1,length(colIndices)), rowIndices),...
                        kron(colIndices                , ones(1,length(rowIndices))) );
                    minNumOfZeros = min(sum((data2DIdx_2_neuronIdx(Idx2D,:)==0),2));
                    % expand matrix, not enough free entries
                    if(nNeuronsPerSeg > minNumOfZeros)
                        data2DIdx_2_neuronIdx(:,end+1:end+nNeuronsPerSeg-minNumOfZeros) = zeros;
                    end
                    for i=1:length(Idx2D)
                        zeroInd = find(data2DIdx_2_neuronIdx(Idx2D(i),:)==0);
                        data2DIdx_2_neuronIdx(Idx2D(i),zeroInd(1:nNeuronsPerSeg)) = ...
                            neuronIdx(neuronIndices);  
                    end
                end
            end  
        end
        % Get mapping of 2D-data to neuron index ##########################
        function [data2DIdx_2_neuronIdx, ...
                  neuronIdxSize] = ...
            one2many(segment_size_neuron, ...
                     segment_shift_neuron, ...
                     num_neuronSegment_shifts)
            % segment size and shift shall be symetric
            nsSize      = segment_size_neuron;  % segment size in 2D neuron formation 
            nsShift     = segment_shift_neuron;
            num_nsShift = num_neuronSegment_shifts;
            % number of neurons per segment
            nNeuronsPerSeg = nsSize^2;

            data2DIdx_2_neuronIdx = [];

            %                      [s1 s4 s7]
            % scalar -> neuron2D = [s2 s5 s8] segments can also overlap
            %                      [s3 s6 s9]

            numSegments_row_N = num_nsShift + 1;
            numSegments_col_N = numSegments_row_N;
            neuron_L_row = num_nsShift * nsShift + nsSize;
            neuron_L_col = neuron_L_row;
            
            % data2D projects to each segment in neuron index
            numSegments_neuron = numSegments_row_N * numSegments_col_N;
            % calculate number of neurons
            nR              = neuron_L_row;
            nC              = neuron_L_col;
            neuronIdx       = (1:nR*nC)';
            neuronIdxSize   = [nR nC];

            for colShiftNum=0:numSegments_col_N-1

                % data column-indices of segment
                N_colIndices = 1+colShiftNum*nsShift:nsSize+colShiftNum*nsShift;
                for rowShiftNum=0:numSegments_row_N-1

                    % data row indices of segment
                    N_rowIndices = 1+rowShiftNum*nsShift:nsSize+rowShiftNum*nsShift;
                    % calculate neuron projection and store indices
                    neuronIndices = sub2ind([nR nC],...
                        kron(ones(1,length(N_colIndices)), N_rowIndices),...
                        kron(N_colIndices                , ones(1,length(N_rowIndices))) );
                    % calculate linear index
                    data2DIdx_2_neuronIdx = ...
                        [data2DIdx_2_neuronIdx neuronIdx(neuronIndices)']; 
                end
            end  
        end
        
        % Get number of segments from #####################################
        function numOfSeg = getNumOfSegments(vecL, segmentSize, segmentShift)
            numOfSeg = (vecL-segmentSize)/segmentShift + 1;
            if(rem(numOfSeg,1) ~= 0)
                error(['(vecL-segmentSize)/segmentShift = ' ...
                    num2str((vecL-segmentSize)/segmentShift) ...
                    ' is not a whole number!']);
            end
        end
        % Get vector length ###############################################
        function vecL = getVecL(numOfSeg, segmentSize, segmentShift)
            vecL = (numOfSeg - 1)*segmentShift + segmentSize;
        end
    end
    
    
end
