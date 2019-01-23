%% test segment mapper 2D to 2D TODO
datS                 = 4;
data2D               = zeros(datS,datS);
data2D(1:datS^2)     = 1:datS^2;
segment_size_data2D  = datS-2;
segment_shift_data2D = 2;

segment_size_neuron  = segment_size_data2D;
segment_shift_neuron = segment_size_data2D;

[data2DIdx_projectsTo_neuronIdx, neuronIdxSize] = ...
mapMatSegments.one2one( size(data2D), ...
                        segment_size_data2D, segment_shift_data2D, ...
                        segment_size_neuron, segment_shift_neuron);
      
%% test segment mapper 1 to many TODO
datS                 = 1;
data2D               = zeros(datS,datS);
data2D(1:datS^2)     = 1:datS^2;
segment_size_neuron      = 1;
segment_shift_neuron     = 1;
num_neuronSegment_shifts = 1;
[data2DIdx_projectsTo_neuronIdx, neuronIdxSize] = ...
            mapMatSegments.one2many(...
            segment_size_neuron, segment_shift_neuron, num_neuronSegment_shifts);
        
%% test segment mapper 1 to many TODO
segment_size_neuron      = 3;
segment_shift_neuron     = 2;
num_neuronSegment_shifts = 3;

[data2DIdx_projectsTo_neuronIdx, neuronIdxSize] = ...
mapMatSegments.one2many(segment_size_neuron, ...
                                segment_shift_neuron, ...
                                num_neuronSegment_shifts);