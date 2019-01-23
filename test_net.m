clear inp;
inp                 = netconfi.getParams_learn_W();

% #######################################################################
s = 20;
inp.nspike          = 8e5;
inp.ne              = s^2 * 4;
inp.ni              = s^2;
inp.ke              = 50;
inp.ki              = 50;
inp.netSizeE_in     = [90 90];
inp.netSizeI_in     = [45 45];
tb = TopoBuilder.mpic_ei2ei(inp.ne,inp.ke,inp.ni,inp.ki,[],[]);
inp.postidx_precomp = tb.postidxMat;
inp.use_synWin_adaption = true;
inp.use_ext_input       = true;
inp.use_synW_adaption   = true;
% ######################################################################

inputType = {'steps'};
inpSelect = inputType{1};
switch inpSelect
    case 'steps' % clear pattern is seen in adaption
        make_input_free_period = 0;
        inputSteps            = 2*80+make_input_free_period;
        input_time            = (1:inputSteps)/10;
        if(make_input_free_period)
            input_time(end)   = inputSteps/10*2;
        end
        endN                  = floor((inp.ne+inp.ni)/20);
        v                     = 1:endN;
        inputMat              = mod(0:inputSteps*endN-1,endN*8)+make_input_free_period;
        inputMatForPlot       = reshape(inputMat, endN(end),inputSteps);
        input_nID             = mat2cell(inputMat, ...
                                         1, endN*ones(1,inputSteps));
        inp.kin               = 1;
end

use_time_sig = 0;
if(use_time_sig)
    T = 0.2;
    tVec = 0:0.001:input_time(end);
    sinSig = abs(sin(2*pi/T * (2*pi/ (T/4) - tVec)));
    nid = 200;
    clock_logical = rand(nid,length(tVec)) < sinSig/4;
    clock_time = [];
    clock_idx  = [];
    for i=1:200
        clock_time = [clock_time tVec(clock_logical(i,:) == 1)];
        clock_idx  = [clock_idx i*ones(1,sum(clock_logical(i,:) == 1))+800];
    end
    [clock_time, sortI] = sort(clock_time);
    clock_idx           = clock_idx(sortI);

    [inp.ext_input_time, sortI] = sort([input_time clock_time]);
    inp.ext_input_nID           = [input_nID mat2cell(clock_idx,1,ones(size(clock_idx)))];
    inp.ext_input_nID           = inp.ext_input_nID(sortI);
else
    inp.ext_input_time = sort(input_time);
    inp.ext_input_nID  = input_nID;
end
inp.ext_input_syn           = ones(size(inp.ext_input_nID));

tic;
%profile on
[rate,sidx,stimes,debugValue] = spikingnet_del_add(inp);
%profile off
%profile viewer
NetOutPresenter.show_infos(rate,sidx,stimes,toc,inp,debugValue);

%% plot results
figure
np = NetOutPresenter;
np = np.addSubplot(stimes,sidx,'.');
if(use_time_sig)
np = np.addSubplot(clock_time,clock_idx,'.r',1);
end
np = np.addSubplot(input_time,inputMatForPlot','.r',1);
%% plot movie
np.plotMovie(0.8,3);