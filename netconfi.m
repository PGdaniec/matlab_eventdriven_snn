classdef netconfi
    %NETCONFI Summary of this class goes here
    %   Detailed explanation goes here
   
    
    methods (Static)    
        function inp = getParams_learn_W_Win()
            inp = netconfi.getParams_exc80inh20_strongSyn();
            inp.use_synW_adaption    = true;    % use algorithm for syn adaption
            inp.use_synWin_adaption  = true;
            inp.use_ext_input        = true;
        end
        
        function inp = getParams_learn_W()
            inp = netconfi.getParams_exc80inh20_strongSyn();
            inp.use_synW_adaption = true;    % use algorithm for syn adaption
            inp.use_ext_input     = true;
        end
        
        function inp = getParams_exc80inh20_strongSyn()
            inp = netconfi.getParams_exc80inh20();
            inp.j0              = 1;         % synaptic base weight
            inp.I0              = 1e-2;      % excitation current
        end
        
        function inp = getParams_exc80inh20()
            % connection topo ###########
            inp.nspike          = 10^5;      % number of spikes
            s = 20;
            inp.ne              = s^2 * 4;   % number of exc neurons
            inp.ni              = s^2;       % number of inh neurons
            inp.ke              = 50;        % number of exc out synapses
            inp.ki              = 50;        % number of inh out synapses
            inp.netSizeE_in     = [90 90];
            inp.netSizeI_in     = [45 45];
                         
            inp.tau_exc         = 0.02;      % time const excitory neuron
            inp.tau_inh         = 0.002;     % time constant inhibitory neuron         
          
            inp.refactory_time    = 0.002;   % after a spike, neuron is not responsive
            inp.use_refactoryTime = true;
            
            % initial condition seed ####
            inp.seedic              = 1;     

            % inputs ####################
            inp.use_synWin_adaption = false; % use input weight adaption
            inp.use_ext_input       = false; % use external input
            
            inp.ext_input_time    = [];      % input time of spike                          [t1 t2 ... tend]
            inp.ext_input_nID     = {};      % target neuron for input spike: {[1 2 n] ... [any number of neurons to target]}
            inp.ext_input_syn     = [];      % target neuron's input synapse of input spike [s1 s2 ... send]
            inp.kin               = 0;       % number of input synapses per neuron
            
            % internal weights ##########
            inp.use_synW_adaption = false;   % use algorithm for syn adaption   
            change              = 1e-1;
            inp.synW_inc_scale  = change;  % increase scale factor for synaptic adaption
            inp.synW_dec_scale  = change;  % decrease scale factor for synaptic adaption
            % limit in update formula, voltage range is [-1 0],
            % exc.limit = 1 -> one spike to fire neuron at reset state
            inp.synW_inc_limit_inh  = 0.5;   % max limit for inhibition
            inp.synW_inc_limit_exc  = 0.1;   % max limit for excitation
            % weights
            inp.j0              = 1e-4;      % synaptic base weight
            inp.I0              = 1e-6;      % excitation current
            inp.synW_exc_f      = -0.2;      % scaling of syn weight
            inp.synW_ext_f      = -0.6;      % scaling of external input synaps
            inp.synW_inh_f      =  1;        % scaling of syn weight 
            
            % connection matrix #########
            tb = TopoBuilder.mpic_ei2ei(inp.ne,inp.ke,inp.ni,inp.ki,[],[]);
            inp.postidx_precomp = tb.postidxMat;
        end             
    end
end