function [rate,sidx,stimes,debugValue] = spikingnet(inp)
    % features:
    % - inputs
    % - refactory time
    % - weight adaption
    % - back-shift technique
    
    debug_max_deltaPhi          = 0;
    debug_min_deltaPhi          = 1;
    debug_num_INC_synAdaptions  = 0;
    debug_num_DEC_synAdaptions  = 0;
    debug_num_refactoryUsed     = 0;
    
    n         = inp.ne + inp.ni;
    ne        = inp.ne;
    ni        = inp.ni;
    nspike    = inp.nspike;
    ke        = inp.ke;
    ki        = inp.ke;
    j0        = inp.j0;
    tau_exc   = inp.tau_exc;
    tau_inh   = inp.tau_inh;
    I0        = inp.I0;
    
    use_ext_input     = inp.use_ext_input;
    if(use_ext_input)   
        ext_input_syn  = inp.ext_input_syn;      % input over this synapse
    end
    
    ext_input_time = inp.ext_input_time;     % input at this time
    ext_input_nID  = inp.ext_input_nID;      % input to this neuron
    [ei_r,~] = size(ext_input_nID{1});
    if(ei_r > 1)
        error('Cells in inp.ext_input_nID must be row-vectors!');
    end
    
    ext_input_k    = 1;
    ext_input_max  = length(ext_input_time);
    kin            = inp.kin;                % input synapse number
     
    use_synW_adaption = inp.use_synW_adaption;
    use_refactoryTime = inp.use_refactoryTime;
    % additive change
    synW_inc_scale    = inp.synW_inc_scale; 
    synW_dec_scale    = inp.synW_dec_scale;
    assert(synW_inc_scale >= 0,'synW_inc_scale must be >= 0');
    assert(synW_dec_scale >= 0,'synW_dec_scale must be >= 0');
    
    synW_inc_limit_inh    = inp.synW_inc_limit_inh;
    synW_inc_limit_exc    = inp.synW_inc_limit_exc;
    use_synWin_adaption   = inp.use_synWin_adaption;
    
    seedic    = inp.seedic;
    % forward phase shift to simplify implementation
    
    % parameter for PTC curve, shortcut   
    if(use_ext_input)
        Iext_exc    = sqrt(max(ki,ke)+kin)*I0;
    else
        Iext_exc    = sqrt(max(ki,ke))*I0;
    end
    
    T_free_exc  = tau_exc*log(1+1/Iext_exc);   % free period seperate for excitation
    % Due to simplicity T_free_inh is the same as for excitory neurons.
    % Actually that means they get a lower external current than the
    % excitory ones, as the time constant is smaller. The calculated
    % Iext_inh is:
    %     T_free_inh = T_free_exc
    %     T_free_inh / tau_inh = ln(1 + 1/Iext_inh) 
    %     -> Iext_inh = 1 / (exp(T_free_inh/tau_inh) - 1)
    
    % T_free_inh  = tau_inh*log(1+1/Iext); % free period seperate for inhibition
    T_free      = T_free_exc;              % free period
    T_free_inh  = T_free_exc;
    Iext_inh = 1 / (exp(T_free_inh/tau_inh) - 1);
    
    % synaptic weight, positive is excitory
    if(use_ext_input)
        synW_val    = -j0/sqrt(max(ki,ke)+kin);
    else
        synW_val    = -j0/sqrt(max(ki,ke));
    end
           
    synW_exc_f  = inp.synW_exc_f;          % scaling of syn weight
    synW_inh_f  = inp.synW_inh_f;          % scaling of syn weight
    synW_ext_f  = inp.synW_ext_f;          % external syn scaling
    
    % generate connections
    postidx_precomp = inp.postidx_precomp;

    % contains outputSynapse x fromNeuronIdx, 
    % negative entries for inh. and pos. for exc.
    synW_mat = [ones(ke,ne)*synW_val.*synW_exc_f ...
                ones(ki,ni)*synW_val.*synW_inh_f]; 
    debugValue.synW_mat_start = synW_mat;               

    % inputSynapse x toNeuronIdx
    synWin_mat = ones(kin,n)*-j0/sqrt(kin+max(ki,ke)).*synW_ext_f;
    debugValue.synWin_mat_start = synWin_mat;
    
    V_R         = -1;                      % voltage reset
    V_T         =  0;                      % voltage treshold
    
    phi_R       = 0.;                      % reset    
    phi_T       = 1.;                      % threshold
    
    rng(seedic);
    % initialize/preallocate spike raster and vector of receiving neurons
    spikeidx              = zeros(nspike,1); % initialize time
    spikephis             = zeros(nspike,1); % spike raster
    % 
    if(use_synW_adaption || use_synWin_adaption)                         
    % col:  idx of receiving neuron
    % row1: phishift of sender neuron 
    % row2: idx of sender neuron, that spiked
    % row3: idx of sender neurons synapse, that spiked
        lastSpikeRX = zeros(3,n);   
    end  
    if(use_synWin_adaption)
        % col:  idx of receiving neuron
        % row1: phishift of sender syn ? ToDo
        % row2: idx of sender synapse, that spiked
        lastInputSpikeRX = zeros(2,n);
    end    
    if(use_refactoryTime || use_synW_adaption || use_synWin_adaption)
        lastSpikeTimeOfNeuron = -ones(n,1);
    end

    % constants for selection
    atTime          = 1;
    fromNeuron      = 2;
    fromSyn         = 3;
    fromInputSyn    = 2;
    phishift        = 0;
    % time equivalient phi values
    deltaT_synW_change   = 0.05;
    deltaPhi_synW_change = abs(phi_T-phi_R)*deltaT_synW_change/T_free;
    if(use_refactoryTime)
        refactory_phi        = inp.refactory_time/T_free;
    end
    % generate heap for neurons
    phi = rand(n,1);

    % main loop ###########################################################
    % #####################################################################
    for s = 1 : nspike
        
        % find next spiking neuron
        j = find(phi == max(phi),1,'first');
        phimax = phi(j);

        % calculate next spike time, 
        % shift global treshold instead of all phis
        dphi = phi_T - phimax - phishift;  
        
        % remember forward shift
        phishift = phishift + dphi; 
        
        if(use_ext_input && ...
           ext_input_k     <= ext_input_max && ...
           phishift*T_free >  ext_input_time(ext_input_k) ...
          )
            % use external input time
            phishift = ext_input_time(ext_input_k) / T_free;
            % use excitory input
            ptc_fun(ext_input_nID{ext_input_k}, ext_input_syn(ext_input_k), true); 
            
            % next external data set
            ext_input_k   = ext_input_k + 1;
        else
            % use normal internal spiking
            postidx = postidx_precomp(j,:);
            
            % evaluate phase transition curve
            ptc_fun(postidx,[] , false);

            % reset spiking neuron, use global shift for proper reset
            % phi(j) contains -phase of last spike for neuron j.
            % Later if phi(j) is stored, the elapsed time since last spike
            % is phi(j)+current_phishift
            phi(j)             = phi_R - phishift;
        end
        
        % store spike raster
        spikephis(s)  = phishift;   % store spiketimes
        spikeidx(s)   = j;          % store spiking neuron
    end
    
    rate    = nspike / (phishift*T_free) / n;
    sidx    = spikeidx;
    stimes  = spikephis*T_free;       
    
    debugValue.debug_max_deltaPhi         = debug_max_deltaPhi;
    debugValue.debug_min_deltaPhi         = debug_min_deltaPhi;
    debugValue.debug_num_INC_synAdaptions = debug_num_INC_synAdaptions;
    debugValue.debug_num_DEC_synAdaptions = debug_num_DEC_synAdaptions;
    debugValue.debug_num_refactoryUsed    = debug_num_refactoryUsed;
    debugValue.synW_mat_end = synW_mat;
    debugValue.synWin_mat   = synWin_mat;

    debugValue.postidx_precomp = postidx_precomp;

    % define phase transition curve ####################################
    function ptc_fun(postidx, extInputSyn_idx, isExtInput) 
        syn_idx = 1;
        for i_l = postidx
            if(i_l == 0)
                % rest is filled up with zeros, if i/e have different num
                % of connections
                break;
            end
            
            % slow excitory and fast inhib.
            if(i_l > ne)
                tau     = tau_inh;
                Iext    = Iext_inh;
            else
                tau     = tau_exc;
                Iext    = Iext_exc;
            end

            % remember values
            if(use_synW_adaption && ~isExtInput)
                lastSpikeRX(atTime,i_l)     = phishift;
                lastSpikeRX(fromNeuron,i_l) = j;
                lastSpikeRX(fromSyn,i_l)    = syn_idx;
                lastSpikeTimeOfNeuron(j)    = phishift;
            end
            if(use_synWin_adaption && isExtInput)
                % row:  idx of receiving neuron
                % col1: phishift of sender syn
                % col2: idx of sender synapse, that spiked
                lastInputSpikeRX(atTime,i_l)      = phishift;
                lastInputSpikeRX(fromInputSyn,i_l)= extInputSyn_idx;
            end

            % if is in refactory period
            if(use_refactoryTime && phishift-lastSpikeTimeOfNeuron(i_l) < refactory_phi)
                % do not accept synapic input spikes, phi(i_l) shall stay
                % the same until phi(i_l) == refactory_phi;
                debug_num_refactoryUsed = debug_num_refactoryUsed + 1;
            else
                last_phi = phi(i_l);
                % for positive input avoid V > 0, this would do a spike
                % restore proper phi, as phi is shifted by global shift
                % range is [-1 0]
                % synW < 0 -> inhibition, synW > 0 -> excitation
                
                % do not use weight matrix for external inputs
                if(~isExtInput)
                    V_phi_afterSpike = ...
                        max(min(Iext - (Iext + 1)* ...
                        exp(-T_free/tau* (phi(i_l)+phishift) ) + ...
                        synW_mat(syn_idx,j), V_T),V_R);
                else
                    % depends on j index, as here we have output synapses of j
                    % synW > 0 is excitory, <0 is inhibitory 
                    V_phi_afterSpike = ...
                        max(min(Iext - (Iext + 1)* ...
                        exp(-T_free/tau* (phi(i_l)+phishift) ) + ...
                        synWin_mat(extInputSyn_idx,i_l), V_T),V_R);
                end
                % range is [0 1]
                phi(i_l) =  tau/T_free * ...
                    log((Iext - V_R) / ...
                        (Iext - V_phi_afterSpike)...
                       )-phishift; % again apply global shift to store right value
                debug_max_deltaPhi = max(debug_max_deltaPhi,phi(i_l)-last_phi);
                debug_min_deltaPhi = min(debug_min_deltaPhi,phi(i_l)-last_phi);
            end
            % TODO check for correctness
            if(use_synW_adaption && ~isExtInput)
                % adapt output syn of neuron j, DEcrease weight ###########
                % there was one and it is less thena 0.1 phi ago
                % when to change
                % lastSpikeTimeOfNeuron(i_l) = phishift at last spike time of i_l
                % lastSpikeTimeOfNeuron >= 0 is a valid value
                if(lastSpikeTimeOfNeuron(i_l)>=0 && ...
                        phishift-lastSpikeTimeOfNeuron(i_l) < deltaPhi_synW_change)
                    
                    % % additive change, avoid signflip for inh/exc     
                    if(j>ne)
                        % for inh. shall be negative or 0
                        synW_mat(syn_idx,j) = min(  synW_mat(syn_idx,j) ...
                                                  + synW_dec_scale,0);    
                    else
                        % for exc. shall be positive or 0
                        synW_mat(syn_idx,j) = max(  synW_mat(syn_idx,j) ...
                                                  - synW_dec_scale,0);    
                    end
                    debug_num_DEC_synAdaptions  = debug_num_DEC_synAdaptions + 1;
                end                
            end
            if(use_synWin_adaption && isExtInput)
                % adapt input syn of neuron j, DEcrease weight ############
                % if target spiked before source -> reduce weight
                % the target is neuron i_l, the source is the input synapse
                % extInputSyn_idx. Use allways exc. weights (>= 0)!
                if(lastSpikeTimeOfNeuron(i_l)>=0 && ...
                        phishift-lastSpikeTimeOfNeuron(i_l) < deltaPhi_synW_change)
                    % synapse to neuron i_l
                    if(sign(synWin_mat(extInputSyn_idx,i_l)) == -1)
                        % for inh. shall be negative or 0
                        synWin_mat(extInputSyn_idx,i_l) = ...
                        	min(  synWin_mat(extInputSyn_idx,i_l) ...
                            	+ synW_dec_scale,0);
                    else
                        % for exc. shall be positive or 0
                        synWin_mat(extInputSyn_idx,i_l) = ...
                        	max(  synWin_mat(extInputSyn_idx,i_l) ...
                            	- synW_dec_scale,0);
                    end
                end
            end

            % increment synW index
            syn_idx = syn_idx+1;
        end
        if(use_synW_adaption && ~isExtInput)
            % adapt syn that target neuron j, INcrease weight #############
            % find presynaptic neuron of j
            n_idx   = lastSpikeRX(fromNeuron,j);
            n_syn_idx = lastSpikeRX(fromSyn,j);
            % there was one and it is less then a 0.1 phi ago
            if(n_idx ~= 0 && phishift-lastSpikeRX(atTime,j) < deltaPhi_synW_change)
                % -synW_inc_limit_inh <--- synW_mat ---> +synW_inc_limit_exc              
                if(j>ne)
                    % for inh. shall be negative or 0
                    synW_mat(n_syn_idx,n_idx) = max(  synW_mat(n_syn_idx,n_idx) ...
                                                    - synW_inc_scale,-synW_inc_limit_inh);    
                else
                    % for exc. shall be positive or 0
                    synW_mat(n_syn_idx,n_idx) = min(  synW_mat(n_syn_idx,n_idx) ...
                                                    + synW_inc_scale,+synW_inc_limit_exc);    
                end         
                debug_num_INC_synAdaptions  = debug_num_INC_synAdaptions + 1;
            end
        end
        if(use_synWin_adaption && isExtInput)
            % adapt input syn of neuron j, INcrease weight ############
        	% TODO
            input_syn_idx   = lastInputSpikeRX(fromInputSyn,j);
            if(input_syn_idx ~= 0 && phishift-lastInputSpikeRX(atTime,j) < deltaPhi_synW_change)
                % -synW_inc_limit_inh <--- synWin_mat ---> +synW_inc_limit_exc    
                if(sign(synWin_mat(input_syn_idx,j)) == -1)
                        % for inh. shall be negative or 0
                    synWin_mat(input_syn_idx,j) = max(  synWin_mat(input_syn_idx,j) ...
                                                      - synW_inc_scale,-synW_inc_limit_inh);    
                else
                    % for exc. shall be positive or 0
                    synWin_mat(input_syn_idx,j) = min(  synWin_mat(input_syn_idx,j) ...
                                                      + synW_inc_scale,+synW_inc_limit_exc);    
                end       
            end
        end
    end    
end
