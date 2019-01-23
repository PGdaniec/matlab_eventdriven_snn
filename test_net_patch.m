[tPatch, sPatch] = get_Input.artificial_rand_patches(5000,"poisson");

%% each patch pixel maps to one neuron in network
tic
inp                   = netconfi.getParams_learn_W_Win();

inp.kin               = 25;
inp.ne                = 15^2;
inp.ni                = 15^2;     
inp.ke                = 9;
inp.ki                = 16;
inp.nspike            = 10^6;


inp.ext_input_time    = tPatch/10+1;
inp.ext_input_nID     = cell(size(sPatch));
inp.ext_input_nID(:)  = {1:inp.ne};   
inp.ext_input_syn     = sPatch;
inp.use_synW_adaption    = true;
inp.use_synWin_adaption  = true;


inp.refactory_time    = 0.002;
inp.j0                = 1e-3;
inp.I0                = 1e-5;

inp.synW_exc_f      = -0.1;
inp.synW_ext_f      = -0.1 + rand(inp.kin,inp.ne+inp.ni)/20-1/40;
inp.synW_inh_f      =  abs(inp.synW_exc_f) * 9/16;        
inp.synW_inc_limit_inh  = 0.1;   
inp.synW_inc_limit_exc  = 0.1 * 9/16;   
            
change              = 0.1;
inp.synW_inc_scale  = change;
inp.synW_dec_scale  = change;

tb = TopoBuilder.mpic_localCompetition(inp.ne,inp.ni,[sqrt(inp.ne) sqrt(inp.ne)],[sqrt(inp.ni) sqrt(inp.ni)]);
inp.postidx_precomp = tb.postidxMat;

[rate,sidx,stimes,debugValue] = spikingnet_del_add(inp);
NetOutPresenter.show_infos(rate,sidx,stimes,toc,inp,debugValue)

r = 1;
c = 1;
sizeNe = [sqrt(inp.ne) sqrt(inp.ne)];
synWin2D = zeros(sizeNe(1)*5,sizeNe(2)*5);
for i=1:prod(sizeNe)
    synWin2D(1+(r-1)*5:r*5,1+(c-1)*5:c*5) = reshape(debugValue.synWin_mat(:,i),[5 5]);
    if(mod(i,sizeNe(1)) == 0)
        r = 1;
        c = c+1;
    else
        r = r+1;
    end
end

imagesc(synWin2D)
%%
np = NetOutPresenter;
np = np.addSubplot(stimes,sidx,'.');
np = np.addSubplot(inp.ext_input_time,sPatch,'.r',1);
%%
close all
np.plotMovie(1,3)

%%
np.showStats(inp,sidx,debugValue);