classdef TopoBuilder
    %SYNWBUILDER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        postidxMat;
    end
    properties (Access = private)
        ni;
        ki;
        ne;
        ke;        
        nMax_exc;    
        netSizeE;
        netSizeI;
    end
    
    methods(Static)
        % Make post index connections (mpic)
        % Make close neighbor excitation and far inhibition
        function tb = mpic_localCompetition(ne,ni,netSizeE_in,netSizeI_in)
            ndistExc = [1 1];
            ke   = 8*1 + 1;
            
            ndistInh = [2 2];
            ki   = 8*2;
            from = 1:ne;
            to   = from; 
            
            tb = TopoBuilder(ne,ke,ni,ki,netSizeE_in,netSizeI_in);  
            tb = tb.connect_e2e_fromToNeighbor(from,to,ndistExc);
            tb = tb.connect_e2i_fromTo        (from,1:ni);
            tb = tb.connect_i2e_fromToNeighbor(1:ni,from,ndistInh);
            
        end
        % Make random connections
        function tb = mpic_ei2ei(ne,ke,ni,ki,~,~)
           tb = TopoBuilder(ne,ke,ni,ki); 
           ke2e = round(ne/(ne+ni) * ke);
           tb = tb.connect_e2e_rand(ke2e);
           tb = tb.connect_e2i_rand(ke-ke2e);
           
           ki2i = round(ni/(ne+ni) * ki);
           tb = tb.connect_i2i_rand(ki2i);
           tb = tb.connect_i2e_rand(ki-ki2i);
        end
        function tb = mpic_e2ei_i2e(ne,ke,ni,ki,~,~)
           tb = TopoBuilder(ne,ke,ni,ki); 
           ke2e = round(ne/(ne+ni) * ke);
           tb = tb.connect_e2e_rand(ke2e);
           tb = tb.connect_e2i_rand(ke-ke2e);
           
           tb = tb.connect_i2e_rand(ki);
        end
    end
    
    methods
        % #################################################################
        function obj = TopoBuilder(ne,ke,ni,ki,netSizeE_in,netSizeI_in)  
            if(nargin < 6)
                netSizeI_in = [];
                if(nargin < 5)
                    netSizeE_in = [];
                end
            end

            assert(isnumeric(ne));
            assert(isnumeric(ke));
            assert(isnumeric(ni));
            assert(isnumeric(ki));
            if(~isempty(netSizeE_in))
                 assert(prod(netSizeE_in) == ne);
            end
            if(~isempty(netSizeI_in))
                 assert(prod(netSizeE_in) == ni);
            end         
            
            obj.ni          = ni;
            obj.ki          = ki;
            obj.ne          = ne;
            obj.ke          = ke;
            obj.nMax_exc    = ne;
            % n x k: post index for neuron n to neuron k
            % if inh/exc neuron has less connections k, then the entry will
            % stay zero
            obj.postidxMat = zeros(ni+ne,max(ki,ke));
            % the nMatSize is a virtual matrix size, if the neurons are
            % connected in 2 dimensions. For example 20 Neurons might be
            % interconnected in a matrix with size of 4 x 5 
            obj.netSizeE = netSizeE_in;
            obj.netSizeI = netSizeI_in;
        end
        
        function inp = update_inp(obj,inp)
            inp.ni = obj.ni;
            inp.ne = obj.ne;
            inp.ki = obj.ki;
            inp.ke = obj.ke;
            inp.netSizeE_in = obj.netSizeE;
            inp.netSizeI_in = obj.netSizeI;
        end
        % #################################################################
        function [] = plotMat_2D(obj)
            [from2D,to2D] = obj.postidxMat_2D();
            figure;
            to2Dinh = to2D(to2D(:,3)==1,:);
            to2Dexc = to2D(to2D(:,3)==0,:);
            from2Dinh = from2D(from2D(:,3)==1,:);
            from2Dexc = from2D(from2D(:,3)==0,:);

            subplot(1,2,1);
            plot(to2Dexc(:,2),to2Dexc(:,1),'*b'); hold on;
            plot(from2Dexc(:,2),from2Dexc(:,1),'*r');
            xlim([0 obj.netSizeE(2)+1]);
            ylim([0 obj.netSizeE(1)+1]);
            title('excitation matrix');

            subplot(1,2,2);
            plot(to2Dinh(:,2),to2Dinh(:,1),'*b'); hold on;
            plot(from2Dinh(:,2),from2Dinh(:,1),'*r');
            xlim([0 obj.netSizeI(2)+1]);
            ylim([0 obj.netSizeI(1)+1]);
            title('inhibition matrix'); 
        end
        
        % #################################################################
        function [from2D,to2D] = postidxMat_2D(obj)
            n = obj.ni + obj.ne;
            % these are the 2d representations of the inhibitor matrix and
            % excitator matrix indices. "isInh" means if 1, this is an
            % inhibition index, 0 is an excitation index
            % [r1 c1 isInh]
            % [rn cn isInh]
            from2D = [];
            to2D   = [];
            IS_INH = 1;
            IS_EXC = 0;
            for i=1:n
                postidx = obj.postidxMat(i,:);
                % remove empty places
                postidx = postidx(postidx ~= 0); 
                if(isempty(postidx))
                   continue; 
                end
                % inhibitor index
                [postRi, postCi] = ind2sub(obj.netSizeE, ...
                postidx(postidx > prod(obj.netSizeE))-prod(obj.netSizeE));
                % excitator index
                [postRe, postCe] = ind2sub(obj.netSizeE, ...
                postidx(postidx <= prod(obj.netSizeE)));
                
                to2D = [to2D; ...
                    postRe' postCe' ones(length(postRe),1)*IS_EXC];
                to2D = [to2D; ...
                    postRi' postCi' ones(length(postRi),1)*IS_INH];
                
                if(i<=prod(obj.netSizeE))
                    % for excitators
                    [fromR, fromC] = ind2sub(obj.netSizeE, i);
                    source = IS_EXC;
                else
                    % for inhibitors
                    [fromR, fromC] = ind2sub(obj.netSizeE, i-prod(obj.netSizeE));
                    source = IS_INH;
                end
                OnesVec = ones(length(postRe),1);
                from2D = [from2D; ...
                    fromR*OnesVec fromC*OnesVec source*OnesVec];  
                OnesVec = ones(length(postRi),1);
                from2D = [from2D; ...
                    fromR*OnesVec fromC*OnesVec source*OnesVec];
            end
        end
        
        % #################################################################
        %tb = tb.connect_e2i('rand',numOfOutConPerN);
        function obj = connect_e2i_rand(obj, numOfOutConPerN)
            obj = obj.connect(@obj.set_postidxMat_e2i, "rand", numOfOutConPerN);
        end
        function obj = connect_i2e_rand(obj, numOfOutConPerN)
            obj = obj.connect(@obj.set_postidxMat_i2e, "rand", numOfOutConPerN);
        end
        function obj = connect_e2e_rand(obj, numOfOutConPerN)
            obj = obj.connect(@obj.set_postidxMat_e2e, "rand", numOfOutConPerN);
        end
        function obj = connect_i2i_rand(obj, numOfOutConPerN)
            obj = obj.connect(@obj.set_postidxMat_i2i, "rand", numOfOutConPerN);
        end
        % #############
        function obj = connect_e2i_fromTo(obj,from,to)
            obj = obj.connect(@obj.set_postidxMat_e2i, "from", from, "to", to);
        end
        function obj = connect_i2e_fromTo(obj,from,to)
            obj = obj.connect(@obj.set_postidxMat_i2e, "from", from, "to", to);
        end
        function obj = connect_e2e_fromTo(obj,from,to)
            obj = obj.connect(@obj.set_postidxMat_e2e, "from", from, "to", to);
        end
        function obj = connect_i2i_fromTo(obj,from,to)
            obj = obj.connect(@obj.set_postidxMat_i2i, "from", from, "to", to);
        end
        % #############
        function obj = connect_e2i_fromToNeighbor(obj,from,to,nDistance)
            obj = obj.set_postidxMat_neighbor(@obj.connect_e2i_fromTo,from,to,nDistance);
        end
        function obj = connect_i2e_fromToNeighbor(obj,from,to,nDistance)
            obj = obj.set_postidxMat_neighbor(@obj.connect_i2e_fromTo,from,to,nDistance);
        end
        function obj = connect_e2e_fromToNeighbor(obj,from,to,nDistance)
            obj = obj.set_postidxMat_neighbor(@obj.connect_e2e_fromTo,from,to,nDistance);
        end
        function obj = connect_i2i_fromToNeighbor(obj,from,to,nDistance)
            obj = obj.set_postidxMat_neighbor(@obj.connect_i2i_fromTo,from,to,nDistance);
        end
    end
    % #####################################################################
    % #####################################################################
    
    methods(Access = protected)
        
        function obj = set_postidxMat_neighbor(obj,conFunc,from,to,nDistance)
        % n neurons can be seen as 2D arrangement of sqrt(n) x sqrt(n)
        % example of excitators:
        %             1   6  11  16  21
        %             2   7  12  17  22
        %             3   8  13  18  23 
        %             4   9  14  19  24
        %             5  10  15  20  25
        %
        % For local inhibition/excitation each excitation neuron excites
        % it's near neighbors and some inhibitors to inhibit wider
        % neighborhood for example
        %             i   i   i   i   i     A: active neuron
        %             i   e   e   e   i     e: excited neighbors
        %             i   e   A   e   i     i: indirect inhibition
        %             i   e   e   e   i
        %             i   i   i   i   i
        % nDistance = [start end]
        
            % check which function handle used
            fInfo = functions(conFunc);
            switch fInfo.function
                case '@(varargin)obj.connect_e2e_fromTo(varargin{:})'
                    targetNetSize   = obj.netSizeE;
                case '@(varargin)obj.connect_i2i_fromTo(varargin{:})'
                    targetNetSize   = obj.netSizeI;
                case '@(varargin)obj.connect_e2i_fromTo(varargin{:})'
                    targetNetSize   = obj.netSizeI;
                case '@(varargin)obj.connect_i2e_fromTo(varargin{:})'
                    targetNetSize   = obj.netSizeE;
                otherwise
                    error('TopoBuilder:set_postidxMat_neighbor wrong function handle as argument!');
            end
            % get 2D indices the get 2D neighbors
            assert(~isempty(targetNetSize),...
                'TopoBuilder:set_postidxMat_neighbor targetNetSize (netSizeE or netSizeI) not given!')

            [rt,ct] = ind2sub(targetNetSize, to);
            L       = length(rt);
            d       = nDistance(2);
            d1      = nDistance(1);
            count   = 0;
            % max number of neighbors
            toNeighbor = zeros(1,sum(nDistance(1):nDistance(2))*8)*L;
            from_mod   = toNeighbor;
            for i=1:L
                % get all neighbors of "to"-point
                rnMin = 1;
                rnMax = targetNetSize(1);
                cnMin = 1;
                cnMax = targetNetSize(2);
                for rn=rt(i)-d:rt(i)+d                    
                    for cn=ct(i)-d:ct(i)+d 
                        if(rn > rt(i)-d1 && rn < rt(i)+d1 && ...
                           cn > ct(i)-d1 && cn < ct(i)+d1    )
                            continue;
                        end
                        toIdx = sub2ind(targetNetSize,mod(rn-1, rnMax)+1,...
                                                      mod(cn-1, cnMax)+1);
                        % is inner point, do not take TODO ERROR IN MODULO CASE!
                        if(to(i) == toIdx)
                            continue;
                        end
                        count = count + 1;
                        toNeighbor(count) = toIdx;
                        from_mod(count)   = from(i);
                    end
                end
            end
            from_mod(from_mod == 0)     = [];
            toNeighbor(toNeighbor == 0) = [];
            % inhibitors are in the upper end, use offset
            obj = conFunc(from_mod,toNeighbor);
        end
        
        % #################################################################
        function obj = connect(obj, set_postidxMat_func, varargin)
            for i=1:2:length(varargin)
                switch varargin{i}
                    case "from"
                        fromVec = varargin{i+1};
                    case "to"
                        toVec = varargin{i+1};    
                        obj = set_postidxMat_func([],fromVec,toVec);
                    case "rand"
                        numOfOutConPerN = varargin{i+1};
                        if(isempty(numOfOutConPerN))
                            obj = set_postidxMat_func();
                        else
                            obj = set_postidxMat_func(numOfOutConPerN);
                        end                        
                    otherwise
                        error(['Argument ' varargin{i} ...
                            ' not known! Use "from" | "rand"!']);
                end
            end            
        end
        
        % #################################################################
        % calculate post index matrix for e2i. The indices are linear
        % indices and must be 1...ne for excitation neurons and 1...ni for
        % inhibition neurons
        function obj = set_postidxMat_e2i(obj,kset,from,to)
            % unpack local for speedup
            nil = obj.ni;
            nel = obj.ne;
            kel = obj.ke;
            kil = obj.ki;
            k   = max(kel,kil); 
            if(nargin < 2)
                kset   = kel;
            end
            if(nargin < 3)
                from = 1:nel;
                to   = [];
            end
            
            postidxMatl = obj.postidxMat;
            count = 0;
            for j=from
                count = count + 1;
                % rule for excitators ####################
                if(isempty(to))
                    postidx = randperm(nil,kset)+nel;
                else
                    % set one target for each from-entry
                    kset    = 1;
                    postidx = to(count) + nel;
                end
                
                % ########################################
                % store connections in empty (0) places
                numFilled = sum(postidxMatl(j,:) ~= 0);
                sel = [false(1,numFilled) true(1,kset) false(1,k-kset-numFilled)];
                postidxMatl(j,sel) = postidx;
            end         
            obj.postidxMat = postidxMatl;
        end
        
        % #################################################################
        % calculate post index matrix
        function obj = set_postidxMat_i2e(obj,kset,from,to)
            % unpack local for speedup
            nil = obj.ni;
            nel = obj.ne;
            kel = obj.ke;
            kil = obj.ki;
            k   = max(kel,kil); 
            if(nargin < 2)
                kset   = kil;
            end
            if(nargin < 3)
                from = 1:nil;
                to   = [];
            end
            
            postidxMatl = obj.postidxMat;
            count = 0;
            for j=from
                count = count + 1;
                % rule for excitators #################### 
                if(isempty(to))
                    postidx = randperm(nel,kset);
                else
                    % set one target for each from-entry
                    kset    = 1;
                    postidx = to(count);
                end
                % ########################################
                % store connections in empty (0) places
                numFilled = sum(postidxMatl(j+nel,:) ~= 0);
                sel = [false(1,numFilled) true(1,kset) false(1,k-kset-numFilled)];
                postidxMatl(j+nel,sel) = postidx;
            end         
            obj.postidxMat = postidxMatl;
        end
        
        % #################################################################
        % calculate post index matrix
        function obj = set_postidxMat_e2e(obj,kset,from,to)
            % unpack local for speedup
            nel = obj.ne;
            kel = obj.ke;
            kil = obj.ki;
            k   = max(kel,kil); 
            if(nargin < 2)
                kset   = kel;
            end
            if(nargin < 3)
                from = 1:nel;
                to   = [];
            end
            
            postidxMatl = obj.postidxMat;
            count = 0;
            for j=from
                count = count + 1;
                % rule for excitators #################### 
                if(isempty(to))
                    postidx = randperm(nel-1,kset);
                else
                    kset    = 1;
                    assert(j ~= to(count), ['TopoBuilder:set_postidxMat_e2e: ' ...
                        'To and from neuron is not allowed to be the same!']);
                    postidx = to(count);
                end
                % ########################################
                if(isempty(to))
                    for i = 1:kset % avoid selfconnections
                        if(postidx(i) >= j)
                            postidx(i) = postidx(i) + 1;
                        end
                    end
                end
                % store connections in empty (0) places
                numFilled = sum(postidxMatl(j,:) ~= 0);
                sel = [false(1,numFilled) true(1,kset) false(1,k-kset-numFilled)];
                postidxMatl(j,sel) = postidx;
                
            end         
            obj.postidxMat = postidxMatl;
        end
        
        % #################################################################
        % calculate post index matrix
        function obj = set_postidxMat_i2i(obj,kset,from,to)           
            % unpack local for speedup
            nel = obj.ne;
            nil = obj.ni;
            kel = obj.ke;
            kil = obj.ki;
            k   = max(kel,kil); 
            if(nargin < 2)
                kset   = kil;
            end
            if(nargin < 3)
                from = 1:nil;
                to   = [];
            end
            
            postidxMatl = obj.postidxMat;
            count = 0;
            for j=from+nel
                count = count + 1;
                % rule for excitators ####################  
                if(isempty(to))
                    postidx = randperm(nil-1,kset)+nel;
                else
                    kset    = 1;
                    assert(j ~= to(count)+nel, ['TopoBuilder:set_postidxMat_e2e: ' ...
                        'To and from neuron is not allowed to be the same!']);
                    postidx = to(count)+nel;
                end
                % ########################################
                if(isempty(to))
                    for i = 1:kset % avoid selfconnections
                        if(postidx(i) >= j)
                            postidx(i) = postidx(i) + 1;
                        end
                    end
                end
                % store connections in empty (0) places
                numFilled = sum(postidxMatl(j,:) ~= 0);
                sel = [false(1,numFilled) true(1,kset) false(1,k-kset-numFilled)];
                postidxMatl(j,sel) = postidx;
            end         
            obj.postidxMat = postidxMatl;
        end    
    end
end

