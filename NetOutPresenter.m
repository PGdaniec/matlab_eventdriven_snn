classdef NetOutPresenter
    %NETOUTPRESENTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n_subplots;
        ah;
        plotCopy;
        plotPos;
    end
    
    methods
        function obj = NetOutPresenter()
            obj.n_subplots = 0;
            obj.ah         = [];
            obj.plotCopy   = {};
            obj.plotPos    = [];
        end
        
        function obj = addSubplot(obj,t,d,pltStr,plotPos)
            St = length(t);
            if(isvector(d))
                Sd = length(d);
            else
                Sd = size(d,1);
            end
            if(St ~=Sd)
               error('NetOutPresenter: t and d must have same number of col!'); 
            end
            if(nargin < 5)
                obj.plotPos(end+1) = obj.n_subplots + 1;
                obj.n_subplots     = obj.n_subplots + 1;
            else
                if(plotPos > obj.n_subplots)
                    obj.n_subplots     = obj.n_subplots + 1;
                    obj.plotPos(end+1) = obj.n_subplots;
                else
                    obj.plotPos(end+1) = plotPos;
                end
            end

            if(isempty(obj.plotCopy))
                ph = plot(t, d, pltStr);
                obj.plotCopy{1} = copy(ph);
                obj.plotPos(1)  = 1;
                ylim([0 max(max(d))+1]); 
            else
                minEndTime  = inf;
                maxSpikeIdx = 0;
                for i=1:length(obj.plotCopy)
                    minEndTime  = min(minEndTime, obj.plotCopy{i}.XData(end));
                    maxSpikeIdx = max(maxSpikeIdx,max(max(obj.plotCopy{i}.YData)));
                end
                minEndTime  = min(minEndTime, t(end));
                maxSpikeIdx = max(maxSpikeIdx,max(max(d)));
                    
                for i=1:length(obj.plotCopy)+1
                    obj.ah(end+1) = subplot(obj.n_subplots,1,obj.plotPos(i));
                    if(i == length(obj.plotCopy)+1)
                        ph = plot(t,d,pltStr);
                        obj.plotCopy{end+1} = copy(ph);
                    else
                        plot(obj.plotCopy{i}.XData,...
                             obj.plotCopy{i}.YData,...
                             obj.plotCopy{i}.Marker, 'color',obj.plotCopy{i}.Color);
                        
                    end
                    hold on;
                    xlim([0 minEndTime]);
                    ylim([0 maxSpikeIdx]);
                end
                linkaxes(obj.ah,'x'); 
            end
        end
        % ##############################
        function [] = plotMovie(obj, frameL_s, playL_s)
            figure('units','normalized','outerposition',[0 0 1 1]);
            minEndTime  = inf;
            maxSpikeIdx = 0;
            for i=1:length(obj.plotCopy)
                for l=1:length(obj.plotCopy{i})
                    minEndTime  = min(minEndTime, obj.plotCopy{i}(l).XData(end));
                    maxSpikeIdx = max(maxSpikeIdx,max(max(obj.plotCopy{i}(l).YData)));
                end
            end
            thisah = [];
            nf = floor(minEndTime / frameL_s);
            for f = 2:nf
                for i=1:length(obj.plotCopy)
                    thisah(i) = subplot(obj.n_subplots,1,obj.plotPos(i)); %#ok
                    for l=1:length(obj.plotCopy{i})
                        sel = obj.plotCopy{i}(l).XData < frameL_s *  f   & ...
                              obj.plotCopy{i}(l).XData > frameL_s * (f-1)      ;
                        plot(obj.plotCopy{i}(l).XData(sel)-frameL_s*(f-1),...
                             obj.plotCopy{i}(l).YData(sel),...
                             obj.plotCopy{i}(l).Marker, 'color',obj.plotCopy{i}(l).Color);
                        hold on;
                    end
                    xlim([0 frameL_s]);
                    ylim([0 maxSpikeIdx]);
                end
                hold off;
                linkaxes(thisah,'x');  
                pause(playL_s/(nf-1));
            end
        end
    end
    % ####################################################################
    methods(Static)
        function [] = show_infos(rate,sidx,stimes,elapT,inp,debugValue)
            numberOfspikesVec = zeros(inp.ne+inp.ni,1);
            for i=1:inp.ne+inp.ni
                s = sum(sidx == i);
                numberOfspikesVec(i) = s;
            end
            numOfNeuronsWithNoSpike = sum(numberOfspikesVec == 0);
               disp(' ');
               disp('##############################');
            fprintf('rate             = %.2f Hz\n',rate);
            fprintf('Sim Time         = %.2f s\n',elapT);
            fprintf('end Time         = %.2f s\n',stimes(end));
            fprintf('max delta phi    = %f\n',debugValue.debug_max_deltaPhi);
            fprintf('min delta phi    = %f\n',debugValue.debug_min_deltaPhi);
            fprintf('syn + mod #      = %d\n',debugValue.debug_num_INC_synAdaptions);
            fprintf('syn - mod #      = %d\n',debugValue.debug_num_DEC_synAdaptions);
            fprintf('n. not resp. #   = %d\n',debugValue.debug_num_refactoryUsed);
            fprintf('n. not spiked. # = %d\n',numOfNeuronsWithNoSpike); 
        end
        % ###################################
        function [] = showStats(inp, sidx, debugValue)
            if(isfield(inp,'use_synW_adaption') && inp.use_synW_adaption == true)
                figure(3);
                subplot(1,2,1);
                h1 = histogram((reshape(debugValue.synW_mat_end(1:inp.ke,1:inp.ne), ...
                                   inp.ne*inp.ke,1)));
                set(gca, 'yscale','log');
                title(['excitory weights | start = ' num2str(debugValue.synW_mat_start(1,1))]);
                subplot(1,2,2);
                h2 = histogram(reshape(debugValue.synW_mat_end(1:inp.ki,inp.ne+1:inp.ni+inp.ne), ...
                                   (inp.ni)*inp.ki,1));
                set(gca, 'yscale','log');                   
                title(['inhibitory weights | start = ' num2str(debugValue.synW_mat_start(1,end))]);

                figure(4);
                subplot(1,2,1)
                imagesc(log10(reshape(debugValue.synW_mat_end(1:inp.ke,1:inp.ne),...
                                     inp.ne,inp.ke)));
                w_exc = debugValue.synW_mat_end(1:inp.ke,1:inp.ne);
                fprintf('mean exc.weight  = %f\n',mean(w_exc(:)));
                title('excitory log10(synWeight)');
                colorbar;
                subplot(1,2,2)
                imagesc(log10(-reshape(debugValue.synW_mat_end(1:inp.ki,inp.ne+1:inp.ni+inp.ne),...
                                      (inp.ni),inp.ki)));
                w_inh = debugValue.synW_mat_end(1:inp.ki,inp.ne+1:inp.ni+inp.ne);
                fprintf('mean inh.weight  = %f\n',mean(w_inh(:)));
                title('inhibitory log10(synWeight)');
                colorbar;
            end

            figure(5)
            numberOfspikesVec = zeros(inp.ne+inp.ni,1);
            for i=1:inp.ne+inp.ni
                s = sum(sidx == i);
                numberOfspikesVec(i) = s;
            end
            subplot(1,2,1)
            plot(numberOfspikesVec(1:inp.ne,1),'*')
            xlabel('number of exc spikes');
            ylabel('count');
            subplot(1,2,2)
            plot(numberOfspikesVec(inp.ne+1:inp.ni+inp.ne,1),'*')
            xlabel('number of inh spikes');
            ylabel('count');
        end

    end
    
end

