classdef Dat2Dconvolver
    %IMGREADER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        
        function [imTransRed] = get_dat2Dconv(dat2D, kernelSize,type_str)
            if(strcmp(type_str,'centerON'))
                s = 1;
            elseif(strcmp(type_str,'centerOFF'))
                s = -1;
            else
                error('type_str must be ''centerON'' or ''centerOFF''');
            end
            kernel_2D  = s*Dat2Dconvolver.makecircularSymKernel_cON(kernelSize);
            convShift  = round(kernelSize/3);

            imTrans    = conv2(dat2D,kernel_2D,'valid'); 
            imTransRed = imTrans(1:convShift:end,1:convShift:end);
        end

        function yMat = makecircularSymKernel_cON(kernelSize)
            x = linspace(-pi,pi,kernelSize);
            yMat  = ((cos(x')+1).*(cos(x)+1)/2)-1;
            yMean = mean(mean(yMat));
            yMat  = yMat - yMean;
        end
        
        function [] = showTestPic()
            im  = imread('cameraman.tif');
            figure;

            kSv = 3:2:20;
            subplot(1,length(kSv)+1,1);
            imagesc(im);

            for i=1:length(kSv)
                kernelSize   = kSv(i);
                subplot(1,length(kSv)+1,i+1);
                [imTransRed] = Dat2Dconvolver.get_dat2Dconv(im, kernelSize, 'centerOFF');
                image(imTransRed./max(max(imTransRed)) * 265);      
            end            
        end

    end
    
end

