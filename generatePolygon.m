%% Modifies the generated void so as to eliminate any 
% extremely-elongated shapes. Considers the alpha parameter. 
%
function [xdatapoints,ydatapoints] = generatePolygon(ctrX, ctrY, aveRadius, irregularity, spikeyness, numVerts, angle, alpha)
%clc;
%clf;
%clear all;
%counter = 1;
%alpha = 1.5;
%for iiii = 0:12
%    for jjjj = 0:10
%        ctrX = 0+80*iiii;
%        ctrY = 0-80*jjjj;
%        if jjjj ~= 10
%            aveRadius = 2*(10-jjjj);
%        else
%            aveRadius = 2;
%        end
%        irregularity = 0;
%        spikeyness = 0.6;
%angle = unifrnd(0, 2*pi);
%        angle = 180;
%        numVerts = 6;
[xdatapoints,ydatapoints] = generatePolygonFunction(ctrX, ctrY, aveRadius, irregularity, spikeyness, numVerts, angle, alpha);
while polyarea(xdatapoints,ydatapoints) < ((aveRadius^2)*numVerts*sin(2*pi/numVerts)/2) || polyarea(xdatapoints,ydatapoints) > ((aveRadius^2)*numVerts*sin(2*pi/numVerts)/2)*1.1
    [xdatapoints,ydatapoints] = generatePolygonFunction(ctrX, ctrY, aveRadius, irregularity, spikeyness, numVerts, angle, alpha);
end
%         plot(xdatapoints,ydatapoints,'LineWidth',1.5);
%         areas(counter) = polyarea(xdatapoints,ydatapoints);
%         counter = counter +1 ;
%         %k = boundary(xdatapoints,ydatapoints,0);
%         hold on;
%         set(gca,'visible','off')
end


