function PlotExtractedDataPoints(A)
%% Function that plots the extracted data points and connect them with lines
hConnections = [];
    for i=1:size(A)
        plot ([A(i,1)],[A(i,2)],'r.');
        hold on;
        plot ([A(i,3)],[A(i,4)],'b.');
        hold on;
        hConnections = [hConnections ; line([A(i,1);A(i,3)], [A(i,2);A(i,4)])];
    end
end 