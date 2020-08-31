function [xdatapoints,ydatapoints] = generatePolygonFunction(ctrX, ctrY, aveRadius, irregularity, spikeyness, numVerts, angle, alpha)
        %{
        Start with the centre of the polygon at ctrX, ctrY,
        then creates the polygon by sampling points on a circle around the centre.
        Randon noise is added by varying the angular spacing between sequential points,
        and by varying the radial distance of each point from the centre.

        Params:
        ctrX, ctrY - coordinates of the "centre" of the polygon
        aveRadius - in px, the average radius of this polygon, this roughly controls how large the polygon is, really only useful for order of magnitude.
        irregularity - [0,1] indicating how much variance there is in the angular spacing of vertices. [0,1] will map to [0, 2pi/numberOfVerts]
        spikeyness - [0,1] indicating how much variance there is in each vertex from the circle of radius aveRadius. [0,1] will map to [0, aveRadius]
        numVerts - self-explanatory

        Returns a list of vertices, in CCW order.

        Website: https://stackoverflow.com/questions/8997099/algorithm-to-generate-random-2d-polygon
        %}


        %irregularity = clip( irregularity, 0,1 ) * 2*pi/ numVerts;
        %spikeyness = clip( spikeyness, 0,1 ) * aveRadius;

        % generate n angle steps
        angleSteps = [];
        lower = (2*pi / numVerts) - irregularity;
        upper = (2*pi / numVerts) + irregularity;
        sum = 0;
        for i =1:numVerts
            angleSteps(i) = unifrnd(lower, upper);
            sum = sum + angleSteps(i);
        end
        % normalize the steps so that point 0 and point n+1 are the same
        k = sum / (2*pi);
        for i =1:numVerts
            angleSteps(i) = angleSteps(i) / k;
        end
        % now generate the points
        points = [];
        for i =1:numVerts
            r_i = clip( normrnd(aveRadius, spikeyness*aveRadius), 0, alpha*aveRadius);
            x = ctrX + r_i* cos(angle);
            y = ctrY + r_i* sin(angle);
            points(i,:)= [(x),(y)];
            angle = angle + angleSteps(i);
        end
        points = [points;points(1,:)];
        xdatapoints = points(:,1);
        ydatapoints = points(:,2);

end



