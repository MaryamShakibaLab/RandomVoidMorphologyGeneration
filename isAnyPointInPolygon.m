%% Check if the generated void intersects with an aggregate. 
% The output of the function is either "True", indicating that 
% the void does intersect with the aggregates and should therefore 
% be neglected, or "False" indicating that the void does not intersect 
% with any aggregate and should therefore be run through 
% the second void intersection check.
function out = isAnyPointInPolygon(circlex_coordinates,circley_coordinates,PS)

out = false;

circle_vertices_polyshape = polyshape(circlex_coordinates,circley_coordinates);

for j=1:length(PS) %length of PS is equal to the number of polygones
    intersection = intersect(circle_vertices_polyshape,PS{1,j});
    if isempty(intersection.Vertices) == 0
        out = true;
        break;
    end
end


end