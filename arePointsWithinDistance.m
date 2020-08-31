function out = arePointsWithinDistance(x_vector,y_vector,offset)
out = true;
% plot(x_vector,y_vector);
for j=1:size(x_vector,1)-1
    for k=1:size(x_vector,1)-1
        if j ~= k
            if ( ( abs(x_vector(j,1) - x_vector(k,1))  <= offset) || (abs(x_vector(j,1) - x_vector(k,1))  <= - offset) ) && (( abs(y_vector(j,1) - y_vector(k,1)) <= offset) || (abs(y_vector(j,1) - y_vector(k,1)) <= - offset ))
                out = false;
                break;
            end
        end
    end
end

end