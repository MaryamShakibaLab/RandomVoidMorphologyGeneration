function M = FindPolygons(A)
%% Create an M matrix that contains all the x and y coordinates arranged in a way
%  such that all polygones are distinguished from one another.
% A is the array that contains the randomnly sorted x y coordinates
% extracted from Autocad.
% Begin by storing the first row of array A into M. Check for
% the next line segment that shares endpoints with the first row.
% Once this line segment has been identified, add it to M, then set the
% coordinates pertaining to that row to be equal to 0 in the original
% Matrix A. Repeat the process until all line segments that are connected
% to each other are identified. Store these line segments into M such that
% each polygone identified is distinguished from the previous one.
%
%% initiate parameters
M=zeros(size(A,1),size(A,2));%initiate the Matrix M

%% Store the first row of A into M
for i=1:4
    M(1,i) = A(1,i);
    A(1,i) = 0;
end

count=1; %counter used to enumerate M's rows
cagg=1;
rum=1; %counter used to skip over the rows whose coordinates are 0 in A

%%
%% Find M
for i=1:size(A,1)
    while 1
        if (M(cagg,1)~=M(rum,3) || M(cagg,2)~=M(rum,4))
            for i1=1:size(A,1)
                if A(i1,4)~=0
                    %cagg=cagg+1;
                    if (M(rum,3) == A(i1,1) && M(rum,4) == A(i1,2))
                        M(rum+1,1) = A(i1,1);
                        M(rum+1,2) = A(i1,2);
                        M(rum+1,3) = A(i1,3);
                        M(rum+1,4) = A(i1,4);
                        A(i1,1)=0;
                        A(i1,2)=0;
                        A(i1,3)=0;
                        A(i1,4)=0;
                        count=count+1;
                        rum=rum+1;
                    elseif (M(rum,3) == A(i1,3) && M(rum,4) == A(i1,4))
                        M(rum+1,1) = A(i1,3);
                        M(rum+1,2) = A(i1,4);
                        M(rum+1,3) = A(i1,1);
                        M(rum+1,4) = A(i1,2);
                        A(i1,1)=0;
                        A(i1,2)=0;
                        A(i1,3)=0;
                        A(i1,4)=0;
                        count=count+1;
                        rum=rum+1;
                    end
                end
            end
        else
            break;
        end
    end
    cagg=count+1;
    rum=count+1;
    if rum ==size(A,1)
        break;
    end
    for i2=1:size(A,1)
        if A(i2,4)~=0
            if cagg~=40
                cagg=cagg+1;
                rum=rum+1;
                count=count+1;
            end
            M(cagg,1) = A(i2,1);
            M(cagg,2) = A(i2,2);
            M(cagg,3) = A(i2,3);
            M(cagg,4) = A(i2,4);
            A(i2,1)=0;
            A(i2,2)=0;
            A(i2,3)=0;
            A(i2,4)=0;
            break;
        end
    end
end
end





