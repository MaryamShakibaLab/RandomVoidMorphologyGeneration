%% FindPolygons stores the coordinates connecting each particle
% in order in a Matrix M
function M = FindPolygons(A)

M=zeros(size(A,1),size(A,2));

for i=1:4
    M(1,i) = A(1,i);
    A(1,i) = 0;
end

count=1;
cagg=1;
rum=1;

%%

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





