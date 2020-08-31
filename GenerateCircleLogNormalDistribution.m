%% GenerateCircleLogNormalDistribution Generates the random circular 
%voids using a log-normal distribution of reported void sizes
%
function [counter_of_air_voids, air_void_content] = GenerateCircleLogNormalDistribution(M,nstep,target)
%% Initiate Vectors

count3=1;
Radius=[];
AreaAV=[];


X = zeros(1,2*(size(M,1))); %Create an empty X row vector to contain
%the starting and ending  X coordinates of all the polygones
%combined which will be extracted from matrix M.
%Example M = [1 2 3 4; 5 6 7 8], then X = [1 3 5 7]
Y = zeros(1,2*(size(M,1))); %Create an empty Y row vector to contain
%the starting and ending  Y coordinates of all the polygones
%combined which will be extracted from matrix M.
%Example M = [1 2 3 4; 5 6 7 8], then Y = [2 4 6 8]

PolyXcoord = []; %Create an empty row vector that will contain only the
% X coordinates of one polygone at a time
PolyYcoord = []; %Create an empty row vector that will contain only the
% Y coordinates of one polygone at a time

%Create the X and Y row vectors
%inititate the X and Y row vectors using the first row of M (i.e. first
%line segment)
% X(1,1) = M(1,1);
% X(1,2) = M(1,3);
% Y(1,1) = M(1,2);
% Y(1,2) = M(1,4);
for j = 1:size(M,1)
    if j~=size(M,1)
        if M(j,1) ~= M(j+1,3) || M(j,2) ~= M(j+1,4)
            X(2*j-1) = M(j,1);
            X(2*j) = M(j,3);
            Y(2*j-1) = M(j,2);
            Y(2*j) = M(j,4);
        end
    else
        X(2*j-1) = M(j,1);
        X(2*j) = M(j,3);
        Y(2*j-1) = M(j,2);
        Y(2*j) = M(j,4);
    end
end
%%
%count the number of polygones that the matrix M has, to use for the next
%section when we will need to know how many polygones to loop for to plot.

counterofpoylgones = 0; %initiate the counter to start from 0;
p = 1; %initiate another counter that skips
%over to next polygone once the previous polygone has been identified
count=1;
while 1
    if  count>=size(M,1)
        break;
    else
        if M(count,1) == M(p,3) && M(count,2) == M(p,4)
            counterofpoylgones = counterofpoylgones+1;
            count=p+1;
        else
            p=p+1;
        end
    end
end

%% Split the X and Y row vectors into other row vectors that contain coordinates of just
% one polygone (i.e. the PolyXcoord and PolyYcoord arrays).

nextcoord = 1; %counter used to store the next coordinate into Xcoord and Ycoord
numpoly =1; %counter used to count the number of coordinate points each polygone has
startcoordpoly = 1; %counter used to store the starting x and y coordinate of the next polygone
XX=cell(1,counterofpoylgones);
YY=cell(1,counterofpoylgones);
for k =1:counterofpoylgones %k is a dummy looper that indicates when to stop plotting,
    %that is, when all polygones have been plotted
    for i=2:length(X)-1
        PolyXcoord(1,1)=X(1,startcoordpoly);
        PolyYcoord(1,1)=Y(1,startcoordpoly);
        %check if the first coordinates of the X and Y row vectors are
        %equal to any of the next ones, if not, set the next coordinate
        %of the Xcoord row vector to be equal to that entry,
        %if yes, add that entry to the vector and stop stop checking
        if PolyXcoord(1,1) ~= X(1,i+startcoordpoly-1) || PolyYcoord(1,1) ~= Y(1,i+startcoordpoly-1)
            PolyXcoord(1,nextcoord+1) = X(1,i+startcoordpoly-1);
            PolyYcoord(1,nextcoord+1) = Y(1,i+startcoordpoly-1);
            nextcoord = nextcoord+1;
            numpoly = numpoly +1;
        elseif PolyXcoord(1,1) == X(1,i+startcoordpoly-1) || PolyYcoord(1,1) == Y(1,i+startcoordpoly-1)
            PolyXcoord(1,nextcoord+1) = X(1,i+startcoordpoly-1);
            PolyYcoord(1,nextcoord+1) = Y(1,i+startcoordpoly-1);
            nextcoord = nextcoord + 1;
            numpoly = numpoly+1;
            break;
        end
    end
    
    XX(1,k)={PolyXcoord};
    YY(1,k)={PolyYcoord};
    
    %plot(PolyXcoord,PolyYcoord);
    %hold on;
    
    %re-generate the counters
    numpoly = numpoly +1;
    startcoordpoly = numpoly;
    nextcoord = 1; %restart ct to be equal to one
    
    PolyXcoord = []; %begin another Xcoord vector (i.e. next polygone)
    PolyYcoord = []; %begin another Ycoord vector (i.e. next polygone)
end
%% Create a vector PS that contains all the polyshapes
% (i.e. shapes of all polygons)

polygone=polyshape();
PS=cell(1,counterofpoylgones);
% Polyg=cell(1,counterofpoylgones);
% Pcell=cell(1,counterofpoylgones);

for ii=1:length(XX)
    polygone(ii) = polyshape(XX(ii),YY(ii));
    %     (polygone);
    PS(1,count3) = {polygone(ii)};
    count3 = count3+1;
    %     hold on;
end

% fields = {'x','y'};
% for polnum = 1:counterofpoylgones
%     Polyg{1,polnum} = PS{1,polnum}.Vertices;
% end

% for pcc=1:counterofpoylgones
%     Pcell{pcc} = num2cell(Polyg{1,pcc},1);
% end

%P = struct('x',{},'b',{});

% for pst=1:length(Pcell)
%     P(pst) = cell2struct(Pcell{1,pst},fields,2);
% end

% min_d=zeros(1,12432);
% for cop=1:counterofpoylgones
%     for copp=1:counterofpoylgones
%         if cop~=copp
%             min_d(count6) = min_dist_between_two_polygons(P(:,cop),P(:,copp),0);
%             count6 = count6 + 1;
%         end
%     end
% end
%r=min(min_d);
%r = 0.1;
%% Create row vectors containing X and Y coordinates of polygones
% delimited by NaN to use for the inpolygon function
acount = 0;
XXcoord =cell(1,2*counterofpoylgones-1);
YYcoord =cell(1,2*counterofpoylgones-1);
for o=1:2*counterofpoylgones-1
    if mod(o,2) == 1
        XXcoord(1,o) = XX(1,o-acount);
        YYcoord(1,o) = YY(1,o-acount);
        acount = acount + 1;
    elseif mod(o,2) == 0
        XXcoord(1,o) = {NaN};
        YYcoord(1,o) = {NaN};
    end
end
F = cell2mat(XXcoord);
G = cell2mat(YYcoord);
%% Uncomment this section to generate aggregates' properties
% VOP = fopen('Vertices_of_Polygones', 'w'); %Vertices of polygones
% for pp=1:counterofpoylgones
%     for pst=1:size(PS{1,pp}.Vertices)
%         fprintf(VOP,'s.Spot(point=(%3.4f, %3.4f))\n',PS{1,pp}.Vertices(pst,1),PS{1,pp}.Vertices(pst,2));
%     end
% end
% fclose(VOP);
% 
% LCPV = fopen('Line_connector_of_Polygone_Vertices', 'w');  %LCPV = Line connector of polygone vertices
% for pm = 1:size(M)
%     fprintf(LCPV,'s.Line(point1=(%3.4f, %3.4f), point2=(%3.4f, %3.4f))\n',M(pm,1),M(pm,2),M(pm,3),M(pm,4));
% end
% fclose(LCPV);
%% Random generation of air voids
max_x = max(X);
min_x = min(X);
max_y = max(Y);
min_y = min(Y);

total_area = (min_x-max_x)*(min_y-max_y); %total area of microstructure

air_void_content = 0;
counter_of_air_voids=0;

th = 0:pi/nstep:2*pi;

result_circles_to_be_plotted_xcoordinates = [];
result_circles_to_be_plotted_ycoordinates = [];

while 1
    radius = lognrnd(-0.429,0.727,1,1)/2;
   
    rng('shuffle');
    rand_perm_1 = randperm(9,1);
    
    rng('shuffle');
    rand_perm_2 = randperm(9,1);
    
    rng('shuffle');
    n1 = rand(1,1);
    
    rng('shuffle');
    n2 = rand(1,1); 
    
    if  0 <= rand_perm_1 && rand_perm_1 <= 1/n1
        circle_origin_xcoordinates = min_x + (max_x-min_x).*n1.*rand_perm_1;
    else
        circle_origin_xcoordinates = min_x + (max_x-min_x).*n1;
    end
    
    if  0 <= rand_perm_2 && rand_perm_2 <= 1/n2
        circle_origin_ycoordinates = min_y + (max_y-min_y).*n2.*rand_perm_2;
    else
        circle_origin_ycoordinates = min_y + (max_y-min_y).*n2;
    end 
    
    [in,~] = inpolygon(circle_origin_xcoordinates,circle_origin_ycoordinates,F,G);
    outer_circle_origin_xcoordinates = circle_origin_xcoordinates(~in); % xout is a column vector that contains the origins of the air voids that lie outside of the polygones (could have more than just one value, depending on the value of nb)
    outer_circle_origin_ycoordinates = circle_origin_ycoordinates(~in);
    
    if isempty(outer_circle_origin_xcoordinates) == 0 %if xout is not empty (i.e. there are such airvoids)
        
        x_vector = radius * cos(th) + outer_circle_origin_xcoordinates; %(size(xout,1)x 2*nstep) vector containing all the xcoordinates of the airvoid circles found in xout
        y_vector = radius * sin(th) + outer_circle_origin_ycoordinates; %vector containing all the ycoordinates of the airvoid circle found in yout
               
        isAnyPointInPolygonResult = isAnyPointInPolygon(x_vector,y_vector,PS);
        isAnyAirVoidInResultAirVoidsResult = isAnyAirVoidInResultAirVoids(x_vector,y_vector,result_circles_to_be_plotted_xcoordinates,result_circles_to_be_plotted_ycoordinates);       
       
        if isAnyPointInPolygonResult == false && isAnyAirVoidInResultAirVoidsResult == false
            
            result_circles_to_be_plotted_xcoordinates(counter_of_air_voids+1,:) = x_vector;
            result_circles_to_be_plotted_ycoordinates(counter_of_air_voids+1,:) = y_vector;
            %plot(result_circles_to_be_plotted_xcoordinates(counter_of_air_voids+1,:),result_circles_to_be_plotted_ycoordinates(counter_of_air_voids+1,:));
            %hold on;
            
            Radius(counter_of_air_voids+1) = radius;
            AreaAV(counter_of_air_voids+1) = pi*radius^2;
            area_air_voids = sum(AreaAV); %area of the air voids
            air_void_content = (area_air_voids/total_area)*100;
            
            XYCO(counter_of_air_voids+1,1) = outer_circle_origin_xcoordinates;
            XYCO(counter_of_air_voids+1,2) = outer_circle_origin_ycoordinates;
            counter_of_air_voids = counter_of_air_voids + 1;            
        end       
    end    
    if ((air_void_content-target<=0.05) && (air_void_content-target>=-0.05)) || air_void_content > target
        AOC = fopen('Origin_of_Void_Circles', 'w'); % AOC= Airvoid Origin Coordinates
        for pt=1:size(XYCO)
            fprintf(AOC,'s.CircleByCenterPerimeter(center=(%3.4f, %3.4f), point1=(%3.4f, %3.4f))\n',XYCO(pt,1),XYCO(pt,2),XYCO(pt,1)+Radius(pt),XYCO(pt,2));
        end
        fclose(AOC);
        break;
    end
end
end
