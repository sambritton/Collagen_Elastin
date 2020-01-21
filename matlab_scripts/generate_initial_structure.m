function T2=generate_initial_structure(N, coord1, coord2, coord3, ...
    mu, sigma, max_edge_length, ro, edge_density, degreeChoice)

T2 = [];

maxPossibleNbrs = ceil(8*(max_edge_length^3)*ro);%triple the average per square region
allLength = inf(length(coord1),maxPossibleNbrs);%holds lengths
allAngle = inf(length(coord1),maxPossibleNbrs);%angles close to 1 means z alignment, angles far from 1 mean no alignment
allId = zeros(length(coord1), maxPossibleNbrs);%holds id's

for i=1:length(coord1)
    for j=1:length(coord1)
        edge = [coord1(i)-coord1(j), coord2(i)-coord2(j),coord3(i)-coord3(j)];
        edgeLength = norm(edge);
        if (0<edgeLength) && (edgeLength<max_edge_length)
            %then we can place an edge. 
            for row = 1:maxPossibleNbrs
                if (allLength(i,row) == inf)
                    %then we have found a place to put our edge length
                    allLength(i,row) = edgeLength;
                    
                    %-1 is 90 degree, 1 is 0/180
                    allId(i,row) = j;
                    break ;
                end
            end
        end
    end
end

%Now compute the density of nodes in a nbhd of radius maxLength around each point.
degreeNear = [];
for i = 1:length(coord1)
    num_nbrs_i = find(allLength(i,:) < max_edge_length);
    degreeNear(i) = length(num_nbrs_i);

end

%we sort the degree so the highest degree will be at the most dense place
[degreeAll,index] = sort(degreeNear);

%% Set the degree of each to the desired degree.
bins = floor(degreeChoice * N);

difference = N - sum(bins);
%choose three largest bins at random.
for i = 1:difference
    [~, ind] = sort(bins,'descend');
    choice = ind(ceil(3*rand));
    bins(choice) = bins(choice) + 1;
end

counter = 1;
for i = 1:length(degreeChoice)
    for j = 1:bins(i)
        %Set the degree
        degreeAll(counter) = i;%(prefDegree - 1) + i;
        counter = counter + 1;
    end
end

    
%% we construct the edges. Here choose the distribution for the edges. 


allPoints = [coord1;coord2;coord3]';

%choose random starting point

random = ceil(length(coord1)*rand);
startPoint = [coord1(random),coord2(random),coord3(random)];
indexStart = random;

%we will run the edge connecting until temp is empty
%we subtract 1 from the point where an edge is made 
%fill temp with the desired degree of each node computed above.
allEdges = zeros(length(coord1),1);
for i = 1:length(coord1)
    allEdges(index(i)) = degreeAll(index(i));
end
sum_temp = sum(allEdges);



%% Making The Graph 
while (edge_density * sum_temp/2 > length(T2))
    %each end point is weighted
    possiblePoints = nonzeros(allId(indexStart,:));
    probOfStep = [];
    
    for i = 1:length(possiblePoints)   
        length_i = [allLength(indexStart, i)];%lengths associated to id
        %angle_i = 0.5*[allAngle(indexStart, i)];
        %%%%%
        probOfStep(i) = lognpdf(length_i,mu,sigma);%weighted
        
        
    end
       
    %choose the endpoint according to this weight.
    indexChoice = find(probOfStep == max(probOfStep));
    prefEdgeIndex = possiblePoints(indexChoice);
    %[row, col] = find(allId(indexStart,:) == prefEdgeIndex)
    prefEdgeLen = allLength(indexStart, indexChoice);
 
    %run this part until a endpoint is chosen each iteration attempts to
    %place a single edge.
    
    counter = length(T2)+1;
    while length(T2) < counter
        
        %Initial choice test. the simple line worked better
        %width = max(probOfStep) - min(probOfStep);
        %random = rand  - (1-prefEdgeLen) + width;
        %draw = lognpdf(random, mu, sigma);
        draw = lognrnd(mu,sigma);
        
        minProb =(min(abs(probOfStep-draw)));
        %probOfStep

        choice = find(abs(probOfStep-draw) == minProb,1);
        %choice
        indexEnd = possiblePoints(choice);
        endPoint = allPoints(indexEnd,:);
        
%% Duplicates: If there is already an edge, lia == (true, true). 
%while this holds, we check all possible edges and then break when there
%are no more. 
        lia = [0,0];
        for i = 1:size(T2,1)
            if ((indexStart == T2(i,1)) && (indexEnd ==  T2(i,2)))...
                        || ((indexStart == T2(i,2)) && (indexEnd == T2(i,1)))
                lia(1,1) = true;
                lia(1,2) = true;
            end
        end

        %to avoid duplicates 
        while (lia(1,1) == true) && (lia(1,2) == true)%check if the edge is already in T2
            %re-choose the end point
            probOfStep(choice) = 0;
            %random = rand  - (1-prefEdgeLen) + width;
            %draw = lognpdf(random, mu, sigma); %+ min(probOfStep);
            draw = lognrnd(mu,sigma);
            minProb =(min(abs(probOfStep-draw)));
            choice = find((abs(probOfStep-draw) == minProb),1);
            indexEnd = possiblePoints(choice);

            endPoint = allPoints(indexEnd,:);
            
            %make sure you reset the lia
            lia = [0,0];
            for i = 1:length(T2)
                if (((indexStart == T2(i,1)) && (indexEnd ==  T2(i,2)))...
                        || ((indexStart == T2(i,2)) && (indexEnd == T2(i,1))))
                    lia(1,1) = true;
                    lia(1,2) = true;
                end
            end
            
            if any(probOfStep ~= 0)
                continue;
            else
                %if you reach here, all probOfStep will be zero.
                break
            end
        end
        
%% Adding Edges is done on a first needed basis. If allEdges has more room and probOfstep is not all zeros 
% a new edge is created.
         %indexEnd
         %endPoint
        if (allEdges(indexEnd) > 0) &&...
                (any(probOfStep ~= 0) &&...
                (true ))
            T2 = [T2; [indexStart, indexEnd]];
            allEdges(indexStart) = allEdges(indexStart) - 1;%account for created edge by subtracting off
            allEdges(indexEnd) = allEdges(indexEnd) - 1;
            
            %indexStart = indexEnd; %use this for random walk constructor
            %startPoint = endPoint; %use this for random walk
            
            %next point is chosen by mean
            possibleStarts = []; %holds top 
            threshold = mean(allEdges); %notice this is different from the bottom
            for i = 1:length(allEdges)
                if (allEdges(i)> threshold)
                    possibleStarts = [possibleStarts; i];
                end
            end
            
            random = ceil(length(possibleStarts)*rand);
            newIndex = possibleStarts(random);
            startPoint = [coord1(newIndex),coord2(newIndex),coord3(newIndex)];
            indexStart = newIndex;
            
            %now that we have a new index, we need to regenerate the
            %probabilities of possible end points.
            possiblePoints = nonzeros(allId(indexStart,:));
            probOfStep = [];

            for i = 1:length(possiblePoints)   
                length_i = [allLength(indexStart, i)];
                %angle_i = 0.5*[allAngle(indexStart, i) ];
                probOfStep(i) = lognpdf(length_i,mu,sigma);%weighted
            end
        %if
        elseif any(probOfStep ~= 0)
            %if the choice was not chosen, we set its prob to zero
            %so it will not be chosen again.
            probOfStep(choice) = 0;
        else
            %you only reach this point if you cannot form an edge and all
            %probabilities of endpoints are zero.
            %chose a new point from top ~50% of allEdges
            
            possibleStarts = []; %holds top ~50% 
            threshold = mean(allEdges);
            for i = 1:length(allEdges)
                if (allEdges(i)> threshold)
                    possibleStarts = [possibleStarts; i];
                end
            end
            
            random = ceil(length(possibleStarts)*rand);
            newIndex = possibleStarts(random);
            startPoint = [coord1(newIndex),coord2(newIndex),coord3(newIndex)];
            indexStart = newIndex;
            %force the next choice to have an available edge in allEdges
            %vector
            while (allEdges(indexStart) == 0)
                random = ceil(length(possibleStarts)*rand);
                newIndex = possibleStarts(random);
                startPoint = [coord1(newIndex),coord2(newIndex),coord3(newIndex)];
                indexStart = newIndex;
            end
            break
        end
        
    end
    counter = counter+1; %this makes the while loop run until an end is chosen
end 

