function [coord1, coord2, coord3] = length_alignment_matcher(T2, coord1, coord2, coord3, E1, E2, E3, ...
    zLine, logn_samples_collagen, error_angles, error_lengths,...
    max_edge_length, orientation_choice, orientation_bin_centers_sort_reduced);


lengths = [];
currentAngles = [];
for i=1:size(T2,1)
    x1 = coord1(T2(i,1));
    x2 = coord1(T2(i,2));
    y1 = coord2(T2(i,1));
    y2 = coord2(T2(i,2));
    z1 = coord3(T2(i,1));
    z2 = coord3(T2(i,2));
    edgeTemp = [x1 - x2, y1-y2, z1-z2];
    normTemp = norm(edgeTemp);
    lengths(i) = normTemp; 
    currentAngleTemp = 2*abs(dot(zLine, edgeTemp/normTemp))-1;
    currentAngles(i) = (currentAngleTemp*(-45) + 45);
   
end 

binSize = 0.1;
BMIN = 0;
BMAX = max_edge_length+1.0; %safety feature
boundMin = min(min(coord1),0.0);
boundMax = max(coord1);
hist = histogram(logn_samples_collagen,'binwidth',binSize,'BinLimits',[BMIN,BMAX]);
preferredNumEdgesPerBin = hist.Values;
binEdges = hist.BinEdges;


hist1 = histogram(lengths,'binwidth',binSize,'BinLimits',[BMIN,BMAX]);
currentNumEdgesPerBin = hist1.Values;

beginErrorEdges = length(T2);%sum(abs(preferredNumEdges - currentNumEdges));
errorCurrentEdges = beginErrorEdges;

%% we are going to scale the degrees  based on the left and right side
preferredNumAnglesPerBin = orientation_choice*length(T2);

angleBinDist = 5;
anglebinEdges = [orientation_bin_centers_sort_reduced-angleBinDist,90];
TempHist = histogram(discretize(currentAngles,anglebinEdges));
currentNumAnglesPerBin=TempHist.Values;

%we'll updatte lenthsUpdate and anglesUpdate
lengthsUpdate = lengths';
anglesUpdate = currentAngles';
maxAngleError = error_angles*beginErrorEdges;

useOrientation=true;
if (useOrientation)
    errorCurrentAngles = sum(abs(currentNumAnglesPerBin - preferredNumAnglesPerBin));
else 
    errorCurrentAngles=0;
end

final_error_angles = error_angles*beginErrorEdges;
final_error_lengths = error_lengths * beginErrorEdges;
iteration = 0;
while ((errorCurrentEdges> final_error_lengths) || (errorCurrentAngles> final_error_angles))
    pointChoice = 1;
    testrnd = rand;
    iteration = iteration + 1;
    if (iteration > 100000)
        break
    end
    %80% choose maximum error bins, 20% choose random
    if (rand < 0.6)
        maxError = max((currentNumEdgesPerBin - preferredNumEdgesPerBin));
        mostErrorBins = find((currentNumEdgesPerBin - preferredNumEdgesPerBin) > 0.7 * maxError);
        
        msize = numel(mostErrorBins);
        binChoice = mostErrorBins(randperm(msize, 1));
        edgeLenMin = (binChoice-1)*binSize;
        edgeLenMax = (binChoice)*binSize;
        [a,b] = find((lengthsUpdate>edgeLenMin) & (lengthsUpdate<edgeLenMax));
        row = randperm(length(a),1);
        col = randperm(2,1);
        pointChoice = T2(row,col);
    else
        pointChoice = ceil(rand*length(coord1));
    end


    %instead of choosing a random point, choose a point from top 20% error
    %bars.
    edgesAreInRange = 0;
    while (edgesAreInRange == 0)
        %choose perturbation 
        nodeChangeIsInRange = 0;
        while (nodeChangeIsInRange == 0)
            %choose node in range

            changeX = E1/4*(rand - 0.5);
            changeY = E2/4*(rand - 0.5);
            changeZ = E3/4*(rand - 0.5);

            prefOrientationAngle=0.71;%tester
            if (useOrientation)
                if (prefOrientationAngle > 0.7)
                    changeZ = E3 *(rand - 0.5);
                elseif ((prefOrientationAngle < 0.7) && (prefOrientationAngle > 0.5))
                    changeZ = 0.5*E3 *(rand - 0.5);
                else
                    changeZ =0.25*E3 *(rand - 0.5);
                end
            end    

            xTestPos = coord1(pointChoice) + changeX;
            yTestPos = coord2(pointChoice) + changeY;
            zTestPos = coord3(pointChoice) + changeZ;
            distMax = max(max(xTestPos,yTestPos),zTestPos);
            distMin = min(min(xTestPos,yTestPos),zTestPos);
            if ((distMax < boundMax) && (distMin>boundMin))
                nodeChangeIsInRange = 1;
            end
        end

        [x,y] = (find(T2 == pointChoice));

        %get nbrs
        neighbors = zeros(length(x),1);
        if (length(neighbors) == 0)
            pointChoice = ceil(rand*length(coord1));
        end
        for nbr = 1:length(x)
            row = x(nbr);
            col = y(nbr);
            nbrCol = 0;
            if (col == 2)
                nbrCol = 1;
            else
                nbrCol = 2;
            end
            neighbors(nbr) = T2(row,nbrCol);
        end

        degree = length(neighbors);
        lengthCurrent = zeros(degree,1);
        lengthTest = zeros(degree,1);
        angleCurrent = zeros(degree,1);
        angleTest = zeros(degree,1);

        for edge = 1:degree
            %calc edges
            nbr = neighbors(edge);
            currentEdge = [coord1(pointChoice)-coord1(nbr),...
                coord2(pointChoice)-coord2(nbr),...
                coord3(pointChoice)-coord3(nbr)];

            lengthCurrent(edge) = norm(currentEdge);

            testEdge = [xTestPos-coord1(nbr),...
                yTestPos - coord2(nbr),...
                zTestPos - coord3(nbr)];

            lengthTest(edge) = norm(testEdge);

            %0-90
            tempAngleCurrent = 2 * abs(dot(zLine,currentEdge/norm(currentEdge) )) - 1;
            angleCurrent(edge) = (tempAngleCurrent*(-45) + 45);

            tempAngleTest = 2 * abs(dot(zLine,testEdge/norm(testEdge) )) - 1;
            angleTest(edge) = (tempAngleTest*(-45) + 45);

        end

        if (max(lengthTest) < BMAX)
            edgesAreInRange = 1;
        end
    end




    binCurrentLength = discretize(lengthCurrent,binEdges);
    binTestLength = discretize(lengthTest,binEdges);
    binCurrentAngle = discretize(angleCurrent, anglebinEdges);
    binTestAngle = discretize(angleTest, anglebinEdges);

    %impose on test histogram.
    %now we have test lengths vs original lengths. Compare them and see if the
    %histogram is more accurate. If so, accept the movement. Else continue. 

    testCurrentNumEdgesPerBin = currentNumEdgesPerBin;

    errorCurrentEdges = sum(abs(currentNumEdgesPerBin - preferredNumEdgesPerBin));
    for num = 1:degree
        binCurrent = binCurrentLength(num);
        binTest = binTestLength(num);
        if isnan(binTest)
            disp('binTest out of bounds');
        end
        testCurrentNumEdgesPerBin(binTest) = testCurrentNumEdgesPerBin(binTest)+1;
        testCurrentNumEdgesPerBin(binCurrent) = testCurrentNumEdgesPerBin(binCurrent) - 1;
    end    
    errorTestLengths = sum(abs(testCurrentNumEdgesPerBin - preferredNumEdgesPerBin));

    %now do the same for angles
   % TempHist1 = histogram(discretize(anglesUpdate,anglebinEdges));
   % currentNumAnglesPerBin=TempHist1.Values;

    testCurrentNumAnglesPerBin = currentNumAnglesPerBin;
    if (useOrientation)
        errorCurrentAngles = sum(abs(currentNumAnglesPerBin - preferredNumAnglesPerBin));
    else 
        errorCurrentAngles=0;
    end

    for num = 1:degree
        binCurrentTempAngle = binCurrentAngle(num);
        binTestTempAngle = binTestAngle(num);
        if isnan(binTestTempAngle)
            disp('binTest out of bounds');
        end
        testCurrentNumAnglesPerBin(binTestTempAngle) = testCurrentNumAnglesPerBin(binTestTempAngle)+1;
        testCurrentNumAnglesPerBin(binCurrentTempAngle) = testCurrentNumAnglesPerBin(binCurrentTempAngle) - 1;
        if (sum(testCurrentNumAnglesPerBin<0)>0)

        end
    end    
    errorTestAngles = sum(abs(testCurrentNumAnglesPerBin - preferredNumAnglesPerBin));


    %final perturbation acceptence
    condition=(errorTestLengths <= errorCurrentEdges);
    
    cond_1 = ((errorTestLengths <= errorCurrentEdges) || (errorTestLengths < final_error_lengths)) ; %good if error decreases or if less than final
    cond_2 = ((errorTestAngles <= errorCurrentAngles) || (errorCurrentAngles<final_error_angles));
    conditionOrientation=(cond_1 && cond_2);
    
    if (useOrientation)
        condition=conditionOrientation;
    end

    if (condition)        %accept change and update point and histogram
        %update angle count
       % changeInAngles = errorTestAngleTemp - errorCurrentAngleTemp;
       % errorCurrentAngles = errorCurrentAngles - changeInAngles;

        coord1(pointChoice) = xTestPos;
        coord2(pointChoice) = yTestPos;
        coord3(pointChoice) = zTestPos;
        currentNumEdgesPerBin = testCurrentNumEdgesPerBin;    
        currentNumAnglesPerBin = testCurrentNumAnglesPerBin;
        %update lengthsUpdate and anglesUpdate to help choose new points
        %update using lengthTest and angleTest
        [a,b] = find(T2 == pointChoice);
        for i = 1:size(a,1)
            lengthsUpdate(a(i)) = lengthTest(i);
            anglesUpdate(a(i)) = angleTest(i); 
        end
    end
%display('current angle, edge')
%display(errorCurrentAngles);
%display(errorCurrentEdges);

%display('target angle, edge')
%display(final_error_angles);
%display(final_error_lengths);

end
%display(errorCurrentAngles);
%display(errorCurrentEdges);