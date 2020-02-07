

%this file makes a graph in the same way a hemGraph does, but in a
%cylindrical shape. No holes are placed in the network. This file is made
%to generate networks.
% x=linspace(0,3,20);
% max_force = 29;
% min_force = 2.1;
% C = 2.35;
% y = @(var)(-var/(log(1 - (max_force-2.1)/29)));
% 
% fun = @(x,alpha)(2.1 + max_force*abs(1-exp(-x/alpha)))
% fun_lin = @(x,alpha)(2.1 + (max_force-2.1)*abs(x/alpha));
% figure
% hold on
% plot(x,fun(x,y(C/2)))
% plot(x,fun_lin(x,(C)))

function Graph=generate_graph()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
rng shuffle

edge_density = 0.9; %Density of needed edges to be filled

%Note: 0.1 is effectively ~5% when edge_density = .95
error_lengths=0.15; %this variable should be less than 1-edge_density
error_angles=0.025;

collagen_density = 0.01;% Overall density of fibers
elastin_density = (0.3) * collagen_density;

collagen_diameter=0.075;
elastin_diameter=0.05;

%Domain lengths
E1 = 10;%xlen
E2 = 10;%ylen
E3 = 10;%zlen

max_edge_length=10;

%controls preferred edge length from paper
%mu = 2.0, 2.75 makes a huge difference, just plot it online
sigma_collagen = 1.0;
mu_collagen = 2.0;

sigma_elastin = 0.5;
mu_elastin = 0.0;

%prefDegree = 4; %minimal preferred
maxLength = max(E1,max(E2,E3));

zLine = [0,0,1]; %orientation line to align fibers

use_len_dist = true;%make true to use length distribution matching

volume = E1  * E2 * E3;
maxDegree = 8;%unused, just use data

%% Note: lots of short fibers take the same volume as a few long fibers. We need
% to calculate the volume that fibers occupy and match it with the
% samples. From the same distribution. 
logn_samples_collagen = [];
logn_samples_elastin = [];

target_collagen_volume = volume * collagen_density;
target_elastin_volume = volume * elastin_density;

current_collagen_volume=0;
current_elastin_volume=0;
tempMaxLength = maxLength+1;
while (current_collagen_volume < target_collagen_volume) 
    sample = tempMaxLength;
    while (sample >= tempMaxLength)
        sample = lognrnd(mu_collagen,sigma_collagen);
    end
    logn_samples_collagen = [logn_samples_collagen,sample];
        
    volume_of_sample = sample * pi*(collagen_diameter / 2)^2; %assume diameter of fiber is 100 nm = 0.1 microns
    current_collagen_volume = current_collagen_volume + volume_of_sample;
       
end

while (current_elastin_volume < target_elastin_volume) 
    sample = tempMaxLength;
    while (sample >= tempMaxLength)
        sample = lognrnd(mu_elastin,sigma_elastin);
    end
    logn_samples_elastin = [logn_samples_elastin,sample];
        
    volume_of_sample = sample * pi*(elastin_diameter / 2)^2; %assume diameter of fiber is 100 nm = 0.1 microns
    current_elastin_volume = current_elastin_volume + volume_of_sample;
       
end

total_collagen_fibers = length(logn_samples_collagen);
total_elastin_fibers = length(logn_samples_elastin);
total_fibers = total_collagen_fibers + total_elastin_fibers;


%% Degree collagen
%biomaterials data
degree_nodes1 = 388;
degree_nodes2 = 0;
degree_nodes3 = 646;
degree_nodes4 = 128;
degree_nodes5 = 14;
degree_nodes6 = 0;
degree_nodes7 = 0;
degree_nodes8 = 0;

degreeBiomat = [degree_nodes1,degree_nodes2,degree_nodes3,degree_nodes4,degree_nodes5,degree_nodes6,degree_nodes7,degree_nodes8];
degree_elastin = degreeBiomat / sum(degreeBiomat);

%redo for elastin
degree_nodes1 = 15;
degree_nodes2 = 15;
degree_nodes3 = 10;
degree_nodes4 = 1;
degree_nodes5 = 0;
degree_nodes6 = 0;
degree_nodes7 = 0;
degree_nodes8 = 0;

degree = [degree_nodes1,degree_nodes2,degree_nodes3,degree_nodes4,degree_nodes5,degree_nodes6,degree_nodes7,degree_nodes8];
degree_collagen = degree / sum(degree);

%% Orientation
%% orientation of fibers data.
%orientation is calculated in percentages along the following degree bins. 
orientation_bin_centers = [95,105,115,125,135,145,155,165,175,5,15,25,35,45,55,65,75,85];
no_flow_orientation = [0.052981146
0.052903428
0.053681837
0.052198755
0.058907136
0.061035482
0.059325730
0.096802632
0.094803359
0.098547209
0.098090158
0.059370446
0.062973519
0.063871322
0.051331607
0.050746768
0.046034342
0.04639513];

low_flow_orientation = [0.042720163
0.043731627
0.046191588
0.050547744
0.05092216
0.055036953
0.054090592
0.05904476
0.069506447
0.070157402
0.073671379
0.068156344
0.059806939
0.059408796
0.054537537
0.048074856
0.048625248
0.045649883];

high_flow_orientation = [0.031188673
0.032979896
0.029314074
0.027764853
0.028602296
0.034106976
0.047823461
0.072355511
0.133170597
0.143824643
0.109551209
0.062680509
0.054696462
0.049680018
0.038842595
0.038551434
0.034153912
0.03071293];

% instead of using the data, just generate 
%orientation_bin_centers
temp_elastin_sample=abs(normrnd(0,35,10000,1));
temp_elastin_sample = temp_elastin_sample(temp_elastin_sample<85);
figure
histogram(temp_elastin_sample)
elastin_orientation_sort = histc(temp_elastin_sample,sort(orientation_bin_centers));
elastin_orientation_sort=elastin_orientation_sort(1:9)/sum(elastin_orientation_sort(1:9));

temp_collagen_sample=abs(normrnd(0,7.5,10000,1));
temp_collagen_sample = temp_collagen_sample(temp_collagen_sample<85);
figure
histogram(temp_collagen_sample)
collagen_orientation_sort = histc(temp_collagen_sample,sort(orientation_bin_centers));
collagen_orientation_sort=collagen_orientation_sort(1:9)/sum(collagen_orientation_sort(1:9));

%sort data by key values
mapObjNoFlow = containers.Map(orientation_bin_centers,no_flow_orientation);
mapObjLowFlow = containers.Map(orientation_bin_centers,low_flow_orientation);
mapObjHighFlow = containers.Map(orientation_bin_centers,high_flow_orientation);

no_flow_orientationSort = cell2mat(values(mapObjNoFlow));
low_flow_orientationSort = cell2mat(values(mapObjLowFlow));
high_flow_orientationSort = cell2mat(values(mapObjHighFlow));
orientation_bin_centers_sort = cell2mat(keys(mapObjNoFlow));


%since the distribution should be symmetric, add data and reduce
orientation_bin_centers_sort_reduced = [];
no_flow_orientation_sort_reduced = [];
low_flow_orientation_sort_reduced = [];
high_flow_orientation_sort_reduced = [];
for i = 1:9
    orientation_bin_centers_sort_reduced(i) = orientation_bin_centers_sort(i);
    no_flow_orientation_sort_reduced(i) = no_flow_orientationSort(i) + no_flow_orientationSort(19-i);
    low_flow_orientation_sort_reduced(i) = low_flow_orientationSort(i) + low_flow_orientationSort(19-i);
    high_flow_orientation_sort_reduced(i) = high_flow_orientationSort(i) + high_flow_orientationSort(19-i);
    
    elastin_flow_orientation_sort_reduced(i) = elastin_orientation_sort(i);
    collagen_flow_orientation_sort_reduced(i) = collagen_orientation_sort(i);
end

%for now always use true. 
useOrientation = true;
orientation_choice_collagen = collagen_flow_orientation_sort_reduced;
orientation_choice_elastin = elastin_flow_orientation_sort_reduced;

%% End Parameters
%num_edges = (degreeBiomat .* N .* [1,2,3,4,5,6,7,8])

%number of fiber points temporary
N_collagen = 2 * round( total_collagen_fibers ./ sum(degree_collagen.*[1,2,3,4,5,6,7,8]) );
N_elastin = 2 * round( total_elastin_fibers ./ sum(degree_elastin.*[1,2,3,4,5,6,7,8]) );
N_total = N_collagen + N_elastin;
%% generate random coordinates of nodes
%this makes a shape in quadrant 1 of length Len and radius Rad
%the center is (Len/2,Rad,Rad)


%% generate random coordinates of nodes
%this makes a shape in quadrant 1 of length Len and radius Rad
%the center is (Len/2,Rad,Rad)
ro_collagen = N_collagen / volume;
collagen_coord1 = rand(1,N_collagen) * E1;
collagen_coord2 = rand(1,N_collagen) * E2;
collagen_coord3 = rand(1,N_collagen) * E3;


ro_elastin = N_elastin / volume;
elastin_coord1 = rand(1,N_elastin) * E1;
elastin_coord2 = rand(1,N_elastin) * E2;
elastin_coord3 = rand(1,N_elastin) * E3;

T2_collagen = generate_initial_structure(N_collagen, collagen_coord1, collagen_coord2, collagen_coord3,...
                                        mu_collagen,sigma_collagen, max_edge_length,...
                                        ro_collagen, edge_density, degree_collagen);

T2_elastin = generate_initial_structure(N_elastin, elastin_coord1, elastin_coord2, elastin_coord3,...
                                        mu_elastin,sigma_elastin, max_edge_length,...
                                        ro_elastin, edge_density, degree_elastin);
                                    

[collagen_coord1, collagen_coord2, collagen_coord3] = length_alignment_matcher(T2_collagen, collagen_coord1, collagen_coord2, collagen_coord3, E1, E2, E3, ...
    zLine, logn_samples_collagen, error_angles, error_lengths, ...
    max_edge_length, orientation_choice_collagen, orientation_bin_centers_sort_reduced);  
                                 
[elastin_coord1, elastin_coord2, elastin_coord3] = length_alignment_matcher(T2_elastin, elastin_coord1, elastin_coord2, elastin_coord3, E1, E2, E3, ...
    zLine, logn_samples_elastin, error_angles, error_lengths, ...
    max_edge_length, orientation_choice_elastin, orientation_bin_centers_sort_reduced); 


%set(gca, 'FontSize', 60)
figure
hold on
for i = 1:length(T2_collagen)
    pt1 = T2_collagen(i,1);
    
    pt2 = T2_collagen(i,2);
    plot3([collagen_coord1(pt1),collagen_coord1(pt2)],[collagen_coord2(pt1),collagen_coord2(pt2)],[collagen_coord3(pt1),collagen_coord3(pt2)],...
        'Color','b','LineWidth',1.5);
end

for i = 1:length(T2_elastin)
    pt1 = T2_elastin(i,1);
    
    pt2 = T2_elastin(i,2);
    plot3([elastin_coord1(pt1),elastin_coord1(pt2)],[elastin_coord2(pt1),elastin_coord2(pt2)],[elastin_coord3(pt1),elastin_coord3(pt2)],...
        'Color','r', 'LineWidth',0.05);
end





angle_current_collagen = zeros(length(T2_collagen),1);
angle_current_elastin = zeros(length(T2_elastin),1);

for edge = 1:length(T2_collagen)
    %calc edges
    pointChoice = T2_collagen(edge,1);
    nbr = T2_collagen(edge,2);
    currentEdge = [collagen_coord1(pointChoice)-collagen_coord1(nbr),...
        collagen_coord2(pointChoice)-collagen_coord2(nbr),...
        collagen_coord3(pointChoice)-collagen_coord3(nbr)];

    %0-90
    tempAngleCurrent = 2 * abs(dot(zLine,currentEdge/norm(currentEdge) )) - 1;
    angle_current_collagen(edge) = (tempAngleCurrent*(-45) + 45);


end

for edge = 1:length(T2_elastin)
    %calc edges
    pointChoice = T2_elastin(edge,1);
    nbr = T2_elastin(edge,2);
    currentEdge = [elastin_coord1(pointChoice)-elastin_coord1(nbr),...
        elastin_coord2(pointChoice)-elastin_coord2(nbr),...
        elastin_coord3(pointChoice)-elastin_coord3(nbr)];

    %0-90
    tempAngleCurrent = 2 * abs(dot(zLine,currentEdge/norm(currentEdge) )) - 1;
    angle_current_elastin(edge) = (tempAngleCurrent*(-45) + 45);


end


figure
hold on
histogram(temp_elastin_sample,'BinWidth',5,'Normalization','probability')
histogram(angle_current_elastin,'BinWidth',5,'Normalization','probability')

title('Alignment Comparison: Elastin')
xlabel('Degree, d') 
ylabel('Probability, P(d)')
legend({'Experimentally Observed Orientation','Simulated Orientation'},'Location','northeast')
hold off


figure
hold on
histogram(temp_collagen_sample,'BinWidth',5,'Normalization','probability')
histogram(angle_current_collagen,'BinWidth',5,'Normalization','probability')

title('Alignment Comparison: Collagen')
xlabel('Degree, d') 
ylabel('Probability, P(d)')
legend({'Experimentally Observed Orientation','Simulated Orientation'},'Location','northeast')
hold off

%% Possibly connect the two structures??

Fixed = 0;
Graph.node_elastin=[elastin_coord1',elastin_coord2',elastin_coord3'];
Graph.node_collagen=[collagen_coord1',collagen_coord2',collagen_coord3'];
Graph.edge_elastin=T2_elastin;
Graph.edge_collagen=T2_collagen;
Graph.fixed_node={Fixed};

DateString=['Rec_collagen_',num2str(collagen_density),'_elastin_' , num2str(elastin_density),'_', num2str(E1), '_', num2str(E2), '_', num2str(E3)];
if (useOrientation == true)   
    DateString=['Rec_collagen_',num2str(collagen_density),'_elastin_' , num2str(elastin_density),'_', num2str(E1), '_', num2str(E2), '_', num2str(E3),'_oriented',];
end
Graph.String=DateString;
mkdir(pwd,'data_files');
savdir = strcat(pwd,'\data_files' );
save(fullfile(savdir,[DateString,'.mat']),'Graph');
 