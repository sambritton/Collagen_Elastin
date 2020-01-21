function generate_graph_parameters(Graph)

%getinfo about graph

Fixed = 0;
collagen_nodes=Graph.node_collagen;
elastin_nodes=Graph.node_elastin;
collagen_edges=Graph.edge_collagen;
elastin_edges=Graph.edge_elastin;
fixed = Graph.fixed_node;


mkdir('xml_data')
savdir = strcat(pwd,"\xml_data");

N_collagen = length(collagen_nodes);
N_elastin = length(elastin_nodes);

coord1_collagen=collagen_nodes(:,1)';
coord2_collagen=collagen_nodes(:,2)';
coord3_collagen=collagen_nodes(:,3)';

coord1_elastin=elastin_nodes(:,1)';
coord2_elastin=elastin_nodes(:,2)';
coord3_elastin=elastin_nodes(:,3)';

fileID = fopen('Settings.txt','r');
xT1=textscan(fileID,'%s','whitespace','\n');
fclose(fileID);
xT1=xT1{1,1};

xT1=[cell(1,0)
'<?xml version="1.0" encoding="UTF-8"?>'
'<data>'
xT1
'    <nodes>'];

for i=1:N_collagen
    xT1=[xT1;['        <node_collagen>',num2str(coord1_collagen(i)),' ',num2str(coord2_collagen(i)),' ',num2str(coord3_collagen(i)),'</node_collagen>']];
end

for i=1:N_elastin
    xT1=[xT1;['        <node_elastin>',num2str(coord1_elastin(i)),' ',num2str(coord2_elastin(i)),' ',num2str(coord3_elastin(i)),'</node_elastin>']];
end
xT1=[xT1
'    </nodes>'];

xT1=[xT1
'    <edges>'];
for i=1:size(collagen_edges,1)
    xT1=[xT1;['        <edge_collagen>',num2str(collagen_edges(i,1)-1),' ',num2str(collagen_edges(i,2)-1),'</edge_collagen>']];
end

for i=1:size(elastin_edges,1)
    xT1=[xT1;['        <edge_elastin>',num2str(elastin_edges(i,1)-1),' ',num2str(elastin_edges(i,2)-1),'</edge_elastin>']];
end

xT1=[xT1;['    </edges>']];
% [~,mx]=max([Rx,Ry,Rz]);mx=mx(1);
% 10 percent threshold


%For the fixed points we choose the ones on the bottom of the square
xT1=[xT1  
'    <fixed>'];

xT1=[xT1
'    </fixed>' '</data>'];

DateString=Graph.String;

fileID = fopen(fullfile(savdir,['data_',Graph.String,'.xml']),'w');
C=xT1;
[nrows,~] = size(C);
formatSpec = '%s\n';
for row = 1:nrows
    fprintf(fileID,formatSpec,C{row,:});
end


