% Function to create .FOLD files.
% Created by Steve Grey 25/06/2018 - University of Bristol
% Email: steven.grey@bristol.ac.uk
% NOTE: I am aware this is not the most efficient way to do what this
% code does. If any user of this code makes improvements I would appriciate
% access to any improved versions if possible.
close all
clear all
clc

%% Start by defining some of the things about the file we are going to create

% 1 = Miura Tube, 2 = Miura Sheet
Model_type = 2;

% .FOLD Version
file.spec = 1.1;

% Made in MATLAB - Keep the quotes
file.creator = '"MATLAB Code"';

% Author - Keep the quotes
file.author = '"Steve Grey"';

% Class:
% - "singleModel": A single origami model, possibly still in multiple frames to represent crease pattern, folded form, etc
file.classes = '"singleModel"';

% Error checking
if Model_type > 2
    error('Update title code to include extra case(s)')
end

% Give it a title - Depends on what we want to do
switch Model_type    
    case 1
        frame.title = '"Miura-ori Tube"';        
    case 2   
        frame.title = '"Miura-ori Sheet"';      
end

% This code is geared towards only doing the folded form here
frame.classes = '"foldedForm"';

% Again this code is only aimed towards 3D here
frame.attributes = '"3D"';

% NOTE: all units are in mm throughout
frame.unit = '"mm"';

%% Next define some of the things about the file we are going to create

% Error check
if Model_type > 2
    error('Update model generation code to include extra case(s)')
end
% Each different model type needs a different process
switch Model_type
    % Case 1 is the Miura tube
    case 1
        % Side lengths of facet
        a = 30;
        b = 30;
        % Angle of facet - Degrees
        alpha = 60;
        % Angle between facet and mid-plane - Degrees
        theta = 60;
        % Unit cells in each direction (m should always be 1 for a tube)
        m = 1;
        n = 10;
        
        % Define unit cell parameters [Schenk & Guest (2013)]
        H = a*sind(theta)*sind(alpha);
        S = b*cosd(theta)*tand(alpha)/sqrt(1+cosd(theta)^2*tan(alpha)^2);
        L = a*sqrt(1-sind(theta)^2*sind(alpha)^2);
        V = b/sqrt(1+cosd(theta)^2*tan(alpha)^2);
        
        % Define the vertices in a unit cell (From back left up to the top
        % down to the front and finally to the bottom. Then moving on to
        % the middle and finally the right)
        unit_cell_vertices(1,:)  = [0 0 0];
        unit_cell_vertices(2,:)  = [0 L H];
        unit_cell_vertices(3,:)  = [0 2*L 0];
        unit_cell_vertices(4,:)  = [0 L -H];
        
        unit_cell_vertices(5,:)  = [S V 0];
        unit_cell_vertices(6,:)  = [S L+V H];
        unit_cell_vertices(7,:)  = [S 2*L+V 0];
        unit_cell_vertices(8,:)  = [S L+V -H];
        
        unit_cell_vertices(9,:)  = [2*S 0 0];
        unit_cell_vertices(10,:) = [2*S L H];
        unit_cell_vertices(11,:) = [2*S 2*L 0];
        unit_cell_vertices(12,:) = [2*S L -H];
             
        % The vertices of the rest of the tube can be defined by offsetting
        % the unit cell in x
        for ii = 1:n
            vertices.coords(8*(ii-1)+1:8*ii,2:3) = unit_cell_vertices(1:8,2:3);
            vertices.coords(8*(ii-1)+1:8*(ii-1)+4,1) = (ii-1)*2*S;
            vertices.coords(8*(ii-1)+5:8*(ii-1)+8,1) = (ii-1)*2*S + S;
%             vertices.coords(12*(ii-1)+1:12*ii,1) = unit_cell_vertices(:,1) + ones(12,1)*(ii-1)*2*S;  
        end
        vertices.coords(8*ii+1:8*ii+4,2:3) = unit_cell_vertices(1:4,2:3);
        vertices.coords(8*(ii)+1:8*(ii)+4,1) = ii*S*2;
        
        
        % Delete any duplicated vertices
%         vertices.coords = unique(vertices.coords,'rows');
        
        % Define the faces in a unit cell on the top
        unit_cell_faces(1,:) = [1 2 6 5];
        unit_cell_faces(2,:) = [2 3 7 6];
        unit_cell_faces(3,:) = [9 5 6 10];
        unit_cell_faces(4,:) = [6 10 11 7];
        % and the bottom
        unit_cell_faces(5,:) = [1 4 8 5];
        unit_cell_faces(6,:) = [5 9 12 8];
        unit_cell_faces(7,:) = [4 8 7 3];
        unit_cell_faces(8,:) = [12 11 7 8];
        
        % Use the unit cell faces to define the vertices which surround
        % every face in the tube
        for ii = 1:n 
            faces.vertices(8*(ii-1)+1:8*ii,:) = unit_cell_faces + ones(8,4)*8*(ii-1);
        end
    % Case 2 is the Miura sheet    
    case 2
        
        % Side lengths of facet
        a = 30;
        b = 30;
        % Angle of facet - Degrees
        alpha = 60;
        % Angle between facet and mid-plane - Degrees
        theta = 60;
        % Unit cells in each direction
        m = 10;
        n = 10;
        
        % Define unit cell parameters [Schenk & Guest (2013)]
        H = a*sind(theta)*sind(alpha);
        S = b*cosd(theta)*tand(alpha)/sqrt(1+cosd(theta)^2*tan(alpha)^2);
        L = a*sqrt(1-sind(theta)^2*sind(alpha)^2);
        V = b/sqrt(1+cosd(theta)^2*tan(alpha)^2);
        
        % Define the vertices in a unit cell (From back left up to the top
        % down to the front and finally to the bottom. Then moving on to
        % the middle and finally the right)
        unit_cell_vertices(1,:) = [0 0 0];
        unit_cell_vertices(2,:) = [0 L H];
        unit_cell_vertices(3,:) = [0 2*L 0];
        
        unit_cell_vertices(4,:) = [S V 0];
        unit_cell_vertices(5,:) = [S L+V H];
        unit_cell_vertices(6,:) = [S 2*L+V 0];
        
        unit_cell_vertices(7,:) = [2*S 0 0];
        unit_cell_vertices(8,:) = [2*S L H];
        unit_cell_vertices(9,:) = [2*S 2*L 0];
        
        count = 0;
        for jj = 1:m
            for ii = 1:n
                count = count + 1;
                vertices.coords((count-1)*9+1:(count-1)*9+9,1) = unit_cell_vertices(:,1) + 2*S*(ii-1);
                vertices.coords((count-1)*9+1:(count-1)*9+9,2) = unit_cell_vertices(:,2) + 2*L*(jj-1);
                vertices.coords((count-1)*9+1:(count-1)*9+9,3) = unit_cell_vertices(:,3);   
            end
        end

        vertices.coords = remove_extra_vertices(vertices.coords);
%         vertices.coords = unique(vertices.coords,'rows','stable');
        
        % Define the faces in a unit cell on the top
        unit_cell_faces(1,:) = [1 2 5 4];
        unit_cell_faces(2,:) = [4 5 8 7];
        unit_cell_faces(3,:) = [2 3 6 5];
        unit_cell_faces(4,:) = [6 5 8 9];
        
        for ii = 1:n
            for jj = 1:m
                if jj == 1
                    node_nums(ii,jj,:) = [1:9]+6*(ii-1);
                elseif jj == 2
                    node_nums(ii,jj,[2 3]) = node_nums(1,jj-1,[2 3]) + 6*n+2 + (ii-1)*4;
                    node_nums(ii,jj,[5 6]) = node_nums(1,jj-1,[5 6]) + 6*n+1 + (ii-1)*4;
                    node_nums(ii,jj,[8 9]) = node_nums(1,jj-1,[8 9]) + 6*n+0 + (ii-1)*4;
                    node_nums(ii,jj,[1 4 7]) = node_nums(ii,jj-1,[3 6 9]);
                else
                    node_nums(ii,jj,:) = node_nums(ii,jj-1,:) + 4*n+2;
                    node_nums(ii,jj,[1 4 7]) = node_nums(ii,jj-1,[3 6 9]);
                    
                end
            end
        end
        
        % Use the unit cell faces to define the vertices which surround
        % every face in the tube
        count = 0;
        for ii = 1:n
            for jj = 1:m            
                for kk = 1:4
                    faces.vertices(4*count + kk,:) = node_nums(ii,jj,unit_cell_faces(kk,:));
                end
                count = count + 1;
            end
        end
end

%% Plot to check
h = figure;
for ii = 1:size(faces.vertices,1)   
    x = vertices.coords(faces.vertices(ii,:),1);
    y = vertices.coords(faces.vertices(ii,:),2);
    z = vertices.coords(faces.vertices(ii,:),3);   
    patch(x,y,z,'b')
    hold on   
     
end

for ii = 1:length(vertices.coords)
    x = vertices.coords(ii,1);
    y = vertices.coords(ii,2);
    z = vertices.coords(ii,3);
    hold on
    scatter3(x,y,z,30,'r') 
end
% Make it look presentable
view(35,30)
axis equal
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.ZAxis.Visible = 'off';
h.Color = [1 1 1];
set(gcf, 'Position', [680, 49, 1037, 948])
h.PaperSize = [25 25];

%% Writing to the file

% Error check
if Model_type > 2
    error('Update filename to include extra case(s)')
end

% Define a descriptive filename
switch Model_type
    case 1
        filename = sprintf('Output\\Miura_tube_n%i_a%i_b%i.FOLD',n,a,b); 
    case 2
        filename = sprintf('Output\\Miura_sheet_n%i_m%i_a%i_b%i.FOLD',n,m,a,b);
end

% Open the file and write { to start
fid = fopen(filename,'w');
fprintf(fid,'{\n');

% Write the contents of the 'file' structure
file_names = fieldnames(file);
for ii = 1:length(file_names)
    command_str = sprintf('var_contents = num2str(file.%s);',file_names{ii});
    eval(command_str)
    fprintf(fid,'  "file_%s": %s,\n',file_names{ii},var_contents);
end

% Write the contents of the 'frame' structure
frame_names = fieldnames(frame);
for ii = 1:length(frame_names)
    command_str = sprintf('var_contents = num2str(frame.%s);',frame_names{ii});
    eval(command_str)
    fprintf(fid,'  "frame_%s": %s,\n',frame_names{ii},var_contents);
end

% Write the vertex locations
fprintf(fid,'  "vertices_coords": [\n');
for ii = 1:size(vertices.coords,1)-1
    fprintf(fid,'    [%g,%g,%g],\n',vertices.coords(ii,1),vertices.coords(ii,2),vertices.coords(ii,3));
end
fprintf(fid,'    [%g,%g,%g]\n',vertices.coords(ii+1,1),vertices.coords(ii+1,2),vertices.coords(ii+1,3));
fprintf(fid,'  ],\n');

% Write the faces
fprintf(fid,'  "faces_vertices": [\n');
for ii = 1:size(faces.vertices,1)-1
    fprintf(fid,'    [%g,%g,%g,%g],\n',faces.vertices(ii,1),faces.vertices(ii,2),faces.vertices(ii,3),faces.vertices(ii,4));
end
fprintf(fid,'    [%g,%g,%g,%g]\n',faces.vertices(ii+1,1),faces.vertices(ii+1,2),faces.vertices(ii+1,3),faces.vertices(ii+1,4));
fprintf(fid,'  ],\n');

fprintf(fid,'}');
[~] = fclose(fid);