clear
clc

% Fill in XXX where the files are located.

cfdPost_files = dir('XXX'); % files from CFD Post used as a check
TecPlot_files = dir('XXX'); % files from TECPLOT
faces_file = dir('XXX'); % files from CFD Post that contain face connectivity
% 1x 2y 3z 4WSS 5WSSX 6WSSY 7WSSZ 8WSS-DIV 9WSSX_DDX 10WSSX_DDY 11WSSX_DDZ 12WSSY_DDX
% 13WSSY_DDY 14WSSY_DDZ 15WSSZ_DDX 16WSSZ_DDY 17WSSZ_DDZ 18pressure
%% Sort based on CFD Post
% This will sort the first tecplot file according to the CFD post file
% This is important for importing back into CFD Post

tecPlot = readtable('XXX'); % fill in name of TecPlot file. 
tecPlot = table2array(tecPlot);

row = find(isnan(tecPlot(:,1)));

t_size = [];

for ii = 3:3:size(row,1)
    temp = (row(ii)-row(ii-1))-1;
    t_size = [t_size;temp];
end

clear temp
max_array = max(t_size);

tecPlot_files = [];
counter = 1;
for ii = 3:3:size(row,1)
    temp = tecPlot(row(ii-1)+1:row(ii)-1,:);
    %temp = padarray(temp,pads,'replicate','post');
    temp(:,1:3) = temp(:,1:3)*1000;
    tecplot(:,:,counter) = temp;
    counter = counter+1;
end

raw = importdata(cfdPost_files(1).name,',');
cfdPost = raw.data;

%% Tracking and Identifying Section
% This will identify a node, and track its position throughout cycle
f = waitbar(0,'1','Name','Tracking Initial Nodes',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);
set(f,'Units','Pixels','Position',[1500 100 360 125]);


m = size(tecplot,1);
tt = size(tecplot,3);
tracked = zeros(size(tecplot(:,:,1)));
tracked(:,1:3) = cfdPost(:,2:4);
nodes = [0:max_array-1]';
tracked = cat(2,nodes,tracked);
tStart = tic;
tecplot_t = tecplot(:,:,1);
for n = 1:m
    if getappdata(f,'canceling')
        break
    end
    check = abs(tracked(n,2:4)-tecplot_t(:,1:3));
    residual = sum(check,2);
    [M,I] = min(residual);
    tracked(n,2:19) = tecplot_t(I,:);        
    tecplot_t(I,:) = []; 
    waitbar(n/m, f, sprintf('Progress: %4.1f %% (i = %i out of %i)',(n/m)*100,n,m));
end
delete(f)
%%
f = waitbar(0,'1','Name','Tracking Nodes',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);
set(f,'Units','Pixels','Position',[1500 100 360 125]);

for t = 2:tt
        tecplot_t = tecplot(:,:,t);
    for n=1:m
        tic
        if getappdata(f,'canceling')
            break
        end
        check = abs(tracked(n,2:4,t-1)-tecplot_t(:,1:3));
        residual = sum(check,2);
        [M,I] = min(residual);
        tracked(n,1,t) = tracked(n,1,1);
        tracked(n,2:19,t) = tecplot_t(I,:);        
        tecplot_t(I,:) = [];
        waitbar(n/m, f, sprintf('Current t: %i out of %i\nProgress: %4.1f %% (i = %i out of %i)',t,tt,(n/m)*100,n,m));
        T = toc;
    end
    tEnd(t-1) = toc(tStart);
end

delete(f)

%% Time-Average WSS Calculations
% This section will calculate TSM and OSI

% Calculates the sum of the WSS magnitudes for all time steps
sumWSS = sum(tracked(:,5,:),3);
% Calculate the sum of the vector components of WSS for all times step 
sumWSSx = sum(tracked(:,6,:),3);
sumWSSy = sum(tracked(:,7,:),3);
sumWSSz = sum(tracked(:,8,:),3);
% Calculates the sum of the Pressure
sumPress = sum(tracked(:,19,:),3);
% tot_time = tStep(end);
% Calculates the time-averaged WSS magnitude
n = size(tecplot,3);
TSM = sumWSS/n;
% Calculates the time-average of each WSS vector component 
AvWSSx = sumWSSx/n;
AvWSSy = sumWSSy/n;
AvWSSz = sumWSSz/n;
APress = sumPress/n;
AvWSS = [AvWSSx AvWSSy AvWSSz];
% Calculates the magnitude of the time-averaged WSS vector
AWSSV = sqrt(AvWSSx.^2 + AvWSSy.^2 + AvWSSz.^2);
% Determines OSI
OSI = .5*(1-(AWSSV./TSM));

%% WSS Topo Skeleton
% Calculates the time-sum of the WSS divergence and Jacobian magnitudes
sumWSSdiv = sum(tracked(:,9,:),3);
sumJacWSS = sum(tracked(:,10:18,:),3);

% Time-averaged quantities
AWSSdiv = sumWSSdiv/n;
AJacWSS = sumJacWSS/n;

%% Topological shear variation index
% This will calcualte the TSVI for later use (for more info, see
% "Morbiducci et al. - 2020 - Wall Shear Stress Topological Skeleton
% Independently Predicts Long-Term Restenosis After Carotid Bifurcation
% Endartecrectomy"

temp_TSVI = (tracked(:,9,:)-AWSSdiv).^2;
int_TSVI= sum(temp_TSVI(:,1,:),3)/n;
TSVI = sqrt(int_TSVI);

%% Point Cloud Information
pointcld = [tracked(:,1:4,1) TSM OSI AWSSdiv APress TSVI];
sprintf = ('The solution inforation is contained in "pointcld".The order of rows are as follows: 1) node number, 2) x coordinate, 3) y cooridnate, 4) z coordinate, 5) TSM, 6) OSI, 7) WSS divergence, 8) Pressure, and 9) TSVI')

%% Finding Fixed Points
% Need info on face connectivity
% Match face connectivity info to identify triangles
face = importdata(faces_file.name,','); % import face data
face = sort(face,2);
counter = 1;
for ii = 1:length(face) % matching face connectivity to node info to create triangles
    i = face(ii,1);j = face(ii,2);k = face(ii,3);
    for jj = 1:size(pointcld,1)
        if pointcld(jj,1) == i
            vertices = cat(2,pointcld(jj,1:4),AvWSS(jj,:),AJacWSS(jj,:));
        end
        if pointcld(jj,1) == j
            vertices2 = cat(2,pointcld(jj,1:4),AvWSS(jj,:),AJacWSS(jj,:));
            vertices = [vertices;vertices2];
        end
        if pointcld(jj,1) == k
            vertices2 = cat(2,pointcld(jj,1:4),AvWSS(jj,:),AJacWSS(jj,:));
            vertices = [vertices;vertices2];
        end
    end
    %% Generate unit normal vector to triangle surface
    % Use cross product
    v1 = [vertices(2,2) - vertices(1,2); vertices(2,3) - vertices(1,3); vertices(2,4) - vertices(1,4)]; 
    v2 = [vertices(3,2) - vertices(1,2); vertices(3,3) - vertices(1,3); vertices(3,4) - vertices(1,4)];
    cross_v = cross(v1,v2);
    mag_v = norm(cross_v);
    unit_v = abs(cross_v/mag_v);

    %% Set up Poincare index aglorithm

    mat1 = [vertices(1,5:7);vertices(2,5:7);unit_v'];
    mat2 = [vertices(2,5:7);vertices(3,5:7);unit_v'];
    mat3 = [vertices(3,5:7);vertices(1,5:7);unit_v'];
    cond1 = det(mat1);
    if abs(cond1)<1e-3 % here, a tolerance is specified to keep only the relevant fixed points
        cond1 = 0;
    end
    cond2 = det(mat2);
    if abs(cond2)<1e-3
        cond2 = 0;
    end
    cond3 = det(mat3);
    if abs(cond3)<1e-3
        cond3 = 0;
    end
    cond = [cond1 cond2 cond3];
    if cond1<0&cond2<0&cond3<0 || cond1>0&cond2>0&cond3>0
        temp_fixedPoints = sign(cond1);
        temp_fixedPoints(1,2:4) = mean(vertices(:,2:4),1);
        temp_fixedPoints2(1,5:13) = mean(vertices(:,8:16),1); % this will find the location/values of the interior point
        jac = temp_fixedPoints2(1,5:13); %information containing the averaged Jacobian matrix
        jac =reshape(jac(:),3,[]).';
        nature = eig(jac)';
        if temp_fixedPoints(1,1) < 0
            fixedPoints(counter,:) = cat(2,nature,temp_fixedPoints);
            counter = counter+1;
            continue
        elseif all(nature<0) || all(nature>0) & isreal(nature)
            fixedPoints(counter,:) = cat(2,nature,temp_fixedPoints);
            counter = counter+1;
            continue
        elseif nature(1)>0&nature(2)>0&nature(3)>0
            fixedPoints(counter,:) = cat(2,nature,temp_fixedPoints);
            counter = counter+1;
            continue
        elseif nature(1)<0&nature(2)<0&nature(3)<0
            fixedPoints(counter,:) = cat(2,nature,temp_fixedPoints);
            counter = counter+1;
            continue
        end
        clear temp_fixedPoints
        clear temp_fixedPoints2
    end
end