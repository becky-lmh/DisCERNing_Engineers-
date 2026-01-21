%% Simulate the detector setup
% Extracts the pion data and reconstructs the tracks to find the best
% configuration

%% Import data
% The table has data including collision number and secondary particle 
collision_tab = importfile2tab("C:\Users\black\OneDrive\matlab\10000_collisions7T.txt", [1, Inf]);
collision_tab450 = importfile2tab("C:\Users\black\OneDrive\matlab\10000_collisions450G.txt", [1, Inf]);
% If printing enabled, it will print the figures
printing = false;

%% Plot the secondary particles
% Print count of seconary particles
plotSecondary(collision_tab,printing);

%% Extract the charged, long-lived particles travelling in direction of beam
% Print number of pions
pionCharge = extractPion(collision_tab);
pion = categorical({'pi+','pi-'});
[~,~,ic] = unique(table2array(pionCharge(:,"Secondary")));
a_counts = accumarray(ic,1);
if printing == true
%     figure()
%     bar(pion,a_counts);
%     grid on
%     set(gca,'FontSize',10)
%     xlabel("Pions")
%     ylabel("Number of particles produced")
    %title("Number of charged pions produced in 7TeV gas-beam collision")
end

pionCharge450 = extractPion(collision_tab450);

%% Plot number of pions per vertex
printVertex(pionCharge, printing);

%% Plot momentum distrubution of pions
% printMom(pionCharge,printing);
%printMom(pionCharge450,0);

if printing == true
    figure()
    histogram(table2array(pionCharge(:,6)),20);
    set(gca,'yscale','log')
    grid on
    xlabel("Momentum / MeV")
    ylabel("Number of charged pions")
    %title("Momentum histogram for charged pions")
    
%     hold on
%     histogram(table2array(pionCharge450(:,6)),25);
%     set(gca,'yscale','log')
%     grid on
%     xlabel("Momentum / MeV")
%     ylabel("Number of charged pions")
end
%% Take pions from collision and reconstruct the track at a random gaussian position
[C,~,ic] = unique(table2array(pionCharge(:,"collisionNum")));
total_collisions = max(table2array(pionCharge(:,"collisionNum")));                             % Selected in Geant4 code
count_start = 1;
count_end = 0;
[particles_length, ~] = size(table2array(pionCharge(:,1)));
xPos = zeros(particles_length,2);
yPos = zeros(particles_length,2);
zPos = zeros(particles_length,2);
skip = 0;
massEnergy = pionCharge.mass_energy_MeV;

% assume position has bivariate gaussian distribution of 0.5 mm
mu = [0 0];
variance = 0.005;
%Sigma = eye(2)*variance^2;
Sigma = [0.005, 0; 0, 0.001];
rng('default') % For reproducibility
r = mvnrnd(mu,Sigma,particles_length);
xPos(:,1) = r(:,1);
yPos(:,1) = r(:,2);

if printing == true
%     figure()
%     grid on
%     c = linspace(1,10,length(xPos(:,1)));
%     scatter(xPos(:,1),yPos(:,1),[],c,'x')
%     axis equal
%     set(gca,'FontSize',10)
%     xlabel("x (m)")
%     ylabel("y (m)")

    figure()
    set(gca,'FontSize',10)
    histogram2(xPos(:,1),yPos(:,1))
end


%%
% Assume a gaussian distribution of z with mean=-0.2mm with width 0.2 mm
 zPos(:,1) = normrnd(-0.0002,0.0002,[particles_length,1]);
% zPos(:,1) = normrnd(-0.5,0.5,[particles_length,1]);

for i = 1:total_collisions
    try
        collision_num = C(i);
    catch
        %fprintf("No pions at collision %.0f, skipped \n",i)
        skip = skip +1;
        continue
    end

    vertex = pionCharge(ismember(pionCharge.collisionNum,collision_num),:);
    [vertices, ~] = size(vertex); 
    xDirection = table2array(vertex(:,3));
    yDirection = table2array(vertex(:,4));
    zDirection = table2array(vertex(:,5));
    
    distance = ones(vertices,1)*1.5;
    count_end = count_end + vertices;

    xPos(count_start:count_end,2) = xPos(count_start:count_end,1) + xDirection.*distance;
    yPos(count_start:count_end,2) = yPos(count_start:count_end,1) + yDirection.*distance;
    zPos(count_start:count_end,2) = zPos(count_start:count_end,1) + zDirection.*distance;
    count_start = count_start + vertices;
end

%% 450 GeV
[C,ia,ic] = unique(table2array(pionCharge450(:,"collisionNum")));
total_collisions = max(table2array(pionCharge450(:,"collisionNum")));                             % Selected in Geant4 code
count_start = 1;
count_end = 0;
[particles_length, ~] = size(table2array(pionCharge450(:,1)));
xPos450 = zeros(particles_length,2);
yPos450 = zeros(particles_length,2);
zPos450 = zeros(particles_length,2);
skip = 0;
massEnergy450 = pionCharge450.mass_energy_MeV;

% assume position has bivariate gaussian distribution of 0.5 mm
rng('default') % For reproducibility
r = mvnrnd(mu,Sigma,particles_length);
xPos450(:,1) = r(:,1);
yPos450(:,1) = r(:,2);

% Assume a gaussian distribution of z with mean=-0.2mm with width 0.2 mm
 zPos450(:,1) = normrnd(-0.0002,0.0002,[particles_length,1]);
% zPos(:,1) = normrnd(-0.5,0.5,[particles_length,1]);

for i = 1:total_collisions
    try
        collision_num = C(i);
    catch
        %fprintf("No pions at collision %.0f, skipped \n",i)
        skip = skip +1;
        continue
    end

    vertex = pionCharge450(ismember(pionCharge450.collisionNum,collision_num),:);
    [vertices, ~] = size(vertex); 
    xDirection = table2array(vertex(:,3));
    yDirection = table2array(vertex(:,4));
    zDirection = table2array(vertex(:,5));
    
    distance = ones(vertices,1)*1.5;
    count_end = count_end + vertices;

    xPos450(count_start:count_end,2) = xPos450(count_start:count_end,1) + xDirection.*distance;
    yPos450(count_start:count_end,2) = yPos450(count_start:count_end,1) + yDirection.*distance;
    zPos450(count_start:count_end,2) = zPos450(count_start:count_end,1) + zDirection.*distance;
    count_start = count_start + vertices;
end

%% Show graph of constructed tracks 
printTracks(xPos,yPos,zPos,printing);

%% Find the pions positions at the 3 detectors 
% Pions in beam pipe will not be detected by SciFi
% Beam radius - 22.5mm
beamRad = 0.0225;
sciFiRad = 0.132;

total_collisions = max(table2array(pionCharge(:,1)));
[particles_length, ~] = size(table2array(pionCharge(:,1)));

inc = 0.05;
totalSpaces = 1.0/inc;
newXPos = zeros(particles_length,totalSpaces+1);
newYPos = zeros(particles_length,totalSpaces+1);
newZPos = zeros(particles_length,totalSpaces+1);
newXPos(:,1) = xPos(:,1);
newYPos(:,1) = yPos(:,1);
newZPos(:,1) = zPos(:,1);

for i = 2:totalSpaces+1
    detPos = inc*(i-1) * ones(particles_length,1);
    % Track line equation is [x0,y0,z0] + t.n where n is vector 
    % Solve for t and sub back in to find values of x and y at detector
    t = (detPos-zPos(:,1))./table2array(pionCharge(:,"p_z"));
    newXPos(:,i) = xPos(:,1) + t.*table2array(pionCharge(:,"p_x"));
    newYPos(:,i) = yPos(:,1) + t.*table2array(pionCharge(:,"p_y"));
    newZPos(:,i) = detPos;
end

% % Remove pions in undetectable beam region
radialPos = (newXPos.*newXPos + newYPos.*newYPos).^0.5;
newXPos(radialPos < beamRad) = -inf;
newYPos(radialPos < beamRad) = -inf;
newXPos(radialPos > sciFiRad) = -inf;
newYPos(radialPos > sciFiRad) = -inf;

allValues = newXPos==-inf;
allValues = allValues == zeros(size(allValues));
numPions = sum(allValues==1);
[a, b] = max(numPions);
maxPionLength = (b-1)*inc;

% 450 GeV
total_collisions = max(table2array(pionCharge450(:,1)));
[particles_length, ~] = size(table2array(pionCharge450(:,1)));

inc = 0.05;
totalSpaces = 1.0/inc;
newXPos450 = zeros(particles_length,totalSpaces+1);
newYPos450 = zeros(particles_length,totalSpaces+1);
newZPos450 = zeros(particles_length,totalSpaces+1);
newXPos450(:,1) = xPos450(:,1);
newYPos450(:,1) = yPos450(:,1);
newZPos450(:,1) = zPos450(:,1);

for i = 2:totalSpaces+1
    detPos = inc*(i-1) * ones(particles_length,1);
    % Track line equation is [x0,y0,z0] + t.n where n is vector 
    % Solve for t and sub back in to find values of x and y at detector
    t = (detPos-zPos450(:,1))./table2array(pionCharge450(:,"p_z"));
    newXPos450(:,i) = xPos450(:,1) + t.*table2array(pionCharge450(:,"p_x"));
    newYPos450(:,i) = yPos450(:,1) + t.*table2array(pionCharge450(:,"p_y"));
    newZPos450(:,i) = detPos;
end

% % Remove pions in undetectable beam region
radialPos = (newXPos450.*newXPos450 + newYPos450.*newYPos450).^0.5;
newXPos450(radialPos < beamRad) = -inf;
newYPos450(radialPos < beamRad) = -inf;
newXPos450(radialPos > sciFiRad) = -inf;
newYPos450(radialPos > sciFiRad) = -inf;

allValues = newXPos450==-inf;
allValues = allValues == zeros(size(allValues));
numPions450 = sum(allValues==1);
[a, b] = max(numPions);
maxPionLength = (b-1)*inc;

%% Show graph of pions at each detector

printNumPions(newZPos,numPions, inc, printing);
printNumPions(newZPos450,numPions450, inc, printing);

% Maximum point is at 0.4 m so this is the position of first detector
%% Add smearing to the x and y positions
% The SciFi screens have resolution of 40um
[smearedX, smearedY] = smear(newXPos,newYPos);

% Remove undetectable pions
% Remove data from pions before 0.45 m
detectorX = smearedX(:,10:end);
detectorY = smearedY(:,10:end);
detectorZ = newZPos(1,10:end);
massEnergy_ = massEnergy;

% Remove pions that are not detected at 0.45 m
remove1 = detectorX(:,1)==-inf;
detectorX(remove1,:)=[];
detectorY(remove1,:)=[];
massEnergy_(remove1)=[];
remove2 = detectorY(:,1)==-inf;
detectorX(remove2,:)=[];
detectorY(remove2,:)=[];
massEnergy_(remove2)=[];

% 450 GeV
[smearedX450, smearedY450] = smear(newXPos450,newYPos450);

% Remove undetectable pions
% Remove data from pions before 0.45 m
detectorX450 = smearedX450(:,10:end);
detectorY450 = smearedY450(:,10:end);
detectorZ450 = newZPos450(1,10:end);
massEnergy_450 = massEnergy450;

% Remove pions that are not detected at 0.45 m
remove1 = detectorX450(:,1)==-inf;
detectorX450(remove1,:)=[];
detectorY450(remove1,:)=[];
massEnergy_450(remove1)=[];
remove2 = detectorY450(:,1)==-inf;
detectorX450(remove2,:)=[];
detectorY450(remove2,:)=[];
massEnergy_450(remove2)=[];

%% Use the position at the 3 detectors to reconstruct the vertex
% Assume first detector at 0.45 m
% Add MS error as random gaussian noise sigma(mm) = 2.85e-4 x 
% Do this for all possible choices and see the close my simulated data is
% to the true data's beam position (0,0) and beam variance 

[numPionsLeft, positions] = size(detectorX);
meanX1 = zeros(66,1);
meanY1 = zeros(66,1);
stdX1 = zeros(66,1);
stdY1 = zeros(66,1);
z2 = zeros(66,1);
z3 = zeros(66,1);

meanX2 = zeros(66,1);
meanY2 = zeros(66,1);
stdX2 = zeros(66,1);
stdY2 = zeros(66,1);

step = 1;
z1=0.4;
clusterDis = 0.000030;

% Using relativistic energy equation to calculate velocity m/s
velocity = sqrt(1 - (139.57039./massEnergy).^2).*3e8;
time = (1./velocity)*detectorZ;
for i = 1:(positions-1)
    for j = 2:positions
        if i==j || i>j
            continue
        else
            % Extract data for 3 detectors
            curDetX = zeros(numPionsLeft,3);
            curDetY = zeros(numPionsLeft,3);
            curDetX(:,1) = detectorX(:,1);
            curDetY(:,1) = detectorY(:,1);
            curDetX(:,2)=detectorX(:,i);
            curDetY(:,2)=detectorY(:,i);
            curDetX(:,3)=detectorX(:,j);
            curDetY(:,3)=detectorY(:,j);
            time = (1./velocity)*detectorZ;
            
            z2(step) = detectorZ(i);
            z3(step) = detectorZ(j);
    
            % Remove pions that aren't detected in all 3 detectors
            remove = curDetX(:,2)==-inf;
            curDetX(remove,:)=[];
            curDetY(remove,:)=[];
            time(remove,:)=[];
            remove = curDetX(:,3)==-inf;
            curDetX(remove,:)=[];
            curDetY(remove,:)=[];
            time(remove,:)=[];
            remove = curDetY(:,2)==-inf;
            curDetX(remove,:)=[];
            curDetY(remove,:)=[];
            time(remove,:)=[];
            remove = curDetY(:,3)==-inf;
            curDetX(remove,:)=[];
            curDetY(remove,:)=[];
            time(remove,:)=[];
%     
%             % Add shifts from multiple scattering
%             MS1 = 0.00028571428 * z1;
%             MS2 = 0.00028571428 * z2;
%             MS3 = 0.00028571428 * z3;
%             MS1_ = normrnd(1,1,size(curDetX(:,1))).*MS1;
%             MS2_ = normrnd(1,1,size(curDetX(:,1))).*MS2;
%             MS3_ = normrnd(1,1,size(curDetX(:,1))).*MS3;
%             MSX = [MS1_ MS2_ MS3_];
%     
%             MS1_ = normrnd(1,1,size(curDetX(:,1))).*MS1;
%             MS2_ = normrnd(1,1,size(curDetX(:,1))).*MS2;
%             MS3_ = normrnd(1,1,size(curDetX(:,1))).*MS3;
%             MSY = [MS1_ MS2_ MS3_];
%     
%             curDetX = curDetX + MSX;
%             curDetY = curDetY + MSY;
    
            % Linear Regression
            dt = [time(:,i)-time(:,1),time(:,j)-time(:,1)];
            curDetZ = [z1, z2(step), z3(step)];
           
            [PVx, PVy] = linearReg(curDetX, curDetY, curDetZ);
            figure()
            c = linspace(1,10,length(PVy));
            scatter(PVx,PVy,[],c,'x')
            axis equal
            set(gca,'FontSize',10)
            xlabel("Detector Distance 2 (m)")
            ylabel("Detector Distance 3 (m)")
            meanX1(step,1)=mean(PVx);
            meanY1(step,1)=mean(PVy);
            stdX1(step,1)=std(PVx);
            stdY1(step,1)=std(PVy);
            
            
%             [PVx, PVy] = kalmanFilter(curDetX, curDetY, curDetZ, dt);
%             random_x = PVx(randperm(size(PVx, 1)), :);
%             random_y = PVy(randperm(size(PVy, 1)), :);
% 
%             [clustersCentroids,~,~] = clusterXYpoints([random_x(1:1000,:), random_y(1:1000,:)],clusterDis,2,'centroid');
%             IP = clustersCentroids(:,1).*clustersCentroids(:,2);
%             theta = atan(clustersCentroids(:,2)./clustersCentroids(:,1));
% %             figure()
% %             scatter(theta,IP,'green');
% %             hold on
% %             x0 = 0.0000001;
% %             y0 = 0.000000004;
% %             y = x0.*sin(theta) - y0.*cos(theta);
% %             scatter(theta,y,'blue')
% %             ylabel("Impact Parameter (m)")
% %             xlabel("Theta")
% %             set(gca,'FontSize',10)
% %             grid on
% %             hold off
%             meanX1(step,1)=mean(clustersCentroids(:,1));
%             meanY1(step,1)=mean(clustersCentroids(:,2));
%             stdX1(step,1)=std(clustersCentroids(:,1));
%             stdY1(step,1)=std(clustersCentroids(:,2));
% 
%             [clustersCentroids,~,~] = clusterXYpoints([random_x(end-1000:end,:), random_y(end-1000:end,:)],clusterDis,2,'centroid');
%             meanX2(step,1)=mean(clustersCentroids(:,1));
%             meanY2(step,1)=mean(clustersCentroids(:,2));
%             stdX2(step,1)=std(clustersCentroids(:,1));
%             stdY2(step,1)=std(clustersCentroids(:,2));
            step = step + 1;
        end
    end
end

%%
% 450 GeV
[numPionsLeft, positions] = size(detectorX450);
meanX1450 = zeros(66,1);
meanY1450 = zeros(66,1);
stdX1450 = zeros(66,1);
stdY1450 = zeros(66,1);

meanX2450 = zeros(66,1);
meanY2450 = zeros(66,1);
stdX2450 = zeros(66,1);
stdY2450 = zeros(66,1);

step = 1;
z1450=0.4;

% Using relativistic energy equation to calculate velocity m/s
velocity450 = sqrt(1 - (139.57039./massEnergy450).^2).*3e8;
time450 = (1./velocity450)*detectorZ450;
for i = 1:(positions-1)
    for j = 2:positions
        if i==j || i>j
            continue
        else
            % Extract data for 3 detectors
            curDetX450 = zeros(numPionsLeft,3);
            curDetY450 = zeros(numPionsLeft,3);
            curDetX450(:,1) = detectorX450(:,1);
            curDetY450(:,1) = detectorY450(:,1);
            curDetX450(:,2)=detectorX450(:,i);
            curDetY450(:,2)=detectorY450(:,i);
            curDetX450(:,3)=detectorX450(:,j);
            curDetY450(:,3)=detectorY450(:,j);
            time450 = (1./velocity450)*detectorZ450;
            
            z2450 = detectorZ450(i);
            z3450 = detectorZ450(j);
    
            % Remove pions that aren't detected in all 3 detectors
            remove = curDetX450(:,2)==-inf;
            curDetX450(remove,:)=[];
            curDetY450(remove,:)=[];
            time450(remove,:)=[];
            remove = curDetX450(:,3)==-inf;
            curDetX450(remove,:)=[];
            curDetY450(remove,:)=[];
            time450(remove,:)=[];
            remove = curDetY450(:,2)==-inf;
            curDetX450(remove,:)=[];
            curDetY450(remove,:)=[];
            time450(remove,:)=[];
            remove = curDetY450(:,3)==-inf;
            curDetX450(remove,:)=[];
            curDetY450(remove,:)=[];
            time450(remove,:)=[];
    
%             % Add shifts from multiple scattering
%             MS1 = 0.00028571428 * z1450;
%             MS2 = 0.00028571428 * z2450;
%             MS3 = 0.00028571428 * z3450;
%             MS1_ = normrnd(1,1,size(curDetX450(:,1))).*MS1;
%             MS2_ = normrnd(1,1,size(curDetX450(:,1))).*MS2;
%             MS3_ = normrnd(1,1,size(curDetX450(:,1))).*MS3;
%             MSX = [MS1_ MS2_ MS3_];
%     
%             MS1_ = normrnd(1,1,size(curDetX450(:,1))).*MS1;
%             MS2_ = normrnd(1,1,size(curDetX450(:,1))).*MS2;
%             MS3_ = normrnd(1,1,size(curDetX450(:,1))).*MS3;
%             MSY = [MS1_ MS2_ MS3_];
%     
%             curDetX450 = curDetX450 + MSX;
%             curDetY450 = curDetY450 + MSY;
    
            % Linear regression
            dt450 = [time450(:,i)-time450(:,1),time450(:,j)-time450(:,1)];
            curDetZ450 = [z1450, z2450, z3450];


            [PVx450, PVy450] = linearReg(curDetX, curDetY, curDetZ);
%             [PVx450, PVy450] = kalmanFilter(curDetX450, curDetY450, curDetZ450, dt450);
%             random_x450 = PVx450(randperm(size(PVx450, 1)), :);
%             random_y450 = PVy450(randperm(size(PVy450, 1)), :);
% 
%             [clustersCentroids450,~,~] = clusterXYpoints([random_x450(1:1000,:), random_y450(1:1000,:)],clusterDis,2,'centroid');
%             meanX1450(step,1)=mean(clustersCentroids450(:,1));
%             meanY1450(step,1)=mean(clustersCentroids450(:,2));
%             stdX1450(step,1)=std(clustersCentroids450(:,1));
%             stdY1450(step,1)=std(clustersCentroids450(:,2));
% 
%             [clustersCentroids,~,~] = clusterXYpoints([random_x450(end-1000:end,:), random_y450(end-1000:end,:)],clusterDis,2,'centroid');
%             meanX2450(step,1)=mean(clustersCentroids(:,1));
%             meanY2450(step,1)=mean(clustersCentroids(:,2));
%             stdX2450(step,1)=std(clustersCentroids(:,1));
%             stdY2450(step,1)=std(clustersCentroids(:,2));
            meanX1450(step,1)=mean(PVx450);
            meanY1450(step,1)=mean(PVy450);
            stdX1450(step,1)=std(PVx450);
            stdY1450(step,1)=std(PVy450);
            step = step + 1;
        end
    end
end
%% Find precision of beam sigma_hit
z2real = z2;
z3real = z3;

sigmaHitX = abs(meanX1 - meanX2);
sigmaHitY = abs(meanY1 - meanY2);
sigmaStdX = abs(stdX1-stdX2);
sigmaStdY = abs(stdY1 - stdY2);

T = sigmaHitY+sigmaHitX + sigmaStdY + sigmaStdX;
%T = sigmaHitY+sigmaHitX;
% T=smoothdata(T);
% T=smoothdata(T);
% T=smoothdata(T);
% figure()
% plotMesh(z2,z3,T)
% % title("Resolution of Primary Vertex (7 TeV)")
% xlabel("Detector Distance 2 (m)")
% ylabel("Detector Distance 2 (m)")
% set(gca,'FontSize',10)
% title(colorbar,"Pre")

% % 3D surface plot
figure()
T1=delaunay(z2,z3);
hObj = trisurf(T1,z2,z3,T);
colorbar
colormap jet;
brighten(0.5);
xlabel("Detector Distance 2 (m)")
ylabel("Detector Distance 2 (m)")
set(gca,'FontSize',10)
title(colorbar,"Pre")

%% 450 GeV
sigmaHitX = abs(meanX1450 - meanX2450);
sigmaHitY = abs(meanY1450 - meanY2450);
sigmaStdX = abs(stdX1450 - stdX2450);
sigmaStdY = abs(stdY1450 - stdY2450);
T450 = sigmaHitY+sigmaHitX + sigmaStdY + sigmaStdX;
% T450 = sigmaHitY+sigmaHitX;
% T450 = smoothdata(T450);
% T450 = smoothdata(T450);
% T450 = smoothdata(T450);
% figure()
% plotMesh(z2,z3,T450)
% %title("Resolution of Primary Vertex (7 TeV)")
% xlabel("Detector Distance 2 (m)")
% ylabel("Detector Distance 2 (m)")
% set(gca,'FontSize',10)
% title(colorbar,'Pre')

% % 3D surface plot
figure()
T1=delaunay(z2,z3);
hObj = trisurf(T1,z2,z3,T450);
colorbar
colormap jet;
brighten(0.5);
xlabel("Detector Distance 2 (m)")
ylabel("Detector Distance 2 (m)")
set(gca,'FontSize',10)
title(colorbar,"Pre")

%% Extrapolation error

% extrapolation error
% sigma_ext = rescale((z2.*z2 + z3.*z3).*0.0001./((z3-z2).*(z3-z2)));
sigma_ext = ( ((z2.*z2 + z3.*z3).*T.^2) ./ ((z3-z2).*(z3-z2)) ) .^0.5;
sigma_ext450 = ((z2.*z2 + z3.*z3).*T450.^2./((z3-z2).*(z3-z2))).^0.5;
% figure()
% plotMesh(z2,z3,sigma_ext)
% xlabel("Detector Distance 2 (m)")
% ylabel("Detector Distance 2 (m)")
% set(gca,'FontSize',10)
% title(colorbar, 'Ext')

% % 3D surface plot
figure()
T1=delaunay(z2,z3);
hObj = trisurf(T1,z2,z3,sigma_ext);
xlabel('Detector Distance 2 (m)')
ylabel('Detector Distance 3 (m)')
zlabel('M')
colorbar
colormap jet;
set(gca,'FontSize',10)
title(colorbar, 'Ext')
brighten(0.5);

% % 3D surface plot
figure()
T1=delaunay(z2,z3);
hObj = trisurf(T1,z2,z3,sigma_ext450);
xlabel('Detector Distance 2 (m)')
ylabel('Detector Distance 3 (m)')
zlabel('M')
colorbar
colormap jet;
set(gca,'FontSize',10)
title(colorbar, 'Ext')
brighten(0.5);

% figure()
% plotMesh(z2,z3,sigma_ext450)
% xlabel("Detector Distance 2 (m)")
% ylabel("Detector Distance 2 (m)")
% set(gca,'FontSize',10)
% title(colorbar, 'Ext')


% 
% figure()
% T=delaunay(z2,z3);
% hObj = trisurf(T,z2,z3,sigma_ext450);
% xlabel('Detector Distance 2 (m)')
% ylabel('Detector Distance 3 (m)')
% zlabel('M')
% colorbar
% colormap jet;
% title(colorbar, 'Extrapolation Error')
% brighten(0.5);

% figure()
% T=delaunay(z2,z3);
% hObj = trisurf(T,z2,z3,sigma_ext);
% xlabel('Detector Distance 2 (m)')
% ylabel('Detector Distance 3 (m)')
% zlabel('Extrapolation Error')
% colorbar
% colormap jet;
% brighten(0.5);

%% Find Accuracy of beam
realMean = 0;
realStd = 5e-4;

% find error
meanXError1 = abs(meanX1 - realMean);
meanYError1 = abs(meanY1 - realMean);
stdXError1 = abs(stdX1 - realStd);
stdYError1 = abs(stdY1 - realStd);
meanError = normalize(1./(meanYError1+meanXError1),'range');
stdError = normalize(1./(stdXError1 + stdYError1),'range');
score1 = (meanError + stdError);

% figure()
% plot(score1)
% grid on
% %title("Accuracy of each detector configuration (7TeV)")
% xlabel("Detector Configuration")
% ylabel("Score")

%Remove outlier
% z2 = detectorDistances.z2;
% z3 = detectorDistances.z3;
% z2(29)=[];
% z3(29)=[];
% score1(29)=[];
% z2(25)=[];
% z3(25)=[];
% score1(25)=[];
% z2(23)=[];
% z3(23)=[];
% score1(23)=[];

% % Heat map
% figure()
% score1 = smoothdata(score1);
% score1 = smoothdata(score1);
% plotMesh(z2,z3,score1)
% xlabel('Detector Distance 2 (m)');
% ylabel('Detector Distance 3 (m)');
% %title('Merit Index')
% brighten(-0.5);
% set(gca,'FontSize',10)
% title(colorbar,'Accuracy')

% % 3D surface plot
figure()
T1=delaunay(z2,z3);
hObj = trisurf(T1,z2,z3,score1);
xlabel('Detector Distance 2 (m)')
ylabel('Detector Distance 3 (m)')
zlabel('M')
colorbar
colormap jet;
set(gca,'FontSize',10)
title(colorbar,'Accuracy')
brighten(0.5);

% 450 GeV
% find error
meanXError1 = abs(meanX1450 - realMean);
meanYError1 = abs(meanY1450 - realMean);
stdXError1 = abs(stdX1450 - realStd);
stdYError1 = abs(stdY1450 - realStd);
meanError450 = normalize(1./(meanYError1+meanXError1),'range');
stdError450 = normalize(1./(stdXError1 + stdYError1),'range');
score1450 = (meanError450 + stdError450);

% % 3D surface plot
figure()
T1=delaunay(z2,z3);
hObj = trisurf(T1,z2,z3,score1450);
xlabel('Detector Distance 2 (m)')
ylabel('Detector Distance 3 (m)')
zlabel('M')
colorbar
colormap jet;
set(gca,'FontSize',10)
title(colorbar,'Accuracy')
brighten(0.5);

% figure()
% plot(score1450)
% set(gca,'FontSize',10)
% grid on
% %title("Accuracy of each detector configuration (7TeV)")
% xlabel("Detector Configuration")
% ylabel("Score")

% Remove outlier
% z2 = detectorDistances.z2;
% z3 = detectorDistances.z3;
% z2(28)=[];
% z3(28)=[];
% score1450(28)=[];
% z2(1)=[];
% z3(1)=[];
% score1450(1)=[];

% % Heat map
% figure()
% score1450 = smoothdata(score1450);
% plotMesh(z2,z3,score1450)
% xlabel('Detector Distance 2 (m)');
% ylabel('Detector Distance 3 (m)');
% %title('Merit Index')
% set(gca,'FontSize',10)
% brighten(-0.5);
% title(colorbar,'Accuracy')

%% Total Merit Index
sigma_IP = normalize(1./sigma_ext,'range');
sigma_IP450 = normalize(1./sigma_ext450,'range');
precise = normalize(1./T,'range');
precise450 = normalize(1./T450,'range');
merit = score1 + sigma_IP + precise;
merit450 = score1450 + sigma_IP450 + precise450;

% Remove outlier
% z2 = z2real;
% z3 = z3real;
% z2(30)=[];
% z3(30)=[];
% merit(30)=[];
% z2(26)=[];
% z3(26)=[];
% merit(26)=[];
% % z2(1)=[];
% % z3(1)=[];
% % merit(1)=[];

% figure()
% plot(merit)
% xlabel("Detector Configuration")
% ylabel("Merit Index")
% set(gca,'FontSize',10)

figure()
merit = smoothdata(merit);
merit = smoothdata(merit);
merit = smoothdata(merit);
T1=delaunay(z2,z3);
hObj = trisurf(T1,z2,z3,merit);
xlabel('Detector Distance 2 (m)')
ylabel('Detector Distance 3 (m)')
colorbar
colormap jet;
brighten(0.5);
title(colorbar,'M')
set(gca,'FontSize',10)

% 450 GeV
%Remove outlier
% z2 = z2real;
% z3 = z3real;
% % z2(10)=[];
% % z3(10)=[];
% % merit450(10)=[];
% % z2(21)=[];
% % z3(21)=[];
% % merit450(21)=[];

% figure()
% plot(merit450)
% xlabel("Detector Configuration")
% ylabel("Merit Index")
% set(gca,'FontSize',10)

figure()
merit450 = smoothdata(merit450);
merit450 = smoothdata(merit450);
% merit450 = smoothdata(merit450);
T1=delaunay(z2,z3);
hObj = trisurf(T1,z2,z3,merit450);
xlabel('Detector Distance 2 (m)')
ylabel('Detector Distance 3 (m)')
colorbar
title(colorbar,'M')
colormap jet;
brighten(0.5);
set(gca,'FontSize',10)

% figure()
% plotMesh(z2,z3,merit450)
% xlabel('Detector Distance 2 (m)');
% ylabel('Detector Distance 3 (m)');
% title(colorbar,'M')
% set(gca,'FontSize',10)
%% Cool Plots
% % Rainbox line plot
% figure()
% hold on
% [Z2,Z3]=meshgrid(z2,z3);
% Z=griddata(z2,z3,score1,Z2,Z3);
% contour(Z2,Z3,Z)

% % Heat map
% figure()
% plotMesh(z2,z3,score1)
% xlabel('Detector Distance 2 (m)');
% ylabel('Detector Distance 3 (m)');
% %title('Merit Index')
% brighten(-0.5);

% % 3D surface plot
% figure()
% T=delaunay(z2,z3);
% hObj = trisurf(T,z2,z3,score1);
% xlabel('Detector Distance 2 (m)')
% ylabel('Detector Distance 3 (m)')
% zlabel('M')
% title('Merit Index')
% colorbar
% colormap jet;
% brighten(0.5);

% % 3D scatter plot
% figure()
% scatter3(z2,z3,score1)
% %scatter3(D(1:70,1),D(1:70,2),D(1:70,4))
% 
% % Colour Blob plot
% figure()
% s = scatter(z2,z3,[],score1,'filled','SizeData',100);
% colorbar
% xlabel('Detector Distance 2 (m)')
% ylabel('Detector Distance 2 (m)')
% colormap jet;
% brighten(0.5);
%% Redo it but start making an imperfect beam and see how my chosen dimensions manage
% 

