%% Imperfect beams
% Varies the proton beam to see how reliable my detector design is


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
    figure()
    bar(pion,a_counts);
    grid on
    xlabel("Pions")
    ylabel("Number of particles produced")
    title("Number of charged pions produced in 7TeV gas-beam collision")
end

pionCharge450 = extractPion(collision_tab450);

%% Plot number of pions per vertex
printVertex(pionCharge, printing);

%% Plot momentum distrubution of pions
printMom(pionCharge,printing);
printMom(pionCharge450,0);
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
%variance = 0.0005;
Sigma = eye(2)*variance^2;
rng('default') % For reproducibility
r = mvnrnd(mu,Sigma,particles_length);
xPos(:,1) = r(:,1);
yPos(:,1) = r(:,2);

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


% 450 GeV
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
mu = [0 0];
Sigma = eye(2)*variance^2;
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

% Maximum point is at 0.5 m so this is the position of first detector
%% Add smearing to the x and y positions
% The SciFi screens have resolution of 100um
[smearedX, smearedY] = smear(newXPos,newYPos);

% Remove undetectable pions
% Remove data from pions before 0.5 m
detectorX = smearedX(:,11:end);
detectorY = smearedY(:,11:end);
detectorZ = newZPos(1,11:end);
massEnergy_ = massEnergy;

% Remove pions that are not detected at 0.5 m
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
% Remove data from pions before 0.5 m
detectorX450 = smearedX450(:,11:end);
detectorY450 = smearedY450(:,11:end);
detectorZ450 = newZPos450(1,11:end);
massEnergy_450 = massEnergy450;

% Remove pions that are not detected at 0.5 m
remove1 = detectorX450(:,1)==-inf;
detectorX450(remove1,:)=[];
detectorY450(remove1,:)=[];
massEnergy_450(remove1)=[];
remove2 = detectorY450(:,1)==-inf;
detectorX450(remove2,:)=[];
detectorY450(remove2,:)=[];
massEnergy_450(remove2)=[];

%% Use the position at the 3 detectors to reconstruct the vertex
% Assume first detector at 0.5 m
% Add MS error as random gaussian noise sigma(mm) = 2.85e-4 x 
% Do this for all possible choices and see the close my simulated data is
% to the true data's beam position (0,0) and beam variance 

[numPionsLeft, positions] = size(detectorX);
meanX1 = zeros(50,1);
meanY1 = zeros(50,1);
stdX1 = zeros(50,1);
stdY1 = zeros(50,1);

meanX2 = zeros(50,1);
meanY2 = zeros(50,1);
stdX2 = zeros(50,1);
stdY2 = zeros(50,1);

step = 1;
z1=0.5;

% Using relativistic energy equation to calculate velocity m/s
velocity = sqrt(1 - (139.57039./massEnergy).^2).*3e8;
time = (1./velocity)*detectorZ;
for i = 3:(positions-1)
    for j = 5:positions
        if i==j
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
            
            z2 = detectorZ(i);
            z3 = detectorZ(j);
    
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
    
            % Add shifts from multiple scattering
            MS1 = 0.00028571428 * z1;
            MS2 = 0.00028571428 * z2;
            MS3 = 0.00028571428 * z3;
            MS1_ = normrnd(1,1,size(curDetX(:,1))).*MS1;
            MS2_ = normrnd(1,1,size(curDetX(:,1))).*MS2;
            MS3_ = normrnd(1,1,size(curDetX(:,1))).*MS3;
            MSX = [MS1_ MS2_ MS3_];
    
            MS1_ = normrnd(1,1,size(curDetX(:,1))).*MS1;
            MS2_ = normrnd(1,1,size(curDetX(:,1))).*MS2;
            MS3_ = normrnd(1,1,size(curDetX(:,1))).*MS3;
            MSY = [MS1_ MS2_ MS3_];
    
            curDetX = curDetX + MSX;
            curDetY = curDetY + MSY;
    
            % Kalman filter
            dt = [time(:,i)-time(:,1),time(:,j)-time(:,1)];
            curDetZ = [z1, z2, z3];
            
            [PVx, PVy] = kalmanFilter(curDetX, curDetY, curDetZ, dt);
            random_x = PVx(randperm(size(PVx, 1)), :);
            random_y = PVy(randperm(size(PVy, 1)), :);

            [clustersCentroids,~,~] = clusterXYpoints([random_x(1:1000,:), random_y(1:1000,:)],0.000033,2,'centroid');
            meanX1(step,1)=mean(clustersCentroids(:,1));
            meanY1(step,1)=mean(clustersCentroids(:,1));
            stdX1(step,1)=std(clustersCentroids(:,1));
            stdY1(step,1)=std(clustersCentroids(:,1));

            [clustersCentroids,~,~] = clusterXYpoints([random_x(end-1000:end,:), random_y(end-1000:end,:)],0.000033,2,'centroid');
            meanX2(step,1)=mean(clustersCentroids(:,1));
            meanY2(step,1)=mean(clustersCentroids(:,1));
            stdX2(step,1)=std(clustersCentroids(:,1));
            stdY2(step,1)=std(clustersCentroids(:,1));
            step = step + 1;
        end
    end
end

% 450 GeV
[numPionsLeft, positions] = size(detectorX450);
meanX1450 = zeros(50,1);
meanY1450 = zeros(50,1);
stdX1450 = zeros(50,1);
stdY1450 = zeros(50,1);

meanX2450 = zeros(50,1);
meanY2450 = zeros(50,1);
stdX2450 = zeros(50,1);
stdY2450 = zeros(50,1);

step = 1;
z1450=0.5;

% Using relativistic energy equation to calculate velocity m/s
velocity450 = sqrt(1 - (139.57039./massEnergy450).^2).*3e8;
time450 = (1./velocity450)*detectorZ450;
for i = 3:(positions-1)
    for j = 5:positions
        if i==j
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
    
            % Add shifts from multiple scattering
            MS1 = 0.00028571428 * z1450;
            MS2 = 0.00028571428 * z2450;
            MS3 = 0.00028571428 * z3450;
            MS1_ = normrnd(1,1,size(curDetX450(:,1))).*MS1;
            MS2_ = normrnd(1,1,size(curDetX450(:,1))).*MS2;
            MS3_ = normrnd(1,1,size(curDetX450(:,1))).*MS3;
            MSX = [MS1_ MS2_ MS3_];
    
            MS1_ = normrnd(1,1,size(curDetX450(:,1))).*MS1;
            MS2_ = normrnd(1,1,size(curDetX450(:,1))).*MS2;
            MS3_ = normrnd(1,1,size(curDetX450(:,1))).*MS3;
            MSY = [MS1_ MS2_ MS3_];
    
            curDetX450 = curDetX450 + MSX;
            curDetY450 = curDetY450 + MSY;
    
            % Kalman filter
            dt450 = [time450(:,i)-time450(:,1),time450(:,j)-time450(:,1)];
            curDetZ450 = [z1450, z2450, z3450];
            
            [PVx450, PVy450] = kalmanFilter(curDetX450, curDetY450, curDetZ450, dt450);
            random_x450 = PVx450(randperm(size(PVx450, 1)), :);
            random_y450 = PVy450(randperm(size(PVy450, 1)), :);

            [clustersCentroids450,~,~] = clusterXYpoints([random_x450(1:1000,:), random_y450(1:1000,:)],0.000033,2,'centroid');
            meanX1450(step,1)=mean(clustersCentroids450(:,1));
            meanY1450(step,1)=mean(clustersCentroids450(:,1));
            stdX1450(step,1)=std(clustersCentroids450(:,1));
            stdY1450(step,1)=std(clustersCentroids450(:,1));

            [clustersCentroids,~,~] = clusterXYpoints([random_x450(end-1000:end,:), random_y450(end-1000:end,:)],0.000033,2,'centroid');
            meanX2450(step,1)=mean(clustersCentroids(:,1));
            meanY2450(step,1)=mean(clustersCentroids(:,1));
            stdX2450(step,1)=std(clustersCentroids(:,1));
            stdY2450(step,1)=std(clustersCentroids(:,1));
            step = step + 1;
        end
    end
end

%% Compare simulated values with real vales
% real mean = [0,0]
% real std = 5e-4

sigmaHitX = abs(meanX1 - meanX2);
sigmaHitY = abs(meanY1 - meanY2);

figure()
plot(sigmaHitY)
grid on
title("Resolution of Primary Vertex (7 TeV)")
xlabel("Config Number")
ylabel("Resolution")

realMean = 0;
realStd = 5e-4;

% weight the values and get a score for each configuatiom/
meanXError1 = meanX1;
meanYError1 = meanY1;
stdXError1 = abs(stdX1 - realStd)/1000;
stdYError1 = abs(stdY1 - realStd)/1000;

% 450 Gev
sigmaHitX450 = abs(meanX1450 - meanX2450);
sigmaHitY450 = abs(meanY1450 - meanY2450);

figure()
plot(sigmaHitY450)
grid on
title("Resolution of Primary Vertex (450 GeV)")
xlabel("Config Number")
ylabel("Resolution")

realMean = 0;
realStd = variance;

% weight the values and get a score for each configuatiom/
meanXError1450 = meanX1450;
meanYError1450 = meanY1450;
stdXError1450 = abs(stdX1450 - realStd)/1000;
stdYError1450 = abs(stdY1450 - realStd)/1000;

%% Add in extrapolation error
detectorDistances = importfile("C:\Users\black\OneDrive\matlab\detectorDistances.txt");
z2 = detectorDistances.z2;
z3 = detectorDistances.z3;

% extrapolation error
% sigma_ext = rescale((z2.*z2 + z3.*z3).*0.0001./((z3-z2).*(z3-z2)));
sigma_ext = (z2.*z2 + z3.*z3).*0.001./((z3-z2).*(z3-z2));
score1 = 1./(meanXError1 + meanYError1 + stdXError1/2 + stdYError1/2 + sigma_ext);
figure()
plot(score1)
grid on
title("Accuracy of each detector configuration (7TeV)")
xlabel("Detector Configuration")
ylabel("Score")

% 450 GeV
score1450 = 1./(meanXError1450 + meanYError1450 + stdXError1450/2 + stdYError1450/2 + sigma_ext);
figure()
plot(score1450)
grid on
title("Accuracy of each detector configuration (450 GeV)")
xlabel("Detector Configuration")
ylabel("Score")