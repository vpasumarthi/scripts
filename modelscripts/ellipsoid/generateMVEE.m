% Frequently modified parameters
materialName = 'BVO';
speciesType = 'hole';
numSpecies = 2;
numTrajRecorded = 1.00E+02;
tFinal = 1.00E-05;
timeInterval = 1.00E-09;
maxLim = 1000; % Hem-2e: 1400, BVO-2e: 700, BVO-2h: 1000
xlimits = [-maxLim, maxLim];
ylimits = [-maxLim, maxLim];
zlimits = [-maxLim, maxLim];
numFrames = 100;

% Not so frequently modified parameters
inputFileName = 'unwrappedTraj.dat';
bohr2ang = 0.529177249;
positionArray = dlmread('unwrappedTraj.dat') * bohr2ang;
numPathStepsPerTraj = round(tFinal / timeInterval) + 1;
positionArraySize = size(positionArray);
nSpecies = positionArraySize(2) / 3;
dataArray = zeros(numPathStepsPerTraj, numTrajRecorded * nSpecies, 3);
numStepsPerFrame = round((numPathStepsPerTraj - 1) / numFrames);
tol = 0.001;
if numSpecies > 1
    speciesTail = 's';
else
    speciesTail = '';
end

for trajIndex = 0:numTrajRecorded-1
    headStart = trajIndex * numPathStepsPerTraj;
    for step =0:numPathStepsPerTraj-1
        stepPosition = positionArray(headStart + step + 1, :);
        for speciesIndex = 0:nSpecies-1
            dataArray(...
                step + 1, trajIndex * nSpecies + speciesIndex + 1, :) = ...
                stepPosition(speciesIndex * 3 + 1: (speciesIndex + 1) * 3);
        end
    end
end

F(numFrames) = struct('cdata',[],'colormap',[]);
outputFileName = strcat(materialName, '_', num2str(numSpecies), ...
                        speciesType, speciesTail, '.avi');
v = VideoWriter(outputFileName);
open(v);
index = 1;
figure('visible', 'off');
xlabelstr = sprintf('x (%c)', 197);
ylabelstr = sprintf('y (%c)', 197);
zlabelstr = sprintf('z (%c)', 197);
for step = 0:numPathStepsPerTraj-1
    if mod(step, numStepsPerFrame) == 0 && step ~= 0
        Pext = dataArray(step + 1, :, :);
        P = reshape(Pext, numTrajRecorded * nSpecies, 3)';
        [A , C] = MinVolEllipse(P, tol);
        plot3(P(1,:),P(2,:),P(3,:),'*')
        hold on
        Ellipse_plot(A,C)
        hold off
        xlabel(xlabelstr)
        ylabel(ylabelstr)
        zlabel(zlabelstr)
        figTitle = ['Diffusion of ', num2str(numSpecies), ' ', ...
                    speciesType, speciesTail, ...
                    ' depicted over ', num2str(numTrajRecorded), ...
                    ' traj in ms-BVO'];
        title(figTitle)
        xlim(xlimits)
        ylim(ylimits)
        zlim(zlimits)
        F(index) = getframe(gcf);
        writeVideo(v,F(index));
        index = index + 1;
    end
end
close(v);