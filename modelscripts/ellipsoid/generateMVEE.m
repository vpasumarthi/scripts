% Frequently modified parameters
materialName = 'ms-BVO';
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
plotEllipsoid = 1;
plotPrincipalAxes = 1;

% Not so frequently modified parameters
inputFileName = 'unwrappedTraj.dat';
bohr2ang = 0.529177249;
positionArray = dlmread('unwrappedTraj.dat') * bohr2ang;
numPathStepsPerTraj = round(tFinal / timeInterval) + 1;
positionArraySize = size(positionArray);
nSpecies = positionArraySize(2) / 3;
posDataArray = zeros(numPathStepsPerTraj, numTrajRecorded * nSpecies, 3);
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
            posDataArray(...
                step + 1, trajIndex * nSpecies + speciesIndex + 1, :) = ...
                stepPosition(speciesIndex * 3 + 1: (speciesIndex + 1) * 3);
        end
    end
end

if plotEllipsoid
    FrameStruct(numFrames) = struct('cdata',[],'colormap',[]);
    outputVideoFileName = strcat(materialName, '_', num2str(numSpecies),...
                                 speciesType, speciesTail, '.avi');
    v = VideoWriter(outputVideoFileName);
    open(v);
end
if plotPrincipalAxes
    semiAxesLengths = zeros(numFrames, 3);
    cartesianSemiAxesLengths = zeros(numFrames, 3);
end
frameIndex = 1;
figure('visible', 'off');
xlabelstr = sprintf('x (%c)', 197);
ylabelstr = sprintf('y (%c)', 197);
zlabelstr = sprintf('z (%c)', 197);
for step = 0:numPathStepsPerTraj-1
    if mod(step, numStepsPerFrame) == 0 && step ~= 0
        stepPosData = reshape(posDataArray(step + 1, :, :), ...
                              numTrajRecorded * nSpecies, 3)';
        [ellipseMatrix , center] = MinVolEllipse(stepPosData, tol);
        if plotEllipsoid
            plot3(stepPosData(1,:),stepPosData(2,:),stepPosData(3,:),'*')
            hold on
            Ellipse_plot(ellipseMatrix,center)
            hold off
            xlabel(xlabelstr)
            ylabel(ylabelstr)
            zlabel(zlabelstr)
            figTitle = ['Diffusion of ', num2str(numSpecies), ' ', ...
                        speciesType, speciesTail, ...
                        ' depicted over ', num2str(numTrajRecorded), ...
                        ' traj in ', materialName];
            title(figTitle)
            xlim(xlimits)
            ylim(ylimits)
            zlim(zlimits)
            FrameStruct(frameIndex) = getframe(gcf);
            writeVideo(v,FrameStruct(frameIndex));
        end
        if plotPrincipalAxes
           [eigVec, eigValMatrix] = eig(ellipseMatrix);
           eigVal = diag(eigValMatrix);
           semiAxesLengths(frameIndex, :) = eigVal.^-0.5';
           cartesianSemiAxesLengths(frameIndex, :) = ...
                            abs(eigVec * semiAxesLengths(frameIndex, :)');
        end
        frameIndex = frameIndex + 1;
    end
end
if plotEllipsoid
    close(v);
end
if plotPrincipalAxes
    figure('visible', 'off');
    plot(semiAxesLengths)
    xlabel('Frame Number')
    ylabel(sprintf('Magnitude of semi-axes (%c)', 197));
    title('Time evolution of ellipsoid shape')
    legend('a', 'b', 'c', 'Location', 'southeast')
    saveas(gcf, 'EllipsoidShape.png')
    figure('visible', 'off');
    plot(cartesianSemiAxesLengths)
    xlabel('Frame Number')
    ylabel(sprintf('Magnitude of cartesian projected semi-axes (%c)', 197));
    title('Time evolution of cartesian projected ellipsoid shape')
    legend('a', 'b', 'c', 'Location', 'southeast')
    saveas(gcf, 'CartesianEllipsoidShape.png')
end