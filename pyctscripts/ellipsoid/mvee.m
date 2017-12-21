function mvee(materialName, speciesType, numSpecies, numTrajRecorded, ...
              tFinal, timeInterval, numFrames, plotEllipsoid, ...
              plotPrincipalAxes, average_ellipsoid, inputFileName, ...
              bohr2ang, tol)

positionArray = dlmread(inputFileName) * bohr2ang;
numPathStepsPerTraj = round(tFinal / timeInterval) + 1;
positionArraySize = size(positionArray);
nSpecies = positionArraySize(2) / 3;
posDataArray = zeros(numPathStepsPerTraj, numTrajRecorded * nSpecies, 3);
numStepsPerFrame = round((numPathStepsPerTraj - 1) / numFrames);

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
finalPosArray = reshape(posDataArray(step + 1, :, :), ...
                        numTrajRecorded * nSpecies, 3);
minPosLimit = min(min(finalPosArray, [], 1));
maxPosLimit = max(max(finalPosArray, [], 1));
posLimits = max(abs(minPosLimit), abs(maxPosLimit));
numDigits = ceil(log10(posLimits));
boundLimits = round(posLimits / 10^(numDigits - 1)) * 10^(numDigits - 1);
xlimits = [-boundLimits, boundLimits];
ylimits = [-boundLimits, boundLimits];
zlimits = [-boundLimits, boundLimits];
for step = 0:numPathStepsPerTraj-1
    if mod(step, numStepsPerFrame) == 0 && step ~= 0
        stepPosData = reshape(posDataArray(step + 1, :, :), ...
                              numTrajRecorded * nSpecies, 3)';
        if average_ellipsoid
            sumEllipseMatrix = zeros(3);
            sumCenter = zeros(3, 1);
            for trajIndex = 0:numTrajRecorded-1
                trajStepPosData = stepPosData(...
                                :, trajIndex * nSpecies + (1:nSpecies));
                [trajEllipseMatrix , trajCenter] = MinVolEllipse(...
                                                    trajStepPosData, tol);
                sumEllipseMatrix = sumEllipseMatrix + trajEllipseMatrix;
                sumCenter = sumCenter + trajCenter;
            end
            ellipseMatrix = sumEllipseMatrix / numTrajRecorded;
            center = sumCenter / numTrajRecorded;
        else
            [ellipseMatrix , center] = MinVolEllipse(stepPosData, tol);
        end
        if plotEllipsoid
%             plot3(stepPosData(1,:),stepPosData(2,:),stepPosData(3,:),'*')
%             hold on
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

end