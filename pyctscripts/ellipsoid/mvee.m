function mvee(materialName, speciesType, numSpecies, numTrajRecorded, ...
    tFinal, timeInterval, numFrames, plotPosData, ellipsoidConstruct, ...
    inputFileName, bohr2ang, nDim, tol)

global speciesTail

% Construct array of positional data
positionArray = dlmread(inputFileName) * bohr2ang;
numPathStepsPerTraj = round(tFinal / timeInterval) + 1;
posDataArray = zeros(numPathStepsPerTraj, numTrajRecorded * numSpecies, ...
    nDim);
for trajIndex = 0:numTrajRecorded-1
    headStart = trajIndex * numPathStepsPerTraj;
    for step =0:numPathStepsPerTraj-1
        stepPosition = positionArray(headStart + step + 1, :);
        for speciesIndex = 0:numSpecies-1
            posDataArray(step + 1, ...
                trajIndex * numSpecies + speciesIndex + 1, :) = ...
                stepPosition(speciesIndex * nDim + 1: ...
                (speciesIndex + 1) * nDim);
        end
    end
end

% Create a video file
if numSpecies > 1
    speciesTail = 's';
else
    speciesTail = '';
end
outputVideoFileName = strcat(materialName, '_', num2str(numSpecies), ...
    speciesType, speciesTail, '.avi');
videoFile = VideoWriter(outputVideoFileName);
open(videoFile);

% Compute axes limits for frame
finalPosArray = reshape(...
    posDataArray(step + 1, :, :), numTrajRecorded * numSpecies, nDim)';
axesLimits = computeFrameLimits(finalPosArray, plotPosData, ...
    ellipsoidConstruct, nDim, numTrajRecorded, numSpecies, tol);

% Generate video frame from ellipsoid analysis
frameIndex = 1;
numStepsPerFrame = round((numPathStepsPerTraj - 1) / numFrames);
semiAxesLengths = zeros(numFrames, nDim);
cartesianSemiAxesLengths = zeros(numFrames, nDim);
for step = 0:numPathStepsPerTraj-1
    if mod(step, numStepsPerFrame) == 0 && step ~= 0
        stepPosData = reshape(posDataArray(step + 1, :, :), ...
            numTrajRecorded * numSpecies, nDim)';
        if ellipsoidConstruct
            sumEllipseMatrix = zeros(nDim);
            sumCenter = zeros(nDim, 1);
            sumSemiAxesLengths = zeros(1, nDim);
            sumCartesianSemiAxesLengths = zeros(1, nDim);
            for trajIndex = 0:numTrajRecorded-1
                trajStepPosData = stepPosData(...
                    :, trajIndex * numSpecies + (1:numSpecies));
                [trajEllipseMatrix , trajCenter] = MinVolEllipse(...
                    trajStepPosData, tol);
                sumEllipseMatrix = sumEllipseMatrix + trajEllipseMatrix;
                sumCenter = sumCenter + trajCenter;
                [trajSemiAxesLengths, trajCartesianSemiAxesLengths] = ...
                    axesLengths(trajEllipseMatrix);
                sumSemiAxesLengths = sumSemiAxesLengths + ...
                    trajSemiAxesLengths;
                sumCartesianSemiAxesLengths = ...
                    sumCartesianSemiAxesLengths ...
                    + trajCartesianSemiAxesLengths;
            end
            ellipseMatrix = sumEllipseMatrix / numTrajRecorded;
            center = sumCenter / numTrajRecorded;
            semiAxesLengths(frameIndex, :) = sumSemiAxesLengths / ...
                numTrajRecorded;
            cartesianSemiAxesLengths(frameIndex, :) = ...
                sumCartesianSemiAxesLengths / numTrajRecorded;
        else
            [ellipseMatrix , center] = MinVolEllipse(stepPosData, tol);
            [semiAxesLengths(frameIndex, :), ...
                cartesianSemiAxesLengths(frameIndex, :)] = ...
                axesLengths(ellipseMatrix);
        end
        videoFrame = generateVideoFrame(...
            numSpecies, speciesType, numTrajRecorded, materialName, ...
            ellipseMatrix, center, plotPosData, stepPosData, axesLimits);
        writeVideo(videoFile, videoFrame);
        frameIndex = frameIndex + 1;
    end
end
close(videoFile);

plotTimeEvolutionSeries(semiAxesLengths, cartesianSemiAxesLengths)
end

function [semiAxesLengths, cartesianSemiAxesLengths] = ...
    axesLengths(ellipseMatrix)
nDim = length(ellipseMatrix);
semiAxesLengths = zeros(1, nDim);
cartesianSemiAxesLengths = zeros(1, nDim);
[eigVec, eigValMatrix] = eig(ellipseMatrix);
eigVal = diag(eigValMatrix);
semiAxesLengths(1, :) = eigVal.^-0.5;
cartesianSemiAxesLengths(1, :) = abs(eigVec * semiAxesLengths');
end

function F = generateVideoFrame(...
    numSpecies, speciesType, numTrajRecorded, materialName, ...
    ellipseMatrix, center, plotPosData, stepPosData, boundLimits)

global speciesTail

figure('visible', 'off');
if plotPosData
    plot3(stepPosData(1,:), stepPosData(2,:), stepPosData(3,:), '*')
    hold on
end
Ellipse_plot(ellipseMatrix, center)
hold off

xlabel(sprintf('x (%c)', 197))
ylabel(sprintf('y (%c)', 197))
zlabel(sprintf('z (%c)', 197))
figTitle = ['Diffusion of ', num2str(numSpecies), ' ', speciesType, ...
    speciesTail, ' depicted over ', num2str(numTrajRecorded), ...
    ' traj in ', materialName];
title(figTitle)
xlim(boundLimits(1, :))
ylim(boundLimits(2, :))
zlim(boundLimits(3, :))
F = getframe(gcf);
end

function boundLimits = computeFrameLimits(finalPosArray, plotPosData, ...
    ellipsoidConstruct, nDim, numTrajRecorded, numSpecies, tol)

if plotPosData
    minPosLimit = min(min(finalPosArray, [], 2));
    maxPosLimit = max(max(finalPosArray, [], 2));
    posLimits = max(abs(minPosLimit), abs(maxPosLimit));
    numDigits = ceil(log10(posLimits));
    boundLimitValue = ceil(posLimits / 10^(numDigits - 1)) ...
        * 10^(numDigits - 1);
    boundLimits = ones(nDim, 2) .* [-1, 1] * boundLimitValue;
else
    if ellipsoidConstruct
        sumEllipseMatrix = zeros(nDim);
        sumCenter = zeros(nDim, 1);
        for trajIndex = 0:numTrajRecorded-1
            trajStepPosData = finalPosArray(...
                :, trajIndex * numSpecies + (1:numSpecies));
            [trajEllipseMatrix , trajCenter] = MinVolEllipse(...
                trajStepPosData, tol);
            sumEllipseMatrix = sumEllipseMatrix + trajEllipseMatrix;
            sumCenter = sumCenter + trajCenter;
        end
        ellipseMatrix = sumEllipseMatrix / numTrajRecorded;
        center = sumCenter / numTrajRecorded;
    else
        [ellipseMatrix , center] = MinVolEllipse(finalPosArray, tol);
    end
    [~, cartesianSemiAxesLengths] = axesLengths(ellipseMatrix);
    maxSemiAxesLength = max(cartesianSemiAxesLengths);
    posLimits = [center - maxSemiAxesLength, center + maxSemiAxesLength];
    numDigits = ceil(log10(abs(posLimits)));
    boundLimits = sign(posLimits) .* ceil(abs(posLimits) ...
        ./ 10.^(numDigits - 1)) .* 10.^(numDigits - 1);
end

end

function plotTimeEvolutionSeries(semiAxesLengths, cartesianSemiAxesLengths)

% Plot time evolution of ellipsoid shape
figure('visible', 'off');
plot(semiAxesLengths)
xlabel('Frame Number')
ylabel(sprintf('Magnitude of semi-axes (%c)', 197));
title('Time evolution of ellipsoid shape')
legend('a', 'b', 'c', 'Location', 'southeast')
saveas(gcf, 'EllipsoidShape.png')

% Plot time evolution of bounding box limits
figure('visible', 'off');
plot(cartesianSemiAxesLengths)
xlabel('Frame Number')
ylabel(sprintf('Magnitude of cartesian projected semi-axes (%c)', 197));
title('Time evolution of cartesian projected ellipsoid shape')
legend('a', 'b', 'c', 'Location', 'southeast')
saveas(gcf, 'CartesianEllipsoidShape.png')

% Plot time evolution of ab/c anisotropy
figure('visible', 'off');
abPlaneBoundingLimits = sum(cartesianSemiAxesLengths(:, 1:2).^2, 2).^0.5;
cDirBoundingLimits = cartesianSemiAxesLengths(:, 3);
degreeOfAnisotropy = abPlaneBoundingLimits ./ cDirBoundingLimits;
plot(degreeOfAnisotropy)
xlabel('Frame Number')
ylabel('Degree of anisotropy')
title('Time evolution of anisotropy in ab-plane vs. c-direction')
saveas(gcf, 'DegreeOfAnistropy.png')

end
