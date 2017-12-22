function mvee(materialName, speciesType, numSpecies, numTrajRecorded, ...
              tFinal, timeInterval, numFrames, plotEllipsoid, ...
              plotPosData, plotPrincipalAxes, average_ellipsoid, ...
              inputFileName, bohr2ang, nDim, tol)

positionArray = dlmread(inputFileName) * bohr2ang;
numPathStepsPerTraj = round(tFinal / timeInterval) + 1;
positionArraySize = size(positionArray);
nSpecies = positionArraySize(2) / nDim;
posDataArray = zeros(numPathStepsPerTraj, numTrajRecorded * nSpecies, ...
                     nDim);
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
                stepPosition(speciesIndex * nDim + 1: ...
                             (speciesIndex + 1) * nDim);
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
frameIndex = 1;
figure('visible', 'off');
xlabelstr = sprintf('x (%c)', 197);
ylabelstr = sprintf('y (%c)', 197);
zlabelstr = sprintf('z (%c)', 197);
finalPosArray = reshape(posDataArray(step + 1, :, :), ...
                        numTrajRecorded * nSpecies, nDim)';
if plotPosData
    minPosLimit = min(min(finalPosArray, [], 2));
    maxPosLimit = max(max(finalPosArray, [], 2));
    posLimits = max(abs(minPosLimit), abs(maxPosLimit));
    numDigits = ceil(log10(posLimits));
    boundLimits = round(posLimits / 10^(numDigits - 1)) ...
                        * 10^(numDigits - 1);
    xlimits = [-boundLimits, boundLimits];
    ylimits = [-boundLimits, boundLimits];
    zlimits = [-boundLimits, boundLimits];
else
    if average_ellipsoid
        sumEllipseMatrix = zeros(nDim);
        sumCenter = zeros(nDim, 1);
        for trajIndex = 0:numTrajRecorded-1
            trajStepPosData = finalPosArray(...
                            :, trajIndex * nSpecies + (1:nSpecies));
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
    [semiAxesLengths, cartesianSemiAxesLengths] = ...
                                                axesLengths(ellipseMatrix);
    maxSemiAxesLength = max(cartesianSemiAxesLengths);
    posLimits = [center - maxSemiAxesLength, center + maxSemiAxesLength];
    numDigits = ceil(log10(abs(posLimits)));
    boundLimits = sign(posLimits) .* ...
                  ceil(abs(posLimits) ./ 10.^(numDigits - 1)) ...
                       .* 10.^(numDigits - 1);
    xlimits = boundLimits(1, :);
    ylimits = boundLimits(2, :);
    zlimits = boundLimits(3, :);
end

if plotPrincipalAxes
    semiAxesLengths = zeros(numFrames, nDim);
    cartesianSemiAxesLengths = zeros(numFrames, nDim);
end

for step = 0:numPathStepsPerTraj-1
    if mod(step, numStepsPerFrame) == 0 && step ~= 0
        stepPosData = reshape(posDataArray(step + 1, :, :), ...
                              numTrajRecorded * nSpecies, nDim)';
        if average_ellipsoid
            sumEllipseMatrix = zeros(nDim);
            sumCenter = zeros(nDim, 1);
            sumSemiAxesLengths = zeros(1, nDim);
            sumCartesianSemiAxesLengths = zeros(1, nDim);
            for trajIndex = 0:numTrajRecorded-1
                trajStepPosData = stepPosData(...
                                :, trajIndex * nSpecies + (1:nSpecies));
                [trajEllipseMatrix , trajCenter] = MinVolEllipse(...
                                                    trajStepPosData, tol);
                sumEllipseMatrix = sumEllipseMatrix + trajEllipseMatrix;
                sumCenter = sumCenter + trajCenter;
                if plotPrincipalAxes
                   [trajSemiAxesLengths, ...
                       trajCartesianSemiAxesLengths] = ...
                                            axesLengths(trajEllipseMatrix);
                   sumSemiAxesLengths = sumSemiAxesLengths + ...
                                                    trajSemiAxesLengths;
                   sumCartesianSemiAxesLengths = ...
                                       sumCartesianSemiAxesLengths + ...
                                       trajCartesianSemiAxesLengths;
                end
            end
            ellipseMatrix = sumEllipseMatrix / numTrajRecorded;
            center = sumCenter / numTrajRecorded;
            if plotPrincipalAxes
                semiAxesLengths(frameIndex, :) = sumSemiAxesLengths / ...
                                                        numTrajRecorded;
                cartesianSemiAxesLengths(frameIndex, :) = ...
                            sumCartesianSemiAxesLengths / numTrajRecorded;
            end
        else
            [ellipseMatrix , center] = MinVolEllipse(stepPosData, tol);
            if plotPrincipalAxes
               [semiAxesLengths(frameIndex, :), ...
                cartesianSemiAxesLengths(frameIndex, :)] = ...
                                                axesLengths(ellipseMatrix);
            end
        end
        if plotEllipsoid
            if plotPosData
                plot3(stepPosData(1,:), stepPosData(2,:), ...
                      stepPosData(3,:), '*')
                hold on
            end
            Ellipse_plot(ellipseMatrix, center)
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
            writeVideo(v, FrameStruct(frameIndex));
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
    ylabel(sprintf('Magnitude of cartesian projected semi-axes (%c)', ...
                                                                    197));
    title('Time evolution of cartesian projected ellipsoid shape')
    legend('a', 'b', 'c', 'Location', 'southeast')
    saveas(gcf, 'CartesianEllipsoidShape.png')
end

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