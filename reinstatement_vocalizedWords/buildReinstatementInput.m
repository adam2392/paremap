function [reinstatement_inputone, reinstatement_inputtwo] = buildReinstatementInput(featureMatOne, featureMatTwo, WITHIN)
    % initialize result matrices
    numEventsOne = size(featureMatOne, 1);
    numEventsTwo = size(featureMatTwo, 1);
    
    % BUILD REINSTATEMENT MAT FROM 
    if ~WITHIN
        numComparisons = size(featureMatOne, 1) * size(featureMatTwo, 1);
        reinstatement_inputone = zeros(numComparisons, size(featureMatOne, 2), size(featureMatOne, 3));
        reinstatement_inputtwo = zeros(numComparisons, size(featureMatTwo, 2), size(featureMatTwo, 3));
    
        for i=1:size(featureMatOne, 1),
            reinstatement_inputone((i-1)*numEventsTwo + 1 : (i)*numEventsTwo, :, :) = ...
                repmat(featureMatOne(i,:, :), numEventsTwo, 1, 1);
            reinstatement_inputtwo((i-1)*numEventsTwo + 1:(i)*numEventsTwo, :, :) = ...
                featureMatTwo;
        end
    else
        numComparisons = (size(featureMatOne, 1) * (size(featureMatTwo, 1) - 1))/2;
        reinstatement_inputone = zeros(numComparisons, size(featureMatOne, 2), size(featureMatOne, 3));
        reinstatement_inputtwo = zeros(numComparisons, size(featureMatTwo, 2), size(featureMatTwo, 3));
        
        for i=1:size(featureMatOne, 1) - 1,
            % intelligently go through and define ranges, so function is
            % fast
            if i==1
                rnge = 1:size(featureMatOne, 1) - 1;
                prevrnge = max(size(rnge));
            else
                rnge = prevrnge+1:sum([size(featureMatOne, 1) - i: size(featureMatOne, 1) - 1]);
                prevrnge = sum([size(featureMatOne, 1) - i: size(featureMatOne, 1) - 1]);
            end
            
            reinstatement_inputone(rnge, :, :) = ...
                repmat(featureMatOne(i,:, :), size(featureMatOne, 1) - i, 1, 1);
            reinstatement_inputtwo(rnge, :, :) = ...
                featureMatTwo(i+1:end,:,:);
        end
    end
    
    reinstatement_inputone = permute(reinstatement_inputone, [1 3 2]);
    reinstatement_inputtwo = permute(reinstatement_inputtwo, [1 3 2]);
end