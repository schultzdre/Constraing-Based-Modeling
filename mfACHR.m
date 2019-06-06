%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function performs the matrix-form ACHR algorithm
% 
% INPUTS
% model - reduced model. That is, all upper and lower bounds are the same
% as the feasible optima. Warmup points must be included as 'warmupPts'
% field. Other required fields are ub, lb and S.
% nSteps - steps of the algorithm for which points distribution will be 
% saved. The number of steps performed will be max(nSteps). Samples are
% saved to a text file. In parallel mode these need to ne divisible by 50.
% file - String. Name of text file for which points will be saved to
% 
% OPTIONAL INPUTS
% parallel - boolean defining whether to parform the sampling in dual core
% parallelization. The function does not support parallelization in more
% than two cores at this time. Default false
% waitwindow - boolean on whether to display progress window
% 
% OUTPUTS
% model - model containing final points distriubtion. Warmup points are
% deleted, but other fields are not.
% runTimeSteps - Time taken to complete each step in nSteps. In parallel 
% mode these need to ne divisible by 50.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,runTimeSteps] = mfACHR(model,nSteps,...
    waitwindow,file,pickup)
    % Define optional parameters
    if nargin < 3 || isempty(waitwindow)
        waitwindow = true;
    end
    if nargin < 4 || isempty(file)
        file = 'mfACHR';
    end
    if nargin < 5 || isempty(pickup)
        pickup = 0;
    end
    
    %Define function constants
    t0 = clock;
    nSamples = size(model.warmupPts,2);
    maxMinTol = 1e-9;
    Nor = null(full(model.S));
    stepstpf = max(nSteps);
    runTimeSteps = zeros(length(nSteps),1);
    
    %Read in warmup points
    if pickup == 0
        points = model.warmupPts;
    else
        load([file '_' num2str(pickup) '.mat'])
        fprintf('Picking up at %d steps\n\n',pickup)
    end
    model = rmfield(model,'warmupPts');
    
    %Initialize waiting bar
    if waitwindow
        h = waitbar(0,'Sampling');
    end
    
    %Sample nSteps
    for i = (pickup+1):stepstpf
        %Find center point
        center = mean(points,2);
        %Get random directions and normalize them
        dir = bsxfun(@minus,points(:,randperm(nSamples)),center);
        dir = bsxfun(@rdivide,dir,sqrt(sum(dir.^2))); 
        %Find distance of each point to the boundaries
        dub = bsxfun(@plus,-points,model.ub);
        dlb = bsxfun(@minus,points,model.lb);
        %Find steps from boundary in each direction
        dub = dub ./ dir;
        dlb = -dlb ./ dir;
        %Find forward and backwards step limits
        Fstep = 1./max(1./[dub;dlb]);
        Bstep = -1./(max(1./(-[dub;dlb])));
        %Calculate maximum stepsize
        Dstep = (Fstep-Bstep).*rand(1,nSamples) + Bstep;
        %Ignore columns that cannot move enough
        Dstep((abs(Fstep) < maxMinTol & abs(Bstep) < maxMinTol) | (Bstep > Fstep)) = 0;
        %Move points
        points = points + bsxfun(@times,dir,Dstep);

        %Account for numerical approximation error every 50 steps
        if mod(i,50) == 0
            if max(max(abs(model.S*points))) > 1e-09
                points = Nor*(Nor'*points);
            end
            if waitwindow
                waitbar(i/stepstpf);
            end
        end

        if ismember(i,nSteps)
            t1 = clock();
            fprintf([file ' ' num2str(i) ' steps performed\n'])
            %dlmwrite([file '_' num2str(i) '.txt'],points)
            save([file '_' num2str(i) '.mat'],'points')
            temptime = etime(t1, t0);
            runTimeSteps(nSteps==i) = temptime;
            fprintf([num2str(temptime) ' seconds\n\n'])
        end
    end
    if waitwindow
        close(h)
    end
    
    points = Nor*(Nor'*points);
    model.points = points;
end















