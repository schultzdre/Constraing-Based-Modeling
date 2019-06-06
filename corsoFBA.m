%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flux = corsoFBA(model,onstr,perc,costas)
% 
% corsoFBA minimizes the overall flux through the entire metabolism
% according to costs assigned to each reaction. The function first
% optimizes the objective function then constrains the objective to a given
% percentage. The model is then simulated again while minimizing the costs
% given. The function has been tested for maximizing and minimizing a given
% reaction, but will not work with multiple objectives.
% 
% This method is described in detail in the following publication:
%   Schultz, A., & Qutub, A. A. (2015). Predicting internal cell fluxes at 
%   sub-optimal growth. BMC systems biology, 9(1), 18.
% 
% INPUTS:
%   model - metabolic reconstruction to be simulated.
%   onstr - string defining whether the objective function should be
%       maximized or optimize. Options are 'max' and 'min'.
%   perc - percentage from optima at which to contrain the objective.
%       numeric value from 0 to 100.
%   costas - cost assigned to each reactions. The options are a numeric
%       value (all reactions will have the same cost), a vector of length
%       model.rxns (reaction model.rxns{i} will have a cost of costas(i)
%       both forward and backwards) or a vector of length twice that of
%       model.rxns (model.rxns{i} will have a cost of costas(i) in the
%       forward direction and a cost of costas(i+length(model.rxns)) in the
%       backwards direction.
% OUTPUTS:
%   flux - Struct of fields:
%       x - flux distribution with same length as model.rxns
%       y - Dual. Same as in optimizeCbModel
%       f - objective value
%       fm - cost assciate with the flux distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function flux = corsoFBA(model,onstr,perc,costas)
%determine objective function
flux1 = optimizeCbModel(model,onstr);
if abs(flux1.f) < 1e-6
%     warning('FBA problem infeasible')
    flux.f = [];
    flux.x = zeros(length(model.rxns),1);
    return
end

%relax results to avoid computational error
flux1.fm = flux1.x(find(model.c))*(perc/100) - 1e-6;   %Lower bound
flux1.f = flux1.x(find(model.c))*(perc/100) + 1e-6;   %Upper bound
%save original model
model1 = model;

%See if cost is of right length
if length(costas) == 1
    costas = ones(length(model.rxns),1);
end
if ~iscolumn(costas)
    costas = costas';
end
if length(costas)==length(model.rxns)
    costas = [costas; costas];
elseif length(costas) ~= 2*length(model.rxns)
    fprintf('Invalid length of costs\n');
    flux = [];
    return
end

%find internal reactions that are actively reversible
orlen = length(model.rxns);
leng = find(model.lb<0 & model.ub>=0); 

%Tailor model
model.S = [model.S -model.S(:,leng); sparse(zeros(1,orlen+length(leng)))];
model.mets{length(model.mets)+1} = 'pseudomet';
model.b = zeros(length(model.mets),1);
model.c = zeros(orlen+length(leng),1);
model.ub = [model.ub; -model.lb(leng)];
model.lb = zeros(orlen+length(leng),1);
model.S(end,:) = [costas(1:orlen); costas(orlen+leng)];
model.rxns = cat(1,model.rxns,strcat(model.rxns(leng),'added'));

%add reaction for pseudomet consumption
model.rxns = cat(1,model.rxns,'EX_pseudomet');
model.ub = [model.ub; 1e20];
model.lb = [model.lb; 0];
temp = zeros(length(model.mets),1);
temp(end) = -1;
model.S = [model.S temp];
model.c = [model.c; 1];

%change bounds on original optimized reaction
t = find(model1.c);
for k = 1:length(t)
    model = changeRxnBounds(model,model1.rxns(t(k)),flux1.f(k),'u');
    model = changeRxnBounds(model,model1.rxns(t(k)),flux1.fm(k),'l');
    if findRxnIDs(model,[model1.rxns{t(k)} 'added']) ~= 0
        model = changeRxnBounds(model,[model1.rxns{t(k)} 'added'],...
            0,'b');
    end
end

%perform FBA
flux2 = optimizeCbModel(model,'min');

%If solution is infeasible, widen the bounds for ATPcons and try again
% count = 0;
% while (flux2.stat) ~= 1 && (count~=100)
%     fprintf('WARNING: bounds widened. Solution might be infeasible\n');
%     flux1.fm = flux1.fm - 1e-4;
%     flux1.f = flux1.f + 1e-4;
%     model = changeRxnBounds(model,model1.rxns{find(model1.c)},flux1.f,'u');
%     model = changeRxnBounds(model,model1.rxns{find(model1.c)},flux1.fm,'l');
%     flux2 = optimizeCbModel(model,'min');
%     count = count+1;
% end

flux.x = flux2.x(1:orlen);
flux.x(leng) = flux.x(leng) - flux2.x((orlen+1):(end-1));
flux.x(abs(flux.x) < 1e-8) = 0;

if isfield(flux1,'y')
    flux.y = flux1.y;
end
if isfield(flux1,'f')
    flux.f = flux1.f*(perc/100);
end
if isfield(flux2,'f')
    flux.fm = flux2.f;
end
end

