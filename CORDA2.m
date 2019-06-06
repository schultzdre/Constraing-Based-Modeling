%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [tissue_model, rescue, HCtoMC, HCtoNC, MCtoNC] = CORDA2(model,metTests,...
%   HC,MC,NC,MCxNCthresh,constraint,constrainby,om,ci)
% CORDA calculates tissue specific reconstruction by determining reaction
% dependency using a cost association. If you wish to use default
% parameters (when possible) set them as empty
% INPUTS:
%   model - general model from which the tissue specific model will be
%       calculated
%   metTests - metabolic tests to be included in the reconstruction. This
%       argument should be a cell array of strings of size nx5, where n is 
%       the number of metabolic tests to be performed. The columns must be
%       as follows:
%       1 - String. Name of the metabolic test to be performed. These
%       should be different than reaction names in the model.
%       2 - Cell array of strings. Name of exchange reactions to be
%       included during the metabolic test. The exchange bound of all other
%       reactions will be set to zero. If this array is empty exchange
%       fluxe will be left as they are.
%       3 - nx2 numerical array. First columns defines respective lower 
%       bounds of exchange reactions defined previously, and second column
%       upper bounds.
%       4 - Cell array if strings. Extra reactions to be added to the
%       model during the test.
%       5 - Numeric. Value to constraint objective by. Upper and lower
%       bounds will be set to this value on the objective function. If this
%       value is empty or zero the standard constraint parameter will be
%       used with constrainby = 'val'.
%       6 - String. Reaction being tested. This is either a model reaction
%       or a string definind the reaction being tested.
%   HC - High confidence reactions. Reactions to be included in the model. Cell
%       array of strings.
%   MC - Medium confidence reactions to be included in the model if they do 
%       not depend on too many NC reactions to carry a flux. Cell array of 
%       strings.
%   NC - Negatice confidence reactions not to be included in the model. 
%       These reactions will be included in the tissue model only if they 
%       are necessary for the flux of HC reactions or for the flux of MCtoNC 
%       or more MC reactions. Cell array of strings.
% OPTIONAL INPUTS
%   MCxNCthresh - Define the threshold to include NC reactions based on MC
%       dependency. NC reactions will be included in the tissue model only
%       if they are deemed relevant for the flux of MCtoNC or more MC
%       reactions (or necessary for HC reactions). Forward and backwards
%       rates of reactions count separately. Default 2.
%   constraint - constraint value. Numerical. Default 1. Negative
%       percentage values will be used as positive.
%   constrainby - type of constraint used when defining reaction dependency. 
%       String. Either 'perc' or 'val'. 'perc' constrains reactions based 
%       on percentage from their optimal value. The reactions is first 
%       optimized and then held at the percentage from optima defined by 
%       the constraint value. 'val' constrains the reactions based on 
%       numerical values. Forward fluxes will be held numerically at the 
%       constraint value and negative fluxes at -constraint. If
%       constraintby is 'val' and a flux cannot reach the value constraint,
%       then that flux will be set to its maximum value. Default 'val'
%   om - cost assigned to reactions while determining reaction dependency.
%       In step 1, MC reactions will have a cost of sqrt(om). In steps 1
%       and 2 NC reactions, as well as OT reactions in step 3, will be
%       assigned a cost of om. Default 1e+04.
%   ci - cost increase. Numeric. After each dependency assessment, the cost
%       of high cost reactions will be increased by a factor of ci. For
%       instance, if ci = 0.01, after each dependency assessment the cost of
%       reactions used will be multiplied by 1.01. This increase is only
%       applied once every time a reaction is tested, so reaction costs do
%       not grow by more than a factor of ci. Default 0.01.
% OUTPUTS:
%   tissue_model - tissue specific model produced (will be a subset of the
%       input model).
%   rescue - nx2 cell array. First colum is composed of strings of the
%       names of MC reactions not included in the model. The second column
%       is composed of cell arrays of the NC reactions implicated as being
%       necesary for the flux of the deleted MC reaction. This output will
%       help in manually curating the model if desired.
%   HCtoMC - dataset. Describes reaction dependencies between HC and MC
%       reactions. That is, if HCtoMC (i,j) = 1, MC reaction j was 
%       implicated as associated with HC reaction i. If HCtoMC(i,j) = 0
%       then the reactions are not associated.
%   HCtoNC - dataset. Similat to HCtoMC but defined the association between
%       HC and NC reactions defined in step one.
%   MCtoNC - dataset. Similat to HCtoMC but defined the association between
%       MC and NC reactions defined in step two.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tissue_model, rescue, HCtoMC, HCtoNC, MCtoNC] = CORDA2(model,metTests,...
  HC,MC,NC,MCxNCthresh,constraint,constrainby,om,ci,printww)
%% Initialize waitbar and optional inputs
%Optional Inputs
if (nargin < 6) || isempty(MCxNCthresh)
    MCxNCthresh = 2;
end
if (nargin < 7) || isempty(constraint)
    constraint = 1;
end
if (nargin < 8) || isempty(constrainby) || strcmp(constrainby,'val')
    constrainby = 1;
else
    constrainby = 2;
end
if (nargin < 9) || isempty(om)
    om = 1e+04;
end
if (nargin < 10) || isempty(ci)
    ci = 0.01;
end
if (nargin < 11) || isempty(printww)
    printww = true;
end

if printww; h = waitbar(0,'initializing waitbar'); end

%Define positive flux thresholf
fluxThreshold = 1e-7;

%Make reaction groups column vectors
if isrow(HC)
    HC = HC';
end
if isrow(MC)
    MC = MC';
end
if isrow(NC)
    NC = NC';
end

%% Add metabolic tests if they are not added yet
if ~isempty(metTests)
    HC = unique(cat(1,HC,metTests(:,1)));
end

%% Ensure model passes all metabolic tests
modelOrg = model;
fails = [];
for tid = 1:size(metTests,1)
    if printww; waitbar(tid/size(metTests,1),h,'Assessing metabolic tests'); end
    model = modelOrg;
    model = addMetabolicTestTrial(model,metTests,tid,constraint,constrainby);
    % Test objective
    flux = optimizeCbModel(model,'max');
    if flux.stat ~= 1
        fails = [fails tid];
    end
end
model = modelOrg;
if ~isempty(fails)
    mess = 'Metabolic test(s) ';
    for j = 1:length(fails)
        mess = [mess num2str(fails(j)) ', '];
    end
    mess(end-1:end) = [];
    mess = [mess ' did not pass'];
    error(mess);
end

%% Decompose model and add cost metabolite
%Build vector of HC (2) MC(1) and NC (-1)
rxnTypes = zeros(length(model.rxns),1);
rxnTypes(ismember(model.rxns,HC)) = 3;
rxnTypes(ismember(model.rxns,MC)) = 2;
rxnTypes(ismember(model.rxns,NC)) = 1;

%find internal reactions that are actively reversible
orlen = length(model.rxns);
leng = find(model.lb < 0 & model.ub >= 0); 

%Tailor model
model.S = [model.S -model.S(:,leng); sparse(zeros(1,orlen+length(leng)))];
model.mets{length(model.mets)+1} = 'pseudomet';
model.b = zeros(length(model.mets),1);
model.c = zeros(orlen+length(leng),1);
model.ub = [model.ub; -model.lb(leng); 1e20];
model.lb = zeros(orlen+length(leng)+1,1);
model.rev = zeros(orlen+length(leng)+1,1);
model.rxns = cat(1,model.rxns,strcat(model.rxns(leng),'_CORDA_rev_rxn'));
model.csense = repmat('E',1,length(model.mets));

%Update reaction groups
tmpNames = strcat(model.rxns(leng),'_CORDA_rev_rxn');
tmpTypes = rxnTypes(leng);
HC = cat(1,HC,tmpNames(tmpTypes == 3));
MC = cat(1,MC,tmpNames(tmpTypes == 2));
NC = cat(1,NC,tmpNames(tmpTypes == 1));
OT = model.rxns(~ismember(model.rxns,cat(1,HC,MC,NC)));

%add reaction for pseudomet consumption
model.rxns = cat(1,model.rxns,'EX_pseudomet');
temp = zeros(length(model.mets),1);
temp(end) = -10;
model.S = [model.S temp];
model.c = [model.c; 1];
model = changeObjective(model,'EX_pseudomet',1);

%% Step 1
%Find ID of MC and NC reactions to associate costs
MCids = findRxnIDs(model,MC);
NCids = findRxnIDs(model,NC);
HCids = findRxnIDs(model,HC);
HCids = HCids(HCids ~= 0); % Remove metabolic tests from IDs
MCpres = false(length(MC),1); % Keep track of reactions used
NCpres = false(length(NC),1);
HCpres = false(length(HCids),1);
%Associate costs
cost = zeros(1,length(model.rxns)) + 1e-03;
cost(MCids) = sqrt(om);
cost(NCids) = om;
cost(end) = -50;
%Save modified model
modelMod = model;
%Other variables for HC
HCtodel = false(length(HC),1); %HCs that are unable to carry flux
HCtoMC = zeros(length(HC),length(MC)); %Dependecy of HCs on MCs
HCtoNC = zeros(length(HC),length(NC));

fprintf('\nStep 1\n')
for j = 1:length(HC)
    if printww; waitbar(j/length(HC),h,'Step 1 - finding support reactions'); end
    % Use modified model and base cost at beginning of the step
    model = modelMod;
    %Close reverse of the reaction (if available)
    if ismember([HC{j} '_CORDA_rev_rxn'],model.rxns) %Close reverse
        model = changeRxnBounds(model,[HC{j} '_CORDA_rev_rxn'],0,'b');
    elseif length(HC{j})>14 && ismember(HC{j}(1:end-14),model.rxns) %Close original if testing reverse
        model = changeRxnBounds(model,HC{j}(1:end-14),0,'b');
    end

    %If reaction is a test set it up
    rxid = findRxnIDs(model,HC{j});
    if rxid == 0
        model = addMetabolicTest(model,metTests,HC,j,constraint,constrainby);
    %Set it up otherwise
    else
        modelt = changeObjective(model,HC{j},1);
        flux = optimizeCbModel(modelt,'max');
        if constrainby == 1
            if flux.f >= constraint
                model = changeRxnBounds(model,HC{j},constraint,'b');
            else
                model = changeRxnBounds(model,HC{j},flux.f,'b');
            end
        else
            model = changeRxnBounds(model,HC{j},0.01*constraint*flux.f,'b');
        end
    end
    
    %Test see which MC and NC reactions are needed
    model.S(end,1:length(cost)) = cost;
    flux = optimizeCbModel(model,'min');
    if flux.stat == 1 && flux.x(findRxnIDs(model,HC{j})) >= fluxThreshold
        MCpres(flux.x(MCids) > fluxThreshold) = true;
        NCpres(flux.x(NCids) > fluxThreshold) = true;
        HCpres(flux.x(HCids) > fluxThreshold) = true;
        prevHCtoMC = HCtoMC(j,:);
        prevHCtoNC = HCtoNC(j,:);
        HCtoMC(j,flux.x(MCids) > fluxThreshold) = 1;
        HCtoNC(j,flux.x(NCids) > fluxThreshold) = 1;
        while ~(all(prevHCtoMC == HCtoMC(j,:)) && ...
                all(prevHCtoNC == HCtoNC(j,:)))
            model.S(end,MCids(prevHCtoMC ~= HCtoMC(j,:))) = ...
            model.S(end,MCids(prevHCtoMC ~= HCtoMC(j,:)))*(1+ci);
            model.S(end,NCids(prevHCtoNC ~= HCtoNC(j,:))) = ...
            model.S(end,NCids(prevHCtoNC ~= HCtoNC(j,:)))*(1+ci);
            flux = optimizeCbModel(model,'min');
            MCpres(flux.x(MCids) > fluxThreshold) = true;
            NCpres(flux.x(NCids) > fluxThreshold) = true;
            HCpres(flux.x(HCids) > fluxThreshold) = true;
            prevHCtoMC = HCtoMC(j,:);
            prevHCtoNC = HCtoNC(j,:);
            HCtoMC(j,flux.x(MCids) > fluxThreshold) = 1;
            HCtoNC(j,flux.x(NCids) > fluxThreshold) = 1;
        end
    elseif flux.stat ~= 1 || flux.x(findRxnIDs(model,HC{j})) < fluxThreshold
        HCtodel(j) = true;
    end
end
model = modelMod;

%Tailor reaction sets
HCtodel(ismember(HC,model.rxns(HCids(HCpres)))) = false;
fprintf([num2str(length(unique(strrep(HC(HCtodel),'_CORDA_rev_rxn','')))) ' blocked reactions removed from HC\n'])
HC(HCtodel) = [];
HCtoMC(HCtodel,:) = [];
HCtoNC(HCtodel,:) = [];
HCtoMC = mat2dataset(HCtoMC,'varnames',MC,'obsnames',HC);
HCtoNC = mat2dataset(HCtoNC,'varnames',NC,'obsnames',HC);

fprintf([num2str(length(find(MCpres))) ' partial MC reactions added to HC\n'])
HC = cat(1,HC,MC(MCpres));
MC(MCpres) = [];
fprintf([num2str(length(find(NCpres))) ' partial NC reactions added to HC\n'])
HC = cat(1,HC,NC(NCpres));
NC(NCpres) = [];

%% Step 2
%% Step 2.1
%Initialize variables
NCids = findRxnIDs(model,NC);
MCxNC = zeros(length(MC),length(NC));
MCtodel = false(length(MC),1);

%Assign costs
cost = zeros(1,length(model.rxns));
cost(NCids) = om;
cost(end) = -50;

fprintf('\nStep 2.1\n')
for j = 1:length(MC)
    if printww; waitbar(j/length(MC),h,'Step 2.1 - Checking MC and NC co-occurence'); end
    model = modelMod;
    
    %Close reverse of the reaction (if available)
    if ismember([MC{j} '_CORDA_rev_rxn'],model.rxns)
        model = changeRxnBounds(model,[MC{j} '_CORDA_rev_rxn'],0,'b');
    elseif length(MC{j})>14 && ismember(MC{j}(1:end-14),model.rxns)
        model = changeRxnBounds(model,MC{j}(1:end-14),0,'b');
    end
    
    %Set it up
    modelt = changeObjective(model,MC{j},1);
    flux = optimizeCbModel(modelt,'max');
    if constrainby == 1
        if flux.f >= constraint
            model = changeRxnBounds(model,MC{j},constraint,'b');
        else
            model = changeRxnBounds(model,MC{j},flux.f,'b');
        end
    else
        model = changeRxnBounds(model,MC{j},0.01*constraint*flux.f,'b');
    end
    
    %Check flux
    model.S(end,:) = cost;
    flux = optimizeCbModel(model,'min');
    if flux.stat ~= 1 || flux.x(findRxnIDs(model,MC{j})) < fluxThreshold
        MCtodel(j) = true;
    else
        while true
            prevMCxNC = MCxNC(j,:);
            MCxNC(j,flux.x(NCids) > fluxThreshold) = 1;
            if all(prevMCxNC == MCxNC(j,:))
                break
            end
            model.S(end,NCids(prevMCxNC ~= MCxNC(j,:))) = ...
            model.S(end,NCids(prevMCxNC ~= MCxNC(j,:)))*(1+ci);
            flux = optimizeCbModel(model,'min');
        end
    end
end
model = modelMod;
fprintf([num2str(length(unique(strrep(MC(MCtodel),'_CORDA_rev_rxn',''))))...
    ' reactions were at least partially deleted from MC\n'])
MCxNC(MCtodel,:) = [];
MC(MCtodel) = [];
MCtoNC = mat2dataset(MCxNC,'varnames',NC,'obsnames',MC);
%% Step 2.2
%Add high occuring NCs to MC
t = sum(MCxNC);
ind = NC(t >= MCxNCthresh);
fprintf([num2str(length(ind)) ' partial reactions from NC are added to HC\n'])
%Fix occurence matrix
MC = cat(1,MC,ind);
MCxNC = [MCxNC;zeros(length(ind),length(NC))];
MCxNC(:,ismember(NC,ind)) = [];
NC(ismember(NC,ind)) = [];

%See which reactions from MC are no longer feasible
model = changeRxnBounds(model,NC,0,'b');
modelMod2 = model;
clear modelMod
MCtodel = false(length(MC),1);
res1 = {};
res2 = {};

fprintf('\nStep 2.2\n')
for j = 1:length(MC)
    if printww; waitbar(j/length(MC),h,'Step 2.2 - Checking MC feasibility'); end
    model = changeObjective(modelMod2,MC{j},1);
    
    %Close reverse of the reaction (if available)
    if ismember([MC{j} '_CORDA_rev_rxn'],model.rxns)
        model = changeRxnBounds(model,[MC{j} '_CORDA_rev_rxn'],0,'b');
    elseif length(MC{j})>14 && ismember(MC{j}(1:end-14),model.rxns)
        model = changeRxnBounds(model,MC{j}(1:end-14),0,'b');
    end
    
    %Check flux
    flux = optimizeCbModel(model,'max');
    if flux.f < fluxThreshold
        MCtodel(j) = true;
    end
    
    if MCtodel(j)
        fprintf([MC{j} ' was deleted. Dependent on: '])
        res1 = cat(1,res1,MC{j});
        tmp = find(MCxNC(j,:));
        if isempty(tmp)
            % NC reaction added can be dependent on other NC reactions that
            % appear less than MCtoNC times, and have thus been blocked.
            fprintf('Undefined')
            res2 = cat(1,res2,' ');
        else
            for k = 1:length(tmp)
                fprintf(NC{tmp(k)})
                if k ~= length(tmp)
                    fprintf(', ')
                end
                if k == 1
                    res2 = cat(1,res2,NC{tmp(k)});
                else
                    res2{length(res2)} = strcat(res2{length(res2)},',',NC{tmp(k)});
                end
            end
        end
        fprintf('\n')
    end
end
model = modelMod2;
fprintf([num2str(length(unique(strrep(MC(MCtodel),'_CORDA_rev_rxn','')))) ' reactions deleted from MC\n'])
MC(MCtodel) = [];
HC = cat(1,HC,MC);
rescue = cat(2,res1,res2);

%% Step 3
%Block reactions not in HC or OT
model = changeRxnBounds(model,...
    model.rxns(~ismember(model.rxns,cat(1,HC,OT,'EX_pseudomet'))),0,'b');
modelMod3 = model;
clear modelMod2
%Define variables
OTids = findRxnIDs(model,OT);
HCxOT = zeros(length(HC),length(OT));
%Define cost
cost = zeros(1,length(model.rxns));
cost(OTids) = om;
cost(end) = -50;

%Parse through HC
fprintf('\nStep 3\n')
for j = 1:length(HC)
    if printww; waitbar(j/length(HC),h,'Step 3 - Define remaining reactions'); end
    model = modelMod3;
    
    %Close reverse of the reaction (if available)
    if ismember([HC{j} '_CORDA_rev_rxn'],model.rxns)
        model = changeRxnBounds(model,[HC{j} '_CORDA_rev_rxn'],0,'b');
    elseif length(HC{j})>14 && ismember(HC{j}(1:end-14),model.rxns)
        model = changeRxnBounds(model,HC{j}(1:end-14),0,'b');
    end
    
    %If reaction is a test set it up
    rxid = findRxnIDs(model,HC{j});
    if rxid == 0
        model = addMetabolicTest(model,metTests,HC,j,constraint,constrainby);
    %Set it up otherwise
    else
        modelt = changeObjective(model,HC{j},1);
        flux = optimizeCbModel(modelt,'max');
        if constrainby == 1
            if flux.f >= constraint
                model = changeRxnBounds(model,HC{j},constraint,'b');
            else
                model = changeRxnBounds(model,HC{j},flux.f,'b');
            end
        else
            model = changeRxnBounds(model,HC{j},0.01*constraint*flux.f,'b');
        end
    end
    
    %optimize HC 
    model.S(end,1:length(cost)) = cost;
    flux = optimizeCbModel(model,'min');
    %Some reactions can become active and be utilized during metabolic
    %tests, but alone, without the appropriate sinks and/or sources, be
    %innactive. These are accounted for here.
    if flux.x(findRxnIDs(model,HC{j})) < fluxThreshold
        fprintf([HC{j} ' cannot carry flux alone\n'])
        continue
    end
    while true
        prevHCxOT = HCxOT(j,:);
        HCxOT(j,flux.x(OTids) > fluxThreshold) = 1;
        if all(prevHCxOT == HCxOT(j,:))
            break
        end
        model.S(end,OTids(prevHCxOT ~= HCxOT(j,:))) = ...
        model.S(end,OTids(prevHCxOT ~= HCxOT(j,:)))*(1+ci);
        flux = optimizeCbModel(model,'min');
    end
end
fprintf([num2str(length(find(sum(HCxOT)))) ' partial reactions added to HC for final model\n'])
HC = cat(1,HC,OT(sum(HCxOT) ~= 0));
if printww; close(h); end

%% Genrate model and revert it
% Pick out reversible reactions and irreversible ones
x = regexp(HC,'CORDA_rev_rxn');
HCrev = HC(~cellfun(@isempty,x));
HCrev = strrep(HCrev,'_CORDA_rev_rxn','');
HCrev = HCrev(ismember(HCrev,modelOrg.rxns));
HC = HC(cellfun(@isempty,x));
HC = HC(ismember(HC,modelOrg.rxns));
HCall = unique(cat(1,HC,HCrev));

% Generate model
tissue_model = removeRxns(modelOrg,modelOrg.rxns(~ismember(modelOrg.rxns,HCall)));

% Set bounds and reversibility
tmp = findRxnIDs(tissue_model,HCall);
tissue_model.rev(tmp) = 0;
tmp = findRxnIDs(tissue_model,HCrev);
tissue_model.rev(tmp) = 1;
tissue_model = changeRxnBounds(tissue_model,HCall,0,'b');
tissue_model = changeRxnBounds(tissue_model,HC,1000,'u');
tissue_model = changeRxnBounds(tissue_model,HCrev,-1000,'l');
end

% Slightly modified version of addReaction from the COBRA toolbox. Even
% if the reaction does not exist in the model rxnIDexists returns the
% new reaction's ID. Also, the warning if the reaction already exists is
% suppressed.
function [model,rxnIDexists] = addReactionNW(model,rxnName,metaboliteList,stoichCoeffList,revFlag,lowerBound,upperBound,objCoeff,subSystem,grRule,geneNameList,systNameList,checkDuplicate)
%addReaction Add a reaction to the model or modify an existing reaction
%
% model = addReaction(model,rxnName,metaboliteList,stoichCoeffList,revFlag,lowerBound,upperBound,objCoeff,subSystem,grRule,geneNameList,systNameList,checkDuplicate)
% model = addReaction(model,rxnName,rxnFormula)
%
%INPUTS
% model             COBRA model structure
% rxnName           Reaction name abbreviation (i.e. 'ACALD')
%                   (Note: can also be a cell array {'abbr','name'}
% metaboliteList    Cell array of metabolite names or alternatively the
%                   reaction formula for the reaction
% stoichCoeffList   List of stoichiometric coefficients (reactants -ve,
%                   products +ve), empty if reaction formula is provided
%
%OPTIONAL INPUTS
% revFlag           Reversibility flag (Default = true)
% lowerBound        Lower bound (Default = 0 or -vMax)
% upperBound        Upper bound (Default = vMax)
% objCoeff          Objective coefficient (Default = 0)
% subSystem         Subsystem (Default = '')
% grRule            Gene-reaction rule in boolean format (and/or allowed)
%                   (Default = '');
% geneNameList      List of gene names (used only for translation from 
%                   common gene names to systematic gene names)
% systNameList      List of systematic names
% checkDuplicate    Check S matrix too see if a duplicate reaction is
%                   already in the model (Deafult true)
%
%OUTPUTS
% model             COBRA model structure with new reaction
% rxnIDexists       Empty if the reaction did not exist previously, or if
%                   checkDuplicate is false. Otherwise it contains the ID
%                   of an identical reaction already present in the model.
%
% Examples:
%
% 1) Add a new irreversible reaction using the formula approach
%
%    model = addReaction(model,'newRxn1','A -> B + 2 C')
%
% 2) Add a the same reaction using the list approach
%
%    model = addReaction(model,'newRxn1',{'A','B','C'},[-1 1 2],false);
%
% Markus Herrgard 1/12/07
%
% Modified the check to see if duplicate reaction already is in model by
% using S matrix coefficients to be able to handle larger matricies
% Richard Que 11/13/2008


parseFormulaFlag = false;
rxnIDexists = [];

if iscell(rxnName)&&length(rxnName)>1
    rxnNameFull = rxnName{2};
    rxnName = rxnName{1};
end

% Figure out if reaction already exists
nRxns = length(model.rxns);
if (sum(strcmp(rxnName,model.rxns)) > 0)
    warning('Reaction with the same name already exists in the model');
    [tmp,rxnID] = ismember(rxnName,model.rxns);
    oldRxnFlag = true;
else
    rxnID = nRxns+1;
    oldRxnFlag = false;
end

% Figure out what input format is used
if (nargin < 4)
    if (~iscell(metaboliteList))
        parseFormulaFlag = true;
    else
        error('Missing stoichiometry information');
    end
else
    if isempty(stoichCoeffList)
        parseFormulaFlag = true;
    else
        if (length(metaboliteList) ~= length(stoichCoeffList))
            error('Incorrect number of stoichiometric coefficients provided');
        end
    end
end

% Reversibility
if (nargin < 5 | isempty(revFlag))
    if (oldRxnFlag)
        revFlag = model.rev(rxnID);
    else
        revFlag = true;
    end
end

% Parse formula
if (parseFormulaFlag)
    rxnFormula = metaboliteList;
    [metaboliteList,stoichCoeffList,revFlag] = parseRxnFormula(rxnFormula);
end

% Missing arguments
if (nargin < 6 | isempty(lowerBound))
    if (oldRxnFlag)
        lowerBound = model.lb(rxnID);
    else
        if (revFlag)
            lowerBound = min(model.lb);
            if isempty(lowerBound)
                lowerBound=-1000;
            end
        else
            lowerBound = 0;
        end
    end
end
if (nargin < 7 | isempty(upperBound))
    if (oldRxnFlag)
        upperBound = model.ub(rxnID);
    else
        upperBound = max(model.ub);
        if isempty(upperBound)
            upperBound=1000;
        end
    end
end
if (nargin < 8 | isempty(objCoeff))
    if (oldRxnFlag)
        objCoeff = model.c(rxnID);
    else
        objCoeff = 0;
    end
end
if (nargin < 9)
    if (oldRxnFlag) && (isfield(model,'subSystems'))
        subSystem = model.subSystems{rxnID};
    else
        subSystem = '';
    end
end
if  isempty(subSystem)
    if (oldRxnFlag) && (isfield(model,'subSystems'))
        subSystem = model.subSystems{rxnID};
    else
        subSystem = '';
    end
end
if (nargin < 10) && (isfield(model,'grRules'))
    if (oldRxnFlag)
        grRule = model.grRules{rxnID};
    else
        grRule = '';
    end
end

if (~exist('checkDuplicate','var'))
   checkDuplicate=true;
end

nMets = length(model.mets);
Scolumn = sparse(nMets,1);

modelOrig = model;

% Update model fields
model.rxns{rxnID,1} = rxnName;
if (revFlag)
    model.rev(rxnID,1) = 1;
else
    model.rev(rxnID,1) = 0;
end
model.lb(rxnID,1) = lowerBound;
model.ub(rxnID,1) = upperBound;
model.c(rxnID,1) = objCoeff;

if (isfield(model,'rxnNames'))
    if exist('rxnNameFull','var')
        model.rxnNames{rxnID,1} = rxnNameFull;
    else
        model.rxnNames{rxnID,1} = model.rxns{rxnID};
    end
end
if (isfield(model,'subSystems'))
    model.subSystems{rxnID,1} = subSystem;
end
if isfield(model,'rxnNotes')
    model.rxnNotes{rxnID,1} = '';
end
if isfield(model,'confidenceScores')
    model.confidenceScores{rxnID,1} = '';
end
if isfield(model,'rxnReferences')
    model.rxnReferences{rxnID,1} = '';
end
if isfield(model,'rxnECNumbers')
    model.rxnECNumbers{rxnID,1} = '';
end


% Figure out which metabolites are already in the model
[isInModel,metID] = ismember(metaboliteList,model.mets);

nNewMets = sum(~isInModel);

% Construct S-matrix column
newMetsCoefs=zeros(0);
for i = 1:length(metaboliteList)
    if (isInModel(i))
        Scolumn(metID(i),1) = stoichCoeffList(i);
    else
        warning(['Metabolite ' metaboliteList{i} ' not in model - added to the model']);
        Scolumn(end+1,1) = stoichCoeffList(i);
        model.mets{end+1,1} = metaboliteList{i};
        newMetsCoefs(end+1) = stoichCoeffList(i);
        if (isfield(model,'metNames'))      %Prompts to add missing info if desired
            model.metNames{end+1,1} = regexprep(metaboliteList{i},'(\[.+\]) | (\(.+\))','') ;
            warning(['Metabolite name for ' metaboliteList{i} ' set to ' model.metNames{end}]);
%             model.metNames(end) = cellstr(input('Enter complete metabolite name, if available:', 's'));
        end
        if (isfield(model,'metFormulas'))
            model.metFormulas{end+1,1} = '';
            warning(['Metabolite formula for ' metaboliteList{i} ' set to ''''']);
%             model.metFormulas(end) = cellstr(input('Enter metabolite chemical formula, if available:', 's'));
        end
        if isfield(model,'metChEBIID')
            model.metChEBIID{end+1,1} = '';
        end
        if isfield(model,'metKEGGID')
            model.metKEGGID{end+1,1} = '';
        end
        if isfield(model,'metPubChemID')
            model.metPubChemID{end+1,1} = '';
        end
        if isfield(model,'metInChIString')
            model.metInChIString{end+1,1} = '';
        end
        if isfield(model,'metCharge')
            model.metCharge(end+1,1) = 0;
        end
    end
end

%printLabeledData(model.mets,Scolumn,1);

if isfield(model,'b')
    model.b = [model.b;zeros(length(model.mets)-length(model.b),1)];
end

% if ~oldRxnFlag, model.rxnGeneMat(rxnID,:)=0; end

if (isfield(model,'genes'))
    if (nargin < 11)
        model = changeGeneAssociation(model,rxnName,grRule);
    else
        fprintf('In addReaction, the class of systNameList is %s', class(systNameList));
        model = changeGeneAssociation(model,rxnName,grRule,geneNameList,systNameList);
    end
end

% Figure out if the new reaction already exists
rxnInModel=false;
if (nNewMets > 0) && isempty(find(newMetsCoefs == 0, 1))
    Stmp = [model.S;sparse(nNewMets,nRxns)];
else
    Stmp = model.S;
    if (checkDuplicate)
       if size(Stmp,2)<6000
           tmpSel = all(repmat((Scolumn),1,size(Stmp,2)) == (Stmp));
           rxnIDexists = full(find(tmpSel));
           if (~isempty(rxnIDexists))
               rxnIDexists=rxnIDexists(1);
               rxnInModel = true;
           end
       else
           for i=1:size(Stmp,2)
               if(Scolumn==Stmp(:,i))
                   rxnInModel=true;
                   rxnIDexists=i;
                   break
               end
           end
       end
    end
end

if (rxnInModel)
%     warning(['Model already has the same reaction you tried to add: ' modelOrig.rxns{rxnIDexists}]);
    model = modelOrig;
else
    if (oldRxnFlag)
        model.S = Stmp;
        model.S(:,rxnID) = Scolumn;
    else
        model.S = [Stmp Scolumn];
    end
    rxnIDexists = findRxnIDs(model,rxnName);
%     printRxnFormula(model,rxnName);
end
end

% Include metabolic text in the model and set appropriate bounds according
% to metTests.
function model = addMetabolicTest(model,metTests,HC,j,constraint,constrainby)
    tid = find(ismember(metTests(:,1),HC{j}));
    % Fix exchange reactions
    if ~isempty(metTests{tid,2})
        % Shut off all other exchanges
        model = changeRxnBounds(model,excOrg,0,'b');
        for k = 1:length(metTests{tid,2})
            % If reaction is in the model
            if ismember(metTests{tid,2}{k},model.rxns)
                % Set lower bound
                if metTests{tid,3}(k,1) >= 0 %If strictly positive
                    model = changeRxnBounds(model,metTests{tid,2}{k},metTests{tid,3}(k,1),'l');
                else
                    model = changeRxnBounds(model,[metTests{tid,2}{k} '_CORDA_rev_rxn'],-metTests{tid,3}(k,1),'u');
                end
                % Set upper bound
                if metTests{tid,3}(k,2) < 0 %If stricty negative
                    model = changeRxnBounds(model,[metTests{tid,2}{k} '_CORDA_rev_rxn'],-metTests{tid,3}(k,2),'l');
                else
                    model = changeRxnBounds(model,metTests{tid,2}{k},metTests{tid,3}(k,2),'u');
                end
            % If not in the model, add it
            else
                % If reaction is reversible, add it both ways
                if ~isempty(regexp(metTests{tid,2}{k},'<=>','ONCE')) || ~isempty(regexp(metTests{tid,2}{k},'<->','ONCE'))
                    metTests{tid,2}{k} = strrep(metTests{tid,2}{k},'<->','<=>');
                    rxnf = metTests{tid,2}{k};
                    rxnb = strsplit(rxnf,'<=>');
                    if length(rxnb)==1
                        if strcmp(rxnf(1:3),'<=>')
                            rxnb = [rxnb{1} ' <=> '];
                        else
                            rxnb = [' <=> ' rxnb{1}];
                        end
                    else
                        rxnb = [rxnb{2} ' <=> ' rxnb{1}];
                    end
                    model = addReactionNW(model,['testRxn' num2str(k) '_f'],rxnf);
                    model = changeRxnBounds(model,['testRxn' num2str(k) '_f'],metTests{tid,3}(k,2),'u');
                    model = addReactionNW(model,['testRxn' num2str(k) '_b'],rxnb);
                    model = changeRxnBounds(model,['testRxn' num2str(k) '_b'],-metTests{tid,3}(k,1),'u');
                % Otherwise, just add it
                else
                    model = addReactionNW(model,['testRxn' num2str(k)],metTests{tid,2}{k});
                    model = changeRxnBounds(model,['testRxn' num2str(k)],metTests{tid,3}(k,1),'l');
                    model = changeRxnBounds(model,['testRxn' num2str(k)],metTests{tid,3}(k,2),'u');
                end
            end
        end
    end
    % Add extra reactions
    if ~isempty(metTests{tid,4})
        % If only one reactino and not a cell array
        if ischar(metTests{tid,4})
            % If reversible decompose
            if ~isempty(regexp(metTests{tid,4},'<=>','ONCE')) || ~isempty(regexp(metTests{tid,4},'<->','ONCE'))
                metTests{tid,4} = strrep(metTests{tid,4},'<->','<=>');
                rxnf = metTests{tid,4};
                rxnb = strsplit(rxnf,'<=>');
                if length(rxnb)==1
                    if strcmp(rxnf(1:3),'<=>')
                        rxnb = [rxnb{1} ' <=> '];
                    else
                        rxnb = [' <=> ' rxnb{1}];
                    end
                else
                    rxnb = [rxnb{2} ' <=> ' rxnb{1}];
                end
                model = addReactionNW(model,'newRxn_f',rxnf);
                model = addReactionNW(model,'newRxn_b',rxnb);
            % Otherwise
            else
                model = addReactionNW(model,'newRxn_f',metTests{tid,4});
            end
        else
            for k = 1:length(metTests{tid,4})
                % If reversible decompose
                if ~isempty(regexp(metTests{tid,4}{k},'<=>','ONCE')) || ~isempty(regexp(metTests{tid,4}{k},'<->','ONCE'))
                    metTests{tid,4}{k} = strrep(metTests{tid,4}{k},'<->','<=>');
                    rxnf = metTests{tid,4}{k};
                    rxnb = strsplit(rxnf,'<=>');
                    if length(rxnb)==1
                        if strcmp(rxnf(1:3),'<=>')
                            rxnb = [rxnb{1} ' <=> '];
                        else
                            rxnb = [' <=> ' rxnb{1}];
                        end
                    else
                        rxnb = [rxnb{2} ' <=> ' rxnb{1}];
                    end
                    model = addReactionNW(model,'newRxn_f',rxnf);
                    model = addReactionNW(model,'newRxn_b',rxnb);
                % Otherwise
                else
                    model = addReactionNW(model,'newRxn_f',metTests{tid,4}{k});
                end
            end
        end
    end
    % Add objective
    if findRxnIDs(model,metTests{tid,6}) == 0
        [model, tmp] = addReactionNW(model,metTests{tid,1},metTests{tid,6});
    else
        tmp = findRxnIDs(model,metTests{tid,6});
        model.rxns{tmp} = metTests{tid,1};
    end
    % Set objective
    if isempty(metTests{tid,5}) || metTests{tid,5} == 0
        if constrainby == 1
            model = changeRxnBounds(model,model.rxns{tmp},constraint,'b');
        else
            modelt = changeObjective(model,model.rxns{tmp},1);
            flux = optimizeCbModel(modelt,'max');
            model = changeRxnBounds(model,model.rxns{tmp},0.01*constraint*flux.f,'b');
        end
    else
        model = changeRxnBounds(model,model.rxns{tmp},metTests{tid,5},'b');
    end
end

% Include metabolic test for initial assessment of metabolic tests.
function model = addMetabolicTestTrial(model,metTests,tid,constraint,constrainby)
    % Fix exchange reactions
    if ~isempty(metTests{tid,2})
        model = changeRxnBounds(model,model.rxns(findExcRxns(model)),0,'b');
        for k = 1:length(metTests{tid,2})
            if ismember(metTests{tid,2}{k},model.rxns)
                model = changeRxnBounds(model,metTests{tid,2}{k},metTests{tid,3}(k,1),'l');
                model = changeRxnBounds(model,metTests{tid,2}{k},metTests{tid,3}(k,2),'u');
            else
                model = addReactionNW(model,['testRxn' num2str(k)],metTests{tid,2}{k});
                model = changeRxnBounds(model,['testRxn' num2str(k)],metTests{tid,3}(k,1),'l');
                model = changeRxnBounds(model,['testRxn' num2str(k)],metTests{tid,3}(k,2),'u');
            end
        end
    end
    % Add extra reactions
    if ~isempty(metTests{tid,4})
        if ischar(metTests{tid,4})
            [model, tmp] = addReactionNW(model,['newRxn' num2str(k)],metTests{tid,4});
            if ~isempty(regexp(metTests{tid,4},'<=>','ONCE')) || ~isempty(regexp(metTests{tid,4},'<->','ONCE'))
                model = changeRxnBounds(model,model.rxns{tmp},-1000,'l');
                model.rev(tmp) = 1;
            end
        else
            for k = 1:length(metTests{tid,4})
                [model, tmp] = addReactionNW(model,['newRxn' num2str(k)],metTests{tid,4}{k});
                if ~isempty(regexp(metTests{tid,4},'<=>','ONCE')) || ~isempty(regexp(metTests{tid,4}{k},'<->','ONCE'))
                    model = changeRxnBounds(model,model.rxns{tmp},-1000,'l');
                    model.rev(tmp) = 1;
                end
            end
        end
    end
    % Add objective
    if findRxnIDs(model,metTests{tid,6}) == 0
        [model, tmp] = addReactionNW(model,metTests{tid,1},metTests{tid,6});
    else
        tmp = findRxnIDs(model,metTests{tid,6});
    end
    % Set objective
    if isempty(metTests{tid,5}) || metTests{tid,5} == 0
        if constrainby == 1
            model = changeRxnBounds(model,model.rxns{tmp},constraint,'b');
        else
            modelt = changeObjective(model,model.rxns{tmp},1);
            flux = optimizeCbModel(modelt,'max');
            model = changeRxnBounds(model,model.rxns{tmp},0.01*constraint*flux.f,'b');
        end
    else
        model = changeRxnBounds(model,model.rxns{tmp},metTests{tid,5},'b');
    end
    model = changeObjective(model,model.rxns{tmp},1);
end