function testSammi(testnum)
if nargin < 1
    testnum = 1;
end

%Get COBRA directory
global CBTDIR;

% Load models to test
load([CBTDIR '/test/models/mat/ecoli_core_model.mat'])
x = regexp(model.mets,'\[(.)\]$','tokens');
x = [x{:}]; x = [x{:}];
model.compartment = x';

% load([CBTDIR '/test/models/mat/iJO1366.mat'])
% model = iJO1366; clear iJO1366;

% load([CBTDIR '/test/models/mat/Recon2.v04.mat'])
% model = modelR204; clear modelR204;

%Set seed for consistency
rng(734)

switch testnum
    case 0
        sammi(model)
    case 1
        %Get unique subsystems
        ss = unique(model.subSystems);
        %Make subsystem struct
        for i = 1:length(ss)
            dat(i).name = ss{i};
            dat(i).rxns = model.rxns(ismember(model.subSystems,ss{i}));
            if rand > 0.5; continue; end
            dat(i).flux = 10*randn(length(dat(i).rxns),1);
            dat(i).flux(randsample(length(dat(i).flux),floor(0.2*length(dat(i).flux)))) = NaN;
        end
        %Plot
        sammi(model,dat);
    case 2
        %Define number of conditions
        n = 5;
        %Make reaction table with random data
        rxntbl = randn(length(model.rxns),n);
        rxntbl(randsample(length(model.rxns)*n,floor(n*length(model.rxns)/10))) = NaN;
        rxntbl = array2table(rxntbl,'VariableNames',sprintfc('condition_%d',1:n),...
            'RowNames',model.rxns);
        %Make metabolites table with random data
        mettbl = randn(length(model.mets),n);
        mettbl(randsample(length(model.mets)*n,floor(0.5*length(model.mets)))) = NaN;
        mettbl = array2table(mettbl,'VariableNames',sprintfc('condition_%d',1:n),...
            'RowNames',model.mets);
        %Make struct
        dat(1).type = {'rxns' 'color'};
        dat(1).data = rxntbl;
        dat(2).type = {'rxns' 'size'};
        dat(2).data = rxntbl;
        dat(3).type = {'mets' 'color'};
        dat(3).data = mettbl;
        dat(4).type = {'mets' 'size'};
        dat(4).data = mettbl;
        dat(5).type = {'links' 'size'};
        dat(5).data = rxntbl;
        %Define secondaries
        secondaries = {'^h\[.\]$','^h20\[.\]$','^o2\[.\]$','^co2\[.\]$',...
            '^atp\[.\]$','^adp\[.\]$','^pi\[.\]$',...
            '^nadh\[.\]$','^nadph\[.\]$','^nad\[.\]$','^nadp\[.\]$'};
        %Plot dividing up by subsystems
        sammi(model,'subSystems',dat,secondaries)
        %sammi(model,'compartment',[])
    case 3
        options.jscode = 'zoom.transform(gMain, d3.zoomIdentity.translate(-1149,-863).scale(2.64));';
        sammi(model,[CBTDIR '/external/visualization/sammif/demo.json'],[],[],options)
    case 4
        rxns = model.rxns(ismember(model.subSystems,'Glycolysis/Gluconeogenesis'));
        sammi(model,rxns)
end
end





