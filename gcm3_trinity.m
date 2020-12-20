%% Fit a model using Trinity
% Cleanup first
clear all
close all

proj_id = 'gcm3';

%% First, enter the data
% Make a structure with the data (note that the name of the field needs to
% match the name of the variable in the JAGS model)

dataName   =  'Filtration Position';
dataStrc   =  getfield(load('data/Kruschke1993'), 'd');

condIdx    =  find(contains(dataStrc.conditionNames, dataName));

data.nP    =  dataStrc.nParticipants;
data.nS    =  dataStrc.nStimuli;
data.nT    =  dataStrc.nTrials;
data.nY    =  dataStrc.nParticipants * dataStrc.nTrials;

data.ptsX  =  zeros(1, data.nY);
data.ptsY  =  zeros(1, data.nY);
data.ctgA  =  zeros(1, data.nY);
data.ctgB  =  zeros(1, data.nY);

data.ptsX(1, 1:dataStrc.nStimuli)  =  dataStrc.points(:,1);
data.ptsY(1, 1:dataStrc.nStimuli)  =  dataStrc.points(:,2);
data.ctgA(1, 1:dataStrc.nStimuli)  =  dataStrc.categoryStructure(:, condIdx) == 1;
data.ctgB(1, 1:dataStrc.nStimuli)  =  dataStrc.categoryStructure(:, condIdx) == 2;

i = 0;

for partIdx = 1:dataStrc.nParticipants
    for stimIdx = 1:dataStrc.nStimuli
        i = i + 1;
        data.pe(i)   =  partIdx;
        data.st(i)   =  stimIdx;
        data.resp(i) =  dataStrc.y(stimIdx, condIdx, partIdx);
    end
end

dataMatrix = [data.pe' data.st' data.resp' data.ptsX' data.ptsY' data.ctgA' data.ctgB'];
csvwrite('/tmp/gcm3-data.csv',dataMatrix)

!echo 'pe,st,resp,ptsX,ptsY,ctgA,ctgB' > gcm3.csv
!cat /tmp/gcm3-data.csv >> gcm3.csv
!rm /tmp/gcm3-data.csv

%% Now, make all inputs that Trinity needs
% Write the model into a variable (cell variable)
model = {
    '# Generalized Context Model'
    'model{'
    ''
    '  lambda ~ dgamma(2, 1)'
    ''
    '  for (i in 1:nP) {'
    '      omega[i] ~ dunif(0, 1)'
    '      beta[i]  ~ dunif(0, 1)'
    ''
    '      for (j in 1:nS) {'
    '          for (k in 1:nS) {'
    '              s[i,j,k] = exp( -lambda * ('
    '                               omega[i]  * abs(ptsX[j] - ptsX[k]) +'
    '                          (1 - omega[i]) * abs(ptsY[j] - ptsY[k])'
    '                        ) )'
    '          }'
    ''
    '          inProdA[i,j] = inprod(s[i,j, ], ctgA[1:nS])'
    '          inProdB[i,j] = inprod(s[i,j, ], ctgB[1:nS])'
    '          sTimesA[i,j] = s[i,j,j] * ctgA[j]'
    '          sTimesB[i,j] = s[i,j,j] * ctgB[j]'
    ''
    '          diffA[i,j] =       beta[i]   * ( inProdA[i,j] - sTimesA[i,j] )'
    '          diffB[i,j] = ( 1 - beta[i] ) * ( inProdB[i,j] - sTimesB[i,j] )'
    ''
    '      }'
    ''
    '  }'
    ''
    '  for (t in 1:nY){'
    ''
    '      resp[t] ~ dbin(theta[t], nT)'
    ''
    '      theta[t] = diffA[pe[t],st[t]] / ( diffA[pe[t],st[t]] + diffB[pe[t],st[t]] )'
    ''
    '  }'
    ''
    '}'
    };

% List all the parameters of interest (cell variable)
params = {
    'beta', 'omega', 'lambda'
    };

% Write a function that generates a structure with one random value for
% each parameter in a field
generator = @()struct(...
    'omega' , rand(data.nP, 1), ...
    'beta'  , rand(data.nP, 1), ...
    'lambda', rand);

% Tell Trinity which engine to use
engine = 'jags';


%% Run Trinity with the CALLBAYES() function
tic
[stats, chains, diagnostics, info] = callbayes(engine, ...
    'model'          ,     model , ...
    'data'           ,      data , ...
    'outputname'     , 'samples' , ...
    'init'           , generator , ...
    'datafilename'   ,   proj_id , ...
    'initfilename'   ,   proj_id , ...
    'scriptfilename' ,   proj_id , ...
    'logfilename'    ,   proj_id , ...
    'nchains'        ,        3  , ...
    'nburnin'        ,      500  , ...
    'nsamples'       ,      500  , ...
    'monitorparams'  ,    params , ...
    'thin'           ,        1  , ...
    'workingdir'     ,    ['/tmp/' proj_id]  , ...
    'verbosity'      ,        0  , ...
    'saveoutput'     ,     true  , ...
    'parallel'       ,  isunix() , ...
    'modules'        ,  {'dic'}  );

fprintf('%s took %f seconds!\n', upper(engine), toc)


%% Inspect the results
% First, inspect the convergence of each parameter
% disp('Convergence statistics:')
% grtable(chains, 1.05)

% Now check some basic descriptive statistics averaged over all chains
% disp('Descriptive statistics for all chains:')
codatable(chains, 'deviance')
