%% FDTR dataset generation using COMSOL
% This script generates the high-fidelity FDTR dataset used for POD-NNROM training.
%
% Main workflow:
% 1) define the baseline material, geometric, and optical parameters;
% 2) specify which parameters are varied and their sampling ranges;
% 3) generate Latin hypercube samples in the prescribed parameter space;
% 4) build and solve the corresponding COMSOL model for each sample;
% 5) extract the probe-weighted complex temperature response and convert it to phase;
% 6) save the signals, sampled parameters, and dataset structure for later ROM training.
%
% Main outputs:
%   Outputs/FDTR_LHS_dataset.mat     % full dataset for ROM training
%   Outputs/FDTR_LHS_params.txt      % sampled parameter log
%   Outputs/FDTR_LHS_signals/        % signal files for individual samples
%

clc; clear; close all;
tic

%% ---------------- User inputs ----------------
freq = logspace(log10(100), log10(30e6), 100).';   % Hz

kz = [100 0.1 1.38];%W/(m*K)
kr = [100 0.1 1.38];
C  = [2.44 0.01 1.65];%MJ/m^3*K
h  = [80 1 300e3];%nm
w0 = 1.13;   % pump radius (1/e^2), um
w1 = 1.36;   % probe radius (1/e^2), um
d_trans = h(1);   % transducer thickness, m
d_sub   = h(3);  % substrate thickness, m  (finite approx semi-infinite)
rMax    = 30;   % computational domain radius, um
% --------  pattern radius  --------
Rpat = 2;      % patterned FDTR (disk radius) um
Pamp = 0.05;        % W

%% -------- Dataset / training set controls --------
Nsamp = 1;                 % Number of training samples to generate
rng(1);                      % Random seed for reproducibility

saveEachTxt = true;          % Save each sample as an individual txt file
filePrefix  = 'FDTR';        % Prefix for per-sample output files

% Isotropy flag: iso_flags(i)=1 forces kr(i)=kz(i) for that layer
iso_flags = [0 1 0];
kr(iso_flags==1) = kz(iso_flags==1);   
ParaChange.kz   = [1 1 1];
ParaChange.kr   = [1 0 1];     % kr(i) is ignored when isotropy forces kr(i)=kz(i)
ParaChange.C    = [1 0 1];
ParaChange.h    = [1 0 0];     
ParaChange.w0   = 1;
ParaChange.w1   = 1;
ParaChange.rMax = 0;
ParaChange.Rpat = 1;

% ---- Sampling ranges in the original engineering units [min; max] ----
Range.kz = [30 0.05 0.1;   200 0.2 10];   % 2x3  [min; max]
Range.kr = [30 0 0.1;   200 0 10];
Range.C  = [2.2 0 0.3;   2.6 0 3];
Range.h  = [50 0 0;   150 0 0];

Range.w0   = [0.8; 2];  % 2x1
Range.w1   = [0.8; 2];
Range.rMax = [0; 0];
Range.Rpat = [1.5; 5];


% ---- Sampling mode: ''lin'' = uniform, ''log'' = log-uniform (min>0 required) ----
RandMode.kz   = 'lin';
RandMode.kr   = 'lin';
RandMode.C    = 'lin';
RandMode.h    = 'lin';
RandMode.w0   = 'lin';
RandMode.w1   = 'lin';
RandMode.rMax = 'lin';
RandMode.Rpat = 'lin';

G    = kz(2)/h(2);       % W/m^2/K

% Options
opts.meshSize = 2;                 % 1(extremely fine) ... 9(coarse)
opts.saveMph  = false;
opts.normalizePumpPower = false;   % if Rpat not >> w0, set true to keep total pump power = Pamp within disk
opts.normalizeProbeW    = false;   % if Rpat comparable to w1, set true to renormalize probe weighting on [0,Rpat]

% ---------- Output folders ----------
outputsDir = fullfile(pwd, 'Outputs');
modelsDir  = fullfile(pwd, 'Models');

if ~exist(outputsDir,'dir'), mkdir(outputsDir); end
if ~exist(modelsDir,'dir'),  mkdir(modelsDir);  end

opts.modelsDir = modelsDir;   % give subfunction a place to save mph

outTxt = fullfile(outputsDir, 'FDTR_COMSOL_phase.txt');

%% ---------------- Materials ----------------
mat.trans.rho = 2700;       % Density is only used by COMSOL material input,No impact on signal
mat.trans.k_r = kr(1);
mat.trans.k_z = kz(1);
mat.trans.Cv  = C(1);

mat.sub.rho = 2200;
mat.sub.k_r = kr(3);
mat.sub.k_z = kz(3);
mat.sub.Cv  = C(3);

%% ---------------- Build + run ----------------
failLog = fullfile(outputsDir, 'FDTR_failed_samples.txt');


if ~exist(failLog,'file')
    fid = fopen(failLog,'w');
    fprintf(fid, ['idx\tmessage\tkz1\tkz2\tkz3\tkr1\tkr2\tkr3\t' ...
                  'C1\tC2\tC3\th1\th2\th3\tw0\tw1\trMax\tRpat\tG\n']);
    fclose(fid);
end


%% -------- LHS dataset + logging (NEW) --------
datasetFile = fullfile(outputsDir, 'FDTR_LHS_dataset.mat');     % Main cached dataset file
paramTxt    = fullfile(outputsDir, 'FDTR_LHS_params.txt');      % Parameter log for all sampled cases
sigDir      = fullfile(outputsDir, 'FDTR_LHS_signals');         % Per-sample signal txt files

useDatasetIfExists = true;      % Reuse an existing dataset when the configuration matches
exportTxtWhenLoad  = false;     % Re-export txt files after loading a cached dataset

if ~exist(sigDir,'dir'), mkdir(sigDir); end

% ---- baseline ----
base.kz   = kz;        % W/mK
base.kr   = kr;        % W/mK
base.C    = C;         % MJ/m^3K   
base.h    = h;         % nm       
base.w1   = w1;        % um
base.w0   = w0;        % um
base.rMax = rMax;      % um
base.Rpat = Rpat;      % um

% ---- Configuration signature used to check dataset compatibility ----
cfg = struct();
cfg.Nsamp      = Nsamp;
cfg.seed       = 1;                
cfg.freq       = freq;
cfg.ParaChange = ParaChange;
cfg.iso_flags  = iso_flags;
cfg.Range      = Range;
cfg.RandMode   = RandMode;
cfg.filePrefix = filePrefix;

% ---- Load the cached dataset directly if the configuration matches ----
if useDatasetIfExists && exist(datasetFile,'file')==2
    tmp = load(datasetFile,'D');
    D = tmp.D;

    if cfg_equal(D.cfg, cfg)
        fprintf('✅ Found existing dataset: %s\n', datasetFile);
        fprintf('   -> Skip COMSOL runs, using stored results.\n');

        % Optional: re-export signal txt files after loading the cached dataset
        if exportTxtWhenLoad && saveEachTxt
            for n = 1:D.Nsamp
                if ~D.okMask(n), continue; end
                outTxt_n = fullfile(sigDir, sprintf('%s_%05d.txt', filePrefix, n));
                writematrix([D.freq(:), D.phi_deg_all(:,n), D.amp_all(:,n)], outTxt_n, 'Delimiter','\t');
            end
            fprintf('   -> Re-exported signal txt into %s\n', sigDir);
        end

        
        toc; 
        return;
    else
        fprintf('⚠️ Dataset exists but config mismatch -> regenerate dataset.\n');
    end
end


varList = build_var_list(ParaChange, iso_flags, Range, RandMode);

if isempty(varList)
    error('ParaChange all values ​= 0: No parameters have changed, making it impossible to create an LHS dataset. At least one parameter must be set to 1.');
end

% ---- LHS sampling  ----
[X_raw, X_SI, varNames] = lhs_sample_all(Nsamp, varList, cfg.seed);



Nf = numel(freq);
phi_deg_all = nan(Nf, Nsamp);
amp_all     = nan(Nf, Nsamp);
okMask      = false(Nsamp,1);
msgCell     = strings(Nsamp,1);
tFull       = nan(Nsamp,1);
G_all       = nan(Nsamp,1);

fprintf('=== Start COMSOL dataset generation (LHS) ===\n');
fprintf('Nsamp = %d, Nvar = %d\n', Nsamp, numel(varList));


paramTxt = fullfile(outputsDir, 'FDTR_LHS_params.txt');

[paramFid, paramTxtUsed] = open_param_log(paramTxt);
cleanupObj = onCleanup(@() safe_fclose(paramFid)); 

write_param_header_fid(paramFid, varList);        
fprintf('Params log -> %s\n', paramTxtUsed);

for n = 1:Nsamp

    
    p = apply_lhs_row_to_params(base, X_raw(n,:), varList, iso_flags);

  
    kz_n = p.kz;                 % W/mK
    kr_n = p.kr;                 % W/mK
    C_SI = p.C * 1e6;            % MJ/m^3K -> J/m^3K
    h_SI = p.h * 1e-9;           % nm -> m
    w1_SI = p.w1 * 1e-6;         % um -> m
    w0_SI = p.w0 * 1e-6;                
    rMax_SI = p.rMax * 1e-6;     % um -> m
    Rpat_SI = p.Rpat * 1e-6;     % um -> m

    
    if Rpat_SI >= rMax_SI*(1-1e-12)
        Rpat_SI = rMax_SI*(1-1e-12);
    end
    if Rpat_SI <= 0, Rpat_SI = 1e-12; end
    if rMax_SI <= 0
        rMax_SI = 1e-9;
        Rpat_SI = min(Rpat_SI, rMax_SI*(1-1e-12));
    end

    
    d_trans = h_SI(1);
    d_sub   = h_SI(3);
    G       = kz_n(2) / h_SI(2);    % W/m^2/K
    G_all(n)= G;

    mat.trans.k_r = kr_n(1);
    mat.trans.k_z = kz_n(1);
    mat.trans.Cv  = C_SI(1);

    mat.sub.k_r = kr_n(3);
    mat.sub.k_z = kz_n(3);
    mat.sub.Cv  = C_SI(3);

try
    % Run COMSOL and extract the phase signal
    t0 = tic;
    model    = build_fdtr_model_Fig2(freq, w0_SI, w1_SI, d_trans, d_sub, rMax_SI, Rpat_SI, Pamp, mat, G, opts);
    bTopMeas = sel_top_transducer(model, d_trans, Rpat_SI);
    [phi_deg, amp, ~] = eval_phase_from_model(model, w1_SI, bTopMeas, opts.normalizeProbeW, Rpat_SI);
    tFull(n) = toc(t0);

    phi_deg_all(:,n) = phi_deg(:);
    amp_all(:,n)     = amp(:);
    okMask(n)        = true;
    msgCell(n)       = "";

   
    if saveEachTxt
        outTxt_n = fullfile(sigDir, sprintf('%s_%05d.txt', filePrefix, n));
        writematrix([freq(:), phi_deg(:), amp(:)], outTxt_n, 'Delimiter','\t');
    end

  
    write_param_line_fid(paramFid, n, true, "", tFull(n), X_SI(n,:), '%.4e');

    fprintf('✅ %d/%d done (t=%.2fs)\n', n, Nsamp, tFull(n));

catch ME
    okMask(n)  = false;
    msgCell(n) = string(strrep(ME.message, sprintf('\n'), ' '));

    % failLog
    fidFail = fopen(failLog, 'a');
    fprintf(fidFail, '%d\t%s\t', n, char(msgCell(n)));
    fprintf(fidFail, '%.15g\t%.15g\t%.15g\t', kz_n(1), kz_n(2), kz_n(3));
    fprintf(fidFail, '%.15g\t%.15g\t%.15g\t', kr_n(1), kr_n(2), kr_n(3));
    fprintf(fidFail, '%.15g\t%.15g\t%.15g\t', C_SI(1),  C_SI(2),  C_SI(3));
    fprintf(fidFail, '%.15g\t%.15g\t%.15g\t', h_SI(1),  h_SI(2),  h_SI(3));
    fprintf(fidFail, '%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n', w0_SI, w1_SI, rMax_SI, Rpat_SI, G);
    fclose(fidFail);

    write_param_line_fid(paramFid, n, false, ME.message, NaN, X_SI(n,:), '%.4e');

    fprintf('❌ Sample %d failed: %s\n', n, char(msgCell(n)));
    continue;
end


if mod(n,10)==0
    safe_fclose(paramFid);
    paramFid = fopen(paramTxtUsed, 'a');
    if paramFid < 0
        error('Cannot reopen params log: %s (close Excel/Notepad++, check Baidu sync lock)', paramTxtUsed);
    end
end
end

% ---- Save the dataset for future reuse ----
D = struct();
D.cfg         = cfg;
D.Nsamp       = Nsamp;
D.freq        = freq;
D.varList     = varList;
D.X_raw       = X_raw;          
D.X_SI        = X_SI;           
D.okMask      = okMask;
D.msgCell     = msgCell;
D.tFull       = tFull;
D.G_all       = G_all;
D.phi_deg_all = phi_deg_all;
D.amp_all     = amp_all;
D.sigDir      = sigDir;
D.paramTxt    = paramTxt;

save(datasetFile, 'D', '-v7.3');
fprintf('✅ Dataset saved: %s (success=%d/%d)\n', datasetFile, nnz(okMask), Nsamp);

toc



%% ========================= SUBFUNCTIONS =========================

function model = build_fdtr_model_Fig2(freq, w0, w1, d1, d2, rMax, Rpat, P, mat, G, opts)
% Build the 2D axisymmetric COMSOL model for one FDTR case.
    import com.comsol.model.*
    import com.comsol.model.util.*

    ModelUtil.clear;
    model = ModelUtil.create('Model');
    if isfield(opts,'modelsDir') && ~isempty(opts.modelsDir)
        model.modelPath(opts.modelsDir);
    else
        model.modelPath(pwd);
    end
    model.label('FDTR_Fig2_conventional_or_patterned.mph');

    % -------- Parameters --------
    model.param.set('d1',   sprintf('%.15g[m]', d1));
    model.param.set('d2',   sprintf('%.15g[m]', d2));
    model.param.set('rMax', sprintf('%.15g[m]', rMax));
    model.param.set('Rpat', sprintf('%.15g[m]', Rpat));
    model.param.set('wo',   sprintf('%.15g[m]', w0));
    model.param.set('w1',   sprintf('%.15g[m]', w1));
    model.param.set('P',    sprintf('%.15g[W]', P));
    model.param.set('G',    sprintf('%.15g[W/(m^2*K)]', G));
    model.param.set('TBR',  '1/G');

    Cp_trans = mat.trans.Cv / mat.trans.rho;
    Cp_sub   = mat.sub.Cv   / mat.sub.rho;

    model.param.set('ktr',   sprintf('%.15g[W/(m*K)]', mat.trans.k_r));
    model.param.set('ktz',   sprintf('%.15g[W/(m*K)]', mat.trans.k_z));
    model.param.set('ksr',   sprintf('%.15g[W/(m*K)]', mat.sub.k_r));
    model.param.set('ksz',   sprintf('%.15g[W/(m*K)]', mat.sub.k_z));
    model.param.set('rho_t', sprintf('%.15g[kg/m^3]',  mat.trans.rho));
    model.param.set('rho_s', sprintf('%.15g[kg/m^3]',  mat.sub.rho));
    model.param.set('Cp_t',  sprintf('%.15g[J/(kg*K)]', Cp_trans));
    model.param.set('Cp_s',  sprintf('%.15g[J/(kg*K)]', Cp_sub));

    % -------- Geometry (Fig2b style; Fig2a is just Rpat=rMax) --------
    model.component.create('comp1', true);
    model.component('comp1').geom.create('geom1', 2);
    model.component('comp1').geom('geom1').axisymmetric(true);
    g = model.component('comp1').geom('geom1');

    g.create('sub', 'Rectangle');        % substrate: 0..rMax, -d2..0
    g.feature('sub').set('pos',  {'0', '-d2'});
    g.feature('sub').set('size', {'rMax', 'd2'});

    g.create('tr', 'Rectangle');         % transducer disk: 0..Rpat, 0..d1
    g.feature('tr').set('pos',  {'0', '0'});
    g.feature('tr').set('size', {'Rpat', 'd1'});

    g.run;

    % -------- Materials --------
    model.component('comp1').material.create('matS', 'Common');
    model.component('comp1').material.create('matT', 'Common');

    domT = sel_domain_transducer(model, Rpat, d1);
    domS = sel_domain_substrate(model, rMax, d2);

    domT = unique(domT(:)).';
    domS = unique(domS(:)).';

    if isempty(domT)
        debugAll = mphselectbox(model,'geom1', [0, rMax; -d2, d1], 'domain');
        debugAll = unique(debugAll(:)).';
        error('Transducer domain selection failed. All domain candidates: %s', mat2str(debugAll));
    end
    if isempty(domS), error('Substrate domain selection failed.'); end

    
    ov = intersect(domT, domS);
    if ~isempty(ov)
        error('Material domain overlap detected: %s', mat2str(ov));
    end

    fprintf('Domain IDs: trans=%s, sub=%s\n', mat2str(domT), mat2str(domS));

    model.component('comp1').material('matT').selection.set(domT);
    model.component('comp1').material('matS').selection.set(domS);

    model.component('comp1').material('matS').propertyGroup('def').set( ...
        'thermalconductivity', {'ksr','0','0','0','ksr','0','0','0','ksz'} );
    model.component('comp1').material('matS').propertyGroup('def').set('density', 'rho_s');
    model.component('comp1').material('matS').propertyGroup('def').set('heatcapacity', 'Cp_s');
    model.component('comp1').material('matS').propertyGroup('def').set('thermalconductivity_symmetry','0');

    model.component('comp1').material('matT').propertyGroup('def').set( ...
        'thermalconductivity', {'ktr','0','0','0','ktr','0','0','0','ktz'} );
    model.component('comp1').material('matT').propertyGroup('def').set('density', 'rho_t');
    model.component('comp1').material('matT').propertyGroup('def').set('heatcapacity', 'Cp_t');
    model.component('comp1').material('matT').propertyGroup('def').set('thermalconductivity_symmetry','0');

    % -------- Physics --------
    model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');
    ht = model.component('comp1').physics('ht');

    % ---- Pump boundary heating on TOP of transducer (0..Rpat at z=d1) ----
    ht.create('bhs1', 'BoundaryHeatSource', 1);
    bTopTr = sel_top_transducer(model, d1, Rpat);
    ht.feature('bhs1').selection.set(bTopTr);

    if isfield(opts,'normalizePumpPower') && opts.normalizePumpPower
        ht.feature('bhs1').set('Qb', 'P*(2/(pi*wo^2))*exp(-2*r^2/wo^2)/(1-exp(-2*Rpat^2/wo^2))');
    else
        ht.feature('bhs1').set('Qb', 'P*(2/(pi*wo^2))*exp(-2*r^2/wo^2)');
    end
    ht.feature('bhs1').set('harmonicPerturbation', true);

    % ---- Thermal contact ONLY under transducer footprint (z=0, 0<r<Rpat) ----
    ht.create('tc1', 'ThermalContact', 1);
    bInt = sel_interface_under_transducer(model, Rpat);
    fprintf('bInt IDs = %s\n', mat2str(bInt));
    ht.feature('tc1').selection.set(bInt);
    ht.feature('tc1').set('ContactModel','EquThinLayer');
    ht.feature('tc1').set('Req','TBR');

    % ---- Outer exposed sample top (patterned case): adiabatic (default is insulation)
    bTopSampleOuter = sel_top_sample_outer(model, rMax, Rpat);
    if ~isempty(bTopSampleOuter)
        ht.create('tiTop','ThermalInsulation',1);
        ht.feature('tiTop').selection.set(bTopSampleOuter);
    end
   
    % ---- Transducer side wall (patterned case): adiabatic ----
    if Rpat < rMax*(1-1e-12)
        bSideTr = sel_transducer_side_boundary(model, Rpat, d1);
        if ~isempty(bSideTr)
            ht.create('tiSideTr','ThermalInsulation',1);
            ht.feature('tiSideTr').selection.set(bSideTr);
        end
    end

    
    % ---- Far boundaries: constant temperature (Fig2 blue) ----
    bRight = sel_far_right_boundary(model, rMax, d2, d1);
    ht.create('T_right','TemperatureBoundary',1);
    ht.feature('T_right').selection.set(bRight);
    ht.feature('T_right').set('T0', 0);

    bBottom = sel_bottom_boundary(model, d2, rMax);
    ht.create('T_bottom','TemperatureBoundary',1);
    ht.feature('T_bottom').selection.set(bBottom);
    ht.feature('T_bottom').set('T0', 0);

    ht.feature('init1').set('Tinit', 0);

    % -------- Mesh --------
    model.component('comp1').mesh.create('mesh1');
    if isfield(opts,'meshSize')
        model.component('comp1').mesh('mesh1').autoMeshSize(opts.meshSize);
    else
        model.component('comp1').mesh('mesh1').autoMeshSize(2);
    end
    model.component('comp1').mesh('mesh1').run;

    % -------- Study --------
    model.study.create('std1');
    model.study('std1').create('frlin', 'Frequencylinearized');
    F = num2str(freq(:).');     
    model.study('std1').feature('frlin').set('plist', F);

    
    try
        model.study('std1').feature('frlin').set('ngen', 5);
    catch
    end

    model.study('std1').run;

    if isfield(opts,'saveMph') && opts.saveMph
        ts = datestr(now,'yyyymmdd_HHMMSS');
        if isfield(opts,'modelsDir') && ~isempty(opts.modelsDir)
            model.save(fullfile(opts.modelsDir, ['FDTR_' ts '.mph']));
        else
            model.save(fullfile(pwd, ['FDTR_' ts '.mph']));
        end
    end
end

function [phi_deg, amp, dbg] = eval_phase_from_model(model, w1, bTopAll, normalizeProbeW, Rpat)
    T = mpheval(model, 'T', 'edim','boundary', 'selection', bTopAll);
    r = T.p(1,:);
    Ts = T.d1;            

    [rS, idx] = sort(r);
    Ts = Ts(:, idx);

    W = (4./w1.^2).*exp(-2.*rS.^2./w1.^2).*rS;

    denom = 1.0;
    if nargin >= 4 && normalizeProbeW
        denom = trapz(rS, W);
        if denom <= 0
            denom = 1.0;
        end
    end

    Nf = size(Ts,1);
    Tavg = zeros(Nf,1);
    for m = 1:Nf
        Tavg(m) = trapz(rS, W .* Ts(m,:)) / denom;
    end

    phi_deg = atan2(imag(Tavg), real(Tavg)) * 180/pi;
    amp     = abs(Tavg);

    dbg.r = rS; dbg.W = W; dbg.denom = denom; dbg.Rpat = Rpat;
end

%% --------- Selection helpers (robust, avoid grabbing axis by epsr) ---------
function b = sel_top_transducer(model, d1, Rpat)
    % Robust top boundary selection: scan epsz windows (avoid ultra-tiny boxes)
    epsr = max(1e-12, 1e-6*Rpat);
    epsz_list = [5e-10, 1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8];  % meters
    b = [];
    for epsz = epsz_list
        cand = mphselectbox(model, 'geom1', [0, Rpat+epsr; d1-epsz, d1+epsz], 'boundary');
        cand = unique(cand(:)).';
        if ~isempty(cand)
            b = cand;
            return;
        end
    end
    debugCand = mphselectbox(model, 'geom1', [0, Rpat+epsr; d1-1e-7, d1+1e-7], 'boundary');
    debugCand = unique(debugCand(:)).';
    error('Top transducer boundary selection failed. Candidates near z=d1: %s', mat2str(debugCand));
end

function b = sel_interface_under_transducer(model, Rpat)
    epsr = max(1e-12, 1e-6*Rpat);
    epsz_list = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]; 
    b = [];
    for epsz = epsz_list
        cand = mphselectbox(model,'geom1', [0, Rpat+epsr; -epsz, +epsz], 'boundary');
        cand = unique(cand(:)).';
        if ~isempty(cand)
            b = cand;
            return;
        end
    end

    debugCand = mphselectbox(model,'geom1', [0, Rpat+epsr; -1e-4, +1e-4], 'boundary');
    debugCand = unique(debugCand(:)).';
    error('Interface boundary selection failed. Candidates near z=0: %s', mat2str(debugCand));
end

function b = sel_top_sample_outer(model, rMax, Rpat)
    if Rpat >= rMax*(1-1e-12)
        b = []; return; 
    end
    epsz = 1e-12;
    epsr = max(1e-12, rMax*1e-9);
    b = mphselectbox(model, 'geom1', [Rpat+epsr, rMax-epsr; -epsz, +epsz], 'boundary');
    b = unique(b(:)).';
end

function b = sel_far_right_boundary(model, rMax, d2, d1)
    epsr = max(1e-12, rMax*1e-9);
    b = mphselectbox(model, 'geom1', [rMax-epsr, rMax+epsr; -d2, d1], 'boundary');
    b = unique(b(:)).';
    if isempty(b), error('Far-right boundary selection failed.'); end
end


function b = sel_bottom_boundary(model, d2, rMax)
    zB = -d2;
    epsz_list = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4];
    epsr = max(1e-12, 1e-9*rMax);
    b = [];
    for epsz = epsz_list
        cand = mphselectbox(model, 'geom1', ...
            [0, rMax+epsr; zB-epsz, zB+epsz], 'boundary');
        cand = unique(cand(:)).';
        if ~isempty(cand)
            b = cand;
            return;
        end
    end

    debugCand = mphselectbox(model, 'geom1', [0, rMax; zB-1e-3, zB+1e-3], 'boundary');
    debugCand = unique(debugCand(:)).';
    error('Bottom boundary selection failed. Candidates near z=-d2: %s', mat2str(debugCand));
end

function b = sel_transducer_side_boundary(model, Rpat, d1)
    epsr = max(1e-12, Rpat*1e-9);
    epsz = max(1e-12, d1*1e-6);

    b = mphselectbox(model, 'geom1', ...
        [Rpat-epsr, Rpat+epsr; 0+epsz, d1-epsz], 'boundary');

    b = unique(b(:)).';
end

function v = sample_vec_masked(v0, vmin, vmax, mask, mode)
    v = v0;
    for i = 1:numel(v0)
        if mask(i)
            v(i) = sample_in_range(vmin(i), vmax(i), mode);
        end
    end
end

function x = sample_scalar_if(x0, xmin, xmax, doChange, mode)
    if ~doChange
        x = x0; return;
    end
    x = sample_in_range(xmin, xmax, mode);
end

function x = sample_in_range(xmin, xmax, mode)
    if xmin > xmax, tmp=xmin; xmin=xmax; xmax=tmp; end
    if strcmpi(mode,'log')
        if xmin <= 0 || xmax <= 0
            error('Log sampling needs xmin>0, xmax>0. Got xmin=%.3g, xmax=%.3g', xmin, xmax);
        end
        a = log10(xmin); b = log10(xmax);
        x = 10^(a + (b-a)*rand());
    else
        x = xmin + (xmax-xmin)*rand();
    end
end

function tag = build_var_tag(ParaChange, iso_flags, kz, kr, C, h, w0, w1, rMax, Rpat)
    tags = strings(0);

    for i = 1:numel(kz)
        if ParaChange.kz(i)
            tags(end+1) = "kz"+i+"_"+safe_num(kz(i));
        end
    end
    for i = 1:numel(kr)

        if ParaChange.kr(i) && ~(iso_flags(i)==1)
            tags(end+1) = "kr"+i+"_"+safe_num(kr(i));
        end
    end
    for i = 1:numel(C)
        if ParaChange.C(i)
            tags(end+1) = "C"+i+"_"+safe_num(C(i));
        end
    end
    for i = 1:numel(h)
        if ParaChange.h(i)
            tags(end+1) = "h"+i+"_"+safe_num(h(i));
        end
    end

    if ParaChange.w0,   tags(end+1) = "w0_"+safe_num(w0); end
    if ParaChange.w1,   tags(end+1) = "w1_"+safe_num(w1); end
    if ParaChange.rMax, tags(end+1) = "rMax_"+safe_num(rMax); end
    if ParaChange.Rpat, tags(end+1) = "Rpat_"+safe_num(Rpat); end

    if isempty(tags)
        tag = "base";
    else
        tag = strjoin(tags, "_");
    end
end

function s = safe_num(x)
    s = sprintf('%.2e', x);
end


function domT = sel_domain_transducer(model, Rpat, d1)
    er_list = max(1e-12, Rpat * [1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]);
    ez_list = max(1e-12, d1   * [1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]);

    domT = [];
    for er = er_list
        for ez = ez_list
            box = [0-er, Rpat+er; 0-ez, d1+ez]; 
            cand = mphselectbox(model,'geom1', box, 'domain');
            cand = unique(cand(:)).';
            if ~isempty(cand)
                domT = cand;
                return;
            end
        end
    end
end

function domS = sel_domain_substrate(model, rMax, d2)
    er_list = max(1e-12, rMax * [1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]);
    ez_list = max(1e-12, d2   * [1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]);

    domS = [];
    for er = er_list
        for ez = ez_list
            box = [0-er, rMax+er; -d2-ez, 0+ez]; 
            cand = mphselectbox(model,'geom1', box, 'domain');
            cand = unique(cand(:)).';
            if ~isempty(cand)
                domS = cand;
                return;
            end
        end
    end
end


%% ========================= LHS + DATASET HELPERS =========================

function tf = cfg_equal(A, B)
tf = true;

if ~isequaln(A.Nsamp, B.Nsamp), tf=false; return; end
if ~isequaln(A.seed,  B.seed),  tf=false; return; end
if ~isequaln(A.ParaChange, B.ParaChange), tf=false; return; end
if ~isequaln(A.iso_flags,  B.iso_flags),  tf=false; return; end
if ~isequaln(A.Range,      B.Range),      tf=false; return; end
if ~isequaln(A.RandMode,   B.RandMode),   tf=false; return; end

fa = A.freq(:); fb = B.freq(:);
if numel(fa) ~= numel(fb), tf=false; return; end
if max(abs(log10(fa./fb))) > 1e-12, tf=false; return; end
end

function varList = build_var_list(ParaChange, iso_flags, Range, RandMode)

varList = struct('name',{},'group',{},'idx',{},'min',{},'max',{},'mode',{},'scaleSI',{});

for i = 1:3
    if ParaChange.kz(i)
        varList(end+1) = make_var(sprintf('kz%d',i),'kz',i,Range.kz(1,i),Range.kz(2,i),RandMode.kz,1.0); 
    end
end


for i = 1:3
    if ParaChange.kr(i) && ~(iso_flags(i)==1)
        varList(end+1) = make_var(sprintf('kr%d',i),'kr',i,Range.kr(1,i),Range.kr(2,i),RandMode.kr,1.0); 
    end
end


for i = 1:3
    if ParaChange.C(i)
        varList(end+1) = make_var(sprintf('C%d',i),'C',i,Range.C(1,i),Range.C(2,i),RandMode.C,1e6); 
    end
end

for i = 1:3
    if ParaChange.h(i)
        varList(end+1) = make_var(sprintf('h%d',i),'h',i,Range.h(1,i),Range.h(2,i),RandMode.h,1e-9);
    end
end


if ParaChange.w0
    varList(end+1) = make_var('w0','w0',1,Range.w0(1),Range.w0(2),RandMode.w0,1e-6);
end
if ParaChange.w1
    varList(end+1) = make_var('w1','w1',1,Range.w1(1),Range.w1(2),RandMode.w1,1e-6); 
end
if ParaChange.rMax
    varList(end+1) = make_var('rMax','rMax',1,Range.rMax(1),Range.rMax(2),RandMode.rMax,1e-6); 
end
if ParaChange.Rpat
    varList(end+1) = make_var('Rpat','Rpat',1,Range.Rpat(1),Range.Rpat(2),RandMode.Rpat,1e-6); 
end
end

function v = make_var(name, group, idx, vmin, vmax, mode, scaleSI)
v.name    = name;
v.group   = group;
v.idx     = idx;
v.min     = min(vmin, vmax);
v.max     = max(vmin, vmax);
v.mode    = lower(mode);    
v.scaleSI = scaleSI;
end

function [X_raw, X_SI, varNames] = lhs_sample_all(Nsamp, varList, seed)

rng(seed);

Nvar = numel(varList);
U = lhsdesign_safe(Nsamp, Nvar);

X_raw = zeros(Nsamp, Nvar);
X_SI  = zeros(Nsamp, Nvar);
varNames = cell(1,Nvar);

for j = 1:Nvar
    vmin = varList(j).min;
    vmax = varList(j).max;
    u = U(:,j);

    if strcmpi(varList(j).mode,'log')
        if vmin<=0 || vmax<=0
            error('Log sampling requires min/max>0. %s got [%.3g, %.3g]', varList(j).name, vmin, vmax);
        end
        a = log10(vmin); b = log10(vmax);
        x = 10.^(a + (b-a)*u);
    else
        x = vmin + (vmax-vmin)*u;
    end

    X_raw(:,j) = x;
    X_SI(:,j)  = x * varList(j).scaleSI;
    varNames{j} = varList(j).name;
end
end

function U = lhsdesign_safe(N, P)
if exist('lhsdesign','file') == 2
    U = lhsdesign(N, P);
    return;
end

U = zeros(N,P);
for j = 1:P
    edges = (0:N)'/N;
    r = rand(N,1)/N;
    U(:,j) = edges(1:end-1) + r;
    U(:,j) = U(randperm(N), j);
end
end

function p = apply_lhs_row_to_params(base, xRow_raw, varList, iso_flags)

p = base;

for j = 1:numel(varList)
    v = xRow_raw(j);
    switch lower(varList(j).group)
        case 'kz'
            p.kz(varList(j).idx) = v;
        case 'kr'
            p.kr(varList(j).idx) = v;
        case 'c'
            p.C(varList(j).idx)  = v;
        case 'h'
            p.h(varList(j).idx)  = v;
        case 'w0'
            p.w0 = v;
        case 'w1'
            p.w1 = v;
        case 'rmax'
            p.rMax = v;
        case 'rpat'
            p.Rpat = v;
    end
end


p.kr(iso_flags==1) = p.kz(iso_flags==1);

end

function [fid, usedPath] = open_param_log(paramTxt)

[folder, name, ext] = fileparts(paramTxt);
usedPath = paramTxt;

if ~exist(folder,'dir'), mkdir(folder); end


fid = fopen(usedPath, 'w'); 
if fid > 0
    return;
end


usedPath = fullfile(tempdir, [name '_WORK' ext]);
fid = fopen(usedPath, 'w');
if fid > 0
    warning('Cannot write to %s (sync lock/permission). Writing params to temp: %s', paramTxt, usedPath);
    return;
end

error('Cannot open params log for writing:\n  %s\n  %s\nTry closing Excel/Notepad++ and check folder permission.', ...
      paramTxt, usedPath);
end


function safe_fclose(fid)
if ~isempty(fid) && fid>0
    fclose(fid);
end
end


function write_param_header_fid(fid, varList)
fprintf(fid, '# idx\tstatus\tt_full(s)\t');
for j = 1:numel(varList)
    fprintf(fid, '%s_SI\t', varList(j).name);
end
fprintf(fid, 'message\n');
end


function write_param_line_fid(fid, idx, ok, msg, t_full, xRow_SI, fmt)
if nargin < 7 || isempty(fmt), fmt = '%.4e'; end

st = 'success';
if ~ok, st = 'fail'; end

msgOut = sanitize_one_line(msg);

fprintf(fid, '%d\t%s\t', idx, st);
if isnan(t_full), fprintf(fid, 'NaN\t'); else, fprintf(fid, '%.3f\t', t_full); end

for j = 1:numel(xRow_SI)
    fprintf(fid, [fmt '\t'], xRow_SI(j));
end

fprintf(fid, '%s\n', msgOut);
end


function s = sanitize_one_line(s)
if isempty(s), s = ''; return; end
s = char(s);
s = strrep(s, sprintf('\r'), ' ');
s = strrep(s, sprintf('\n'), ' ');
s = strrep(s, sprintf('\t'), ' ');
end