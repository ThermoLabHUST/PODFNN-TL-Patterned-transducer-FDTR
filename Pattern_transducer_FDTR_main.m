%% FDTR main program: signal generation, sensitivity analysis, and inverse fitting
% This script is the main analysis program for the FDTR forward model and the POD-NNROM.
%
% Main functions:
% 1) generate FDTR signals using either the full-order COMSOL model or the POD-NNROM;
% 2) compute phase sensitivities of selected parameters;
% 3) perform automatic inverse fitting using PSO and optional quasi-Newton refinement;
% 4) compare experimental/synthetic target signals with simulated fitting results.
%
% Available forward-model options:
% - Full-order COMSOL model
% - POD-NNROM surrogate model
%
% Main outputs:
%   Outputs/FDTR_COMSOL_phase.txt    % full-order simulated signal
%   Outputs/FDTR_ROM_phase.txt       % ROM-predicted signal
%   Outputs/FDTR_fit_result.mat      % fitting results structure
%   Outputs/FDTR_fit_result.txt      % fitting summary
%
% Notes:
% - The script uses switch flags to control signal generation, sensitivity analysis,
%   ROM usage, and inverse fitting.
% - The ROM model must be trained in advance before ROM-based prediction or fitting.
% - The current implementation is written for the 3-layer patterned FDTR case used in this study.

clc; clear; close all;
tic

outputsDir = fullfile(pwd, 'Outputs');
modelsDir  = fullfile(pwd, 'Models');

%% ---------------- User inputs ----------------
freq = logspace(log10(100), log10(30e6), 100).';   % Hz
% ===== Main parameters =====
% Units here:
%   kz, kr : W/mK
%   C      : J/m^3K
%   h      : m
%   w0,w1,rMax,Rpat : m

kz = [60 0.1 1.38];       % W/mK
kr = [60 0.1 1.38];       % W/mK
C  = [2.44 0.01 1.65] * 1e6;  % J/m^3K
h  = [120 1 300e3] * 1e-9;     % m
w0 = 1.645e-6;                 % pump (m     
 w1 = 1.175e-6;                %probe (m
Pamp = 0.05;                  % W
d_trans = h(1);               % m
d_sub   = h(3);               % m
rMax    = 30e-6;             % m
Rpat    = 2e-6;               % m Pattern Radius
G    = kz(2) / h(2);          % W/m^2/K

%% ---------------- Options / switches ----------------
runSim     = 1;      % Generate one signal curve and export it to txt
sensplot   = 0;      % Compute phase sensitivities for the selected parameters
sensDelta  = 1.01;   % Relative perturbation factor used in the sensitivity study

USE_ROM_GEN  = 1;     % 0 = full-order COMSOL, 1 = POD-NNROM for signal generation
USE_ROM_SENS = 0;     % 0 = full-order COMSOL, 1 = POD-NNROM for sensitivity evaluation

SAVE_AS_TARGET_FOR_FIT =0;   % Save the generated phase curve as Data process/1.txt

opts.meshSize = 2;                 % 1(extremely fine) ... 9(coarse)
opts.saveMph  = true;
opts.normalizePumpPower = false;
opts.normalizeProbeW    = false;

runFit           = 0;   % Enable automatic inverse fitting
USE_ROM_FIT      = 1;  % 0 = full-order COMSOL, 1 = POD-NNROM during fitting
RUN_QN_AFTER_PSO = 1;   % Run quasi-Newton refinement after PSO
fitRes = [];
% Isotropy flag: 1 means kr(i)=kz(i)
isotrop = [1 1 0];

% Parameter-selection switches for inverse fitting
 xu_flags = [0 0 1 ;   % kr
            0 0 1  ;   % kz
            0 0 0  ;   % C
            0 0 0  ;   % h
            0 0 0  ];  % [w0,w1,Rpat]


%% ---------------- Auto-fit switches ----------------
% Parameter-selection matrix for fitting (5x3)
% Rows:
%   Row 1: kr(1:3)
%   Row 2: kz(1:3)
%   Row 3: C(1:3)
%   Row 4: h(1:3)
%   Row 5: [w0, w1, Rpat]
%
% Notes:
%   1) kr(2) and C(2) are not active in the current Pattern FDTR model
%   2) h(2) is kept for matrix consistency; use it only when you really want to fit G through kz(2)/h(2)
%   3) If isotrop(i)=1, fit kz(i) only and leave kr(i) off


% ========= Fitting bounds entered in engineering units =========
% kr, kz : W/mK
% C      : MJ/m^3K
% h      : nm
% w0,w1,Rpat : um

fitRange.lb.kr = [30,    0.05,  0.1];
fitRange.ub.kr = [200,  0.2, 10];

fitRange.lb.kz = [30,    0.05,   0.1];
fitRange.ub.kz = [200,  0.2,     10];

fitRange.lb.C  = [1.0,  C(2)/1e6, 0.3];
fitRange.ub.C  = [5.0,  C(2)/1e6, 3];

fitRange.lb.h  = [50,   h(2)/1e-9, 1e5];
fitRange.ub.h  = [150,  h(2)/1e-9, 1e6];

fitRange.lb.spot = [0.8, 0.8, 2];
fitRange.ub.spot = [2, 2, 5];

fitOpts = struct();
fitOpts.dataDir           = fullfile(pwd, 'Data process');
fitOpts.outputsDir        = outputsDir;
fitOpts.USE_ROM_FIT       = USE_ROM_FIT;
fitOpts.runQN             = RUN_QN_AFTER_PSO;
fitOpts.fitness_threshold = 0.005;     
fitOpts.display           = 'iter';
fitOpts.useAmp            = false;
fitOpts.outOfBoundPenalty = 1e6;
fitOpts.saveFitMph        = false;
fitOpts.plotFit           = true;

%  ROM_FIT settings
    fitOpts.swarmSize = 200;
    fitOpts.maxIter   = 20;


%% ---------- Output folders ----------

if ~exist(outputsDir,'dir'), mkdir(outputsDir); end
if ~exist(modelsDir,'dir'),  mkdir(modelsDir);  end

opts.modelsDir = modelsDir;

outTxtFull = fullfile(outputsDir, 'FDTR_COMSOL_phase.txt');
outTxtROM  = fullfile(outputsDir, 'FDTR_ROM_phase.txt');

%% ---------------- Materials ----------------
mat.trans.rho = 2700;
mat.trans.k_r = kr(1);
mat.trans.k_z = kz(1);
mat.trans.Cv  = C(1);

mat.sub.rho = 2200;
mat.sub.k_r = kr(3);
mat.sub.k_z = kz(3);
mat.sub.Cv  = C(3);

%% ====== ROM model load (only if needed) ======
rom_model = [];
if  USE_ROM_GEN == 1 || USE_ROM_SENS == 1 || USE_ROM_FIT == 1
     rom_model_file = fullfile(outputsDir, 'FDTR_PODNNROM_model.mat');
    if exist(rom_model_file, 'file') ~= 2
        error('ROM model not found: %s.', rom_model_file);
    end

    Srom = load(rom_model_file, 'rom_model');
    rom_model = Srom.rom_model;
    fprintf('✅ ROM model loaded: %s\n', rom_model_file);

    if numel(freq) ~= numel(rom_model.freq) || ...
            max(abs(log10(freq(:) ./ rom_model.freq(:)))) > 1e-12
        error('Current freq does not match the ROM training frequencies.');
    end
end


%% ---------------- Build + run ----------------
phi0 = [];
amp0 = [];

if runSim == 1
    if USE_ROM_GEN == 0
        % ===== Full-order COMSOL =====
        fprintf('Using full-order COMSOL model for signal generation.\n');

        model = build_fdtr_model_Fig2(freq, w0, w1, d_trans, d_sub, rMax, Rpat, Pamp, mat, G, opts);
        bTopMeas = sel_top_transducer(model, d_trans, Rpat);
        [phi_deg, amp, dbg] = eval_phase_from_model(model, w1, bTopMeas, opts.normalizeProbeW, Rpat); 

        outTxt = outTxtFull;

    else
        % ===== POD-NNROM (phase-only) =====
        fprintf('Using POD-NNROM model for signal generation (phase-only).\n');

        phi_deg = fdtr_rom_run_case_phase(freq, kz, kr, C, h, w0, w1, Rpat, rMax, rom_model);
        amp     = nan(size(phi_deg));   
        dbg     = struct(); 

        outTxt = outTxtROM;
    end

    data = [freq(:), phi_deg(:), amp(:)];
    writematrix(data, outTxt, 'Delimiter','\t');
    fprintf('✅ Saved: %s\n', outTxt);

    figure;
    semilogx(freq, phi_deg, 'o-'); grid on;
    xlabel('Frequency (Hz)');
    ylabel('\phi (deg)');
    if USE_ROM_GEN == 0
        title('FDTR phase from full-order COMSOL');
    else
        title('FDTR phase from POD-NNROM');
    end

    if SAVE_AS_TARGET_FOR_FIT == 1
        dataProcessDir = fullfile(pwd, 'Data process');
        if ~exist(dataProcessDir, 'dir'), mkdir(dataProcessDir); end
        delete(fullfile(dataProcessDir, '*.txt'));

        targetTxt = fullfile(dataProcessDir, '1.txt');
        writematrix([freq(:), phi_deg(:)], targetTxt, 'Delimiter','\t');
        fprintf('✅ Synthetic target saved: %s\n', targetTxt);
    end

    % Reuse the generated baseline curve in the sensitivity analysis
    phi0 = phi_deg;
    amp0 = amp;
end

%% ---------------- Sensitivity ----------------
if sensplot == 1
    FDTR_Sensplots_unified(freq, kz, kr, C, h, w0, w1, Rpat, rMax, Pamp, ...
        mat, opts, outputsDir, sensDelta, phi0, amp0, USE_ROM_SENS, rom_model);
end

%% ---------------- Auto fit ----------------

if runFit == 1
   fitRes = FDTR_auto_fit_unified(freq, kz, kr, C, h, w0, w1, Rpat, rMax, Pamp, ...
    mat, opts, xu_flags, fitRange, isotrop, fitOpts, rom_model);
end

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

    % -------- Geometry --------
    model.component.create('comp1', true);
    model.component('comp1').geom.create('geom1', 2);
    model.component('comp1').geom('geom1').axisymmetric(true);
    g = model.component('comp1').geom('geom1');

    g.create('sub', 'Rectangle');
    g.feature('sub').set('pos',  {'0', '-d2'});
    g.feature('sub').set('size', {'rMax', 'd2'});

    g.create('tr', 'Rectangle');
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

    if isempty(domT), error('Transducer domain selection failed.'); end
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

    ht.create('bhs1', 'BoundaryHeatSource', 1);
    bTopTr = sel_top_transducer(model, d1, Rpat);
    ht.feature('bhs1').selection.set(bTopTr);

    if isfield(opts,'normalizePumpPower') && opts.normalizePumpPower
        ht.feature('bhs1').set('Qb', 'P*(2/(pi*wo^2))*exp(-2*r^2/wo^2)/(1-exp(-2*Rpat^2/wo^2))');
    else
        ht.feature('bhs1').set('Qb', 'P*(2/(pi*wo^2))*exp(-2*r^2/wo^2)');
    end
    ht.feature('bhs1').set('harmonicPerturbation', true);

    ht.create('tc1', 'ThermalContact', 1);
    bInt = sel_interface_under_transducer(model, Rpat);
    fprintf('bInt IDs = %s\n', mat2str(bInt));
    ht.feature('tc1').selection.set(bInt);
    ht.feature('tc1').set('ContactModel','EquThinLayer');
    ht.feature('tc1').set('Req','TBR');

    bTopSampleOuter = sel_top_sample_outer(model, rMax, Rpat);
    if ~isempty(bTopSampleOuter)
        ht.create('tiTop','ThermalInsulation',1);
        ht.feature('tiTop').selection.set(bTopSampleOuter);
    end

    if Rpat < rMax*(1-1e-12)
        bSideTr = sel_transducer_side_boundary(model, Rpat, d1);
        if ~isempty(bSideTr)
            ht.create('tiSideTr','ThermalInsulation',1);
            ht.feature('tiSideTr').selection.set(bSideTr);
        end
    end

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
% Probe-weight the complex surface temperature and convert it to phase/amplitude.
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

    dbg.r = rS;
    dbg.W = W;
    dbg.denom = denom;
    dbg.Rpat = Rpat;
end

%% --------- Selection helpers ---------
function domT = sel_domain_transducer(model, Rpat, d1)
    
    er_list = max(1e-12, Rpat * [1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]);
    ez_list = max(1e-12, d1   * [1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]);

    domT = [];
    for er = er_list
        for ez = ez_list
            box = [0-er, Rpat+er; 0-ez, d1+ez];
            cand = mphselectbox(model, 'geom1', box, 'domain');
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
            cand = mphselectbox(model, 'geom1', box, 'domain');
            cand = unique(cand(:)).';
            if ~isempty(cand)
                domS = cand;
                return;
            end
        end
    end
end

function b = sel_top_transducer(model, d1, Rpat)
    epsr = max(1e-12, 1e-6*Rpat);
    epsz_list = [5e-10, 1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8];

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



function [phi_deg, amp] = fdtr_comsol_run_case(freq, kz, kr, C, h, w0, w1, Rpat, rMax, Pamp, matBase, opts, saveMphFlag)
    d_trans = h(1);
    d_sub   = h(3);
    G = kz(2)/h(2);

    mat = matBase;
    mat.trans.k_r = kr(1);
    mat.trans.k_z = kz(1);
    mat.trans.Cv  = C(1);

    mat.sub.k_r   = kr(3);
    mat.sub.k_z   = kz(3);
    mat.sub.Cv    = C(3);

    opts2 = opts;
    opts2.saveMph = saveMphFlag;

    model = build_fdtr_model_Fig2(freq, w0, w1, d_trans, d_sub, rMax, Rpat, Pamp, mat, G, opts2);
    bTopMeas = sel_top_transducer(model, d_trans, Rpat);
    [phi_deg, amp] = eval_phase_from_model(model, w1, bTopMeas, opts2.normalizeProbeW, Rpat);
end

function phi_deg = fdtr_rom_run_case_phase(freq, kz, kr, C, h, w0, w1, Rpat, rMax, rom_model)
    phi_deg = predict_fdtr_phase_from_main_params_SI(freq, kz, kr, C, h, w0, w1, rMax, Rpat, rom_model);
end

%% ===== ROM helpers =====
function phi_deg = predict_fdtr_phase_from_main_params_SI(freq, kz, kr, C, h, w0, w1, rMax, Rpat, rom_model)
% Convert the current SI-parameter set into the ROM input vector and predict the phase.
% main units:
%   kz/kr : W/mK
%   C     : J/m^3K
%   h     : m
%   w0/w1/rMax/Rpat : m
%
% ROM input:
%   training used X_SI directly

freq = freq(:);
if numel(freq) ~= numel(rom_model.freq) || ...
        max(abs(log10(freq ./ rom_model.freq(:)))) > 1e-12
    error('The current frequency (freq) is inconsistent with the ROM training frequency.');
end

Nvar = numel(rom_model.varList);
x_all_SI = zeros(1, Nvar);

for j = 1:Nvar
    switch lower(rom_model.varList(j).group)
        case 'kz'
            v = kz(rom_model.varList(j).idx);      % W/mK
        case 'kr'
            v = kr(rom_model.varList(j).idx);      % W/mK
        case 'c'
            v = C(rom_model.varList(j).idx);       % J/m^3K
        case 'h'
            v = h(rom_model.varList(j).idx);       % m
        case 'w0'
            v = w0;                                % m
        case 'w1'
            v = w1;                                % m
        case 'rmax'
            v = rMax;                              % m
        case 'rpat'
            v = Rpat;                              % m
        otherwise
            error('Unknown group in rom_model.varList: %s', rom_model.varList(j).group);
    end
    x_all_SI(j) = v;
end

y_rom = predict_fdtr_rom_from_SI_local(x_all_SI, rom_model);

% current ROM is phase-only
phi_deg = y_rom(1:numel(freq));
phi_deg = phi_deg(:);
end

function y_rom = predict_fdtr_rom_from_SI_local(x_all_SI, rom_model)
Nvar = numel(rom_model.varList);
xvec = zeros(Nvar,1);

for j = 1:Nvar
    if strcmpi(rom_model.varList(j).mode, 'log')
        xvec(j) = log10(max(x_all_SI(j), realmin));
    else
        xvec(j) = x_all_SI(j);
    end
end

x_norm = mapminmax('apply', xvec, rom_model.Xsettings);
y_norm = rom_model.net(x_norm);
coeffs = mapminmax('reverse', y_norm, rom_model.Ysettings);
coeffs = coeffs(:);

y_rom = rom_model.POD_basis * coeffs + rom_model.snapshot_mean;
end



%% ========================= AUTO-FIT FUNCTIONS =========================

function fitRes = FDTR_auto_fit_unified(freq_default, kz, kr, C, h, w0, w1, Rpat, rMax, Pamp, ...
    matBase, opts, xu_flags, fitRange, isotrop, fitOpts, rom_model)
% Main inverse-fitting driver: load target curves, run PSO/QN, and save results.
   
    Nlayers = numel(h);
    if Nlayers ~= 3
        error('The current automatic fitting version is written based on a 3-layer PatternFDTR. Please keep the length of h/kr/kz/C to 3.');
    end

    check_xu_flags_fdtr(xu_flags, isotrop);

    targets = load_fdtr_target_files(fitOpts.dataDir);
    nFiles = numel(targets);
    fprintf('✅ Loaded %d target file(s) from: %s\n', nFiles, fitOpts.dataDir);

    
    % kr, kz : W/mK
    % C      : MJ/m^3K
    % h      : nm
    % w0,w1,Rpat : um
    all_params_initial_fit = pack_fdtr_params_fit_units(kr, kz, C, h, w0, w1, Rpat);

    
    lb_all = pack_fdtr_bounds_fit_units(fitRange.lb);
    ub_all = pack_fdtr_bounds_fit_units(fitRange.ub);

    flag_vec = reshape(xu_flags.', 1, []);

    if numel(flag_vec) ~= numel(all_params_initial_fit)
        error('The flattened length of xu_flags is inconsistent with the total length of the parameters.');
    end

    params_to_fit_index = find(flag_vec == 1);
    if isempty(params_to_fit_index)
        error('xu_flags has no parameters to be fitted.');
    end

    params_initial = all_params_initial_fit(params_to_fit_index);
    params_lb      = lb_all(params_to_fit_index);
    params_ub      = ub_all(params_to_fit_index);

    if any(params_initial < params_lb) || any(params_initial > params_ub)
        error('initial values ​​have parameters that fall outside the fitting boundary. Please check fitRange.lb/fitRange.ub.');
    end

    param_names_all = get_fdtr_param_names_fit_units();
    param_names_fit = param_names_all(params_to_fit_index);

    fprintf('\n========== FDTR Auto Fit ==========\n');
    if fitOpts.USE_ROM_FIT == 1
        fprintf('Forward model : POD-NNROM\n');
    else
        fprintf('Forward model : Full-order COMSOL\n');
    end
    fprintf('Fitting-space units:\n');
    fprintf('  kr, kz          : W/mK\n');
    fprintf('  C               : MJ/m^3K\n');
    fprintf('  h               : nm\n');
    fprintf('  w0, w1, Rpat    : um\n');
    fprintf('Target metric : mean phase RMSE (deg)\n');
    fprintf('Parameters to fit:\n');
    for i = 1:numel(params_to_fit_index)
        fprintf('  %-24s init = %.6g,  lb = %.6g,  ub = %.6g\n', ...
            param_names_fit{i}, params_initial(i), params_lb(i), params_ub(i));
    end
    fprintf('===================================\n\n');

    optsFit = opts;
    optsFit.saveMph = false;  

    objfun = @(p) computeResiduals_FDTR( ...
        p, params_to_fit_index, all_params_initial_fit, isotrop, targets, ...
        rMax, Pamp, matBase, optsFit, fitOpts.USE_ROM_FIT, rom_model, ...
        params_lb, params_ub, fitOpts.outOfBoundPenalty);

    % ---------- PSO ----------
    swarm0 = params_lb + rand(fitOpts.swarmSize, numel(params_initial)) .* (params_ub - params_lb);

    options_pso = optimoptions('particleswarm', ...
        'SwarmSize', fitOpts.swarmSize, ...
        'MaxIterations', fitOpts.maxIter, ...
        'InitialSwarmMatrix', swarm0, ...
        'Display', fitOpts.display, ...
        'FunctionTolerance', 1e-12, ...
        'UseParallel', false, ...
        'OutputFcn', @(optimValues, state) stopPSOIfThresholdReached_FDTR( ...
            optimValues, state, fitOpts.fitness_threshold));

    [best_param_fit_pso, best_fitness_pso, ~, output_pso] = particleswarm( ...
        objfun, numel(params_initial), params_lb, params_ub, options_pso);

    fprintf('\nPSO result:\n');
    for i = 1:numel(best_param_fit_pso)
        fprintf('  %-24s = %.10g\n', param_names_fit{i}, best_param_fit_pso(i));
    end
    fprintf('  mean phase RMSE = %.10g deg\n', best_fitness_pso);

    best_param_fit = best_param_fit_pso;
    best_fitness   = best_fitness_pso;
    best_stage     = 'PSO';

    % ---------- Quasi-Newton ----------
    if fitOpts.runQN
        fprintf('\nSwitching to quasi-Newton refinement...\n');

        options_qn = optimoptions('fminunc', ...
            'Algorithm', 'quasi-newton', ...
            'Display', 'iter', ...
            'OptimalityTolerance', 1e-12, ...
            'StepTolerance', 1e-12, ...
            'MaxIterations', 100);

        [best_param_fit_qn, best_fitness_qn] = fminunc(objfun, best_param_fit_pso, options_qn);

        best_param_fit_qn = min(max(best_param_fit_qn, params_lb), params_ub);
        best_fitness_qn   = objfun(best_param_fit_qn);

        fprintf('\nQuasi-Newton result:\n');
        for i = 1:numel(best_param_fit_qn)
            fprintf('  %-24s = %.10g\n', param_names_fit{i}, best_param_fit_qn(i));
        end
        fprintf('  mean phase RMSE = %.10g deg\n', best_fitness_qn);

        if best_fitness_qn < best_fitness_pso
            best_param_fit = best_param_fit_qn;
            best_fitness   = best_fitness_qn;
            best_stage     = 'PSO + QN';
        end
    end

    % ---------- Recover best full parameter set ----------
    all_params_best_fit = all_params_initial_fit;
    all_params_best_fit(params_to_fit_index) = best_param_fit;

    [kr_fit, kz_fit, C_fit_MJ, h_fit_nm, w0_fit_um, w1_fit_um, Rpat_fit_um, ...
        kr_si, kz_si, C_si, h_si, w0_si, w1_si, Rpat_si] = ...
        unpack_fdtr_params_fit_units(all_params_best_fit, isotrop);

    G_fit = kz_si(2) / h_si(2);

    fprintf('\nFinal accepted result (%s):\n', best_stage);
    fprintf('  kr   (W/mK)      = [%g  %g  %g]\n', kr_fit(1), kr_fit(2), kr_fit(3));
    fprintf('  kz   (W/mK)      = [%g  %g  %g]\n', kz_fit(1), kz_fit(2), kz_fit(3));
    fprintf('  C    (MJ/m^3K)   = [%g  %g  %g]\n', C_fit_MJ(1), C_fit_MJ(2), C_fit_MJ(3));
    fprintf('  h    (nm)        = [%g  %g  %g]\n', h_fit_nm(1), h_fit_nm(2), h_fit_nm(3));
    fprintf('  w0   (um)        = %g\n', w0_fit_um);
    fprintf('  w1   (um)        = %g\n', w1_fit_um);
    fprintf('  Rpat (um)        = %g\n', Rpat_fit_um);
    fprintf('  G    (W/m^2/K)   = %.6e\n', G_fit);
    fprintf('  mean phase RMSE  = %.10g deg\n', best_fitness);

    % ---------- Generate best-fit curves ----------
    simCurves = cell(nFiles,1);
    for i = 1:nFiles
        fi = targets(i).freq;

        if fitOpts.USE_ROM_FIT == 1
            sim_phi = fdtr_rom_run_case_phase_fit(fi, kz_si, kr_si, C_si, h_si, ...
                w0_si, w1_si, Rpat_si, rMax, rom_model);
            sim_amp = nan(size(sim_phi));
        else
            [sim_phi, sim_amp] = fdtr_comsol_run_case(fi, kz_si, kr_si, C_si, h_si, ...
                w0_si, w1_si, Rpat_si, rMax, Pamp, matBase, optsFit, false);
        end

        simCurves{i}.freq = fi;
        simCurves{i}.phi  = sim_phi(:);
        simCurves{i}.amp  = sim_amp(:);

        outCurve = fullfile(fitOpts.outputsDir, sprintf('FDTR_fit_curve_%d.txt', i));
        writematrix([fi(:), targets(i).phi(:), sim_phi(:)], outCurve, 'Delimiter','\t');
    end

    % ---------- Save summary ----------
    fitRes = struct();
    fitRes.best_stage          = best_stage;
    fitRes.best_fitness        = best_fitness;
    fitRes.params_to_fit_index = params_to_fit_index;
    fitRes.param_names_fit     = {param_names_fit{:}};
    fitRes.best_param_fit      = best_param_fit;
    fitRes.output_pso          = output_pso;

    fitRes.kr_fit_WmK     = kr_fit;
    fitRes.kz_fit_WmK     = kz_fit;
    fitRes.C_fit_MJm3K    = C_fit_MJ;
    fitRes.h_fit_nm       = h_fit_nm;
    fitRes.w0_fit_um      = w0_fit_um;
    fitRes.w1_fit_um      = w1_fit_um;
    fitRes.Rpat_fit_um    = Rpat_fit_um;
    fitRes.G_fit_Wm2K     = G_fit;

    fitRes.kr_fit_SI      = kr_si;
    fitRes.kz_fit_SI      = kz_si;
    fitRes.C_fit_SI       = C_si;
    fitRes.h_fit_SI       = h_si;
    fitRes.w0_fit_SI      = w0_si;
    fitRes.w1_fit_SI      = w1_si;
    fitRes.Rpat_fit_SI    = Rpat_si;

    fitRes.targets        = targets;
    fitRes.simCurves      = simCurves;

    save(fullfile(fitOpts.outputsDir, 'FDTR_fit_result.mat'), 'fitRes');

    txtout = fullfile(fitOpts.outputsDir, 'FDTR_fit_result.txt');
    fid = fopen(txtout, 'w');
    fprintf(fid, 'FDTR auto-fit summary\n');
    fprintf(fid, 'Forward model: %s\n', ternary_str(fitOpts.USE_ROM_FIT==1, 'POD-NNROM', 'Full-order COMSOL'));
    fprintf(fid, 'Best stage   : %s\n', best_stage);
    fprintf(fid, 'Metric       : mean phase RMSE (deg)\n');
    fprintf(fid, 'Best fitness : %.12g\n\n', best_fitness);

    fprintf(fid, 'Fitted subset (fitting-space units):\n');
    for i = 1:numel(best_param_fit)
        fprintf(fid, '  %-24s = %.12g\n', param_names_fit{i}, best_param_fit(i));
    end

    fprintf(fid, '\nRecovered full parameters:\n');
    fprintf(fid, 'kr   (W/mK)      = [%.12g  %.12g  %.12g]\n', kr_fit(1), kr_fit(2), kr_fit(3));
    fprintf(fid, 'kz   (W/mK)      = [%.12g  %.12g  %.12g]\n', kz_fit(1), kz_fit(2), kz_fit(3));
    fprintf(fid, 'C    (MJ/m^3K)   = [%.12g  %.12g  %.12g]\n', C_fit_MJ(1), C_fit_MJ(2), C_fit_MJ(3));
    fprintf(fid, 'h    (nm)        = [%.12g  %.12g  %.12g]\n', h_fit_nm(1), h_fit_nm(2), h_fit_nm(3));
    fprintf(fid, 'w0   (um)        = %.12g\n', w0_fit_um);
    fprintf(fid, 'w1   (um)        = %.12g\n', w1_fit_um);
    fprintf(fid, 'Rpat (um)        = %.12g\n', Rpat_fit_um);
    fprintf(fid, 'G    (W/m^2/K)   = %.12g\n', G_fit);
    fclose(fid);

    fprintf('✅ Fit summary saved:\n  %s\n  %s\n', ...
        fullfile(fitOpts.outputsDir, 'FDTR_fit_result.mat'), txtout);

    % ---------- Plot ----------
    if isfield(fitOpts, 'plotFit') && fitOpts.plotFit
        figure('Name', 'FDTR Auto Fit');
        hold on;
        for i = 1:nFiles
            semilogx(targets(i).freq, targets(i).phi, 'o', 'DisplayName', sprintf('Exp %d', i));
            semilogx(simCurves{i}.freq, simCurves{i}.phi, '-', 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Fit %d', i));
        end
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('\phi (deg)');
        title(sprintf('FDTR auto-fit (%s), mean phase RMSE = %.4g deg', best_stage, best_fitness));
        legend('Location','best');
    end
end

function F = computeResiduals_FDTR(p, idx_to_fit, all_params_initial_fit, isotrop, targets, ...
    rMax, Pamp, matBase, optsFit, USE_ROM_FIT, rom_model, lb, ub, penaltyBase)

    vioL = max(lb - p, 0);
    vioU = max(p - ub, 0);
    vio  = vioL + vioU;
    if any(vio > 0)
        F = penaltyBase + 1e8 * sum(vio.^2);
        return;
    end

    all_params_fit = all_params_initial_fit;
    all_params_fit(idx_to_fit) = p;

    try
        [~, ~, ~, ~, ~, ~, ~, kr_si, kz_si, C_si, h_si, w0_si, w1_si, Rpat_si] = ...
            unpack_fdtr_params_fit_units(all_params_fit, isotrop);

    
        if any(kr_si <= 0) || any(kz_si <= 0) || any(C_si <= 0) || any(h_si <= 0) || ...
           w0_si <= 0 || w1_si <= 0 || Rpat_si <= 0 || Rpat_si >= rMax
            F = penaltyBase;
            return;
        end

        G = kz_si(2) / h_si(2);
        if ~isfinite(G) || G <= 0
            F = penaltyBase;
            return;
        end

        nFiles = numel(targets);
        rmse_list = nan(nFiles,1);

        for i = 1:nFiles
            fi = targets(i).freq;

            if USE_ROM_FIT == 1
                sim_phi = fdtr_rom_run_case_phase_fit(fi, kz_si, kr_si, C_si, h_si, ...
                    w0_si, w1_si, Rpat_si, rMax, rom_model);
            else
                [sim_phi, ~] = fdtr_comsol_run_case(fi, kz_si, kr_si, C_si, h_si, ...
                    w0_si, w1_si, Rpat_si, rMax, Pamp, matBase, optsFit, false);
            end

            dv = wrap_phase_diff_deg(sim_phi(:) - targets(i).phi(:));
            rmse_list(i) = sqrt(mean(dv.^2));
        end

        F = mean(rmse_list);

        if ~isfinite(F)
            F = penaltyBase;
        end

    catch
        F = penaltyBase;
    end
end

function stop = stopPSOIfThresholdReached_FDTR(optimValues, state, fitness_threshold)
    stop = false;
    if strcmpi(state, 'iter')
        if optimValues.bestfval <= fitness_threshold
            fprintf('Mean phase RMSE reached threshold. Stop PSO.\n');
            stop = true;
        end
    end
end

function targets = load_fdtr_target_files(dataDir)
    files = dir(fullfile(dataDir, '*.txt'));
    if isempty(files)
        error('No .txt files were found in the Data process folder：%s', dataDir);
    end

    fileNumbers = nan(numel(files),1);
    pureNumeric = true;
    for i = 1:numel(files)
        [~, name] = fileparts(files(i).name);
        fileNumbers(i) = str2double(name);
        if isnan(fileNumbers(i))
            pureNumeric = false;
        end
    end

    if pureNumeric
        [~, idx] = sort(fileNumbers);
    else
        [~, idx] = sort({files.name});
    end
    files = files(idx);

    targets = struct('name',{},'freq',{},'phi',{},'amp',{});
    for i = 1:numel(files)
        fp = fullfile(dataDir, files(i).name);
        A = readmatrix(fp);

        if size(A,2) < 2
            error('The %s file must have at least two columns: [freq(Hz), phase(deg)]', files(i).name);
        end

        freq_i = A(:,1);
        phi_i  = A(:,2);

        good = isfinite(freq_i) & isfinite(phi_i);
        freq_i = freq_i(good);
        phi_i  = phi_i(good);

        [freq_i, sidx] = sort(freq_i(:));
        phi_i = phi_i(sidx);

        amp_i = [];
        if size(A,2) >= 3
            amp_raw = A(good,3);
            amp_i = amp_raw(sidx);
        end

        targets(i).name = files(i).name;
        targets(i).freq = freq_i(:);
        targets(i).phi  = phi_i(:);
        targets(i).amp  = amp_i(:);
    end
end

function all_params_fit = pack_fdtr_params_fit_units(kr_si, kz_si, C_si, h_si, w0_si, w1_si, Rpat_si)
    C_fit_MJ = C_si(:).' / 1e6;
    h_fit_nm = h_si(:).' / 1e-9;
    spot_fit = [w0_si, w1_si, Rpat_si] / 1e-6;

    all_params_fit = [kr_si(:).', kz_si(:).', C_fit_MJ, h_fit_nm, spot_fit];
end

function all_bounds_fit = pack_fdtr_bounds_fit_units(S)
    all_bounds_fit = [S.kr(:).', S.kz(:).', S.C(:).', S.h(:).', S.spot(:).'];
end

function [kr_fit, kz_fit, C_fit_MJ, h_fit_nm, w0_fit_um, w1_fit_um, Rpat_fit_um, ...
          kr_si, kz_si, C_si, h_si, w0_si, w1_si, Rpat_si] = ...
          unpack_fdtr_params_fit_units(all_params_fit, isotrop)

    Nlayers = numel(isotrop);

    kr_fit = all_params_fit(1:Nlayers);
    kz_fit = all_params_fit(Nlayers+1:2*Nlayers);
    C_fit_MJ = all_params_fit(2*Nlayers+1:3*Nlayers);
    h_fit_nm = all_params_fit(3*Nlayers+1:4*Nlayers);

    aux = all_params_fit(4*Nlayers+1:4*Nlayers+3);
    w0_fit_um   = aux(1);
    w1_fit_um   = aux(2);
    Rpat_fit_um = aux(3);

    for i = 1:Nlayers
        if isotrop(i) == 1
            kr_fit(i) = kz_fit(i);
        end
    end

   
    kr_si   = kr_fit;
    kz_si   = kz_fit;
    C_si    = C_fit_MJ * 1e6;
    h_si    = h_fit_nm * 1e-9;
    w0_si   = w0_fit_um * 1e-6;
    w1_si   = w1_fit_um * 1e-6;
    Rpat_si = Rpat_fit_um * 1e-6;
end

function names = get_fdtr_param_names_fit_units()
    names = { ...
        'kr(1) [W/mK]', ...
        'kr(2) [unused]', ...
        'kr(3) [W/mK]', ...
        'kz(1) [W/mK]', ...
        'kz(2) [W/mK]', ...
        'kz(3) [W/mK]', ...
        'C(1) [MJ/m^3K]', ...
        'C(2) [unused]', ...
        'C(3) [MJ/m^3K]', ...
        'h(1) [nm]', ...
        'h(2) [nm]', ...
        'h(3) [nm]', ...
        'w0 [um]', ...
        'w1 [um]', ...
        'Rpat [um]'};
end

function check_xu_flags_fdtr(xu_flags, isotrop)
    if ~isequal(size(xu_flags), [5,3])
        error('xu_flags must be a 5×3 matrix.');
    end

    if xu_flags(1,2) == 1
        error('kr(2) is not used in the current PatternFDTR model and cannot be set as a parameter to be fitted.');
    end
    if xu_flags(3,2) == 1
        error('C(2) is not used in the current PatternFDTR model and cannot be set as a parameter to be fitted.');
    end

    for i = 1:3
        if isotrop(i) == 1 && xu_flags(1,i) == 1
            error('When isotrop=1 for layer %d, do not enable the fitting switch for kr(%d), only enable kz(%d).', i, i, i);
        end
    end

    if xu_flags(2,2) == 1 && xu_flags(4,2) == 1
        warning('G = kz(2)/h(2), you have enabled fitting for both kz(2) and h(2), which is usually not discernible. It is recommended to fit only one of them.');
    end
end

function dv = wrap_phase_diff_deg(dphi)
    dv = mod(dphi + 180, 360) - 180;
end

function phi_deg = fdtr_rom_run_case_phase_fit(freq_target, kz, kr, C, h, w0, w1, Rpat, rMax, rom_model)
    freq_target = freq_target(:);
    freq_rom    = rom_model.freq(:);

    phi_full = predict_fdtr_phase_from_main_params_SI(freq_rom, kz, kr, C, h, w0, w1, rMax, Rpat, rom_model);

    if numel(freq_target) == numel(freq_rom)
        ratio_err = max(abs(log10(freq_target ./ freq_rom)));
        if ratio_err <= 1e-12
            phi_deg = phi_full(:);
            return;
        end
    end

    if min(freq_target) < min(freq_rom) || max(freq_target) > max(freq_rom)
        error('The frequency to be fitted exceeds the training frequency range of the ROM.');
    end

    phi_deg = interp1(log10(freq_rom), phi_full(:), log10(freq_target), 'pchip');
end

function s = ternary_str(cond, a, b)
    if cond
        s = a;
    else
        s = b;
    end
end