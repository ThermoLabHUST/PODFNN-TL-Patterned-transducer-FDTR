%% FDTR POD-FNN training
% This script trains a POD-based neural-network reduced-order model
% (POD-FNNROM) for FDTR signal prediction using a precomputed COMSOL dataset.
%
% Main workflow:
% 1) load the FDTR_LHS_dataset.mat dataset;
% 2) build phase-only or phase+amplitude snapshots;
% 3) perform POD/SVD and retain the dominant modes;
% 4) train a feedforward neural network to map input parameters to POD coefficients;
% 5) export POD mean/modes and save the trained ROM model.
%
% Main output:
%   Outputs/FDTR_PODNNROM_model.mat
%
% Optional output:
%   Outputs/POD_exports/   % exported POD mean and retained modes
%
% Notes:
% - The dataset input X_SI is used directly as the ROM input.
% - If varList(j).mode = 'log', log10 scaling is applied before normalization.

clc; clear; close all;

%% ===== paths =====
outputsDir  = fullfile(pwd, 'Outputs');
datasetFile = fullfile(outputsDir, 'FDTR_LHS_dataset.mat');
 modelFile   = fullfile(outputsDir, 'FDTR_PODNNROM_model.mat');

reuseModelIfExists = false;   % true: If the model already exists, load it directly and skip training.
useAmpInSnapshot   = false;    % true: snapshot = [phase; log10(amp)],false:only trains phase
ampTransform       = 'log10';
ampFloor           = 1e-30;

%% ===== POD / NN options =====
podOpt.energyThres = 0.999999;
podOpt.maxModes    = 20;

nnOpt.hidden   = [32 32];
nnOpt.trainFcn = 'trainbr';  
nnOpt.epochs   = 500;

%% ===== load dataset =====
if exist(datasetFile, 'file') ~= 2
    error('Dataset not found: %s', datasetFile);
end

S = load(datasetFile, 'D');
D = S.D;

needFields = {'freq','varList','X_SI','okMask','phi_deg_all','amp_all'};
for k = 1:numel(needFields)
    if ~isfield(D, needFields{k})
        error('Dataset missing fields D.%s', needFields{k});
    end
end

ok = D.okMask(:);
if nnz(ok) < 10
    error('With too few successful samples (%d), training the NN is not recommended.', nnz(ok));
end

X = D.X_SI(ok, :);                  
phi_all = D.phi_deg_all(:, ok);     
amp_all = D.amp_all(:, ok);         
freq = D.freq(:);
varList = D.varList;

if useAmpInSnapshot
    if strcmpi(ampTransform, 'log10')
        snapshots = [phi_all; log10(max(amp_all, ampFloor))];
    else
        snapshots = [phi_all; amp_all];
    end
else
    snapshots = phi_all;
end

fprintf('Loaded dataset:\n');
fprintf('  success samples = %d\n', size(X,1));
fprintf('  Nvar            = %d\n', size(X,2));
fprintf('  Nfreq           = %d\n', numel(freq));
fprintf('  snapshot dim    = %d\n', size(snapshots,1));

%% ===== build POD =====
[POD_basis, singular_values, pod_coeffs, snapshot_mean] = ...
    build_POD_basis_fdtr(snapshots, X, podOpt);

num_modes_current = size(POD_basis, 2);

%% ===== export mean + all POD modes to txt =====
podOutDir = fullfile(outputsDir, 'POD_exports');
if ~exist(podOutDir, 'dir')
    mkdir(podOutDir);
end

Nf = numel(freq);
num_modes = size(POD_basis, 2);


energy_each = singular_values.^2 / sum(singular_values.^2);
energy_cum  = cumsum(energy_each);


energy_each_sel = energy_each(1:num_modes);
energy_cum_sel  = energy_cum(1:num_modes);
sv_sel          = singular_values(1:num_modes);


summaryFile = fullfile(podOutDir, 'POD_summary.txt');
fid = fopen(summaryFile, 'w');
if fid < 0
    error('Unable to create file:%s', summaryFile);
end

fprintf(fid, 'POD export summary\n');
fprintf(fid, 'useAmpInSnapshot = %d\n', useAmpInSnapshot);
fprintf(fid, 'ampTransform     = %s\n', ampTransform);
fprintf(fid, 'Nfreq            = %d\n', Nf);
fprintf(fid, 'num_modes        = %d\n', num_modes);
fprintf(fid, 'energy_threshold = %.8f\n', podOpt.energyThres);
fprintf(fid, 'max_modes        = %d\n\n', podOpt.maxModes);

fprintf(fid, 'mode_index\tsingular_value\tenergy_percent\tcumulative_energy_percent\n');
for m = 1:num_modes
    fprintf(fid, '%d\t%.12e\t%.8f\t%.8f\n', ...
        m, sv_sel(m), 100*energy_each_sel(m), 100*energy_cum_sel(m));
end
fclose(fid);

fprintf('✅ Exported: %s\n', summaryFile);

% ---- current setting: phase only ----
if ~useAmpInSnapshot
    % mean phase
    meanPhaseFile = fullfile(podOutDir, 'POD_mean_phase.txt');
    meanPhaseData = [freq(:), snapshot_mean(:)];
    writematrix(meanPhaseData, meanPhaseFile, 'Delimiter', 'tab');
    fprintf('✅ Exported: %s\n', meanPhaseFile);

    % all phase modes
    for m = 1:num_modes
        modeFile = fullfile(podOutDir, sprintf('POD_mode_%03d_phase.txt', m));
        modeData = [freq(:), POD_basis(:,m)];
        writematrix(modeData, modeFile, 'Delimiter', 'tab');
        fprintf('✅ Exported: %s\n', modeFile);
    end

else
    % ---- phase + amp snapshot case ----
    % snapshot = [phase; aux]
    % aux means log10(amp) or amp, depending on ampTransform
    mean_phase = snapshot_mean(1:Nf);
    mean_aux   = snapshot_mean(Nf+1:end);

    meanPhaseFile = fullfile(podOutDir, 'POD_mean_phase.txt');
    writematrix([freq(:), mean_phase(:)], meanPhaseFile, 'Delimiter', 'tab');
    fprintf('✅ Exported:%s\n', meanPhaseFile);

    if strcmpi(ampTransform, 'log10')
        meanAuxFile = fullfile(podOutDir, 'POD_mean_log10amp.txt');
    else
        meanAuxFile = fullfile(podOutDir, 'POD_mean_amp.txt');
    end
    writematrix([freq(:), mean_aux(:)], meanAuxFile, 'Delimiter', 'tab');
    fprintf('✅ Exported: %s\n', meanAuxFile);

    for m = 1:num_modes
        mode_phase = POD_basis(1:Nf, m);
        mode_aux   = POD_basis(Nf+1:end, m);

        modePhaseFile = fullfile(podOutDir, sprintf('POD_mode_%03d_phase.txt', m));
        writematrix([freq(:), mode_phase(:)], modePhaseFile, 'Delimiter', 'tab');
        fprintf('✅ Exported: %s\n', modePhaseFile);

        if strcmpi(ampTransform, 'log10')
            modeAuxFile = fullfile(podOutDir, sprintf('POD_mode_%03d_log10amp.txt', m));
        else
            modeAuxFile = fullfile(podOutDir, sprintf('POD_mode_%03d_amp.txt', m));
        end
        writematrix([freq(:), mode_aux(:)], modeAuxFile, 'Delimiter', 'tab');
        fprintf('✅ Exported: %s\n', modeAuxFile);
    end
end

% ---- optional: export all phase modes into one matrix file for convenience ----
% column 1 = freq, column 2 = mean, column 3... = modes
if ~useAmpInSnapshot
    allModesPhaseFile = fullfile(podOutDir, 'POD_all_phase_modes_in_one.txt');
    allModesPhaseData = [freq(:), snapshot_mean(:), POD_basis];
    writematrix(allModesPhaseData, allModesPhaseFile, 'Delimiter', 'tab');
    fprintf('✅ Exported:%s\n', allModesPhaseFile);
else
    allModesPhaseFile = fullfile(podOutDir, 'POD_all_phase_modes_in_one.txt');
    allModesPhaseData = [freq(:), snapshot_mean(1:Nf), POD_basis(1:Nf,:)];
    writematrix(allModesPhaseData, allModesPhaseFile, 'Delimiter', 'tab');
    fprintf('✅ Exported:%s\n', allModesPhaseFile);

    if strcmpi(ampTransform, 'log10')
        allModesAuxFile = fullfile(podOutDir, 'POD_all_log10amp_modes_in_one.txt');
    else
        allModesAuxFile = fullfile(podOutDir, 'POD_all_amp_modes_in_one.txt');
    end
    allModesAuxData = [freq(:), snapshot_mean(Nf+1:end), POD_basis(Nf+1:end,:)];
    writematrix(allModesAuxData, allModesAuxFile, 'Delimiter', 'tab');
    fprintf('✅ Exported:%s\n', allModesAuxFile);
end



%% ===== load / train NN =====
if reuseModelIfExists && exist(modelFile, 'file') == 2
    T = load(modelFile, 'rom_model');
    rom_model = T.rom_model;

    if isfield(rom_model,'num_modes') && rom_model.num_modes == num_modes_current
        fprintf('✅ existing PODNN model has been detected and will be loaded directly:%s\n', modelFile);
    else
        fprintf('The existing model has a mismatch in the number of modalities; retraining is required.\n');
        rom_model = train_fdtr_rom_model(X, pod_coeffs, varList, freq, ...
            POD_basis, snapshot_mean, singular_values, ...
            useAmpInSnapshot, ampTransform, ampFloor, nnOpt);
        save(modelFile, 'rom_model');
    end
else
    fprintf('Start training POD-NNROM...\n');
    rom_model = train_fdtr_rom_model(X, pod_coeffs, varList, freq, ...
        POD_basis, snapshot_mean, singular_values, ...
        useAmpInSnapshot, ampTransform, ampFloor, nnOpt);
    save(modelFile, 'rom_model');
end

fprintf('✅ The model has been saved to %s\n', modelFile);

%% ===== training-set reconstruction check =====
Y_pred = zeros(size(snapshots));
for i = 1:size(X,1)
    Y_pred(:,i) = predict_fdtr_rom_from_SI(X(i,:), rom_model);
end

Nf = numel(freq);
if useAmpInSnapshot
    phi_pred = Y_pred(1:Nf,:);
    aux_pred = Y_pred(Nf+1:end,:);
    phi_true = snapshots(1:Nf,:);
    aux_true = snapshots(Nf+1:end,:);

    rmse_phi = mean(sqrt(mean((phi_pred - phi_true).^2, 1)));
    r2_phi   = mean(1 - sum((phi_pred - phi_true).^2,1) ./ sum((phi_true - mean(phi_true,1)).^2,1));

    rmse_aux = mean(sqrt(mean((aux_pred - aux_true).^2, 1)));
    fprintf('\n===== Training-set ROM check =====\n');
    fprintf('Phase mean RMSE      = %.4e deg\n', rmse_phi);
    fprintf('Phase mean R2        = %.4f\n', r2_phi);
    if strcmpi(ampTransform,'log10')
        fprintf('log10(Amp) mean RMSE = %.4e\n', rmse_aux);
    else
        fprintf('Amp mean RMSE        = %.4e\n', rmse_aux);
    end
else
    phi_pred = Y_pred;
    phi_true = snapshots;
    rmse_phi = mean(sqrt(mean((phi_pred - phi_true).^2, 1)));
    r2_phi   = mean(1 - sum((phi_pred - phi_true).^2,1) ./ sum((phi_true - mean(phi_true,1)).^2,1));
    fprintf('\n===== Training-set ROM check =====\n');
    fprintf('Phase mean RMSE = %.4e deg\n', rmse_phi);
    fprintf('Phase mean R2   = %.4f\n', r2_phi);
end

fprintf('✅ FDTR POD-NNROM training finished.\n');


%% ========================= helpers =========================
function [POD_basis, singular_values, pod_coeffs, snapshot_mean] = ...
    build_POD_basis_fdtr(snapshots, X, podOpt)

fprintf('\n========== Building POD basis ==========\n');
snapshot_mean = mean(snapshots, 2);
Xc = snapshots - snapshot_mean;

[U, S, ~] = svd(Xc, 'econ');
singular_values = diag(S);

energy = cumsum(singular_values.^2) / sum(singular_values.^2);
num_modes = find(energy >= podOpt.energyThres, 1);
num_modes = max(num_modes, 2);
num_modes = min(num_modes, podOpt.maxModes);

POD_basis = U(:,1:num_modes);
pod_coeffs = POD_basis' * Xc;

fprintf('Number of retained modes = %d, cumulative energy = %.5f%%\n', num_modes, energy(num_modes)*100);

for m = 1:num_modes
    coeff_m = pod_coeffs(m,:);
    [max_val, idx] = max(abs(coeff_m));
    fprintf('  mode %2d: energy = %.4f%%, |coeff|max = %.3e, sample #%d\n', ...
        m, 100*singular_values(m)^2/sum(singular_values.^2), max_val, idx);

    for j = 1:size(X,2)
        fprintf('    x%-2d = %.4e\n', j, X(idx,j));
    end
end
fprintf('========================================\n\n');
end

function rom_model = train_fdtr_rom_model(X, pod_coeffs, varList, freq, ...
    POD_basis, snapshot_mean, singular_values, ...
    useAmpInSnapshot, ampTransform, ampFloor, nnOpt)

[num_modes, Ns] = size(pod_coeffs);
Nvar = size(X,2);

fprintf('\n========== Training NN ==========\n');
fprintf('Nvar = %d, Ns = %d, Nmodes = %d\n', Nvar, Ns, num_modes);

Xfeat = build_feature_matrix(X, varList);   

[Xn, Xsettings] = mapminmax(Xfeat, 0, 1);
[Yn, Ysettings] = mapminmax(pod_coeffs, 0, 1);

net = feedforwardnet(nnOpt.hidden, nnOpt.trainFcn);
net.trainParam.epochs = nnOpt.epochs;
net.divideParam.trainRatio = 0.8;
net.divideParam.valRatio   = 0.1;
net.divideParam.testRatio  = 0.1;

net = train(net, Xn, Yn);

Yn_pred = net(Xn);
Y_pred  = mapminmax('reverse', Yn_pred, Ysettings);

RMSE = sqrt(mean((pod_coeffs - Y_pred).^2, 2));
R2   = 1 - sum((pod_coeffs - Y_pred).^2, 2) ./ sum((pod_coeffs - mean(pod_coeffs,2)).^2);

for i = 1:num_modes
    fprintf('  mode %2d: RMSE = %.3e, R2 = %.4f\n', i, RMSE(i), R2(i));
end
fprintf('=================================\n\n');

rom_model = struct();
rom_model.net             = net;
rom_model.Xsettings       = Xsettings;
rom_model.Ysettings       = Ysettings;
rom_model.num_modes       = num_modes;
rom_model.varList         = varList;
rom_model.freq            = freq(:);
rom_model.POD_basis       = POD_basis;
rom_model.snapshot_mean   = snapshot_mean;
rom_model.singular_values = singular_values;
rom_model.useAmp          = useAmpInSnapshot;
rom_model.ampTransform    = ampTransform;
rom_model.ampFloor        = ampFloor;
end

function Xfeat = build_feature_matrix(X, varList)
Ns = size(X,1);
Nvar = size(X,2);
Xfeat = zeros(Nvar, Ns);

for j = 1:Nvar
    if strcmpi(varList(j).mode, 'log')
        Xfeat(j,:) = log10(max(X(:,j), realmin))';
    else
        Xfeat(j,:) = X(:,j)';
    end
end
end

function y_rom = predict_fdtr_rom_from_SI(x_row_SI, rom_model)
Nvar = numel(rom_model.varList);
xvec = zeros(Nvar,1);

for j = 1:Nvar
    if strcmpi(rom_model.varList(j).mode, 'log')
        xvec(j) = log10(max(x_row_SI(j), realmin));
    else
        xvec(j) = x_row_SI(j);
    end
end

x_norm = mapminmax('apply', xvec, rom_model.Xsettings);
y_norm = rom_model.net(x_norm);
coeffs = mapminmax('reverse', y_norm, rom_model.Ysettings);
coeffs = coeffs(:);

y_rom = rom_model.POD_basis * coeffs + rom_model.snapshot_mean;
end