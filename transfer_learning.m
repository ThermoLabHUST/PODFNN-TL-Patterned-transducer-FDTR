%% FDTR transfer learning (fixed-POD fine-tuning)
% This script fine-tunes a pretrained FDTR POD-FNNROM on a new target-domain dataset.
%
% Transfer-learning strategy:
% - keep the original POD basis, snapshot mean, and normalization settings fixed;
% - project the new snapshots onto the old POD space;
% - continue training the pretrained neural network using the projected coefficients.
%
% Main workflow:
% 1) load the pretrained source-domain ROM model;
% 2) load the new target-domain dataset;
% 3) project the new snapshots into the original POD basis;
% 4) optionally mix a fraction of old samples with new samples;
% 5) fine-tune the existing network and save the transferred model.
%
% Main output:
%   Outputs/FDTR_PODNNROM_model_TL_transfer_finetune.mat
%
% Notes:
% - The source and target datasets must use the same frequency grid and varList.
clc; clear; close all;

%% ===== paths =====
outputsDir = fullfile(pwd, 'Outputs');

oldModelFile   = fullfile(outputsDir, 'FDTR_PODNNROM_model.mat');
oldDatasetFile = fullfile(outputsDir, 'FDTR_LHS_dataset.mat');          % Original source dataset A
newDatasetFile = fullfile(outputsDir, 'FDTR_LHS_dataset_TL.mat');   % New target dataset B

newModelFile   = fullfile(outputsDir, 'FDTR_PODNNROM_model_TL_transfer_finetune.mat');

%% ===== transfer-learning options =====
transferOpt.useMixedFinetune = true;   % true: Mixed fine-tuning with old and new; false: Fine-tuning only with new.
transferOpt.oldKeepFraction  = 0.1;    % If useMixedFinetune=true, what percentage of old data will be retained for transfer learning
transferOpt.shuffleOldData   = true;

% transfer learning settings
fineTuneOpt.epochs     = 300;      
fineTuneOpt.trainRatio = 0.8;
fineTuneOpt.valRatio   = 0.1;
fineTuneOpt.testRatio  = 0.1;

fineTuneOpt.keepOldTrainFcn = true;
fineTuneOpt.newTrainFcn     = 'trainbr';   

% Diagnostic tests before and after fine-tuning of old models on new data
doPreFineTuneCheck  = true;
doPostFineTuneCheck = true;

%% ===== load old model =====
if exist(oldModelFile, 'file') ~= 2
    error('model files not found: %s', oldModelFile);
end

T = load(oldModelFile, 'rom_model');
if ~isfield(T, 'rom_model')
    error('The variable does not exist in the file %s.', oldModelFile);
end
old_rom_model = T.rom_model;

needModelFields = {'net','Xsettings','Ysettings','num_modes','varList','freq', ...
    'POD_basis','snapshot_mean','singular_values','useAmp','ampTransform','ampFloor'};
for k = 1:numel(needModelFields)
    if ~isfield(old_rom_model, needModelFields{k})
        error('model is missing the field rom_model.%s', needModelFields{k});
    end
end

fprintf('✅  %s\n', oldModelFile);
fprintf('   num_modes = %d\n', old_rom_model.num_modes);
fprintf('   Nfreq     = %d\n', numel(old_rom_model.freq));
fprintf('   Nvar      = %d\n', numel(old_rom_model.varList));

%% ===== load new dataset =====
B = load_fdtr_dataset(newDatasetFile, 'NEW dataset B');

assert_dataset_model_compatible(B, old_rom_model);

fprintf('\n========== NEW dataset summary ==========\n');
fprintf('New success samples = %d\n', B.N);
fprintf('Nvar                = %d\n', size(B.X,2));
fprintf('Nfreq               = %d\n', numel(B.freq));
fprintf('========================================\n\n');

%% ===== build NEW snapshots using OLD model settings =====
snap_new = build_snapshots(B.phi_all, B.amp_all, ...
    old_rom_model.useAmp, old_rom_model.ampTransform, old_rom_model.ampFloor);

coeff_new = old_rom_model.POD_basis' * (snap_new - old_rom_model.snapshot_mean);

fprintf('✅ The new data has been projected into the old POD space.\n');
fprintf('   snapshot dim = %d\n', size(snap_new,1));
fprintf('   coeff size   = [%d x %d]\n', size(coeff_new,1), size(coeff_new,2));

%% ===== optional: old model check on NEW dataset before fine-tuning =====
if doPreFineTuneCheck
    fprintf('\n========== Before fine-tuning: old model on NEW dataset ==========\n');
    Ypred_old_on_new = zeros(size(snap_new));
    for i = 1:B.N
        Ypred_old_on_new(:,i) = predict_fdtr_rom_from_SI(B.X(i,:), old_rom_model);
    end
    report_prediction_metrics(snap_new, Ypred_old_on_new, ...
        old_rom_model.useAmp, old_rom_model.ampTransform, old_rom_model.freq, ...
        'Old model on NEW dataset');
    fprintf('=================================================================\n');
end

%% ===== optionally load OLD dataset for mixed fine-tuning =====
A = [];
snap_old = [];
coeff_old = [];
X_finetune = B.X;
coeff_finetune = coeff_new;
sourceTag = 2 * ones(B.N,1);    % 1=OLD, 2=NEW

if transferOpt.useMixedFinetune
    A = load_fdtr_dataset(oldDatasetFile, 'OLD dataset A');
    assert_dataset_model_compatible(A, old_rom_model);

    snap_old = build_snapshots(A.phi_all, A.amp_all, ...
        old_rom_model.useAmp, old_rom_model.ampTransform, old_rom_model.ampFloor);

    coeff_old = old_rom_model.POD_basis' * (snap_old - old_rom_model.snapshot_mean);

    idxOldKeep = select_old_samples(A.N, transferOpt.oldKeepFraction, transferOpt.shuffleOldData);

    X_finetune = [A.X(idxOldKeep,:); B.X];
    coeff_finetune = [coeff_old(:,idxOldKeep), coeff_new];
    sourceTag = [ones(numel(idxOldKeep),1); 2*ones(B.N,1)];

    fprintf('\n========== Fine-tuning dataset summary ==========\n');
    fprintf('Mode                  = OLD + NEW mixed fine-tuning\n');
    fprintf('Old selected samples   = %d / %d\n', numel(idxOldKeep), A.N);
    fprintf('New samples            = %d\n', B.N);
    fprintf('Total finetune samples = %d\n', size(X_finetune,1));
    fprintf('=================================================\n\n');
else
    fprintf('\n========== Fine-tuning dataset summary ==========\n');
    fprintf('Mode                  = NEW-only fine-tuning\n');
    fprintf('New samples            = %d\n', B.N);
    fprintf('Total finetune samples = %d\n', size(X_finetune,1));
    fprintf('=================================================\n\n');
end

%% ===== build finetune feature matrix =====
Xfeat_finetune = build_feature_matrix(X_finetune, old_rom_model.varList);   % [Nvar x Ns]

Xn_finetune = mapminmax('apply', Xfeat_finetune, old_rom_model.Xsettings);
Yn_finetune = mapminmax('apply', coeff_finetune, old_rom_model.Ysettings);

fprintf('Xn range = [%.3f, %.3f]\n', min(Xn_finetune(:)), max(Xn_finetune(:)));
fprintf('Yn range = [%.3f, %.3f]\n', min(Yn_finetune(:)), max(Yn_finetune(:)));

fprintf('Xn outside [0,1]: %.2f%%\n', 100*mean(Xn_finetune(:)<0 | Xn_finetune(:)>1));
fprintf('Yn outside [0,1]: %.2f%%\n', 100*mean(Yn_finetune(:)<0 | Yn_finetune(:)>1));


%% ===== continue training from old net =====
net_ft = old_rom_model.net;

if fineTuneOpt.keepOldTrainFcn
    fprintf('old trainFcn = %s\n', net_ft.trainFcn);
else
    net_ft.trainFcn = fineTuneOpt.newTrainFcn;
    fprintf('new trainFcn = %s\n', net_ft.trainFcn);
end

net_ft.trainParam.epochs = fineTuneOpt.epochs;
net_ft.divideParam.trainRatio = fineTuneOpt.trainRatio;
net_ft.divideParam.valRatio   = fineTuneOpt.valRatio;
net_ft.divideParam.testRatio  = fineTuneOpt.testRatio;

fprintf('\n========== Start fine-tuning ==========\n');
fprintf('epochs      = %d\n', net_ft.trainParam.epochs);
fprintf('trainFcn    = %s\n', net_ft.trainFcn);
fprintf('train/val/test = %.2f / %.2f / %.2f\n', ...
    net_ft.divideParam.trainRatio, net_ft.divideParam.valRatio, net_ft.divideParam.testRatio);
fprintf('=======================================\n\n');

net_ft = train(net_ft, Xn_finetune, Yn_finetune);

%% ===== build transferred model =====
rom_model = old_rom_model;  
rom_model.net = net_ft;

rom_model.transfer_mode        = 'standard_finetune_fixed_POD';
rom_model.oldModelFile         = oldModelFile;
rom_model.oldDatasetFile       = oldDatasetFile;
rom_model.newDatasetFile       = newDatasetFile;
rom_model.useMixedFinetune     = transferOpt.useMixedFinetune;
rom_model.oldKeepFraction      = transferOpt.oldKeepFraction;
rom_model.finetune_epochs      = fineTuneOpt.epochs;
rom_model.finetune_trainFcn    = net_ft.trainFcn;
rom_model.num_finetune_samples = size(X_finetune,1);
rom_model.num_new_samples      = B.N;

if transferOpt.useMixedFinetune
    rom_model.num_old_samples_used = sum(sourceTag == 1);
else
    rom_model.num_old_samples_used = 0;
end

save(newModelFile, 'rom_model');
fprintf('✅ The fine-tuned model has been saved to %s\n', newModelFile);

%% ===== post-finetune checks =====
if doPostFineTuneCheck
    fprintf('\n========== After fine-tuning: NEW dataset ==========\n');
    Ypred_new_after = zeros(size(snap_new));
    for i = 1:B.N
        Ypred_new_after(:,i) = predict_fdtr_rom_from_SI(B.X(i,:), rom_model);
    end
    report_prediction_metrics(snap_new, Ypred_new_after, ...
        rom_model.useAmp, rom_model.ampTransform, rom_model.freq, ...
        'Fine-tuned model on NEW dataset');
    fprintf('===================================================\n');

    if transferOpt.useMixedFinetune
        fprintf('\n========== After fine-tuning: OLD dataset ==========\n');
        Ypred_old_after = zeros(size(snap_old));
        for i = 1:A.N
            Ypred_old_after(:,i) = predict_fdtr_rom_from_SI(A.X(i,:), rom_model);
        end
        report_prediction_metrics(snap_old, Ypred_old_after, ...
            rom_model.useAmp, rom_model.ampTransform, rom_model.freq, ...
            'Fine-tuned model on OLD dataset');
        fprintf('===================================================\n');
    end
end

fprintf('\n✅ Standard transfer learning (fixed POD fine-tuning) finished.\n');

%% ========================= helpers =========================

function S = load_fdtr_dataset(datasetFile, label)
    if exist(datasetFile, 'file') ~= 2
        error('Dataset not found: %s', datasetFile);
    end

    T = load(datasetFile, 'D');
    if ~isfield(T, 'D')
        error('The variable D does not exist in the file %s.', datasetFile);
    end
    D = T.D;

    needFields = {'freq','varList','X_SI','okMask','phi_deg_all','amp_all'};
    for k = 1:numel(needFields)
        if ~isfield(D, needFields{k})
            error('Dataset %s is missing field D.%s', datasetFile, needFields{k});
        end
    end

    ok = D.okMask(:);
    if nnz(ok) < 2
        error('%s Too few successful samples（%d）。', label, nnz(ok));
    end

    S = struct();
    S.label    = label;
    S.file     = datasetFile;
    S.freq     = D.freq(:);
    S.varList  = D.varList;
    S.X        = D.X_SI(ok, :);
    S.phi_all  = D.phi_deg_all(:, ok);
    S.amp_all  = D.amp_all(:, ok);
    S.N        = size(S.X, 1);

    if size(S.phi_all,2) ~= S.N || size(S.amp_all,2) ~= S.N
        error('The sample size in %s is mismatched for X / phi / amp.', label);
    end
    if size(S.phi_all,1) ~= numel(S.freq) || size(S.amp_all,1) ~= numel(S.freq)
        error('The freq in %s does not match the first dimension of phi/amp.', label);
    end
end

function assert_dataset_model_compatible(S, rom_model)
    if numel(S.freq) ~= numel(rom_model.freq)
        error('The number of frequency points in %s is inconsistent with that in the old model.', S.label);
    end

    if max(abs(S.freq(:) - rom_model.freq(:))) > 1e-12 * max(1, max(abs(rom_model.freq)))
        error('%s is inconsistent with the old model freq, so fixed POD transfer learning cannot be performed directly.', S.label);
    end

    if numel(S.varList) ~= numel(rom_model.varList)
        error('The length of %s is inconsistent with the old model varList.', S.label);
    end

    for j = 1:numel(S.varList)
        nameS = get_struct_field_or_default(S.varList(j), 'name', '');
        nameM = get_struct_field_or_default(rom_model.varList(j), 'name', '');
        modeS = get_struct_field_or_default(S.varList(j), 'mode', '');
        modeM = get_struct_field_or_default(rom_model.varList(j), 'mode', '');

        if ~strcmpi(char(nameS), char(nameM))
            error('The %d parameter name is inconsistent with the old model:dataset=%s, model=%s', ...
                S.label, j, char(nameS), char(nameM));
        end
        if ~strcmpi(char(modeS), char(modeM))
            error('The %s parameter %d, mode, is inconsistent with the old model:dataset=%s, model=%s', ...
                S.label, j, char(modeS), char(modeM));
        end
    end
end

function val = get_struct_field_or_default(s, fieldName, defaultVal)
    if isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = defaultVal;
    end
end

function snapshots = build_snapshots(phi_all, amp_all, useAmpInSnapshot, ampTransform, ampFloor)
    if useAmpInSnapshot
        if strcmpi(ampTransform, 'log10')
            snapshots = [phi_all; log10(max(amp_all, ampFloor))];
        else
            snapshots = [phi_all; amp_all];
        end
    else
        snapshots = phi_all;
    end
end

function idxKeep = select_old_samples(Nold, keepFraction, doShuffle)
    if keepFraction <= 0
        idxKeep = [];
        return;
    end

    keepFraction = min(max(keepFraction, 0), 1);
    Nkeep = max(1, round(Nold * keepFraction));

    idx = (1:Nold)';
    if doShuffle
        idx = idx(randperm(Nold));
    end

    idxKeep = idx(1:Nkeep);
    idxKeep = sort(idxKeep);
end

function Xfeat = build_feature_matrix(X, varList)
  
    Ns = size(X,1);
    Nvar = size(X,2);
    Xfeat = zeros(Nvar, Ns);

    for j = 1:Nvar
        mode_j = get_struct_field_or_default(varList(j), 'mode', 'linear');
        if strcmpi(mode_j, 'log')
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
        mode_j = get_struct_field_or_default(rom_model.varList(j), 'mode', 'linear');
        if strcmpi(mode_j, 'log')
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

function report_prediction_metrics(Ytrue, Ypred, useAmpInSnapshot, ampTransform, freq, label)
    if isempty(Ytrue) || isempty(Ypred)
        fprintf('%s: Empty data, skip.\n', label);
        return;
    end

    Nf = numel(freq);

    if useAmpInSnapshot
        phi_true = Ytrue(1:Nf,:);
        aux_true = Ytrue(Nf+1:end,:);

        phi_pred = Ypred(1:Nf,:);
        aux_pred = Ypred(Nf+1:end,:);

        rmse_phi = mean(sqrt(mean((phi_pred - phi_true).^2, 1)));
        r2_phi   = mean(1 - sum((phi_pred - phi_true).^2,1) ./ sum((phi_true - mean(phi_true,1)).^2,1));

        rmse_aux = mean(sqrt(mean((aux_pred - aux_true).^2, 1)));

        fprintf('%s:\n', label);
        fprintf('  Phase mean RMSE      = %.4e deg\n', rmse_phi);
        fprintf('  Phase mean R2        = %.4f\n', r2_phi);
        if strcmpi(ampTransform,'log10')
            fprintf('  log10(Amp) mean RMSE = %.4e\n', rmse_aux);
        else
            fprintf('  Amp mean RMSE        = %.4e\n', rmse_aux);
        end
    else
        phi_true = Ytrue;
        phi_pred = Ypred;

        rmse_phi = mean(sqrt(mean((phi_pred - phi_true).^2, 1)));
        r2_phi   = mean(1 - sum((phi_pred - phi_true).^2,1) ./ sum((phi_true - mean(phi_true,1)).^2,1));

        fprintf('%s:\n', label);
        fprintf('  Phase mean RMSE = %.4e deg\n', rmse_phi);
        fprintf('  Phase mean R2   = %.4f\n', r2_phi);
    end
end



function [POD_basis, snapshot_mean, singular_values] = build_pod_basis_from_snapshots(snapshots, num_modes)
    if isempty(snapshots) || size(snapshots,2) < 2
        error('The number of snapshots is too small to establish a POD basis.');
    end

    snapshot_mean = mean(snapshots, 2);
    Xc = snapshots - snapshot_mean;

    [U, S, ~] = svd(Xc, 'econ');
    s = diag(S);

    r = min([num_modes, size(U,2), numel(s)]);
    if r < num_modes
        warning('The number of available POD modalities is only %d, which is less than the requested %d.', r, num_modes);
    end

    POD_basis = U(:,1:r);
    singular_values = s;
end