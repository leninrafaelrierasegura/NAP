function results = audit_matlab_conversion()
% Reproduce only the supplied polynomial/residue conversion, without editing
% the originals, executing their Mac-only wrapper, or generating 999 fits.
% First run audit_original_generators.py, then from NAP:
% matlab -batch "addpath('test/brasil'); audit_matlab_conversion"
folder = fullfile(fileparts(mfilename('fullpath')), 'generator_audit');
cases = dir(fullfile(folder, 'm*_beta*'));
rows = [];
for j = 1:numel(cases)
    input = load(fullfile(folder, cases(j).name, 'conversion_inputs.mat'));
    c = fliplr(poly(input.zeros));
    b = fliplr(poly(input.poles));
    [r,p,k] = residue(c,b);
    y = input.y(:);
    reconstructed = sum(r(:).' ./ (y-p(:).'), 2) + k;
    original = input.actual(:);
    % Include the missing gain only for this diagnostic. Production generation
    % should bypass monomial expansion and convert MP poles/residues directly.
    corrected = input.gain * reconstructed;
    root_only_at_one = sum(r ./ (1-p)) + k;
    rows = [rows; double(input.m), input.beta, real(root_only_at_one), ...
        max(abs(imag(p))), max(abs(imag(r))), min(real(r)), real(k), ...
        max(abs(corrected-original))]; %#ok<AGROW>
    save(fullfile(folder, cases(j).name, 'matlab_conversion.mat'), 'r','p','k');
end
results = array2table(rows, 'VariableNames', {'m','beta','root_only_at_one', ...
    'max_imaginary_poles','max_imaginary_residues','min_real_residue', ...
    'constant','gain_corrected_vs_barycentric_max_gap'});
writetable(results, fullfile(folder,'matlab_conversion.csv'));
disp(results)
end
