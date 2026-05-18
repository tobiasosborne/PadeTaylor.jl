% motygin_capture.m
% Cross-validation driver: Motygin HeunG / HeunC vs Mathematica oracle.
%
% Regime A parameters:
%   HeunG: a=2, q=1, alpha=1, beta=2, gamma=3, delta=4
%          Mathematica convention: HeunG[a,q,alpha,beta,gamma,delta,z]
%          Motygin convention:     HeunL(a,q,alpha,beta,gamma,delta,z)
%          Both normalise to HeunG(0)=1.  Argument order matches exactly.
%
%   HeunC: q=1, alpha=1, gamma=2, delta=3, epsilon=-1
%          Mathematica convention: HeunC[q,alpha,gamma,delta,epsilon,z]
%          Motygin convention:     HeunC(q,alpha,gamma,delta,epsilon,z)
%          Both normalise to HeunC(0)=1.  Argument order matches exactly.
%
% BRANCH-CUT NOTE:
%   HeunL0 / HeunL returns NaN for real z >= 1 (lies on the branch cut
%   [1, a] and [a, +inf]).  For z in {1.5, 1.9, 2.5, 3.0} we approach
%   from just above the real axis: z + i*eta, eta = 1e-10.
%   Mathematica's HeunG is defined continuously on the principal sheet
%   from below (standard branch convention), so we also request the value
%   from above (eta > 0) to match Motygin's analytic continuation direction.
%   Both oracles share the same limiting value (the function is real-valued
%   approaching from either side inside [0,1]), so for z in (0,1) no offset
%   is needed.
%
% OUTPUT FORMAT:
%   One line per evaluation:
%   Func|...|z=<real>+<imag>i|motygin_val=<re>+<im>i|motygin_err=<err>|wrnmsg=<msg>
%
% Oleg Motygin's code lives at:
%   external/heun-refs/motygin-Heun/
%   external/heun-refs/motygin-confluent-Heun/

% Add Motygin's code directories to the Octave path.
% We use paths relative to the location of this script.
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fullfile(script_dir, '..', '..', '..');
addpath(fullfile(repo_root, 'external', 'heun-refs', 'motygin-Heun'));
addpath(fullfile(repo_root, 'external', 'heun-refs', 'motygin-confluent-Heun'));

% HeunOpts initialises the global tuning parameters.  Call once.
HeunOpts();

% Small positive imaginary offset to approach real branch cut from above.
eta = 1e-10;

% --- Regime A: HeunG ---
a = 2; q = 1; alpha = 1; beta = 2; gamma = 3; delta_p = 4;
z_grid_G = [0.1, 0.5, 0.9, 1.5, 1.9, 2.5, 3.0];

fprintf('# motygin_capture output\n');
fprintf('# Motygin HeunG convention: HeunL(a,q,alpha,beta,gamma,delta,z), normalised HeunG(0)=1\n');
fprintf('# Motygin HeunC convention: HeunC(q,alpha,gamma,delta,epsilon,z), normalised HeunC(0)=1\n');
fprintf('# Branch-cut offset eta=%g applied to real z >= 1 for HeunG\n', eta);
fprintf('# Branch-cut offset eta=%g applied to real z >= 1 for HeunC\n', eta);
fprintf('#\n');
fprintf('# === Regime A: HeunG ===\n');

for k = 1:length(z_grid_G)
    z = z_grid_G(k);
    z_eval = z;
    branch_note = '';
    if (imag(z) == 0) && (real(z) >= 1)
        z_eval = z + 1i * eta;
        branch_note = sprintf('+%gi_offset', eta);
    end

    [val, dval, err, numb, wrnmsg] = HeunL(a, q, alpha, beta, gamma, delta_p, z_eval);

    % Format value
    re_v = real(val);
    im_v = imag(val);
    if im_v >= 0
        val_str = sprintf('%.15e+%.15ei', re_v, im_v);
    else
        val_str = sprintf('%.15e%.15ei', re_v, im_v);
    end

    % Format z for display (always show real part, omit imaginary in label)
    fprintf('HeunG|func=HeunG|a=2|q=1|alpha=1|beta=2|gamma=3|delta=4|z=%.1f%s|motygin_val=%s|motygin_err=%.3e|nterms=%d|wrnmsg=%s\n', ...
        z, branch_note, val_str, err, numb, wrnmsg);
end

% --- Regime A: HeunC ---
qC = 1; alphaC = 1; gammaC = 2; deltaC = 3; epsilonC = -1;
z_grid_C = [0.1, 0.5, 0.9, 1.5, 1.9, 2.5, 3.0];

fprintf('#\n');
fprintf('# === Regime A: HeunC ===\n');

for k = 1:length(z_grid_C)
    z = z_grid_C(k);
    z_eval = z;
    branch_note = '';
    % HeunC branch cut is [1, +inf) on the real axis.
    if (imag(z) == 0) && (real(z) >= 1)
        z_eval = z + 1i * eta;
        branch_note = sprintf('+%gi_offset', eta);
    end

    [val, dval, err, numb, wrnmsg] = HeunC(qC, alphaC, gammaC, deltaC, epsilonC, z_eval);

    re_v = real(val);
    im_v = imag(val);
    if im_v >= 0
        val_str = sprintf('%.15e+%.15ei', re_v, im_v);
    else
        val_str = sprintf('%.15e%.15ei', re_v, im_v);
    end

    fprintf('HeunC|func=HeunC|q=1|alpha=1|gamma=2|delta=3|epsilon=-1|z=%.1f%s|motygin_val=%s|motygin_err=%.3e|nterms=%d|wrnmsg=%s\n', ...
        z, branch_note, val_str, err, numb, wrnmsg);
end

fprintf('#\n');
fprintf('# === capture complete ===\n');
