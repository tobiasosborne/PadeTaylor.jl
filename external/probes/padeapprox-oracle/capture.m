function capture()
% capture.m — run padeapprox.m on test cases, write oracles.txt as Julia literals.

addpath('../../chebfun');
fprintf('Octave %s; padeapprox.m loaded.\n', version);

fid = fopen('oracles.txt', 'w');
fprintf(fid, '# Pade oracle outputs from padeapprox.m, Octave %s.\n', version);
fprintf(fid, '# Format: Julia-readable literals.\n\n');

% --- Case 2.1.2 — exp(z), (m, n) = (20, 20) ----------------------------------
c = zeros(41, 1); fact = 1;
for k = 0:40
  c(k+1) = 1.0 / fact; fact = fact * (k+1);
end
write_case(fid, 'test_2_1_2_exp_20_20', c, 20, 20, 1e-14);

% --- Case 2.1.3 — log(1.2 - z), (m, n) = (20, 20) ---------------------------
c = zeros(41, 1);
c(1) = log(1.2);
for k = 1:40
  c(k+1) = -1.0 / (k * (1.2)^k);
end
write_case(fid, 'test_2_1_3_log12_20_20', c, 20, 20, 1e-14);

% --- Case 2.1.4 — tan(z^4), (m, n) = (20, 20) -------------------------------
f_tan_z4 = @(z) tan(z.^4);
[r, a, b, mu, nu, poles, ~] = padeapprox(f_tan_z4, 20, 20, 1e-14);
fprintf(fid, '# Case test_2_1_4_tan_z4_20_20 — tan(z^4), (m,n) = (20,20)\n');
write_complex_vec(fid, 'test_2_1_4_tan_z4_20_20_poles', poles);
fprintf(fid, 'test_2_1_4_tan_z4_20_20_mu = %d\n', mu);
fprintf(fid, 'test_2_1_4_tan_z4_20_20_nu = %d\n\n', nu);

N = 2048; r_radius = 1;
z = r_radius * exp(2i*pi*(0:N-1)'/N);
fc = fft(f_tan_z4(z))/N;
tc = 1e-15*norm(fc);
fc(abs(fc) < tc) = 0;
if (norm(imag(fc), inf) < tc)
  fc = real(fc);
end
fc = fc ./ r_radius.^(0:(N-1))';
fc_first = fc(1:41);
write_complex_vec(fid, 'test_2_1_4_tan_z4_20_20_coefs', fc_first);

% --- Case 2.1.5 — 1 + z^2, (m, n) = (1, 1) — defect-1 ill-posed ------------
c = [1; 0; 1];
write_case(fid, 'test_2_1_5_one_plus_zsq_1_1', c, 1, 1, 1e-14);

% --- Case 2.1.6 — noisy 1/(1-z), tol = 1e-5 — should recover (0, 1) --------
randn('state', 42);
n_coefs = 21;
c = ones(n_coefs, 1) + 1e-6 * randn(n_coefs, 1);
write_case(fid, 'test_2_1_6_noisy_geom_10_10', c, 10, 10, 1e-5);

% --- Case 2.1.1 — closed-form check exp(z) (m,n) = (2,2) ------------------
c = zeros(5, 1); fact = 1;
for k = 0:4
  c(k+1) = 1.0 / fact; fact = fact * (k+1);
end
write_case(fid, 'test_2_1_1_exp_2_2', c, 2, 2, 1e-14);

fclose(fid);
fprintf('Wrote oracles.txt.\n');
end

function write_case(fid, name, c, m, n, tol)
  [r, a, b, mu, nu, ~, ~] = padeapprox(c, m, n, tol);
  fprintf(fid, '# Case %s — (m,n) = (%d,%d), tol = %g\n', name, m, n, tol);
  write_real_vec(fid, [name '_a'], a);
  write_real_vec(fid, [name '_b'], b);
  fprintf(fid, '%s_mu = %d\n', name, mu);
  fprintf(fid, '%s_nu = %d\n\n', name, nu);
end

function write_real_vec(fid, name, v)
  fprintf(fid, '%s = [', name);
  for k = 1:length(v)
    if k > 1, fprintf(fid, ', '); end
    fprintf(fid, '%.18e', v(k));
  end
  fprintf(fid, ']\n');
end

function write_complex_vec(fid, name, v)
  fprintf(fid, '%s = ComplexF64[', name);
  for k = 1:length(v)
    if k > 1, fprintf(fid, ', '); end
    fprintf(fid, '%.18e + %.18eim', real(v(k)), imag(v(k)));
  end
  fprintf(fid, ']\n');
end
