function capture()
% capture.m — BVP oracle using DMSUITE Chebyshev spectral method.
%
% Produces oracles.txt with Julia-readable literals for:
%   Group 1: D₁ and D₂ matrix primitives from chebdif.m
%   Group 2: Linear BVP  u'' = u, solution cosh(z)/cosh(1)
%   Group 3: Nonlinear PI BVP on real segment [-18,-14]
%   Group 4: Nonlinear PI BVP on real segment [0.1,0.5] (smooth)
%   Group 5: Barycentric interpolation via chebint.m
%
% Citation: references/bvp_recipe.md §1–5
%           references/markdown/FW2011_painleve_methodology_JCP230/
%               FW2011_painleve_methodology_JCP230.md:186 (Eq 3.2)
%           external/DMSUITE/chebdif.m:38,45-57
%           external/DMSUITE/chebint.m:27-36
%
% Usage:
%   cd external/probes/bvp-oracle
%   octave --quiet --eval "capture()"

addpath('../../DMSUITE');
fprintf('Octave %s; DMSUITE loaded.\n', version);

fid = fopen('oracles.txt', 'w');
fprintf(fid, '# BVP oracle outputs from DMSUITE chebdif/chebint, Octave %s.\n', version);
fprintf(fid, '# Format: Julia-readable literals.\n');
fprintf(fid, '# Citation: references/bvp_recipe.md §1-5\n\n');

% =========================================================================
% GROUP 1 — D₂ matrix sanity (deterministic, no Newton; pure DMSUITE)
% Citation: references/bvp_recipe.md §2-3;
%           external/DMSUITE/chebdif.m:38,45-57
% =========================================================================

% --- test_d2_N4: chebdif(5, 2), N=4, 5 nodes ---
fprintf('Computing Group 1 cases...\n');
[x4, DM4] = chebdif(5, 2);
D1_4 = DM4(:,:,1);
D2_4 = DM4(:,:,2);

fprintf(fid, '# Case test_d2_N4 — chebdif(5,2), N=4, 5 nodes\n');
fprintf(fid, '# Citation: references/bvp_recipe.md §2 (D₁ formula), §3 (D₂=D₁² option a)\n');
fprintf(fid, '#           external/DMSUITE/chebdif.m:38,45-57\n');
write_real_vec(fid, 'test_d2_N4_nodes', x4);
write_real_matrix(fid, 'test_d2_N4_D1', D1_4);
write_real_matrix(fid, 'test_d2_N4_D2', D2_4);
fprintf(fid, '\n');

% --- test_d2_N8: chebdif(9, 2), N=8, 9 nodes ---
[x8, DM8] = chebdif(9, 2);
D1_8 = DM8(:,:,1);
D2_8 = DM8(:,:,2);

fprintf(fid, '# Case test_d2_N8 — chebdif(9,2), N=8, 9 nodes\n');
fprintf(fid, '# Citation: references/bvp_recipe.md §2-3\n');
write_real_vec(fid, 'test_d2_N8_nodes', x8);
write_real_matrix(fid, 'test_d2_N8_D1', D1_8);
write_real_matrix(fid, 'test_d2_N8_D2', D2_8);
fprintf(fid, '\n');

% --- test_d2_N16: chebdif(17, 2), N=16, 17 nodes; partial D2 only ---
[x16, DM16] = chebdif(17, 2);
D2_16 = DM16(:,:,2);

fprintf(fid, '# Case test_d2_N16 — chebdif(17,2), N=16, 17 nodes; partial D2\n');
fprintf(fid, '# Citation: references/bvp_recipe.md §3\n');
write_real_vec(fid, 'test_d2_N16_nodes', x16);
% Emit diagonal and first superdiagonal only (regression check; 17×17 is verbose)
write_real_vec(fid, 'test_d2_N16_D2_diag', diag(D2_16));
write_real_vec(fid, 'test_d2_N16_D2_supdiag', diag(D2_16, 1));
fprintf(fid, '\n');

% =========================================================================
% GROUP 2 — Linear BVP: u'' = u on [-1,1]; u(z) = cosh(z)/cosh(1)
% BCs: u(-1) = 1.0, u(1) = 1.0 (both equal cosh(1)/cosh(1) = 1)
% Exact: u(t) = cosh(t)/cosh(1) for t in [-1,1]
% Citation: references/bvp_recipe.md §4 (Newton iteration, linear case)
% =========================================================================

fprintf('Computing Group 2: Linear BVP cases...\n');

% Linear BVP; z_a = -1, z_b = 1; F(z,u) = u, dF/du = 1
% Scale factor: (1/4)(z_b - z_a)^2 = (1/4)*4 = 1
% Residual:  R = (D2 u)_int - u_int
% Jacobian:  J = D2_int_int - I_int
% BCs: u(t=-1) = u_a = cosh(-1)/cosh(1) = 1, u(t=1) = u_b = 1

% --- test_bvp_lin_N8: N=8, 9-node grid ---
N = 8;
tol = 1e-14;
[t8, DM8_loc] = chebdif(N+1, 2);
D2loc = DM8_loc(:,:,2);

u_bc_a = cosh(-1) / cosh(1);  % t=-1 endpoint
u_bc_b = cosh(1) / cosh(1);   % t=+1 endpoint

u8 = cosh(t8) / cosh(1);  % Exact initial guess (warm start; converges in 0 Newton steps)
% Use a generic initial guess (zeros + BCs) to exercise the Newton path
u8_init = zeros(N+1, 1);
u8_init(1) = u_bc_b;    % t=+1 is index 1 in DMSUITE ordering
u8_init(N+1) = u_bc_a;  % t=-1 is index N+1

int_idx = 2:N;  % Interior node indices (1-based; boundary at 1 and N+1)
D2int = D2loc(int_idx, int_idx);

% h = z_b - z_a = 2; scale = (1/4)*h^2 = 1
scale = 1.0;

u_iter = u8_init;
n_iter_lin8 = 0;
for iter = 1:20
    u_int = u_iter(int_idx);
    R = D2loc(int_idx, :) * u_iter - scale * u_int;  % residual: D2*u - u (linear)
    if norm(R, inf) < tol
        break;  % check BEFORE update so final u is the converged state
    end
    % J = D2_int_int - scale*I
    J = D2int - scale * eye(N-1);
    du = J \ R;
    u_iter(int_idx) = u_int - du;
    n_iter_lin8 = n_iter_lin8 + 1;
end
u_lin8 = u_iter;  % save for Group 5 barycentric eval

u_ref_8 = cosh(t8) / cosh(1);
err_8 = norm(u_iter - u_ref_8, inf);
fprintf('  Lin N=8: %d Newton iters, error = %.3e\n', n_iter_lin8, err_8);

fprintf(fid, '# Case test_bvp_lin_N8 — u''''=u on [-1,1], N=8 (9 nodes)\n');
fprintf(fid, '# Citation: references/bvp_recipe.md §4 (Newton), exact: cosh(t)/cosh(1)\n');
fprintf(fid, 'test_bvp_lin_N8_newton_iters = %d\n', n_iter_lin8);
fprintf(fid, 'test_bvp_lin_N8_inf_err_vs_ref = %.18e\n', err_8);
write_real_vec(fid, 'test_bvp_lin_N8_nodes', t8);
write_real_vec(fid, 'test_bvp_lin_N8_u', u_iter);
write_real_vec(fid, 'test_bvp_lin_N8_u_ref', u_ref_8);
fprintf(fid, '\n');

% --- test_bvp_lin_N16: N=16, 17-node grid ---
N = 16;
% Residual floor for N=16 linear BVP: ~5e-13 (cond(D2_16)~10^3 × eps).
% Solution error vs exact is ~1e-15 (near machine eps: spectral accuracy is excellent).
% Use tol = 1e-12 to exit cleanly after the 1-step linear convergence.
tol = 1e-12;
[t16, DM16_loc] = chebdif(N+1, 2);
D2loc = DM16_loc(:,:,2);

u16_init = zeros(N+1, 1);
u16_init(1) = u_bc_b;
u16_init(N+1) = u_bc_a;

int_idx = 2:N;
D2int = D2loc(int_idx, int_idx);
scale = 1.0;

u_iter = u16_init;  % reuse u_iter for N=16 solve
n_iter_lin16 = 0;
for iter = 1:20
    u_int = u_iter(int_idx);
    R = D2loc(int_idx, :) * u_iter - scale * u_int;
    if norm(R, inf) < tol
        break;  % check before update
    end
    J = D2int - scale * eye(N-1);
    du = J \ R;
    u_iter(int_idx) = u_int - du;
    n_iter_lin16 = n_iter_lin16 + 1;
end

t16_loc = t16;
u_ref_16 = cosh(t16_loc) / cosh(1);
err_16 = norm(u_iter - u_ref_16, inf);
fprintf('  Lin N=16: %d Newton iters, error = %.3e\n', n_iter_lin16, err_16);

fprintf(fid, '# Case test_bvp_lin_N16 — u''''=u on [-1,1], N=16 (17 nodes)\n');
fprintf(fid, '# Citation: references/bvp_recipe.md §4, exact: cosh(t)/cosh(1)\n');
fprintf(fid, 'test_bvp_lin_N16_newton_iters = %d\n', n_iter_lin16);
fprintf(fid, 'test_bvp_lin_N16_inf_err_vs_ref = %.18e\n', err_16);
write_real_vec(fid, 'test_bvp_lin_N16_nodes', t16_loc);
write_real_vec(fid, 'test_bvp_lin_N16_u', u_iter);
write_real_vec(fid, 'test_bvp_lin_N16_u_ref', u_ref_16);
fprintf(fid, '\n');

% =========================================================================
% GROUP 3 — Nonlinear PI BVP on real segment [-18, -14]
% ODE: u'' = 6u^2 + z
% Initial guess: u_init = sqrt(-z/6) (FW 2011 asymptotic, valid at large |z|)
% Boundary conditions: from real Painlevé I on the negative real axis
%   u(-18) ≈ sqrt(18/6) = sqrt(3);  u(-14) ≈ sqrt(14/6)
% These are approximate BCs from the asymptotic. We use exact asymptotic.
% Citation: references/bvp_recipe.md §4;
%           references/markdown/FW2011_painleve_methodology_JCP230/
%               FW2011_painleve_methodology_JCP230.md:186 (Eq 3.2)
% =========================================================================

fprintf('Computing Group 3: Nonlinear PI real segment [-18,-14]...\n');

N = 20;
z_a = -18.0;
z_b = -14.0;
h = z_b - z_a;  % = 4.0

[t20, DM20] = chebdif(N+1, 2);
D1_20 = DM20(:,:,1);
D2_20 = DM20(:,:,2);

% Affine map: z(t) = (z_b+z_a)/2 + (z_b-z_a)/2 * t
% Citation: references/bvp_recipe.md §1 (affine map)
z_j = (z_b + z_a)/2 + (z_b - z_a)/2 * t20;

% BCs from asymptotic: u(z) ≈ sqrt(-z/6) for large |z| on real negative axis
u_bc_a_pi = sqrt(abs(z_a) / 6.0);  % z_a = -18, t = -1 -> index N+1
u_bc_b_pi = sqrt(abs(z_b) / 6.0);  % z_b = -14, t = +1 -> index 1

% Initial guess from asymptotic
u_init = sqrt(abs(z_j) / 6.0);
u_init(1) = u_bc_b_pi;
u_init(N+1) = u_bc_a_pi;

int_idx = 2:N;
scale_sq = (h/2)^2;  % (1/4)*(z_b-z_a)^2; Citation: bvp_recipe.md §4 Eq 3.2

D2_20int = D2_20(int_idx, int_idx);
% Tolerance: the Float64 floor for this problem is ~8e-13 (condition(J)*eps~1e-12).
% Set tol = 1e-12 to exit when Newton is in the noise floor (iter ~4).
% bvp_recipe.md §4: "100·eps(real(T)) is a reasonable default"; here 100*eps ≈ 2e-14
% but the floating-point residual floor for N=20, h=4 is ~1e-12 due to cond(J).
tol_nl = 1e-12;

u_iter = u_init;
n_iter_pi_real = 0;
res_final = inf;
for iter = 1:20
    u_int = u_iter(int_idx);
    z_int = z_j(int_idx);
    % Residual: (D2 u)_int - scale_sq * (6 u_int^2 + z_int)
    % Citation: references/bvp_recipe.md §4 Eq. R formula
    lhs = D2_20(int_idx, :) * u_iter;
    rhs = scale_sq * (6 * u_int.^2 + z_int);
    R = lhs - rhs;
    res_final = norm(R, inf);
    if res_final < tol_nl
        break;  % check BEFORE update; final u is already converged state
    end
    % Jacobian: D2_int_int - scale_sq * diag(12 * u_int)
    % Citation: references/bvp_recipe.md §4 J_PI formula
    J = D2_20int - scale_sq * diag(12 * u_int);
    du = J \ R;
    u_iter(int_idx) = u_int - du;
    n_iter_pi_real = n_iter_pi_real + 1;
end
u_pi20 = u_iter;  % save for Group 5 barycentric eval
fprintf('  PI real N=20: %d Newton iters, residual = %.3e\n', n_iter_pi_real, res_final);

% Derivative at endpoints via D1
Du_20 = D1_20 * u_iter;  % full D1 * u
Du1_endpoints = [Du_20(1); Du_20(N+1)];  % t=+1 and t=-1 endpoints
% Scale back to z-derivative: dz/dt = (z_b - z_a)/2 = h/2
% dU/dz = (dU/dt) * (dt/dz) = (dU/dt) / (h/2)
% Citation: references/bvp_recipe.md §1 (affine map jacobian)
Du1_z_endpoints = Du1_endpoints / (h/2);

fprintf(fid, '# Case test_bvp_pi_real_N20 — PI BVP on [-18,-14], N=20 (21 nodes)\n');
fprintf(fid, '# ODE: u'''' = 6u^2 + z; IC: sqrt(-z/6) asymptotic\n');
fprintf(fid, '# Citation: references/bvp_recipe.md §4;\n');
fprintf(fid, '#           FW2011_painleve_methodology_JCP230.md:186 (Eq 3.2)\n');
fprintf(fid, 'test_bvp_pi_real_N20_z_a = %.18e\n', z_a);
fprintf(fid, 'test_bvp_pi_real_N20_z_b = %.18e\n', z_b);
fprintf(fid, 'test_bvp_pi_real_N20_bc_a = %.18e\n', u_bc_a_pi);
fprintf(fid, 'test_bvp_pi_real_N20_bc_b = %.18e\n', u_bc_b_pi);
fprintf(fid, 'test_bvp_pi_real_N20_newton_iters = %d\n', n_iter_pi_real);
fprintf(fid, 'test_bvp_pi_real_N20_residual_inf = %.18e\n', res_final);
write_real_vec(fid, 'test_bvp_pi_real_N20_nodes_t', t20);
write_real_vec(fid, 'test_bvp_pi_real_N20_nodes_z', z_j);
write_real_vec(fid, 'test_bvp_pi_real_N20_u', u_iter);
% du/dz at endpoints (IVP-matching diagnostic; FW 2011 line 192)
fprintf(fid, '# du/dz at t=+1 (z=z_b=-14) and t=-1 (z=z_a=-18)\n');
fprintf(fid, 'test_bvp_pi_real_N20_dup_at_zb = %.18e\n', Du1_z_endpoints(1));
fprintf(fid, 'test_bvp_pi_real_N20_dup_at_za = %.18e\n', Du1_z_endpoints(2));
fprintf(fid, '\n');

% =========================================================================
% GROUP 4 — Nonlinear PI BVP on real segment [0.1, 0.5]
% ODE: u'' = 6u^2 + z (Painleve I, NOT Weierstrass P)
% BCs from FW 2011 ICs via mpmath.odefun at 40 dps:
%   u(0.1) = 1.281737179308048091e+00  (PI; integrates u''=6u^2+z)
%   u(0.5) = 4.034017771895840099e+00  (PI; NOT the Weierstrass value 4.0044...)
% The problems-oracle (oracles_wolfram.txt u_at_0_5 = 4.0044...) is for
%   u'' = 6u^2 (equianharmonic ℘); that is a DIFFERENT ODE.
% Citation: references/bvp_recipe.md §4;
%           FW2011_painleve_methodology_JCP230.md:186 (u''=6u^2+z, Eq. 3.2)
% =========================================================================

fprintf('Computing Group 4: Nonlinear PI smooth segment [0.1,0.5]...\n');

N = 24;
z_a4 = 0.1;
z_b4 = 0.5;
h4 = z_b4 - z_a4;  % = 0.4

[t24, DM24] = chebdif(N+1, 2);
D1_24 = DM24(:,:,1);
D2_24 = DM24(:,:,2);

z_j4 = (z_b4 + z_a4)/2 + (z_b4 - z_a4)/2 * t24;

% BCs from mpmath.odefun integrating u''=6u^2+z (PI) at 40 dps,
% from FW 2011 ICs u(0)=1.071822516416917, u'(0)=1.710337353176786.
% NOTE: these are PI ICs, NOT Weierstrass P values.
% The problems-oracle values (u_at_0_5=4.00446...) are for u''=6u^2, a DIFFERENT equation.
% Citation: FW2011_painleve_methodology_JCP230.md:186 (u''=6u^2+z, Eq. 3.2)
u_bc_a4 = 1.281737179308048091e+00;  % u(0.1) from PI mpmath.odefun at 40 dps
u_bc_b4 = 4.034017771895840099e+00;  % u(0.5) from PI mpmath.odefun at 40 dps

% Initial guess: linear interpolation + small asymptotic adjustment
u_init4 = u_bc_a4 + (u_bc_b4 - u_bc_a4) * (t24 + 1) / 2;  % linear interp in t
u_init4(1) = u_bc_b4;    % t=+1 -> z=z_b
u_init4(N+1) = u_bc_a4;  % t=-1 -> z=z_a

int_idx4 = 2:N;
scale_sq4 = (h4/2)^2;  % (1/4)*(z_b-z_a)^2
D2_24int = D2_24(int_idx4, int_idx4);
% Float64 floor for N=24, h=0.4: ~1.7e-12 (due to cond(D2_24)~10^4 × eps).
% Use tol = 3e-12 to exit cleanly at the noise floor after ~6 iters.
tol_nl4 = 3e-12;

u_iter4 = u_init4;
n_iter_pi_cmplx = 0;
res_final4 = inf;
for iter = 1:20
    u_int4 = u_iter4(int_idx4);
    z_int4 = z_j4(int_idx4);
    lhs4 = D2_24(int_idx4, :) * u_iter4;
    rhs4 = scale_sq4 * (6 * u_int4.^2 + z_int4);
    R4 = lhs4 - rhs4;
    res_final4 = norm(R4, inf);
    if res_final4 < tol_nl4
        break;  % check BEFORE update
    end
    J4 = D2_24int - scale_sq4 * diag(12 * u_int4);
    du4 = J4 \ R4;
    u_iter4(int_idx4) = u_int4 - du4;
    n_iter_pi_cmplx = n_iter_pi_cmplx + 1;
end
fprintf('  PI smooth N=24: %d Newton iters, residual = %.3e\n', n_iter_pi_cmplx, res_final4);

% du/dz at endpoints
Du_24 = D1_24 * u_iter4;
Du1_24_z = [Du_24(1); Du_24(N+1)] / (h4/2);

fprintf(fid, '# Case test_bvp_pi_complex_N24 — PI BVP on [0.1,0.5], N=24 (25 nodes)\n');
fprintf(fid, '# ODE: u'''' = 6u^2 + z (Painleve I); BCs from FW ICs via mpmath.odefun 40dps:\n');
fprintf(fid, '#   u(0.1)=1.281737179308048091e+00 (PI mpmath.odefun from FW ICs),\n');
fprintf(fid, '#   u(0.5)=4.034017771895840099e+00 (PI mpmath.odefun from FW ICs)\n');
fprintf(fid, '# NOTE: These are DIFFERENT from the problems-oracle Weierstrass P values\n');
fprintf(fid, '#   (u_at_0_5=4.0044...) which are for u''''=6u^2, not Painleve I.\n');
fprintf(fid, '# Citation: references/bvp_recipe.md §4;\n');
fprintf(fid, '#           FW2011_painleve_methodology_JCP230.md:186 (u''''=6u^2+z)\n');
fprintf(fid, 'test_bvp_pi_complex_N24_z_a = %.18e\n', z_a4);
fprintf(fid, 'test_bvp_pi_complex_N24_z_b = %.18e\n', z_b4);
fprintf(fid, 'test_bvp_pi_complex_N24_bc_a = %.18e\n', u_bc_a4);
fprintf(fid, 'test_bvp_pi_complex_N24_bc_b = %.18e\n', u_bc_b4);
fprintf(fid, 'test_bvp_pi_complex_N24_newton_iters = %d\n', n_iter_pi_cmplx);
fprintf(fid, 'test_bvp_pi_complex_N24_residual_inf = %.18e\n', res_final4);
write_real_vec(fid, 'test_bvp_pi_complex_N24_nodes_t', t24);
write_real_vec(fid, 'test_bvp_pi_complex_N24_nodes_z', z_j4);
write_real_vec(fid, 'test_bvp_pi_complex_N24_u', u_iter4);
fprintf(fid, 'test_bvp_pi_complex_N24_dup_at_zb = %.18e\n', Du1_24_z(1));
fprintf(fid, 'test_bvp_pi_complex_N24_dup_at_za = %.18e\n', Du1_24_z(2));
fprintf(fid, '\n');

% =========================================================================
% GROUP 5 — Barycentric interpolation via chebint
% Citation: references/bvp_recipe.md §5 (barycentric formula);
%           external/DMSUITE/chebint.m:27-36
%           BerrutTrefethen2004_barycentric_SIAMReview.md:127-133, 164-171
% =========================================================================

fprintf('Computing Group 5: Barycentric interpolation...\n');

% --- test_baryeval_lin_N8 ---
% Evaluate the N=8 linear-BVP solution at t* = [-0.7, -0.3, 0.0, 0.4, 0.8]
t_star_lin = [-0.7; -0.3; 0.0; 0.4; 0.8];
p_lin = chebint(u_lin8, t_star_lin);  % u_lin8: N=8 solution from Group 2
u_ref_lin_star = cosh(t_star_lin) / cosh(1);

fprintf(fid, '# Case test_baryeval_lin_N8 — barycentric eval of N=8 linear BVP solution\n');
fprintf(fid, '# t* in {-0.7,-0.3,0.0,0.4,0.8}; reference: cosh(t*)/cosh(1)\n');
fprintf(fid, '# Citation: references/bvp_recipe.md §5 (barycentric formula)\n');
fprintf(fid, '#           external/DMSUITE/chebint.m:27-36\n');
write_real_vec(fid, 'test_baryeval_lin_N8_t_star', t_star_lin);
write_real_vec(fid, 'test_baryeval_lin_N8_p', p_lin);
write_real_vec(fid, 'test_baryeval_lin_N8_ref', u_ref_lin_star);
fprintf(fid, '\n');

% --- test_baryeval_pi_real_N20 ---
% From Group 3 (PI real N=20) solution, evaluate at t* in {-0.6,-0.2,0.2,0.6}
% and also emit the corresponding z* values after affine map.
t_star_pi = [-0.6; -0.2; 0.2; 0.6];
p_pi = chebint(u_pi20, t_star_pi);  % u_pi20: N=20 PI solution from Group 3
z_star_pi = (z_b + z_a)/2 + (z_b - z_a)/2 * t_star_pi;  % z* on [-18,-14]

fprintf(fid, '# Case test_baryeval_pi_real_N20 — barycentric eval of PI real N=20 solution\n');
fprintf(fid, '# t* in {-0.6,-0.2,0.2,0.6}; segment [-18,-14]\n');
fprintf(fid, '# Citation: references/bvp_recipe.md §5;\n');
fprintf(fid, '#           FW2011_painleve_methodology_JCP230.md:186-190 (barycentric, §3.2)\n');
write_real_vec(fid, 'test_baryeval_pi_real_N20_t_star', t_star_pi);
write_real_vec(fid, 'test_baryeval_pi_real_N20_z_star', z_star_pi);
write_real_vec(fid, 'test_baryeval_pi_real_N20_p', p_pi);
fprintf(fid, '\n');

fclose(fid);
fprintf('Wrote oracles.txt.\n');
end

% =========================================================================
% Helper functions
% =========================================================================

function write_real_vec(fid, name, v)
  v = v(:);
  fprintf(fid, '%s = [', name);
  for k = 1:length(v)
    if k > 1, fprintf(fid, ', '); end
    fprintf(fid, '%.18e', v(k));
  end
  fprintf(fid, ']\n');
end

function write_real_matrix(fid, name, M)
  [nr, nc] = size(M);
  fprintf(fid, '%s = [', name);
  for i = 1:nr
    if i > 1, fprintf(fid, '; '); end
    for j = 1:nc
      if j > 1, fprintf(fid, ' '); end
      fprintf(fid, '%.18e', M(i,j));
    end
  end
  fprintf(fid, ']\n');
end
