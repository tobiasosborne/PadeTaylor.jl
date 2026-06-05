"""
capture.py — independent sympy oracle verification for test/corpus_vector_bvp_test.jl
(epic padetaylor-p1v0, bead padetaylor-g65t).  Run: python3 capture.py

The corpus_vector_bvp file pins the GENERAL two-point boundary condition
    B_a·y(z_a) + B_b·y(z_b) = g
of `VectorBVP.vector_bvp_solve` (gap #10 — ADR-0023's entire rationale, the
coupling/periodic BC assembly, untested until now: the committed
vector_bvp_test.jl exercises only the Dirichlet SELECTOR B = [1 0;0 0]).

Every block re-derives the oracle from the CLOSED FORM y = (sin z, cos z),
NOT from the Julia routine under test (CLAUDE.md Rule 5 / Law 1).  We verify,
symbolically, that (i) the ODE residual is 0, (ii) the chosen g equals
B_a·y(z_a)+B_b·y(z_b) for the closed form, and (iii) the BC pins a UNIQUE
solution (cases 1,2) or is genuinely RANK-DEFICIENT/underdetermined (case 3).
"""
import sympy as sp

z = sp.symbols('z', real=True)
A, B = sp.symbols('A B')

# Closed form: y = (sin z, cos z).  ODE y' = (y2, -y1).
y1, y2 = sp.sin(z), sp.cos(z)
yexact = sp.Matrix([y1, y2])
rhs = sp.Matrix([y2, -y1])
ode_resid = sp.simplify(sp.Matrix([sp.diff(y1, z), sp.diff(y2, z)]) - rhs)
print("# ODE residual  y' - (y2,-y1)  (must be [0,0]):", ode_resid.T.tolist())


def yat(zv):
    return sp.Matrix([sp.sin(zv), sp.cos(zv)])


# General solution of y' = (y2,-y1):  y1 = A sin z + B cos z,  y2 = y1'.
gy1 = A * sp.sin(z) + B * sp.cos(z)
gy2 = sp.diff(gy1, z)


def uniqueness(Ba, Bb, g, za, zb):
    gza = sp.Matrix([gy1.subs(z, za), gy2.subs(z, za)])
    gzb = sp.Matrix([gy1.subs(z, zb), gy2.subs(z, zb)])
    return sp.solve(list(Ba * gza + Bb * gzb - g), [A, B], dict=True)


# === CASE 1 — vector-antiperi-sincos-bvp: B_a=I, B_b=-I on [0, pi/2] ===
print("\n# === CVB.1  antiperiodic  B_a=I, B_b=-I  on [0, pi/2] ===")
za1, zb1 = sp.Integer(0), sp.pi / 2
Ba1, Bb1 = sp.eye(2), -sp.eye(2)
g1 = sp.simplify(Ba1 * yat(za1) + Bb1 * yat(zb1))
print("#   g = B_a y(0) + B_b y(pi/2) =", g1.T.tolist(), " (catalogue (-1,1))")
print("#   y(pi/4) =", sp.simplify(yat(sp.pi / 4)).T.tolist(), " (= (sqrt2/2, sqrt2/2))")
print("#   [Ba|Bb] rank =", sp.Matrix.hstack(Ba1, Bb1).rank(), "(full rank 2)")
print("#   uniqueness (A,B):", uniqueness(Ba1, Bb1, g1, za1, zb1), "-> UNIQUE A=1,B=0 = (sin,cos)")

# === CASE 2 — vector-endpoint-mixing-bvp: B_a=[[0,0],[0,1]], B_b=[[1,0],[0,0]] on [0,pi] ===
print("\n# === CVB.2  mixing  B_a=[[0,0],[0,1]], B_b=[[1,0],[0,0]]  on [0, pi] ===")
za2, zb2 = sp.Integer(0), sp.pi
Ba2 = sp.Matrix([[0, 0], [0, 1]])
Bb2 = sp.Matrix([[1, 0], [0, 0]])
g2 = sp.simplify(Ba2 * yat(za2) + Bb2 * yat(zb2))
print("#   g = B_a y(0) + B_b y(pi) =", g2.T.tolist(), " (catalogue [0,1])")
print("#   row1 pins y1(pi)=sin(pi)=", sp.sin(zb2), " ; row2 pins y2(0)=cos(0)=", sp.cos(za2))
print("#   y(pi/2) =", sp.simplify(yat(sp.pi / 2)).T.tolist(), " (= (1,0))")
print("#   [Ba|Bb] rank =", sp.Matrix.hstack(Ba2, Bb2).rank(), "(full rank 2)")
print("#   uniqueness (A,B):", uniqueness(Ba2, Bb2, g2, za2, zb2), "-> UNIQUE A=1,B=0 = (sin,cos)")

# === CASE 3 — vector-rank-deficient-bc: rank-1 combined operator on [pi/6, pi/2] ===
print("\n# === CVB.3  rank-deficient  B_a=[[1,0],[2,0]], B_b=0  on [pi/6, pi/2] ===")
za3, zb3 = sp.pi / 6, sp.pi / 2
Ba3 = sp.Matrix([[1, 0], [2, 0]])
Bb3 = sp.zeros(2, 2)
g3 = sp.simplify(Ba3 * yat(za3) + Bb3 * yat(zb3))
M3 = sp.Matrix.hstack(Ba3, Bb3)
print("#   [Ba|Bb] =", M3.tolist())
print("#   [Ba|Bb] rank =", M3.rank(), "(of 2) -> RANK-DEFICIENT (row2 = 2*row1)")
print("#   g (from closed form) =", g3.T.tolist(), " ; consistent (g2==2*g1):",
      sp.simplify(g3[1] - 2 * g3[0]) == 0, " ; non-trivial (g != 0): True")
print("#   uniqueness (A,B):", uniqueness(Ba3, Bb3, g3, za3, zb3),
      "-> FREE PARAMETER => underdetermined => singular Jacobian (Rule-1 throw)")
