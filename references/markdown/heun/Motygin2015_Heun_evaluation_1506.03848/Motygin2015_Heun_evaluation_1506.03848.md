# On evaluation of the Heun functions

Oleg V. Motygin

Institute of Problems in Mechanical Engineering, Russian Academy of Sciences, V.O., Bol'shoy pr., 61, 199178 St. Petersburg, Russia email: o.v.motygin@gmail.com

#### Abstract

In the paper we deal with the Heun functions — solutions of the Heun equation, which is the most general Fuchsian equation of second order with four regular singular points. Despite the increasing interest to the equation and numerous applications of the functions in a wide variety of physical problems, it is only Maple amidst known software packages which is able to evaluate the Heun functions numerically. But the Maple routine is known to be imperfect: even at regular points it may return infinities or end up with no result. Improving the situation is difficult because the code is not publicly available. The purpose of the work is to suggest and develop alternative algorithms for numerical evaluation of the Heun functions. A procedure based on power series expansions and analytic continuation is suggested which allows us to avoid numerical integration of the differential equation and to achieve reasonable efficiency and accuracy. Results of numerical tests are given.

#### 1 Introduction

In the present paper we deal with the Heun functions which are solutions of the equation introduced by Karl Heun in 1889 [3] as a generalization of the hypergeometric equation. Heun's equation is the most general Fuchsian equation of second order with four regular singular points; we refer to [12, 14, 9, 15] for a comprehensive mathematical treatment of the topic. At the same time, it is only in recent years when the equation has become popular in the physics literature. Now the Heun equation appears in many fields of modern physics and it is used to describe a wide variety of physical phenomena — a comprehensive literature on physical applications can be found in [4]. As a good source of information on the current development of the field one should also mention "The Heun project": http://theheunproject.org/.

Despite the increasing interest to the equation, at the moment the only, to author's knowledge, software package able to evaluate the Heun functions numerically is Maple. However, the code is far from being perfect, it is not difficult to encounter problems, when for ordinary parameters, at regular points the routine returns infinities or spends tens seconds with no result. The quality of the code is a serious obstacle for (potentially very numerous) applications of the functions. Improving the situation is virtually impossible because the code is not publicly available.

The purpose of the present work is to develop alternative algorithms for numerical evaluation of the Heun functions. We suggest a procedure based on power series expansions and analytic continuation which allows us to avoid numerical integration of the differential equation and to achieve reasonable efficiency and accuracy. Program code is presented in [11]. Results of numerical tests and comparison with a case when Heun's function reduces to a simple algebraic function are given.

The algorithm is applicable for computation of the multi-valued Heun functions. In the present paper we also define their single-valued counterparts by fixation of branch cuts. (Notably, in most studies of the Heun functions the subject is not discussed; the author has also been unable to find information on branch cuts in Maple's documentation on the Heun function.) For the single-valued functions an improvement of the algorithm for points close to the singular ones is described.

It should be noted that the developed algorithms are not intended to be universal — more or less ordinary parameters are assumed. Surely, numerical problems are expected and special treatment is needed for the cases of merging singular points (see [7]) or large accessory parameter.

### 2 Statement and basic notations

We start by writing Heun's equation in the standard form (see [2])

$$H''(z) + \left(\frac{\gamma}{z} + \frac{\delta}{z - 1} + \frac{\varepsilon}{z - a}\right)H'(z) + \frac{\alpha\beta z - q}{z(z - 1)(z - a)}H(z) = 0 \tag{1}$$

with the Riemann P-symbol:

$$P \left\{ \begin{array}{cccc} 0 & 1 & a & \infty \\ 0 & 0 & 0 & \alpha \,; & z \\ 1-\gamma & 1-\delta & 1-\varepsilon & \beta \end{array} \right\}.$$

The parameter q ∈ C is usually referred to as an accessory or auxiliary parameter and α, β, γ, δ, ε (also belonging to C) are exponent-related parameters connected via the Fuchsian relation

$$\alpha + \beta + 1 = \gamma + \delta + \varepsilon.$$

The equation has four regular singular points located at z = 0, 1, a, ∞. It will be assumed below that a ∈ C, a 6= {0, 1, ∞}. In the notation of the solutions sometimes we will omit some of the parameters so that H(z) = H(a, . . . ; z) = H(a, q, α, β, γ, δ; z).

The Frobenius method can be used to derive local power-series solutions to (1). There are 8 local solutions of equation (1) (two per a singular point). In § 3 we will present the two local power-series solutions in a neighbourhood of the point z = 0. One of them is analytic in a vicinity of zero and if γ is not a nonpositive integer, this solution, normalized to unity, is called the local Heun function (see [12]). It is usually denoted by Hl(a, q, α, β, γ, δ; z). For the second Frobenius local solution we will use the notation Hs(a, q, α, β, γ, δ; z).

When γ is a nonpositive integer, one solution of (1) is analytic but equal to zero at z = 0, whereas the second solution can be normalized to unity but generally is not analytic. It can be an arguable point, but we will denote by Hl(a, q, α, β, γ, δ; z) the normalized solution and by Hs(a, q, α, β, γ, δ; z) the analytic one.

It is known that generally Hl(a, . . . ; z) is a multi-valued function with branch points at z = 1, a and ∞ and, so, to define a single-valued function we should choose branch cuts. In the present work we fix the branch cuts B1<sup>∞</sup> = (1, +∞) and Ba<sup>∞</sup> = (a, e i arg(a)∞) connecting the points 1 and a to ∞, respectively (see fig. 1). The second function Hs(a, . . . ; z) has the same branch cuts B1<sup>∞</sup> and Ba∞. Besides, the function Hs(z) — generally, and the function Hl(z) — for γ ∈ {0} ∪ Z<sup>−</sup> (Z<sup>−</sup> means the set of negative integers), have a branch point at z = 0. So, we will also define the branch cut B0<sup>∞</sup> = (−∞, 0).

It should be noted that the choice of branch cuts in Maple is unclear, but, definitely, differs from that in the present paper. With the present fixation the function definition domain is star-like with respect to zero which is natural as we will define functions Hl(z), Hs(z) by analytic continuation from a vicinity of z = 0 (see § 5). Besides, there is good agreement between branch cuts of the single-valued functions Hl, Hs and their representations near singular points (see § 6).

![](_page_1_Figure_12.jpeg)

Figure 1: Branch cuts.

#### 3 Power series expansion at the point z = 0

Power series expansion of the local Heun function Hl, such that  $Hl(a, q, \alpha, \beta, \gamma, \delta; 0) = 1$ , is well-known since [3] for  $\gamma \notin \{0\} \cup \mathbb{Z}^-$ . We have

$$Hl(a, q, \alpha, \beta, \gamma, \delta; z) = \sum_{n=0}^{\infty} b_n z^n,$$
 (2)

where the coefficients  $b_n$  are subject to the following three-term recurrence:

$$P_n b_n = Q_n b_{n-1} + R_n b_{n-2}. (3)$$

Here

$$P_{n} = an(\gamma - 1 + n),$$

$$Q_{n} = q + (n - 1)[(a + 1)(\gamma + n - 2) + \varepsilon + a\delta],$$

$$R_{n} = -(n - 2 + \alpha)(n - 2 + \beta)$$
(4)

and the initial conditions are

$$b_{-1} = 0, \qquad b_0 = 1.$$

The Heun function Hl(z) is analytic in the circle  $|z| < R_0$ , where  $R_0$  is the distance from zero to the nearest singular point,  $R_0 = \min\{1, |a|\}$ . Then the famous Cauchy's theorem on the expansion of an analytic function into a power series (see e.g. Theorem 16.7 in [10, Part I]) guarantees that the series (2) converges to Hl inside the circle  $|z| < R_0$ . At this, of course, there might be issues with stability of the recurrence process. For the stability it can be useful to write (2) in the form  $Hl(z) = \sum_{n=0}^{\infty} \tilde{b}_n(z)$ , where  $\tilde{b}_n(z) = b_n z^n$  and the recurrence takes the form:  $P_n \tilde{b}_n = z Q_n \tilde{b}_{n-1} + z^2 R_n \tilde{b}_{n-2}$ .

In the case  $\gamma \in \mathbb{Z}$  the local Frobenius solution corresponding to the smaller exponent (0 or  $1 - \gamma$ ) may contain a logarithmic factor (see e.g. [5]). So, for  $\gamma \in \{0\} \cup \mathbb{Z}^-$  we are looking for the solution of (1), satisfying Hl(0) = 1, in the following form:

$$Hl(a,q,\alpha,\beta,\gamma,\delta;z) = \sum_{n=0,\,n\neq n_*}^{\infty} c_n z^n + \log(z) \sum_{n=n_*}^{\infty} s_n z^n, \tag{5}$$

where  $n_* = 1 - \gamma$ . Note that a solution with the sought property Hl(0) = 1 could be found for any  $c_{n_*}$ . We fix  $c_{n_*} = 0$  for definiteness.

Substituting (5) into (1) we find

$$\log(z)\,\mathcal{L}\left(\sum_{n=n_*}^{\infty}s_nz^n\right) + \mathcal{L}\left(\sum_{n=0,\,n\neq n_*}^{\infty}c_nz^n\right) + \hat{\mathcal{L}}\left(\sum_{n=n_*}^{\infty}s_nz^n\right) = 0,\tag{6}$$

where  $\mathcal{L}$  is the operator of the Heun equation such that (1) is written as  $\mathcal{L}H = 0$ . Besides,

$$\left(\hat{\mathscr{L}}\psi\right)(z) = \frac{2}{z}\psi'(z) + \frac{1}{z}\left(\frac{\gamma-1}{z} + \frac{\delta}{z-1} + \frac{\varepsilon}{z-a}\right)\psi(z).$$

Obviously, (6) separates into two equations

$$\mathcal{L}\left(\sum_{n=n_*}^{\infty} s_n z^n\right) = 0,\tag{7}$$

$$\mathscr{L}\left(\sum_{n=0,\,n\neq n_*}^{\infty} c_n z^n\right) + \hat{\mathscr{L}}\left(\sum_{n=n_*}^{\infty} s_n z^n\right) = 0,\tag{8}$$

Now we collect in (6) terms having identical asymptotic nature as  $z \to 0$ . First we find that coefficients  $c_n$  for  $n = 1, ..., n_* - 1$  are submitted to the recurrence (3):  $P_n c_n = Q_n c_{n-1} + R_n c_{n-2}$ , where  $P_n$ ,  $Q_n$ ,  $R_n$  are defined by (4) and the initial conditions are  $c_{-1} = 0$ ,  $c_0 = 1$ .

Further we find from (7) and (8) that the coefficients  $s_n$  for  $n = n_* + 1, n_* + 2, \ldots$  are submitted to the same recurrence relationships (3):  $P_n s_n = Q_n s_{n-1} + R_n s_{n-2}$ , where  $s_{n_*-1} = 0$  and another initial conditions include coefficients  $c_{n_*-1}$ ,  $c_{n_*-2}$ :

$$an_* s_{n_*} = c_{n_*-1} \left[ q - \gamma(\varepsilon + a\delta - a - 1) \right] - c_{n_*-2} \left[ (1 + \gamma)(2 - \delta - \varepsilon) + \alpha\beta \right].$$

At the next step we can define coefficients  $c_n$  for  $n = n_* + 1, n_* + 2, ...$  From (8) we obtain the following relationship:

$$P_n c_n = Q_n c_{n-1} + R_n c_{n-2} + S_n s_n + T_n s_{n-1} + U_n s_{n-2},$$

$$\tag{9}$$

where

$$S_n = a(1 - \gamma - 2n),$$
  $T_n = \varepsilon + a\delta + (a+1)(\gamma + 2n - 3),$   $U_n = 4 - 2n - \alpha - \beta.$ 

In this way for  $\gamma \in \{0\} \cup \mathbb{Z}^-$ , using (5) we obtain a zero-exponent Frobenius solution, equal to unity at z = 0. The second solution locally can be defined as follows (see also (7)):

$$Hs(a, q, \alpha, \beta, \gamma, \delta; z) = \sum_{n=n_*}^{\infty} \check{s}_n z^n, \tag{10}$$

where  $P_n \check{s}_n = Q_n \check{s}_{n-1} + R_n \check{s}_{n-2}$  for  $n > n_*$  and  $\check{s}_{n_*} = 1$ ,  $\check{s}_{n_*-1} = 0$ .

As it was mentioned above,  $c_{n_*}$  in (5) could be arbitrary. In other words, the linear combination

$$Hl(a, q, \alpha, \beta, \gamma, \delta; z) + C Hs(a, q, \alpha, \beta, \gamma, \delta; z)$$

is equal to unity at z = 0 for an arbitrary constant C. In this sense, for  $\gamma \in \{0\} \cup \mathbb{Z}^-$  the choice of solution Hl is non-unique.

Consider now the function Hs(z) for arbitrary  $\gamma$ . We should distinguish two situations:  $\gamma = 1$  and  $\gamma \neq 1$ . For  $\gamma \neq 1$  we can use the following representation (see Table 2 in [9], index  $[0_-][1_+][a_+][\infty_-]$ ):

$$Hs(a, q, \alpha, \beta, \gamma, \delta; z) = z^{1-\gamma} Hl(a, q - (\gamma - 1)(\varepsilon + a\delta), \beta - \gamma + 1, \alpha - \gamma + 1, 2 - \gamma, \delta; z). \tag{11}$$

Notably, the later formula includes (10) as a particular case, justifying our way to introduce Hl and Hs for non-positive integer  $\gamma$ .

For  $\gamma = 1$ , repeating the arguments used to derive representation of Hl in the case  $\gamma \in \{0\} \cup \mathbb{Z}^-$ , we can find the following local representation

$$Hs(a, q, \alpha, \beta, \gamma, \delta; z) = \sum_{n=1}^{\infty} d_n z^n + \log(z) \sum_{n=0}^{\infty} t_n z^n.$$
 (12)

Here  $P_n t_n = Q_n t_{n-1} + R_n t_{n-2}$ , where  $t_{-1} = 0$ ,  $t_0 = 1$  and (cf. (9))

$$P_n d_n = Q_n d_{n-1} + R_n d_{n-2} + S_n t_n + T_n t_{n-1} + U_n t_{n-2}, \quad d_{-1} = d_0 = 0.$$

#### 4 Power series expansion at an arbitrary point

Further we will extend the local Heun functions outside the circle of convergence of the series (2), (5), (12). To avoid numerical integration of the differential equation (1) we will use analytic continuation process and the power series expansion which is derived in this section.

We seek the solution  $H_{(z_0,H_0,H_0')}(a,q,\alpha,\beta,\gamma,\delta;z)$  to the equation (1) satisfying the conditions

$$H(a, q, \alpha, \beta, \gamma, \delta; z_0) = H_0, \quad \frac{\partial}{\partial z} H(a, q, \alpha, \beta, \gamma, \delta; z) \Big|_{z=z_0} = H'_0.$$
 (13)

Here  $z_0$  is an arbitrary point that is assumed not to belong to the set  $\{0, a, 1, \infty\}$ . We will derive power series expansion of the Heun function in the following form:

$$H_{(z_0, H_0, H_0')}(a, q, \alpha, \beta, \gamma, \delta; z) = \sum_{n=0}^{\infty} c_n (z - z_0)^n.$$
(14)

Substituting (14) into (1) and collecting terms at the same power of  $z - z_0$ , we obtain the following 4-term recurrent relationship defining  $c_n$ :

$$\mathcal{P}_n c_n = \mathcal{Q}_n c_{n-1} + \mathcal{R}_n c_{n-2} + \mathcal{S}_n c_{n-3}, \tag{15}$$

where

$$\mathcal{P}_{n} = -n(n-1)z_{0}(z_{0}-1)(z_{0}-a),$$

$$\mathcal{Q}_{n} = (n-1)\Big\{ \big[ \gamma + \delta + \varepsilon + 3(n-2) \big] z_{0}^{2} + \big[ (a+1)(4-2n-\gamma) - \varepsilon - a\delta \big] z_{0} + a(\gamma + n - 2) \Big\},$$

$$\mathcal{R}_{n} = \Big\{ (n-2) \big[ 2(\gamma + \delta + \varepsilon) + 3(n-3) \big] + \alpha\beta \Big\} z_{0} - q - (n-2) \big[ (a+1)(\gamma + n - 3) + \varepsilon + a\delta \big],$$

$$\mathcal{S}_{n} = (n-3) \big( \gamma + \delta + \varepsilon + n - 4 \big) + \alpha\beta.$$

Obviously, the solution (14) satisfies the conditions (13) if the recurrence process starts with the initial conditions

$$c_{-1} = 0,$$
  $c_0 = H_0,$   $c_1 = H'_0.$ 

The series (14) converges inside the circle  $|z - z_0| < R$ , where R is the distance to the nearest singular point,  $R = \min\{|z_0|, |z_0 - 1|, |z_0 - a|\}$ .

#### 5 Basic algorithm

First we consider evaluation of Hl(z) for  $\gamma \notin \{0\} \cup \mathbb{Z}^-$ . We introduce the projection operator  $\mathscr{P}_z^N$  which, being applied to an analytic function, truncates its power series expansion at the point z to the first N terms. Using the expansion (2), we evaluate

$$(\mathscr{P}_0^N H l)(z) = \sum_{n=0}^N b_n z^n, \quad (\mathscr{P}_0^N H l)'(z) = \sum_{n=1}^N n \, b_n z^{n-1},$$
 (16)

as approximation of Hl(z) and Hl'(z) in a vicinity of z=0.

In our algorithm the number N in the representations (16) is not fixed, it will be defined as we proceed with recurrent computation of  $\tilde{b}_n(z) = b_n z^n$  and summation until a termination condition is satisfied. Namely, we stop the process when  $(\mathscr{P}_0^N H l)'(z)$  and  $(\mathscr{P}_0^{N-1} H l)'(z)$  are not distinguishable in the used computer arithmetics and  $|\tilde{b}_n(z)| < \epsilon$ , where  $\epsilon$  is the machine epsilon.

To estimate the quality of the approximation, in view of (1) we compute the value

$$\hat{Hl}(z) = \frac{1}{q - \alpha\beta z} \left\{ z(z - 1)(z - a) \left( \mathscr{P}_0^N Hl \right)''(z) + \left[ \gamma(z - 1)(z - a) + \delta z(z - a) + \varepsilon z(z - 1) \right] \left( \mathscr{P}_0^N Hl \right)'(z) \right\},$$

where  $(\mathscr{P}_0^N H l)''(z) = \sum_{n=2}^N n(n-1) b_n z^{n-2}$ . Then we suppose proximity of

$$r_0(z) = \left| \hat{Hl}(z) - \left( \mathscr{P}_0^N Hl \right)(z) \right| \tag{17}$$

to the true error of the approximation  $|Hl(z) - (\mathcal{P}_0^N Hl)(z)|$ . Note that numerical computation of  $\hat{Hl}(z)$  near the point  $z = z_* = q/(\alpha\beta)$  is unreliable due to essential loss of significance. It can be suggested for a vicinity of  $z_*$ , say, for  $\{z : |q - \alpha\beta z| < 0.01\}$ , to use an estimate based on properties of the series, e.g.

$$\hat{r}_0(z) = \sqrt{N}|\tilde{b}_N(z)| + \epsilon N |(\mathscr{P}_0^N Hl)(z)|. \tag{18}$$

We can write the described algorithm as a function  $\mathcal{H}_0(z)$  which returns 4-tuple

$$\mathcal{H}_0: z \mapsto [f, f', r, N],$$

where N+1 is the number of terms in power series, defined by the termination condition, r is the value computed with (17) or (18),  $f = (\mathscr{P}_0^N Hl)(z)$  and  $f' = (\mathscr{P}_0^N Hl)'(z)$ .

The scheme of computation of Hl(z) in the case  $\gamma \in \{0\} \cup \mathbb{Z}^-$  is analogous, though a bit more involved. We use the expansion (5) and, instead of (16), define the function  $\mathcal{H}_0(z)$  starting from the expression

$$\sum_{n=0, n \neq n_*}^{N} c_n z^n + \log(z) \sum_{n=n_*}^{N} s_n z^n.$$

Assume that  $|z| < \varkappa R_0$ , where  $R_0 = \min\{1, |a|\}$  is the radius of convergence of the series in (2) or (5) and  $\varkappa \in (0,1)$  is some coefficient chosen so that N defined by the termination condition can be expected to be moderate, say  $\varkappa = 0.5$ . Then we can use the numerical algorithm  $\mathcal{H}_0$  for evaluation of the function H(z), its derivative and for estimation of the approximation error.

Consider further the case  $|z| \ge \varkappa R_0$ . First we define an auxiliary algorithm. Let  $z_0$  be an arbitrary point not belonging to  $\{0, 1, a, \infty\}$ . Using (14) we define

$$\left(\mathscr{P}_{z_0}^N H_{(z_0,H_0,H_0')}\right)(z) = \sum_{n=0}^N c_n (z-z_0)^n, \quad \left(\mathscr{P}_{z_0}^N H_{(z_0,H_0,H_0')}\right)'(z) = \sum_{n=1}^N n \, c_n (z-z_0)^{n-1}.$$

as approximations of  $H_{(z_0,H_0,H'_0)}(z)$  and  $H'_{(z_0,H_0,H'_0)}(z)$  for z close to  $z_0$ . Here coefficients  $c_n$  are defined by the recurrence relationships (15). So we proceed with recurrent computation of  $\tilde{c}_n(z) = c_n(z-z_0)^n$  and summation until the termination condition (analogous to that described above) is satisfied.

We compute

$$\begin{split} \hat{H}_{(z_0,H_0,H_0')}(z) &= \frac{1}{q - \alpha\beta z} \bigg\{ z(z-1)(z-a) \big( \mathscr{P}_{z_0}^N H_{(z_0,H_0,H_0')} \big)''(z) \\ &+ \Big[ \gamma(z-1)(z-a) + \delta z(z-a) + \varepsilon z(z-1) \Big] \big( \mathscr{P}_{z_0}^N H_{(z_0,H_0,H_0')} \big)'(z) \bigg\}, \end{split}$$

and the value

$$r_{(z_0,H_0,H_0')}(z) = |\hat{H}_{(z_0,H_0,H_0')}(z) - (\mathscr{P}_{z_0}^N H_{(z_0,H_0,H_0')})(z)|, \tag{19}$$

supposing its proximity to the true error of the approximation  $|H_{(z_0,H_0,H'_0)}(z) - (\mathscr{P}^N_{z_0}H_{(z_0,H_0,H'_0)})(z)|$ . In view of essential loss of significance in computation of  $\hat{H}_{(z_0,H_0,H'_0)}(z)$ , near  $z=z_*$  we define

$$\hat{r}_{(z_0, H_0, H_0')}(z) = \sqrt{N} |\tilde{c}_N(z)| + \epsilon N | (\mathscr{P}_{z_0}^N H_{(z_0, H_0, H_0')})(z) |.$$
(20)

We can write the described algorithm as a function  $\mathcal{H}_{(z_0,H_0,H'_0)}(z)$  which returns 4-tuple

$$\mathcal{H}_{(z_0,H_0,H'_0)}: z \mapsto [f,f',r,N],$$

where N+1 is the number of terms in power series defined by the termination condition, r is the value computed with (19) or (20),  $f = \left(\mathscr{P}_{z_0}^N H_{(z_0,H_0,H_0')}\right)(z)$  and  $f' = \left(\mathscr{P}_{z_0}^N H_{(z_0,H_0,H_0')}\right)'(z)$ .

Now we are in a position to describe the algorithm based on analytic continuation along a path from zero to z. Consider first the simplest version of the algorithm when the path is the line segment (0, z). At the first step, we compute

$$[Hl_1, Hl'_1, r_1, N_1] = \mathcal{H}_0(z_1),$$

where  $z_1 = e^{i \arg(z)} \varkappa R_0$  (see Fig. 2). Further, for p = 1, 2, and so on, we define

$$R_n = \min\{|z_n|, |z_n - 1|, |z_n - a|\},\$$

$$z_{p+1} = \begin{cases} z_p + e^{i \arg(z)} \varkappa R_p & \text{if } |z - z_p| > \varkappa R_p, \\ z & \text{if } |z - z_p| \le \varkappa R_p, \end{cases}$$

and compute

$$[\mathit{Hl}_{p+1},\mathit{Hl}'_{p+1},r_{p+1},N_{p+1}]=\mathcal{H}_{(z_p,\mathit{Hl}_p,\mathit{Hl}'_p)}(z_{p+1}).$$

![](_page_6_Figure_0.jpeg)

Figure 2: Analytic continuation using power series.

The iterations stops when  $z_{p+1} = z$ . Finally, we have  $Hl_{p+1}$ ,  $Hl'_{p+1}$  as approximations of Hl(z) and Hl'(z), respectively. We also compute the values  $r_{\Sigma} = r_1 + \ldots + r_{p+1}$  and  $N_{\Sigma} = N_1 + \ldots + N_{p+1} + p + 1$ . Here  $N_{\Sigma}$  is the total number of power series terms which can be used as a measure of computer load and  $r_{\Sigma}$  may be an indicator of the approximation quality.

The described above algorithm of continuation along a line segment is readily generalized for the case when 0 and z are connected by a polyline  $\Upsilon$ . This gives us a way to compute the multi-valued Heun function. The resulting procedure can be considered as a function  $\mathcal{H}(\Upsilon)$  which returns 4-tuple

$$\mathring{\mathcal{H}}: \Upsilon \mapsto [f, f', r_{\Sigma}, N_{\Sigma}],$$

where f and f' are the resulting approximations of the Heun function at z and its derivative. It should be noted that the function Hs(z) for  $\gamma \neq 1$  is expressed through Hl(z) via (11) and to define  $\mathcal{H}s(\Upsilon)$  through  $\mathcal{H}l(\Upsilon)$  one should only compute the multi-valued function  $z^{-\gamma}$  along  $\Upsilon$ . In the case  $\gamma = 1$  the procedure of analytic continuation described above for Hl can be applied with simple modification—it should start from the expansion (12).

It is easy to observe that the size of the step in the described analytic continuation is small for points of the polyline  $\Upsilon$  close to one of the singular points. This also means an increase of the number of used circular elements in the continuation procedure which, in its turn, may lead to loss of accuracy. The influence of the singular points can be reduced by a choice of the path of continuation.

Let us give some details of application of the above algorithm for computation of the single-valued Heun functions. We will only use simple paths consisting of two or three line segments. We denote  $\zeta_1 = a$  and  $\zeta_2 = 1$  if |a| < 1 and  $\zeta_1 = 1$  and  $\zeta_2 = a$  otherwise. Let  $R_j = \min\{|\zeta_j|, |\zeta_1 - \zeta_2|\}$  (j = 1, 2) and  $\eta_j$  be the closest to  $\zeta_j$  point of the line segment (0, z) (see fig. 3). Note that  $\zeta_j \neq \eta_j$  due to the assumption  $z \notin \mathcal{B}_{1\infty} \cup \mathcal{B}_{a\infty}$ . Let  $d_j = |\zeta_j - \eta_j|$  be the distance from  $\zeta_j$  to (0, z). If  $\eta_j$  is an internal point of (0, z) and  $d_j < R_j/2$ , then we introduce new point

$$\zeta_j = \zeta_j + \exp(i\pi/2 + i\arg(z))\min\{R_j/2, |z - \zeta_j|\}\operatorname{sign}(\operatorname{Im}(z/\zeta_j)).$$

In this way we obtain the sequence of points  $\Pi$ , which is either [0, z] or  $[0, \zeta_1, z]$  or  $[0, \zeta_2, z]$  or  $[0, \zeta_1, \zeta_2, z]$ . Then we define the path  $\Upsilon$  that consequently connects the points of  $\Pi$  by line segments (see an example of three-segment path in fig. 3) and define  $\mathcal{H}(z) = \mathcal{H}(\Upsilon)$  and  $\mathcal{H}(z) = \mathcal{H}(\Upsilon)$ .

![](_page_7_Figure_0.jpeg)

Figure 3: Path from zero to z consisting of three line segments

## 6 Computation of single-valued Heun functions near singular points

In this section we consider computation of single-valued Heun functions. We have already noted that the number of circular elements in the analytic continuation procedure becomes large as z approaches one of the singular points  $\{1, a, \infty\}$ . So, we suggest an improvement of the algorithm; for this we will compute local solutions in a vicinity of the singular point, using results of [9] and writing these solutions in terms of Hl(z) and Hs(z).

Let z be close to 1. We use the following representation of the Heun function:

$$Hl(a, q, \alpha, \beta, \gamma, \delta; z) = C_1 Hl(1 - a, \alpha\beta - q, \alpha, \beta, \delta, \gamma; 1 - z) + C_2 Hs(1 - a, \alpha\beta - q, \alpha, \beta, \delta, \gamma; 1 - z),$$
(21)

where  $C_1$ ,  $C_2$  are some constants, the Heun functions in the right-hand side (the first and the second local solutions in a vicinity of z=1) correspond to the index  $[1_+0_+][a_+][\infty_+]$  in Table 2 [9]. There is a complication when the point 1 belongs to  $\mathcal{B}_{a\infty}$  ( $a \in (0,1)$ ). In this case the coefficients in (21) are generally different for z belonging to the upper and the lower half-spaces:  $C_1^{(\pm)}$ ,  $C_2^{(\pm)}$  for  $\{z:\pm \operatorname{Im}(z)>0\}$ .

Consider branch cuts of  $Hl(1-a,\ldots;1-z)$  and  $Hs(1-a,\ldots;1-z)$ . The points of the branch cut, corresponding to  $\mathcal{B}_{a\infty}$  for  $Hl(a,\ldots;z)$  and  $Hs(a,\ldots;z)$ , are defined by the equation 1-z=s(1-a), where  $s\in(1,+\infty)$ . The set of points z=1+s(a-1) constitutes the ray emanating from a to  $\exp\{i\arg(a-1)\}\infty$  and, thus, it is located outside the circle of radius |a-1| with center at z=1. Analogously, consider the branch cut of  $Hl(1-a,\ldots;1-z)$  and  $Hs(1-a,\ldots;1-z)$  corresponding to the branch cut  $\mathcal{B}_{1\infty}$ . The points of the branch cut are defined by the equation 1-z=s, where  $s\in(1,+\infty)$ , thus, lying outside the circle of radius 1 with center at z=1. The last possible branch cut of  $Hl(1-a,\ldots;1-z)$  and  $Hs(1-a,\ldots;1-z)$  is defined as 1-z=q,  $q\in(-\infty,0)$  and coincides with the branch cut  $\mathcal{B}_{1\infty}$  of the function Hl(z) in the left hand side of (21).

In order to use the representation (21) we should define the coefficients  $C_1$ ,  $C_2$  (or, in the special case,  $C_1^{(\pm)}$ ,  $C_2^{(\pm)}$ ). Since for the general Heun equation an explicit solution to the two-point connection problem is not known, we will find the matching coefficients in the following way. First, we define a matching point  $z = z_1^*$ , preferably being equally close to zero and 1 and distant from z = a. In our algorithm we fix

$$z_1^{\star} = \frac{1}{2} + \frac{\mathrm{i}s}{\sqrt{2}}, \text{ where } s = \begin{cases} \operatorname{sign}(\operatorname{Im}(z)) & \text{if } a \in (0,1), \\ -\operatorname{sign}(\operatorname{Im}(a)) & \text{otherwise.} \end{cases}$$

Further, we apply the algorithms  $\mathcal{H}$  and  $\mathcal{H}s$  described in § 5 and find

$$[f_0, f'_0, r_0, N_0] = \mathcal{H}(a, q, \alpha, \beta, \gamma, \delta, z_1^*),$$
  

$$[f_1, f'_1, r_1, N_1] = \mathcal{H}(1 - a, \alpha\beta - q, \alpha, \beta, \delta, \gamma, 1 - z_1^*),$$
  

$$[f_2, f'_2, r_2, N_2] = \mathcal{H}(1 - a, \alpha\beta - q, \alpha, \beta, \delta, \gamma, 1 - z_1^*).$$

Then we solve the linear system

$$\begin{pmatrix} f_1 & f_2 \\ -f_1' & -f_2' \end{pmatrix} \begin{pmatrix} C_1 \\ C_2 \end{pmatrix} = \begin{pmatrix} f_0 \\ f_0' \end{pmatrix}$$

and obtain the matching coefficients  $C_1$ ,  $C_2$ . It may also be reasonable to keep the found values  $C_1 = C_1(a, q, \alpha, \beta, \gamma, \delta)$ ,  $C_2 = C_2(a, q, \alpha, \beta, \gamma, \delta)$  in computer memory.

On finding  $C_1$ ,  $C_2$  (by computation or in the computer memory) we define

$$\mathcal{H}^{(1)}: z \mapsto [f, f', r, N],$$

where

$$f = C_1 f_1 + C_2 f_2, \quad f' = -C_1 f'_1 - C_2 f'_2, \quad r = |C_1| r_1 + |C_2| r_2, \quad N = N_1 + N_2,$$
$$[f_1, f'_1, r_1, N_1] = \mathcal{H}(1 - a, \alpha\beta - q, \alpha, \beta, \delta, \gamma, 1 - z),$$
$$[f_2, f'_2, r_2, N_2] = \mathcal{H}s(1 - a, \alpha\beta - q, \alpha, \beta, \delta, \gamma, 1 - z).$$

The described scheme can be repeated literally to define  $\mathcal{H}_{s}^{(1)}$  based on the representation

$$Hs(a,q,\alpha,\beta,\gamma,\delta;z) = C_1' Hl(1-a,\alpha\beta-q,\alpha,\beta,\delta,\gamma;1-z) + C_2' Hs(1-a,\alpha\beta-q,\alpha,\beta,\delta,\gamma;1-z), \tag{22}$$

where  $C_1'$  and  $C_2'$  are some coefficients to be found (or two sets  ${C'}_1^{(-)}$ ,  ${C'}_2^{(-)}$  and  ${C'}_1^{(+)}$ ,  ${C'}_2^{(+)}$ —in the special case  $a \in (0,1)$ ).

We note that finding  $C_1$ ,  $C_2$  or  $C'_1$ ,  $C'_2$  demands computation of  $Hl(a, ...; z_1^*)$  or  $Hs(a, ...; z_1^*)$ , and both terms  $Hl(1-a, ...; 1-z_1^*)$  and  $Hs(1-a, ...; 1-z_1^*)$  in the right hand side of (21) or (22). So, if the matching constants are not known, the algorithms  $\mathcal{H}^{(1)}$  and  $\mathcal{H}^{(1)}$  are preferable over  $\mathcal{H}$  and  $\mathcal{H}^{(1)}$  in a sufficiently small vicinity of z=1. In the code used in § 7 the algorithms are applied for  $|z-1| < 0.05 R_1$  ( $R_1 = \min\{1, |a-1|\}$ ) when the matching constants are not known and for  $|z-1| < 0.25 R_1$  when the coefficients are already precomputed.

In a very similar way we consider the case when z is located near a. We will use the following representation of the Heun function:

$$Hl(a, q, \alpha, \beta, \gamma, \delta; z) = D_1 Hl\left(\frac{a-1}{a}, \alpha\beta - \frac{q}{a}, \alpha, \beta, \varepsilon, \gamma; \frac{a-z}{a}\right) + D_2 Hs\left(\frac{a-1}{a}, \alpha\beta - \frac{q}{a}, \alpha, \beta, \varepsilon, \gamma; \frac{a-z}{a}\right), \quad (23)$$

where  $D_1$ ,  $D_2$  are some constants and the two functions in the right-hand side correspond to index  $[a_+0_+1_+][\infty_+]$  in Table 2 [9]. There is a complication when the point a belongs to  $\mathscr{B}_{1\infty}$  or  $\mathscr{B}_{0\infty}$ . In this case we generally have different coefficients in (21) in the upper and the lower half-spaces:  $D_1^{(\pm)}$ ,  $D_2^{(\pm)}$  for  $\{z:\pm \operatorname{Im}(z)>0\}$ .

Consider branch cuts of  $Hl((a-1)/a,\ldots;\zeta)$  and  $Hs((a-1)/a,\ldots;\zeta)$ , where  $\zeta=(a-z)/a$ . The points of the branch cut that emanates from the singularity  $\zeta=(a-1)/a$  are defined in the z-plane by the equation a-z=s(a-1), where  $s\in(1,+\infty)$ . This set of points constitutes the ray emanating from 1 to  $\exp\{i\arg(1-a)\}\infty$  and located outside the circle of the radius |1-a| with center at z=a. The branch cut of  $Hl((a-1)/a,\ldots;\zeta)$  and  $Hs((a-1)/a,\ldots;\zeta)$  emanating from  $\zeta=1$  is a ray going from zero to  $\exp\{i\arg(-a)\}\infty$  (z=a(1-s)); it is located outside the circle of radius |a| with center at z=a. The third possible branch cut of  $Hl((a-1)/a,\ldots;\zeta)$  and  $Hs((a-1)/a,\ldots;\zeta)$  is defined as  $a-z=aq,\ q\in(-\infty,0)$  and coincides with the branch cut  $\mathscr{B}_{a\infty}$  of the function  $Hl(a,\ldots;z)$ .

In order to find the coefficients  $D_1$ ,  $D_2$  (or, in the special case,  $D_1^{(\pm)}$ ,  $D_2^{(\pm)}$ ), first, we define a matching point  $z = z_a^*$ , preferably being equally close to zero and a and distant from z = 1. We fix

$$z_a^{\star} = \frac{a}{2} + \frac{\mathrm{i}s}{\sqrt{2}}, \text{ where } s = \begin{cases} \mathrm{sign}(\mathrm{Im}(z)) & \text{if } a \in \mathscr{B}_{0\infty} \cup \mathscr{B}_{1\infty}, \\ a|a|^{-1} \, \mathrm{sign}(\mathrm{Im}(a)) & \text{otherwise.} \end{cases}$$

Further, we find

$$\begin{split} [f_0,f_0',r_0,N_0] &= \mathcal{H}(a,q,\alpha,\beta,\gamma,\delta,z_a^\star), \\ [f_1,f_1',r_1,N_1] &= \mathcal{H}\Big(\frac{a-1}{a},\alpha\beta-\frac{q}{a},\alpha,\beta,\varepsilon,\gamma,\frac{a-z_a^\star}{a}\Big), \\ [f_2,f_2',r_2,N_2] &= \mathcal{H}\!\!s\Big(\frac{a-1}{a},\alpha\beta-\frac{q}{a},\alpha,\beta,\varepsilon,\gamma,\frac{a-z_a^\star}{a}\Big), \end{split}$$

and  $D_1$ ,  $D_2$  are defined as solution to the linear system

$$\begin{pmatrix} f_1 & f_2 \\ -a^{-1}f_1' & a^{-1}f_2' \end{pmatrix} \begin{pmatrix} D_1 \\ D_2 \end{pmatrix} = \begin{pmatrix} f_0 \\ f_0' \end{pmatrix}.$$

It may also be reasonable to keep the found values  $D_1 = D_1(a, q, \alpha, \beta, \gamma, \delta), D_2 = D_2(a, q, \alpha, \beta, \gamma, \delta)$ in computer memory.

On finding  $D_1$ ,  $D_2$  (by computation or in the computer memory) we define

$$\mathcal{H}^{(a)}: z \mapsto [f, f', r, N],$$

where

$$f = D_{1}f_{1} + D_{2}f_{2}, \quad f' = -D_{1}a^{-1}f'_{1} - D_{2}a^{-1}f'_{2}, \quad r = |D_{1}|r_{1} + |D_{2}|r_{2}, \quad N = N_{1} + N_{2},$$

$$[f_{1}, f'_{1}, r_{1}, N_{1}] = \mathcal{H}\left(\frac{a-1}{a}, \alpha\beta - \frac{q}{a}, \alpha, \beta, \varepsilon, \gamma, \frac{a-z}{a}\right),$$

$$[f_{2}, f'_{2}, r_{2}, N_{2}] = \mathcal{H}s\left(\frac{a-1}{a}, \alpha\beta - \frac{q}{a}, \alpha, \beta, \varepsilon, \gamma, \frac{a-z}{a}\right).$$

The described scheme can be repeated literally to define  $\mathcal{H}^{(a)}$  based on the representation

$$Hs(a, q, \alpha, \beta, \gamma, \delta; z) = D'_1 Hl\left(\frac{a-1}{a}, \alpha\beta - \frac{q}{a}, \alpha, \beta, \varepsilon, \gamma; \frac{a-z}{a}\right) + D'_2 Hs\left(\frac{a-1}{a}, \alpha\beta - \frac{q}{a}, \alpha, \beta, \varepsilon, \gamma; \frac{a-z}{a}\right),$$

where  $D_1'$  and  $D_2'$  (or, in the special case,  $D_1'^{(\pm)}$ ,  $D_2'^{(\pm)}$ ) are some coefficients to be found. We note again that the algorithms  $\mathcal{H}^{(a)}$  and  $\mathcal{H}^{(a)}$  are preferable over  $\mathcal{H}$  and  $\mathcal{H}$ s in a sufficiently small vicinity of z=a. In the code used in §7 the algorithms are applied for  $|z-a|<0.05\,R_a$  $(R_a = \min\{|a|, |1-a|\})$  when the matching constants are not known and for  $|z-a| < 0.25 R_a$  when the coefficients are already precomputed.

Consider now vicinity of the point  $z = \infty$ . The branch cuts  $\mathscr{B}_{0\infty}$ ,  $\mathscr{B}_{1\infty}$  and  $\mathscr{B}_{a\infty}$  split the vicinity of infinity  $V_{\infty} = \{z : |z| > \max(1, |a|)\}$  into three sectors, further denoted as  $S^{(1)}$ ,  $S^{(2)}$ ,  $S^{(3)}$ . Coefficients connecting the function Hl(z) with two local solutions at infinity are defined for each of the sectors separately. Note that there is a special case when a is real, then the number of sectors decreases to two.

We use the following representation of the Heun function near  $z = \infty$  for  $z \in S^{(k)}$  (k = 1, 2, 3):

$$Hl(a, q, \alpha, \beta, \gamma, \delta; z) = E_1^{(k)} z^{-\alpha} Hl\left(\frac{1}{a}, \frac{q + \alpha(\delta - \beta)}{a} + \alpha(\varepsilon - \beta), \alpha, \alpha - \gamma + 1, \alpha - \beta + 1, \delta; \frac{1}{z}\right) + E_2^{(k)} z^{-\alpha} Hs\left(\frac{1}{a}, \frac{q + \alpha(\delta - \beta)}{a} + \alpha(\varepsilon - \beta), \alpha, \alpha - \gamma + 1, \alpha - \beta + 1, \delta; \frac{1}{z}\right), \quad (24)$$

where  $E_1^{(k)}$ ,  $E_2^{(k)}$  are some constants. Further we will omit the superscript for brevity. The two functions in the right-hand side correspond to the index  $[\infty_{+}0_{+}][1_{+}][a_{+}]$  in Table 2 [9].

Obviously, the branch cuts  $\mathscr{B}_{a\infty}$  and  $\mathscr{B}_{1\infty}$  transform for  $Hl(1/a,\ldots,1/z)$  and  $Hs(1/a,\ldots,1/z)$ to line segments connecting zero with z=a and z=1 respectively. These branch cuts are located outside  $V_{\infty}$ . The branch cut  $\mathscr{B}_{0\infty}$  suits both sides of (24).

To find the coefficients  $E_1$ ,  $E_2$ , first, we define a matching point  $z=z_{\infty}^{\star}$ . In our algorithm we fix

$$z_{\infty}^{\star} = C_{\infty} R_{\infty} e^{i\omega_k}$$
, for  $z \in S^{(k)}$ ,  $k = 1, 2, 3$ ,

where  $\omega_k$  is the angle of the mean line of the sector  $S^{(k)}$ ,  $R_{\infty} = \max\{1, |a|\}$  and  $C_{\infty} > 1$  ( $C_{\infty} = 2$  in the computations, presented below). Then, we find

$$[f_0, f'_0, r_0, N_0] = \mathcal{H}(a, q, \alpha, \beta, \gamma, \delta, z_{\infty}^{\star}),$$

$$[f_1, f'_1, r_1, N_1] = \mathcal{H}\left(\frac{1}{a}, \frac{q + \alpha(\delta - \beta)}{a} + \alpha(\varepsilon - \beta), \alpha, \alpha - \gamma + 1, \alpha - \beta + 1, \delta, \frac{1}{z_{\infty}^{\star}}\right),$$

$$[f_2, f'_2, r_2, N_2] = \mathcal{H}\left(\frac{1}{a}, \frac{q + \alpha(\delta - \beta)}{a} + \alpha(\varepsilon - \beta), \alpha, \alpha - \gamma + 1, \alpha - \beta + 1, \delta, \frac{1}{z_{\infty}^{\star}}\right),$$

and the matching coefficients are defined as solution to the linear system

$$\begin{pmatrix} [z_{\infty}^{\star}]^{-\alpha}f_1 & [z_{\infty}^{\star}]^{-\alpha}f_2 \\ -[z_{\infty}^{\star}]^{-\alpha-1}(f_1'/z_{\infty}^{\star} + \alpha f_1) & -[z_{\infty}^{\star}]^{-\alpha-1}(f_2'/z_{\infty}^{\star} + \alpha f_2) \end{pmatrix} \begin{pmatrix} E_1 \\ E_2 \end{pmatrix} = \begin{pmatrix} f_0 \\ f_0' \end{pmatrix}.$$

Finally, we define a function  $\mathcal{H}^{(\infty)}$  returning values of  $f(z) \approx Hl(z), f'(z) \approx Hl'(z)$  for sufficiently large |z|. On finding  $E_1$ ,  $E_2$  (by computation or in the computer memory) we compute

$$\mathcal{H}^{(\infty)}: z \mapsto [f, f', r, N],$$

where

$$f = E_{1}z^{-\alpha}f_{1} + E_{2}z^{-\alpha}f_{2}, \quad f' = -E_{1}z^{-\alpha-1}(f'_{1}/z + \alpha f_{1}) - E_{2}z^{-\alpha-1}(f'_{2}/z + \alpha f_{2}),$$

$$r = |E_{1}z^{-\alpha}|r_{1} + |E_{2}z^{-\alpha}|r_{2}, \quad N = N_{1} + N_{2},$$

$$[f_{1}, f'_{1}, r_{1}, N_{1}] = \mathcal{H}\left(\frac{1}{a}, \frac{q + \alpha(\delta - \beta)}{a} + \alpha(\varepsilon - \beta), \alpha, \alpha - \gamma + 1, \alpha - \beta + 1, \delta, \frac{1}{z}\right),$$

$$[f_{2}, f'_{2}, r_{2}, N_{2}] = \mathcal{H}\left(\frac{1}{a}, \frac{q + \alpha(\delta - \beta)}{a} + \alpha(\varepsilon - \beta), \alpha, \alpha - \gamma + 1, \alpha - \beta + 1, \delta, \frac{1}{z}\right).$$

The scheme can be repeated literally to define  $\mathcal{H}_s^{(\infty)}$  based on the representation for  $z \in S^{(k)}$ 

$$\begin{split} Hs(a,q,\alpha,\beta,\gamma,\delta;z) &= E_1'^{(k)} z^{-\alpha} \, Hl\Big(\frac{1}{a},\frac{q+\alpha(\delta-\beta)}{a} + \alpha(\varepsilon-\beta),\alpha,\alpha-\gamma+1,\alpha-\beta+1,\delta;\frac{1}{z}\Big) \\ &+ E_2'^{(k)} z^{-\alpha} \, Hs\Big(\frac{1}{a},\frac{q+\alpha(\delta-\beta)}{a} + \alpha(\varepsilon-\beta),\alpha,\alpha-\gamma+1,\alpha-\beta+1,\delta;\frac{1}{z}\Big), \end{split}$$

where  $E'_1^{(k)}$  and  $E'_2^{(k)}$  are some coefficients to be found. We note again that the algorithms  $\mathcal{H}^{(\infty)}$  and  $\mathcal{H}^{(\infty)}$  are preferable over  $\mathcal{H}$  and  $\mathcal{H}$  for sufficiently large |z|. In the code used in the following section the algorithms are applied for  $|z| > 5 C_{\infty} R_{\infty}$  when the matching coefficients are not known and for  $|z| > C_{\infty} R_{\infty}$  otherwise.

#### Numerical results 7

In this section we present some results of numerical evaluation of the function Hl. For this we will use the algorithm  $\mathcal{H}$  with the improvements described in the previous section  $(\mathcal{H}^{(1)}, \mathcal{H}^{(a)})$  and  $\mathcal{H}^{(\infty)}$ .

It is known that in some cases Heun functions can be expressed in terms of hypergeometric functions (see [6, 13, 8, 16, 1] and references therein). For a test of the present algorithm we use the formula given in [1]:

$$Hl\left(4, \frac{9}{4}, \frac{3}{2}, \frac{3}{2}, \frac{1}{2}, 2, z\right) = {}_{2}F_{1}\left(\frac{1}{2}, \frac{1}{2}, \frac{1}{2}, 1 - (z - 1)^{2}\left(1 - \frac{z}{4}\right)\right) = h(z) := \frac{2}{\sqrt{4 - z}(1 - z)},$$

where  ${}_{2}F_{1}$  is the ordinary hypergeometric function.

We define

$$\Delta(z) := Hl\left(4, \frac{9}{4}, \frac{3}{2}, \frac{3}{2}, \frac{1}{2}, 2, z\right) - \frac{2}{\sqrt{4-z}(1-z)}$$

![](_page_11_Figure_0.jpeg)

![](_page_11_Figure_1.jpeg)

Figure 5: Values of max 1.5, log<sup>10</sup> N(z) .

and check numerically the identity ∆(z) = 0. Calculations are performed in the numerical computing environment GNU Octave and double precision (64-bit) arithmetics (the machine epsilon is about 2.22 · 10−16).

Figure 4 shows in a semilogarithmic scale results of computations of the relative error

$$\Lambda(z) = \frac{|\Delta(z)|}{1 + |h(z)|} + \frac{|\Delta'(z)|}{1 + |h'(z)|}.$$

on the grid (Re z,Im z) ∈ L(1000, [−20, 20]) × L(1000, [−20, 20]), where L(n, χ) is the set consisting of n linearly spaced in the interval χ values (including interval's end-points). It is easy to note that accuracy loss occurs mainly near the singularity z = 4. The maximum error for Λ(z) found on the grid is  $1.9635 \cdot 10^{-14}$ . So, as one can see, the accuracy of computation of the function and its derivative is rather satisfactory for the used double-float arithmetics.

In Fig. 5 we present (in a semilogarithmic scale) the value N(z) which means the total number of terms in power series used to compute  $Hl\left(4,\frac{9}{4},\frac{3}{2},\frac{3}{2},\frac{1}{2},2,z\right)$ . This value varies from N(0)=1 and can characterize the time of computation. To be more specific, in the Table we present typical times of evaluation in Octave of  $Hl\left(4,\frac{9}{4},\frac{3}{2},\frac{3}{2},\frac{1}{2},2,z\right)$  via the basic algorithm  $\mathcal{H}$  and via the algorithm  $\mathcal{H}$  improved by  $\mathcal{H}^{(1)}$ ,  $\mathcal{H}^{(a)}$  and  $\mathcal{H}^{(\infty)}$ . The computer in use has a 3.4 GHz Intel Core i5 CPU and 8 Gb of RAM.

| z                | Time (sec.), basic algorithm |      | improved algorithm subsequent evaluations |
|------------------|------------------------------|------|-------------------------------------------|
| 20i              | 0.04                         | 0.05 | 0.007                                     |
| $20 + i\epsilon$ | 0.075                        | 0.05 | 0.007                                     |
| -20              | 0.035                        | 0.05 | 0.003                                     |
| 0.99             | 0.03                         | 0.02 | 0.003                                     |
| 4 + 0.01i        | 0.075                        | 0.05 | 0.002                                     |

#### Acknowledgements

The author is indebted to Professor A.M. Ishkhanyan for drawing attention to the subject and useful discussions. The author thanks Professor A.Ya. Kazakov for interest in the results of the present work and useful remarks.

#### References

- [1] T. Birkandan, Physical examples of the Heun-to-hypergeometric reduction, arXiv:1401.0449v2 [math-ph].
- [2] A. Erdélyi, W. Magnus, F. Oberhettinger, F.G. Tricomi, *Higher Transcendental Functions*, Vol. III, McGraw-Hill, New York, 1955.
- [3] K. Heun, Zur Theorie der Riemann'schen Functionen zweiter Ordnung mit vier Verzweigungspunkten, *Math. Ann.*, 1889, **33**, p. 161–179.
- [4] M. Hortaçsu, Heun Functions and their uses in Physics, arXiv:1101.0471.
- [5] E. Kamke, Differentialgleichungen: Lösungsmethoden und Lösungen. Bd. 1: Gewöhnliche Differentialgleichungen, Leipzig: Akad. Verlag, 1944.
- [6] K. Kuiken, Heun's equation and the hypergeometric equation, SIAM J. Math. Anal., 1979, 10, pp. 655–657.
- [7] W. Lay, S. Yu. Slavyanov, Heun's equation with nearby singularities, *Proc. R. Soc. Lond. A*, 1999, 455, pp. 4347–4361.
- [8] R. S. Maier, On reducing the Heun equation to the hypergeometric equation, *J. Differential Equations*, 2005, Vol. 213, Issue 1, p. 171–203.
- [9] R. S. Maier, The 192 solutions of the Heun equation, Mathematics of Computation, 2007, Vol. 76(258), p. 811–843; arXiv:math/0408317.
- [10] A. I. Markushevich, *Theory of Functions of a Complex Variable*, 3 volumes in one, Chelsea Publishing Company, New York, 1977.
- [11] O.V. Motygin, Matlab/Octave code for evaluation of the Heun functions, 2015, online: https://github.com/motygin/Heun\_functions/.

- [12] A. Ronveaux (Ed.), Heun's Differential Equations, Oxford University Press, Oxford, 1995.
- [13] A. V. Shanin, R. V. Craster, Removing false singular points as a method of solving ordinary differential equations, Euro. Jnl of Applied Mathematics, 2002, vol. 13, pp. 617–639.
- [14] S. Yu. Slavyanov, W. Lay, Special Functions, Oxford University Press, Oxford, 2000.
- [15] B. D. Sleeman, V. B. Kuznetzov, Heun functions, In: F.W. J. Olver, D. M. Lozier, R. F. Boisvert, et al., NIST Handbook of Mathematical Functions, Cambridge University Press, 2010.
- [16] G. Valent, Heun functions versus elliptic functions, Difference equations, special functions and orthogonal polynomials, World Sci. Publ., Hackensack, NJ, 2007, pp. 664–686, arXiv:math-ph/0512006.