![](_page_0_Picture_1.jpeg)

Contents lists available at [ScienceDirect](http://www.elsevier.com/locate/physd)

## Physica D

journal homepage: [www.elsevier.com/locate/physd](http://www.elsevier.com/locate/physd)

![](_page_0_Picture_5.jpeg)

# Painlevé IV: A numerical study of the fundamental domain and beyond

![](_page_0_Picture_7.jpeg)

Jonah A. Reeger [a,](#page-0-0)[∗](#page-0-1) , Bengt Fornberg [b](#page-0-2)

- <span id="page-0-0"></span><sup>a</sup> *Air Force Institute of Technology, Department of Mathematics and Statistics, 2950 Hobson Way, Wright-Patterson AFB, OH 45433, USA*
- <span id="page-0-2"></span><sup>b</sup> *University of Colorado, Department of Applied Mathematics, 526 UCB, Boulder, CO 80309, USA*

## h i g h l i g h t s

- First comprehensive study of the solution space to the fourth Painlevé equation.
- Extensions to parameter regimes that were unreachable by previous approaches.
- Previously known solutions found to be 'non-typical' of the general case.
- Solutions featuring singularity-free half-planes explored.

## a r t i c l e i n f o

#### *Article history:* Received 20 August 2013 Received in revised form 8 April 2014 Accepted 9 April 2014 Available online 24 April 2014 Communicated by P.D. Miller

*Keywords:* Painlevé equation Painlevé transcendent PIV Pole field Connection formula

## a b s t r a c t

The six Painlevé equations were introduced over a century ago, motivated by rather theoretical considerations. Over the last several decades, these equations and their solutions, known as the Painlevé transcendents, have been found to play an increasingly central role in numerous areas of mathematical physics. Due to extensive dense pole fields in the complex plane, their numerical evaluation remained challenging until the recent introduction of a fast 'pole field solver' (Fornberg and Weideman, 2011). The fourth Painlevé equation has two free parameters in its coefficients, as well as two free initial conditions. After summarizing key analytical results for PIV, the present study applies this new computational tool to the fundamental domain and a surrounding region of the parameter space. We confirm existing analytic and asymptotic knowledge about the equation and also explore solution regimes which have not been described in the previous literature.

<span id="page-0-3"></span>Published by Elsevier B.V.

## **1. Introduction**

With the increasing presence of the Painlevé equations in the reduction of partial differential equations (PDEs) [\[1–4\]](#page-11-0) and the subjects of combinatorics [\[5–7\]](#page-11-1), orthogonal polynomials [\[8–12\]](#page-11-2), statistical physics [\[13–18\]](#page-12-0), integrable continuous dynamical systems [\[19](#page-12-1)[,20\]](#page-12-2) and quantum physics [\[21–24\]](#page-12-3) a greater understanding of the solution space for each of the six equations is important. A collection of applications specific to the PIV equation is presented in [\[25\]](#page-12-4). In the past, solutions that are pole free along the real axis have proven to be particularly relevant. As a resource for the future, one present goal has been to identify such cases, as well as

The solutions of the six Painlevé equations (PI–PVI) are free from movable branch points, but with the possibility of movable poles or movable isolated essential singularities ([\[26\]](#page-12-5), Section 32.2). This Painlevé property inspired the introduction of a novel numerical approach [\[27\]](#page-12-6) – combining a Padé-based ODE solver [\[28\]](#page-12-7) with a partly randomized integration path strategy – that allowed for the first time rapid numerical solutions of the Painlevé equations over extended regions in the complex plane. It was first used for P<sup>I</sup> [\[27\]](#page-12-6) and later for PII [\[29\]](#page-12-8). It was then applied to the fourth Painlevé equation

$$\frac{d^2}{dz^2}u(z) = \frac{1}{2u(z)} \left(\frac{d}{dz}u(z)\right)^2 + \frac{3}{2}u(z)^3 + 4zu(z)^2 + 2\left(z^2 - \alpha\right)u(z) + \frac{\beta}{u(z)},$$
(1)

those with pole-free sectors in the complex-plane, throughout the PIV equations four-parameter solution space.

<span id="page-0-1"></span><sup>∗</sup> Corresponding author. Tel.: +1 7246893830. *E-mail addresses:* [jonah.reeger@gmail.com,](mailto:jonah.reeger@gmail.com) [jonah.reeger@afit.edu](mailto:jonah.reeger@afit.edu) (J.A. Reeger), [fornberg@colorado.edu](mailto:fornberg@colorado.edu) (B. Fornberg).

in the special case of α = β = 0 [\[30\]](#page-12-9). As in these three previous numerical studies, computational explorations in this paper are limited to solutions *u*(*z*)that are real when *z* is real, although some of the presented theory includes solutions that are not always real on the real axis.

For a small portion of the two-dimensional (α, β)-parameter space there exist examples of solutions expressible as rational functions or in terms of special functions, such as the parabolic cylinder function. These well-documented solutions appear, however, as only isolated points or one-parameter families in the fourdimensional space of parameters and initial conditions (ICs). Much of the present study is focused on the distribution of singularities for solutions to [\(1\).](#page-0-3) These are all first-order poles with residue +1 or −1.

The solutions presented in this paper are parameterized by α, β and the values for *u*(0) and *u* ′ (0). Another way to parameterize the solution space is through four ''Stokes multipliers'', which provides a link between the PIV equation and a certain Riemann–Hilbert problem [\[31–35\]](#page-12-10). This approach is particularly well suited for analytical work such as connection formulas and far-field asymptotics. Distant pole field structures can also be approximated via suitable transformations [\[36,](#page-12-11)[37\]](#page-12-12) (however, the present focus is more on pole-free regions). The two parameterization approaches can be related to each other utilizing, for instance, the software RHPackage [\[38\]](#page-12-13) to solve the Riemann–Hilbert problem (given a set of Stokes multipliers and parameters to define α and β) to determine the corresponding set of values for *u*(0) and *u* ′ (0).

#### *1.1. Organization of the paper*

Section [2](#page-1-0) recalls some available theory, including symmetries in PIV and different solution transformations. Section [3](#page-1-1) discusses closed form solutions of PIV, in particular solutions in terms of rational and elementary special functions and also the asymptotic behaviors presented in the literature. This is followed in Section [4](#page-2-0) by the numerical approach used here to explore the much larger space of solutions for which no closed form solutions are available. Sections [5](#page-4-0) and [6](#page-6-0) describe such explorations of the parameter and solution spaces, first highlighting the ''fundamental domain'' and then extending into inspections of the previously unexplored region of β > 0, for which no instances of closed form solutions or transformations have been reported.

#### <span id="page-1-0"></span>**2. Symmetries and solution hierarchies**

This section describes the known symmetries in the PIV equation and transformations that relate solutions for different parameter choices.

#### <span id="page-1-6"></span>*2.1. Symmetries in the equation*

Let PIV(α, β) be the set of all solutions of [\(1\)](#page-0-3) for the particular α and β. Direct inspection of [\(1\)](#page-0-3) shows that if *u*(*z*) ∈ PIV(α, β), then [\[39\]](#page-12-14)

$$-u(-z) \in P_{IV}(\alpha, \beta),$$
 (2)

$$-iu(-iz) \in P_{IV}(-\alpha, \beta), \text{ and}$$
 (3)

$$iu(iz) \in P_{IV}(-\alpha, \beta).$$
 (4)

Incidentally the first of these symmetries also holds for PIII (for all parameter choices), but never for any of the other Painlevé equations. Due to these symmetries, any solution presented in this paper has at least one other counterpart for the same choice of α and β.

## <span id="page-1-7"></span>*2.2. The Bäcklund and Schlesinger transformations*

The equations PII through PVI have collections of transformations relating solutions for given parameters to those of different choices. For instance, [\[39–42\]](#page-12-14) collectively present sixteen such transformations for PIV. Some of these transformations were not always presented correctly. Updated expressions along with computational verification of their forms can be found in [\[43\]](#page-12-15).

## <span id="page-1-1"></span>**3. Closed form solutions and asymptotic approximations**

Before discussing the closed form solutions presented in the literature, we note again that these at most form two-dimensional manifolds in the four-dimensional solution space. That is, they provide a very limited view of the solution types that are possible.

## *3.1. Rational solutions*

The PIV equation has six different sequences of parameter choices leading to rational solutions expressible in terms of either Generalized Hermite or Generalized Okamoto polynomials [\[44\]](#page-12-16), with two particular choices leading to the only known entire solutions, −2*z* and −(2/3)*z*. The locations of the parameter choices in the (α, β)-plane leading to such solutions are shown in [Fig. 1](#page-2-1) as dark (blue) and light (yellow) hexagrams for Generalized Hermite and Generalized Okamoto polynomials, respectively.

## *3.2. Special function solutions*

In addition to the rational solutions, PIV admits solutions that are described by combinations of special functions; cf. [\[26\]](#page-12-5), chapters 12 and 13. This includes solutions expressible in terms of parabolic cylinder functions, *D*<sup>ν</sup> (ζ ) [\[5](#page-11-1)[,45\]](#page-12-17), and, as discovered more recently, solutions in terms of the confluent hypergeometric function, <sup>1</sup>*F*1(*a*, *b*; ζ ), [\[23,](#page-12-18)[24\]](#page-12-19). In either case, for each of the appropriate choices of α and β there is a one parameter family of solutions that are expressible in terms of these special functions. [Fig. 1](#page-2-1) displays the locations of all such parameter choices as black curves.

Three distinct types of solutions have been proposed in the form of determinants involving parabolic cylinder functions [\[5](#page-11-1)[,45\]](#page-12-17). However, only one of these expressions has been confirmed numerically [\[43\]](#page-12-15).

## <span id="page-1-5"></span>*3.3. Asymptotic approximations*

Beyond the known closed form solutions, it is noted in [\[26\]](#page-12-5), section 32.11, that when β = 0, there are solutions that decay asymptotically along the real axis either as *z* → +∞ or *z* → −∞. These solutions result from assuming that the second derivative term in [\(1\)](#page-0-3) is negligible.

When assuming instead that both the first and second derivative terms in [\(1\)](#page-0-3) are negligible, the method of dominant balance (see, e.g., [\[46\]](#page-12-20), section 3.4) leads to the quartic equation

<span id="page-1-2"></span>
$$\frac{3}{2}\hat{w}(z)^4 + 4z\hat{w}(z)^3 + 2(z^2 - \alpha)\hat{w}(z)^2 + \beta = 0.$$
 (5)

<span id="page-1-4"></span>Each root of [\(5\)](#page-1-2) provides a leading asymptotic term for solutions that are smooth as *z* → ±∞. Any number of further terms then follow by substitution into [\(1\).](#page-0-3) For instance, [\(6\)](#page-1-3) through [\(9\)](#page-2-2) illustrate the first two terms.

<span id="page-1-3"></span>
$$w_{+1}^{+}(z;\alpha,\beta) = \frac{\sqrt{-2\beta}}{2z} + \frac{\alpha\sqrt{-2\beta} + 2\beta}{4z^{3}} + O\left(\frac{1}{z^{5}}\right)$$
 (6)

<span id="page-1-8"></span>
$$w_{-1}^{+}(z;\alpha,\beta) = -\frac{\sqrt{-2\beta}}{2z} + \frac{-\alpha\sqrt{-2\beta} + 2\beta}{4z^{3}} + O\left(\frac{1}{z^{5}}\right)$$
(7)

<span id="page-2-1"></span>![](_page_2_Figure_2.jpeg)

**Fig. 1.** Two views of the Weyl Chambers. The shaded region indicates the real part of the fundamental domain. Both figures show several of the chambers and locations of the rational and special function solutions to PIV (dark hexagrams (blue) represent generalized Hermite type, light hexagrams (yellow) show generalized Okamoto type, and lines (black) show parabolic cylinder and confluent hypergeometric types). (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

$$w_{+1}^{-}(z;\alpha,\beta) = -\frac{2}{3}z + \frac{\alpha}{z} - \frac{2 + 6\alpha^2 + 9\beta}{8z^3} + 0\left(\frac{1}{z^5}\right)$$
(8)

$$w_{-1}^{-}(z;\alpha,\beta) = -2z - \frac{\alpha}{z} + \frac{2 + 6\alpha^2 + \beta}{8z^3} + O\left(\frac{1}{z^5}\right). \tag{9}$$

Solutions matching these behaviors as *z* → +∞ and *z* ∈ R do not oscillate and are free from poles for all *z* ∈ R greater than some finite value *z*<sup>0</sup> ∈ R. No other such solutions were observed in the numerical explorations. An analogous statement could be made as *z* → −∞ by the symmetry [\(2\).](#page-1-4) With the assumption *u*(*z*) ∈ R for *z* ∈ R, only the latter two are available as asymptotic approximations when β > 0. Later in this paper, ICs leading to solutions matching [\(6\)–\(9\)](#page-1-3) as *z* → +∞ will be marked in several figures (described as pole counting diagrams) as shown in [Fig. 4.](#page-4-1)

As with the information presented in this section, the rest of this paper will discuss only the asymptotic behaviors as *z* → +∞, *z* ∈ R, since the symmetry [\(2\)](#page-1-4) makes it clear that there are analogous solutions with similar asymptotic behaviors as *z* → −∞. This is also seen by comparing the left and right frames in [Fig. 3.](#page-4-2)

#### <span id="page-2-3"></span>*3.4. The parameter space and the Weyl chambers*

Based on the various symmetries, solution hierarchies and known closed form solutions, the parameter space of PIV with β ≤ 0 can be described in terms of the so-called Weyl chambers (see e.g., [\[47\]](#page-12-21), section II-A, or [\[39\]](#page-12-14), section 26). These chambers feature a complete regularity in the (α, <sup>√</sup> −2β)-plane.

The significance of the Weyl Chambers, when extended to complexes α and β, is that a single chamber in theory provides all of the information to construct solutions for every arbitrary (α, β) pair ([\[39\]](#page-12-14), section 25). The real part of this chamber is shown as a shaded gray region in [Fig. 1,](#page-2-1) often called the fundamental domain. An explicit description of the fundamental domain can be found in [\[39\]](#page-12-14), Section 26.

Notice that in the case of solutions that are real along the real axis, every parameter choice leading to a rational or special function solution of PIV has β ≤ 0. This is also true for the decaying asymptotic approximation discussed in Section [3.3.](#page-1-5) For this reason, part of this study will be devoted to the unexplored region of β > 0.

## <span id="page-2-0"></span>**4. The numerical method and exploration approach**

<span id="page-2-2"></span>Explorations of the four-dimensional space of parameters and ICs require a fast numerical method and a systematic approach for comparing solutions of different parameter choices. These techniques are discussed here.

The extensive pole fields appearing in these solutions have motivated the development of various solution techniques. However, many of the previous methods were limited in the choice of the parameters in the coefficients by considering special forms of the equation (e.g. Riemann–Hilbert problems [\[31\]](#page-12-10)), constrained to the real axis [\[48–50\]](#page-12-22), or restricted to a small domain around the origin [\[51\]](#page-12-23). The present method extends the 'pole vaulting' idea [\[48\]](#page-12-22) in three fundamental ways: (i) use of a 'pole friendly' ODE integrator [\[28\]](#page-12-7), (ii) not using any rigid choices of diversion paths around a pole, but instead utilizing a freely branching network of paths in the complex plane, and (iii) targeting paths toward whole regions in the complex plane (rather than only toward other real axis locations). Some of the existing numerical methods were surveyed in [\[52\]](#page-12-24).

## *4.1. A brief description of the numerical method*

The numerical scheme introduced in [\[27\]](#page-12-6) can achieve very high orders of accuracy, with minimal loss of accuracy even in the vicinity of poles. This is accomplished utilizing a flexible path selection strategy that can efficiently cover large areas of the complex plane starting from arbitrary ICs and for any choice of α and β. When integrating from one start location to a single end location this scheme uses the following strategy, which will be called pole avoidance:

- 1. Choose the location of the initial condition as the first expansion point.
- 2. Compute the Padé approximation about the expansion point.
- 3. Evaluate the Padé approximation a distance *h* away in each of five directions in a swath directed toward the target point and choose as the next expansion point the one with the smallest solution magnitude.
- 4. Unless the target point has been reached, return to step 2.

<span id="page-3-1"></span>![](_page_3_Figure_2.jpeg)

**Fig. 2.** Left: Pole locations and residues in the solution [\(10\)](#page-3-0) of the PIV equation. Poles of residue +1 are marked with dark (blue) circles and poles of residue −1 are marked with light (yellow) circles. Center: Log base 10 of the relative error in the numerical solution [\(10\)](#page-3-0) of [\(1\)](#page-0-3) computed at machine precision using order 14 Padé approximations (the numerator and denominator both polynomials of order 7) and step size of 0.0625. Right: Log base 10 of the relative error in the numerical solution [\(10\)](#page-3-0) of [\(1\)](#page-0-3) computed at 25 digits of accuracy using order 40 Padé approximations (the numerator and denominator both polynomials of order 20) and step size of 0.25. In both cases the pole field solver is started with a high precision initial condition at *z* = 0.5. Notice the similarity in the error patterns of the numerical solutions. (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

This pole avoidance strategy is effective when finding the solution to an initial value problem (IVP) at a single point. However, if the solution is desired at many different points (for instance, for the visualization of the solution over a region in the complex plane) the method is extended to the pole field solver.

- 1. Set up a coarse grid of target points in the complex plane.
- 2. Select the target points in random order.
- 3. Apply the pole avoidance strategy to reach a predetermined neighborhood of the current target point, starting from the closest point that has already been evaluated. In the first step this is the location of the IC.
- 4. Once all of the coarse grid target points have been accounted for, set up a fine grid of the desired evaluation points.
- 5. Compute from the end of each of the previous paths (using a single Padé expansion) values of nearby fine grid evaluation points.

The majority of the explorations conducted in this paper utilized Padé approximations with the numerator and denominator both polynomials of order 7, resulting in an expansion with truncation error *O*(*h* <sup>14</sup>). The coefficients of these polynomials were computed using the Toeplitz approach as in [\[27\]](#page-12-6). The step size was typically chosen as *h* = 0.0625. In [\[43\]](#page-12-15) it is illustrated that the numerical method, when all computations are performed to machine precision, with this combination of polynomial order and continuation step size can easily produce solutions with relative error on the order of 10−<sup>5</sup> or better, except when computing over large regions with no poles.

The initial value problem is well-posed (and the numerical procedure is correspondingly well-conditioned) when meandering through a pole field, allowing high orders of accuracy to be maintained. However, across smooth areas the initial value problem is exponentially unstable. The center frame of [Fig. 2](#page-3-1) illustrates this result for the following closed form solution of PIV [\[53\]](#page-12-25) whose pole locations and residues are shown in the left frame of the figure. The solution is

$$u(z) = -z - \frac{d}{dz} \ln \left( \frac{W(v_0(z), v_1(z))}{W(v_0(z), v_1(z), v_2(z))} \right)$$

$$v_0(z) = \exp \left( -\frac{1}{2} z^2 \right) {}_{1}F_{1} \left( -\frac{9}{4}, \frac{1}{2}; z^2 \right)$$
(10)

$$v_1(z) = -\frac{9}{\sqrt{2}} \exp\left(-\frac{1}{2}z^2\right) z_1 F_1\left(-\frac{5}{4}, \frac{3}{2}, z^2\right)$$

$$v_2(z) = -\frac{3}{2} \exp\left(-\frac{1}{2}z^2\right) \left[3_1 F_1\left(-\frac{5}{4}, \frac{3}{2}, z^2\right) - 5z^2_1 F_1\left(-\frac{1}{4}, \frac{5}{2}, z^2\right)\right]$$

with W the Wronskian and <sup>1</sup>*F*1(*a*, *b*; ζ ) the confluent hypergeometric function.

To verify the features in the solutions presented in this paper, we utilized the extended precision solver discussed in [\[43\]](#page-12-15) with all computations performed with precision on the order of 10−<sup>25</sup> . Both the step sizes and the orders of the Padé approximations were adjusted to ensure adequate accuracy. Reliance on extended precision solvers provides only a partial remedy to computing solutions over smooth regions. The right frame of [Fig. 2](#page-3-1) shows that accuracy still decreases when considering points increasingly far from a pole field when using an extended precision solver. The loss of accuracy would have been reduced in this example if the paths had been constrained to remain within them. A complete remedy against the loss of accuracy over smooth regions can be obtained by combining the initial value problem within pole fields with boundary value problems that link these near-field solutions to far-field asymptotics.

Further tests are available for checking the accuracy in a solution that has no closed form. For instance, the randomized path selection strategy lends to comparisons of the same solution computed twice (as this would use different paths). Since the present study only considers solutions for which *u*(*z*) ∈ R when *z* ∈ R, differences between the upper (Im(*z*) > 0) and lower (Im(*z*) < 0) half-planes can likewise be used for error estimation.

#### *4.2. Pole counting*

<span id="page-3-0"></span>The pole field solver makes it possible to rapidly obtain solutions for a variety of initial conditions. To explore the differences in solution characteristics for each fixed choice of α and β, the number of poles on either the positive or negative real axis is examined for varying (*u*(0), *u* ′ (0)) ∈ R 2 . This, paired with the asymptotic behavior discussed in Section [3.3,](#page-1-5) allows the characterization of the numerous solution possibilities for each fixed α and β. [Fig. 3](#page-4-2)

<span id="page-4-2"></span>![](_page_4_Figure_2.jpeg)

**Fig. 3.** Number of poles on the positive and negative real axes for  $\alpha = 0$  and  $\beta = 0$ . For a description of the markers and shading see Fig. 4.

<span id="page-4-1"></span>![](_page_4_Figure_4.jpeg)

**Fig. 4.** Legend and gray scale/pattern bar for Figs. 3, 6 and 8–10. The legend shows the markers indicating the ICs that generate the dominant asymptotic behaviors (6)–(9) and closed form solutions. If a marker occurs on a curve, then the dominant behavior or type of closed form solution occurs for all of the ICs along that curve. If a marker is emphasized by containing an "×", then it indicates an isolated IC matching the dominant behavior or the IC generates an isolated rational solution. The gray-scale/pattern bar on the right indicates the number of poles on the positive or negative real axis.

(adapted from [30]) provides a prototypical example in the case of  $\alpha = \beta = 0$ . This figure displays the number of poles on the positive and negative real axes for each choice of initial conditions shown, and each of the frames is denoted a pole counting diagram.

Consider, for now, only the right frame in Fig. 3, since the left is completely analogous due to the symmetries discussed in Section 2.1. Each of the ICs marked by a curve or contained within a shaded region generates a solution with a finite number of poles on the positive real axis. The gray scale/pattern bar in Fig. 4 indicates the exact number of poles for a given initial condition with darker and lighter shades indicating odd and even numbers of poles, respectively. On the other hand, ICs neither contained in a shaded region nor marked by a curve generate solutions with an infinity of poles on the corresponding half (positive/negative) of the real axis.

In this case of  $\alpha = \beta = 0$ , each of the shaded regions in the right half-plane contains ICs that generate solutions with an odd number of poles on the positive real axis, while the u(0), u'(0) values along the isolated curves lead to solutions with an even number. The opposite holds in the left half-plane.

Most of the ICs in the shaded regions generate solutions that oscillate as  $z \to +\infty$ ; however, each initial condition marked by a curve, located at the boundary of a shaded region, or designated by an isolated marker has no oscillations as  $z \to +\infty$ . These solutions

<span id="page-4-3"></span>![](_page_4_Figure_10.jpeg)

Fig. 5. The eight sectors of poles in the solutions of P<sub>IV</sub>.

are precisely those approximated by the roots of the quartic equation (5) as  $z \to +\infty$ . The appropriate root is indicated by the symbols shown in left frame of Fig. 4. When two markers appear along the same curve, those ICs generate solutions matching both behaviors (in separate intervals of the real axis), as shown in, for example, Figs. 7 and 15.

#### <span id="page-4-0"></span>5. Numerical illustrations of the fundamental domain

Solution types occurring for parameter choices in the fundamental domain are discussed in the following sections. These (and subsequent) sections describe some solutions as having adjacent pole-free sectors. This terminology arises from asymptotic (and numerical) knowledge that the poles in the solutions of  $P_{\rm IV}$  align in the eight sectors shown in Fig. 5. Further discussions of these sectors are available in [30].

## 5.1. An exploration of the fundamental domain

In Section 3.4 the fundamental domain was introduced, and it was noted that solutions for all parameter choices in theory can be found by applying the transformations discussed in Section 2.2 to the solutions in this domain. However, the literature describes solutions in this domain only for the cases  $\alpha=\beta=0$  (numerical and asymptotically decaying solutions), ( $\alpha=0$ ,  $\beta=-2/9$ ) (a rational solution), along the line  $\beta=0$  (asymptotically decaying solutions), and along the curve  $\beta=-2(\alpha-1)^2$  (asymptotically

<span id="page-5-0"></span>![](_page_5_Figure_2.jpeg)

**Fig. 6.** Number of poles on the positive real axis for  $\alpha=0$  and  $\beta=-2/9$ ,  $\alpha=0$  and  $\beta=-2$ , and  $\alpha=1$  and  $\beta=0$ . A detailed description of the markers and shading is given in Fig. 4.

decaying, rational and special function solutions). All of these occur on the boundary of the fundamental domain.

#### *5.1.1.* Parameter choices with rational or special function solutions

It should again be noted that, for each of the parameter choices  $(\alpha=0,\beta=-\frac{2}{9})$  and along the curve  $\beta=-2(\alpha-1)^2$ , the closed form or asymptotic solutions only lead to a single solution or a one parameter family of solutions in the u(0) versus u'(0) plane. To gain some insight into arbitrary ICs (in the same manner as Fig. 3) the frames in Fig. 6 show the number of poles appearing on the positive real axis for each of the two remaining vertices of the fundamental domain, as well as the case  $(\alpha=0,\beta=-2/9)$ .

Within the frames of Figs. 3 and 6 it is easy to see that often the ICs of solutions matching (6) through (9) along the real axis as  $z \to +\infty$  appear as the boundaries of regions of initial conditions generating solutions with a finite number of poles on the real axis. In some cases, however, such ICs appear along curves that are not part of these boundaries. In such a case, these curves correspond to the ICs of solutions that decay asymptotically (see Section 3.3) or those of one parameter families of solutions in terms of special functions.

Further study of the last two frames of Fig. 6 highlights a peculiar behavior of the solutions matching (6) through (9) when the  $\alpha$  and  $\beta$  choices occur at the vertex of a Weyl chamber. For these solutions, two or three of the behaviors  $w_{\mu}^+$ ,  $\mu=\pm 1$ , and  $w_{-1}^-$  are present in the same solution (as indicated by two separate markers occurring along the same line), but in different segments of the positive real axis. Take, for instance, the ICs for ( $\alpha=0$ ,  $\beta=-2$ ) indicated by the arrow in the second frame of Fig. 6. Along the curve containing these ICs there are three separate markers. The solutions in a neighborhood of these particular ICs are shown in Fig. 7, illustrating that different dominant asymptotic behaviors can occur in the same solution (again, in different segments of the real axis).

#### 5.1.2. Parameter choices along the boundary $\beta = 0$

When the boundary  $\beta=0$  is considered, the literature generally only describes solutions to  $P_{\rm IV}$  that decay asymptotically as  $z\to +\infty$  [26,54]. Connection formulas are available that indicate some of these solutions also match -2z and -(2/3)z as  $z\to -\infty$  [26,54,55]. Through the symmetries discussed in Section 2.1 there must also be solutions that decay asymptotically as  $z\to -\infty$  and match -2z or -(2/3)z as  $z\to +\infty$ . The frames in Fig. 3 and the right frame in Fig. 6 illustrate the locations of the ICs that generate such solutions. They are those ICs that match  $w_\mu^+$ ,

 $\mu=\pm 1$ , in the case of asymptotic decay and  $w_{-1}^-$  and  $w_{+1}^-$  in the cases of -2z and -(2/3)z, respectively. In all cases along the boundary  $\beta=0$  we find that the behaviors  $w_{\mu}^+$ ,  $\mu=\pm 1$ , are identical in line with (6) and (7).

## 5.1.3. Parameter choices along the boundary $\beta = -2(\alpha - 1)^2$

 $P_{IV}$  has a one-parameter family of solutions expressible in terms of the parabolic cylinder function or confluent hypergeometric function for each choice of  $\alpha$  and  $\beta$  along the boundary described by  $\beta=-2(\alpha-1)^2$ . The leading order asymptotic behavior of these solutions can differ for distinct choices of  $\alpha$  and  $\beta$  along this boundary. That is, these special function solutions match the dominant behavior of  $w_{-1}^- \sim -2z$  in the case of  $(\alpha=0,\ \beta=-2)$  and  $w_{-1}^+$  in the case of  $(\alpha=1,\ \beta=0)$ . The center and right frames of Fig. 6 provide examples of choices of parameters along this boundary.

# 5.1.4. Parameter choices along the boundary $\alpha=0$ and interior to the fundamental domain

Along the boundary  $\alpha=0$  and interior to the fundamental domain solutions matching each of (6) through (9) are generated from distinct ICs. There are now no solutions that match two of these behaviors at the same time, but in different segments of the real axis. As a prototype, the distinctness of these ICs can be witnessed in the leftmost frame of Fig. 6.

## 5.1.5. A note on connection formulae

Consider the left and right frames of Fig. 3, showing the number of poles along the negative and positive real axes, respectively. One finds that a segment of the curve extending from the origin and down to the right in the right frame cuts across the shaded region that extends from the origin up and to the right in the left frame. Along this segment  $P_{\rm IV}$  therefore has solutions that are smooth in both directions. A similar analysis of the pole counting diagrams for any choice of  $\alpha$  when  $\beta=0$  would result in an analogous family of solutions that are smooth in both directions. These appear to be the only examples of solutions that have connection formulae available in the literature [26,54,55].

Examination of Fig. 6 (together with the symmetry (2)) shows that similar comparisons of the number of poles on the positive and negative real axes simultaneously will again identify solutions that are smooth in both directions for regions of ICs near u(0) = u'(0) = 0 in cases where  $\beta$  is negative. It turns out that such regions (sometimes only a curve) appear to exist for all parameter choices within the fundamental domain. Fig. 8 illustrates this for a choice interior to the fundamental domain ( $\alpha = 0.25$ ,  $\beta = -0.125$ ).

<span id="page-6-2"></span>![](_page_6_Figure_2.jpeg)

Fig. 7. Solutions (along the real axis (top) and pole locations and residues (bottom)) with adjacent pole-free sectors for  $\alpha=0$  and  $\beta=-2$ . The initial condition u(0) for each column is shown at the top and, in each column, u'(0)=0 and  $u_0=3.170110354518507$ . The initial condition ( $u(0)=u_0$  and u'(0)=0) is marked with an arrow in the center frame of Fig. 6. In the bottom row of figures poles of residue +1 are marked with dark (blue) circles and poles of residue -1 are marked with light (yellow) circles. (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

<span id="page-6-1"></span>![](_page_6_Figure_4.jpeg)

**Fig. 8.** Number of poles on the negative real axis (left), entire real axis (center), and positive real axis (right) for  $\alpha=0.25$  and  $\beta=-0.125$ .

In the following section, Fig. 10 will show that similar regions will also occur outside the fundamental domain when  $\beta < 0$ , however, with the difference that there now may be a finite number of poles on the real axis in either one or both directions. In contrast, positive choices of  $\beta$  do not seem to produce any such regions of ICs

## <span id="page-6-0"></span>6. Solution patterns outside the fundamental domain

The  $(\alpha, \beta)$  space is far too wide to complete an exhaustive survey of here. Therefore, the rest of this paper focuses on the unexplored space of  $\beta > 0$  and highlights some solution types

that seem to appear for all  $\alpha$  and  $\beta$ . In this section, the terminology "nearly pole-free half-plane" refers to a solution with a finite number of poles for all z such that  $\operatorname{Re}(z) > z_0$  or  $\operatorname{Re}(z) < z_0$  for some finite  $z_0 \in \mathbb{R}$ .

## 6.1. The unexplored space of positive beta

Information about  $P_{IV}$  with  $\beta>0$  is noticeably absent from the literature. For instance, all known closed form solutions occur only when  $\beta$  is nonpositive. Even the Bäcklund and Schlesinger transformations are only applicable to  $\beta$ -values that are nonpositive (assuming that u(z) is real when z is real). Exploration of such

<span id="page-7-1"></span>![](_page_7_Figure_2.jpeg)

<span id="page-7-0"></span>**Fig. 9.** Number of poles on the positive real axis for parameter choices where  $\alpha=0$  and  $\beta>0$ . A detailed description of the markers and shading is given in Fig. 4. Similar diagrams for  $\alpha\neq0$  and larger  $\beta$  are shown in Fig. 10.

![](_page_7_Figure_4.jpeg)

**Fig. 10.** Number of poles on the positive real axis at the six  $(\alpha, \beta)$  locations marked in Fig. 11 (all exterior to the fundamental domain). The initial conditions for solutions asymptotic to  $w_{-1}^-$  in the top middle, bottom left, and bottom right frames occur outside of the domain shown at (u(0) = -4.6822, u'(0) = 20.7787), (u(0) = -10.7942, u'(0) = 120.3759), and (u(0) = 49.4606, u'(0) = -2442.3215), respectively. A detailed description of the markers and shading is given in Fig. 4.

cases and knowledge of the tronqueé like solutions that appear in the  $\alpha=\beta=0$  case suggests that solutions with  $\beta>0$  also feature noteworthy characteristics. For instance, there are further analogues to the solution that is pole free for a half-plane.

The dominant asymptotic behaviors (6) and (7) no longer occur as solutions that are real along the real axis, due to the term  $\sqrt{-2\beta}$ . Therefore, the Fig. 9 is much simpler than its counterparts for  $\beta \leq 0$  with a single IC matching the behavior of  $w_{+1}^{-1} \sim -(2/3)z$  and ICs along the boundaries of regions with finite poles generating solutions that match  $w_{-1}^{-1} \sim -2z$ .

#### 6.2. Parameters larger in magnitude

This section illustrates some  $\alpha$ ,  $\beta$  choices slightly larger in magnitude. When  $\beta>0$  there is little difference in the pole counting diagrams from the choices presented in the earlier figures. However, choices of  $\beta\leq 0$  become far more complicated without indicating the existence of further types of solutions with special characteristics. Even parameter choices in adjacent Weyl chambers generate significantly different behaviors near u(0)=u'(0)=0. Fig. 10 displays some choices for  $\alpha$  and  $\beta$  larger in magnitude.

<span id="page-8-0"></span>![](_page_8_Figure_2.jpeg)

**Fig. 11.** Number of poles on the positive (right) and negative (left) real axis for solutions matching w <sup>+</sup><sup>1</sup> ∼ −(2/3)*<sup>z</sup>* as *<sup>z</sup>* → +∞ and *<sup>z</sup>* ∈ <sup>R</sup> and each α and β. The solid curves indicate the boundaries of the Weyl chambers, while the dashed lines show the boundaries of regions of finite poles on both the positive and negative real axes. Note that in this case β > 0 implies an infinity of poles along R <sup>−</sup>. The circles (red) containing an × indicate those parameters shown in [Fig. 10.](#page-7-0) The changes in shading occur simultaneously in the left and right frames corresponding to a pole moving from one half of the real axis (positive/negative) to the other. (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

<span id="page-8-1"></span>![](_page_8_Figure_4.jpeg)

**Fig. 12.** Zero and pole locations of solutions to [\(1\)](#page-0-3) with various even integer values of α with β = 0. The values of *m* in the titles correspond to the discussion at the beginning of Section [6.3.1.](#page-9-0) Poles of residue +1 are marked with dark (blue) circles and poles of residue −1 are marked with light (yellow) circles. Locations of zeros are marked by ×'s (red). (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

## *6.3. Solutions with a nearly pole-free half-plane*

It was noted in [\[30\]](#page-12-9) that, when α = β = 0, there exists a particular solution with an entire half-plane free of poles. The present study has found that solutions with a nearly pole-free halfplane are not confined to only this special choice of α and β. In fact, evidence suggests that for each α and β there exists at least one such solution, and very likely only one. The likelihood that there is

<span id="page-9-1"></span>![](_page_9_Figure_2.jpeg)

Fig. 13. Zero and pole locations for three solutions to  $P_{IV}$ . The middle figure displays the zeros and poles of the solution matching the bottom left frame of Fig. 12. This figure matches -(2/3)z as  $z \to +\infty$  for the case of  $\alpha = 6$  and  $\beta = 0$ , where  $u(0) \approx -0.15560961$  and  $u'(0) \approx 0.30039611$ . The top and bottom frames occur by perturbing  $\beta$  by +0.1 and -0.1, respectively, while leaving  $\alpha$ , u(0) and u'(0) unchanged. Poles of residue +1 are marked with dark (blue) circles and poles of residue -1 are marked with light (yellow) circles. Locations of zeros are marked by  $\times$ 's (red). (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

only one such solution for each  $\alpha$  and  $\beta$  pair makes this solution a prime candidate for comparing and making connections between all parameter choices.

For each  $\alpha$  and  $\beta$  this special solution type matches the root  $w_{-1}^- \sim -(2/3)z$  as  $z \to +\infty$  and  $z \in \mathbb{R}$ . Knowing this, computing the initial conditions leading to such a solution is a simple matter of solving a boundary value problem (BVP). Applying the familiar methodology of counting poles along the positive, and now negative, real axes allows the identification of further special characteristics of these solutions.

In Fig. 11 the pole counts are shown along the negative and positive real axes (left and right frames, respectively) overlayed with the Weyl chambers marked by solid curves. Also in these frames, dashed lines mark the boundaries of regions in the  $\alpha$  versus  $\beta$  plane where these solutions have only a finite number of poles on the negative real axis. Notice that these dashed curves form a regular structure similar to that of the Weyl chambers, with the parabolas offset by one unit on the  $\alpha$  axis and the horizontal lines occurring at  $\beta$  values where these new parabolas and those from the Weyl chambers intersect.

## <span id="page-9-0"></span>6.3.1. The tops of the parabolas

To begin, consider the parameter choices at the tops of these new parabolas. These occur at  $\alpha=2m$  and  $\beta=0$ ,  $m\in\mathbb{Z}$ . In these cases the poles nearest the origin form very regular patterns. Examples for several different choices of m are shown in Fig. 12. Notice the pole structure near the center of these figures. When m<0 poles of residue +1 align in a structure similar to the roots with a positive real part of the degree m Okamoto I polynomial, while poles of residue -1 appear similar to the roots of the degree m-1 polynomial. On the other hand, when m>0 the poles of residue +1 (likewise, -1) align in a structure similar to all of the roots of the order -10 -11 (likewise, -12) polynomials. Note that the Okamoto I polynomials in this context are singly indexed as in [47] while those in the rational solutions of -12 are doubly indexed generalized Okamoto polynomials as in [44].

Although this study is focused on pole locations for  $P_{IV}$  solutions, Fig. 12 offers a good opportunity to also make some brief comments about their zeros. Based on how u(z) appears in the denominator of the two terms on the right hand side of (1), we can

<span id="page-10-0"></span>![](_page_10_Figure_2.jpeg)

**Fig. 14.** Solutions (pole locations and residues) normal to the parabola  $\beta = -2(\alpha - 2)^2$ . All frames depict the solutions asymptotic to -(2/3)z as  $z \to +\infty$ . The center frames occur directly along the parabolas where  $\alpha = \alpha_0^1 = 1.25$  (top) and  $\alpha = \alpha_0^R = 2.75$  (bottom). The left and right frames in both the top and bottom then depict the solutions along the line normal to the parabola at  $\alpha = \alpha_0$  at  $\alpha_0 \pm 10^{-4}$ . Poles of residue +1 are marked with dark (blue) circles and poles of residue -1 are marked with light (yellow) circles. (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

#### deduce:

- If  $\beta = 0$ , then u(z) = 0 implies that u'(z) = 0. That is, every zero has at least multiplicity two.
- If  $\beta \neq 0$ , every zero is simple.
- If  $\beta > 0$ , u(z) cannot have a zero along the real axis (given our assumption that u(z) is real valued there).

Perturbing  $\beta$  slightly away from zero in the last subplot of Fig. 12 (but keeping  $\alpha$ , u(0), and u'(0) unchanged) causes pole fields to enter from both sides, and leads to the solutions shown in the top and bottom frames of Fig. 13. Although some zero pairs are too tight to be clearly distinguishable as such, the figure nevertheless illustrates all of the observations above.

## 6.3.2. Along the parabolas

When  $\alpha$  and  $\beta$  are taken directly on the parabolas the solutions asymptotic to -(2/3)z are nonoscillatory as  $z \to -\infty$ . An example of this is shown in the center frame of Fig. 14. Now, if  $\alpha$  or  $\beta$  are varied slightly such that the choice of parameters no longer falls on one of the parabolas, these solutions can have either an infinity of poles or oscillate as  $z \to -\infty$ . Examples of this are also shown in the left and right frames of Fig. 14.

#### 6.3.3. When $\beta$ is positive

For  $\beta>0$ , Fig. 11 showed that all of the solutions asymptotic to  $w_{-1}^-\sim -(2/3)z$  as  $z\to +\infty$ ,  $z\in\mathbb{R}$ , have an infinity of poles on the negative real axis. These solutions also do not generally have an entire half-plane free of poles. Instead, numerical evidence points

to a value  $z_0 \in \mathbb{R}$  such that for all z with  $\text{Re}(z) > z_0$  the solution has no poles.

## 6.3.4. Other solutions with a pole-free half-plane

The solutions asymptotic to -(2/3)z as  $z\to +\infty$  are not the only solutions that have a half-plane pole free. There are, of course, the rational solutions. Likewise, there are solutions expressible in terms of parabolic cylinder or confluent hypergeometric functions that also feature a pole-free half-plane. Generally, these other solutions with a pole-free half-plane match different roots of (5) as  $z\to +\infty$  than  $w_{-1}^-\sim -(2/3)z$ .

## 6.4. Solutions with adjacent pole-free sectors

Common among the solutions that match  $w_{\mu}^{\pm}$  as  $z \to +\infty$  and  $z \in \mathbb{R}$  is the absence of poles in sectors adjacent to the positive real axis. The solutions are similar to the tronquée solutions of  $P_{\rm I}$  discussed in [30]. For both  $P_{\rm I}$  and  $P_{\rm IV}$  (with  $\alpha=\beta=0$ ) these solutions are characterized by at least two adjacent pole-free sectors. Fig. 7 first illustrated that for certain choices of  $\alpha$  and  $\beta$  the two or three of the behaviors  $w_{\mu}^{\pm}$ ,  $\mu=\pm 1$ , can occur along the positive real axis in the same solution, but in different segments of the axis. Fig. 15 also displays two solutions where the asymptotic behaviors of  $w_{\mu}^{+}$ ,  $\mu=\pm 1$ , and  $w_{-1}^{-}$  are simultaneously present (along different segments of the real axis) in the case of ( $\alpha=0$ ,  $\beta=-2$ ).

Analyzing Figs. 7 and 15, it can be seen that not only are some of the canonical sectors shown in Fig. 5 free from poles, but there are also regions sandwiched between two rows of poles that are

<span id="page-11-3"></span>![](_page_11_Figure_2.jpeg)

**Fig. 15.** Solution types (along the real axis (top) and pole locations and residues (bottom)) with adjacent pole-free sectors for  $\alpha=0$  and  $\beta=-2$ . In all frames u'(0)=0. The left and right frames both show that these solutions simultaneously match the roots  $w_{\mu}^+$ ,  $\mu=\pm 1$ , and  $w_{-1}^-$ . Poles of residue +1 are marked with dark (blue) circles and poles of residue -1 are marked with light (yellow) circles. Locations of zeros are marked by  $\times$ 's (red). (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

also pole free. The sequence shown in Fig. 7 also suggests that as a critical value of  $u(0) = u_0$  is approached with u'(0) fixed, the width of these sandwiched pole-free regions can be increased.

#### 7. Conclusions

This study of the fourth Painlevé equation started by numerically confirming various previous analytic and asymptotic results. A further exploration of the fundamental domain then identified solutions for general  $(\alpha, \beta)$ -values with noteworthy characteristics, such as numerous families of solutions with adjacent polefree sectors. Also, solutions with a nearly pole-free half-plane were found.

Most of the observations in this study were obtained numerically, leaving analytical considerations of some of the illustrated solution types an open topic. Further numerical tasks of interest include enlarging the explored  $(\alpha, \beta)$ -region, both to values of larger magnitude and (more importantly) to complex values for the parameters. Another extension, requiring only minimal changes to the numerical method, would be to drop the assumption that  $u(z) \in \mathbb{R}$  for  $z \in \mathbb{R}$ .

#### Acknowledgments

Discussions with Peter Clarkson, Harvey Segur, and André Weideman are gratefully acknowledged. J. A. R. was supported by

the Department of Defense. The views expressed in this article are those of the authors and do not reflect the official policy or position of the United States Air Force, Department of Defense, or U.S. Government. B. F. wishes to thank the National Science Foundation (NSF) for support under NSF Grant DMS-0914647.

#### References

- <span id="page-11-0"></span> A.S. Fokas, M.J. Ablowitz, On a unified approach to transformations and elementary solutions of the Painlevé equations, J. Math. Phys. 23 (1982) 2033–2042.
- [2] P.J. Olver, Applications of Lie Groups to Differential Equations, 2nd edn, in: Graduate Texts in Mathematics, vol. 107, Springer-Verlag, New York, 1993.
- [3] R. Wong, H. Zhang, On the connection formulas of the fourth Painlevé transcendent, Anal. Appl. (Singap.) 7 (2009) 419–448.
- [4] R. Wong, H. Zhang, On the connection formulas of the third Painlevé transcendent, Discrete Contin. Dyn. Syst. 23 (2009) 541–560.
- <span id="page-11-1"></span>[5] P.J. Forrester, N.S. Witte, Application of τ-function theory of Painlevé equations to random matrices: PIV PII and the GUE, Comm. Math. Phys. 219 (2001) 357–398
- [6] P.J. Forrester, N.S. Witte, Application of τ-function theory of Painlevé equations to random matrices: PIII, the LUE, JUE, and CUE, Comm. Pure Appl. Math. 55 (2002) 679–727.
- [7] C.A. Tracy, H. Widom, Level-spacing distributions and the Airy kernel, Comm. Math. Phys. 159 (1994) 151–174.
- <span id="page-11-2"></span>[8] A.S. Fokas, A.R. Its, A.V. Kitaev, Discrete Painlevé equations and their appearance in quantum gravity, Comm. Math. Phys. 142 (1991) 313–344.
- [9] G. Freud, On the coefficients in the recursion formulae of orthogonal polynomials, Proc. Roy. Irish Acad. Sect. A 76 (1976) 1–6.
- [10] E. Brézin, C. İtzykson, G. Parisi, J. Zuber, Planar diagrams, Comm. Math. Phys. 59 (1978) 35–51.

- [11] A.S[. Fokas, A.R. Its, X. Zhou, Continuous and discrete Painlevé equations,](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref11) in: D. Levi, P. Winternitz (Eds.), Painlevé transcendents: their asymptotics and physical applications, in: NATO Adv. Sci. Inst. Ser. B Phys., vol. 278, Plenum, New York, 1990, pp. 33–47. Proc. NATO Adv. Res. Workshop, Sainte-Adéle, Canada.
- [12] A.P[. Magnus, Painlevé-type differential equations for the recurrence coeffi](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref12)cients of semi-classical orthogonal polynomials, J. Comput. Appl. Math. (1995).
- <span id="page-12-0"></span>[13] B.[M. McCoy, Spin systems, statistical mechanics and Painlevé functions,](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref13) in: D. Levi, P. Winternitz (Eds.), Painlevé transcendents: their asymptotics and physical applications, in: NATO Adv. Sci. Inst. Ser. B Phys., vol. 278, Plenum, New York, 1992, pp. 377–391. Proc. NATO Adv. Res. Workshop, Sainte-Adéle, Canada, 1990.
- [14] B.[M. McCoy, C.A. Tracy, T.T. Wu, Painlevé functions of the third kind, J. Math.](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref14) Phys. (1977).
- [15] M. [Jimbo, T. Miwa, Y. Môri, M. Sato, Density matrix of an impenetrable Bose](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref15) gas and the fifth Painlevé transcendent, Physica D (1980).
- [16] F.H[.L. Essler, H. Frahm, A.R. Its, V.E. Korepin, Painlevé transcendent describes](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref16) quantum correlation function of the *XXZ* antiferromagnet away from the freefermion point, J. Phys. A (1996).
- [17] E. [Kanzieper, Replica field theories, Painlevé transcendents, and exact](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref17) correlation functions, Phys. Rev. Lett. 89 (2002) 250201-1–250201-4.
- [18] E. [Barouch, B.M. McCoy, T.T. Wu, Zero-field susceptibility of the two](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref18)dimensional Ising model near *t<sup>c</sup>* , Phys. Rev. Lett. 31 (1973) 1409–1411.
- <span id="page-12-1"></span>[19] T. [Bountis, H. Segur, F. Vivaldi, Integral Hamiltonian systems and the Painlevé](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref19) property, Phys. Rev. A 25 (3) (1982) 1257–1264.
- <span id="page-12-2"></span>[20] B. [Grammaticos, A.R. Ramani, V. Papageorgiou, Do integrable mappings have](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref20) the Painlevé property?, Phys. Rev. Lett. 67 (1991) 1825–1828.
- <span id="page-12-3"></span>[21] P. [Di Francesco, P. Ginsparg, J. Zinn-Justin, 2D gravity and random matrices,](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref21) Phys. Rep. 254 (1995) 1–133.
- [22] N. Seiberg, D. Shih, Flux vacua and branes of the minimal superstring, J. High Energy Phys. (2005) Electronic. E-print number: [hep-th/0412315.](http://arxiv.org/hep-th/0412315) 38 pp.
- <span id="page-12-18"></span>[23] D. [Bermúdez, D.J. Fernández, Supersymmetric quantum mechanics and](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref23) Painlevé IV equation, SIGMA 7 (2011).
- <span id="page-12-19"></span>[24] D. [Bermúdez, D.J. Fernández, Non-herminitiaon Hamiltonians and the Painlevé](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref24) IV equation with real parameters, Phys. Lett. A 375 (2011) 2974–2978.
- <span id="page-12-4"></span>[25] A.P[. Basson, P.A. Clarkson, A.C. Hicks, On the application of solutions of the](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref25) fourth Painlevé equation to various physically motivated nonlinear partial differential equations, Adv. Differential Equations 1 (1996) 175–198.
- <span id="page-12-5"></span>[26] F.[W.J. Olver, D.W. Lozier, R.F. Boisvert, C.W. Clark, NIST Handbook of](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref26) Mathematical Functions, Cambridge Univ. Press, 2010.
- <span id="page-12-6"></span>[27] B. [Fornberg, J.A.C. Weideman, A numerical methodology for the Painlevé](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref27) equations, J. Comput. Phys. 230 (2011) 5957–5973.
- <span id="page-12-7"></span>[28] I.M[. Willers, A new integration algorithm for ordinary differential equations](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref28) based on continued fraction approximations, Comm. ACM 17 (1974) 504.
- <span id="page-12-8"></span>[29] B. Fornberg, J. A.C. Weideman, A computational exploration of the second Painlevé equation, Found. Comput. Math. (2013) [http://dx.doi.org/10.1007/s10208-013-9156-x.](http://dx.doi.org/doi:10.1007/s10208-013-9156-x)
- <span id="page-12-9"></span>[30] J.A. [Reeger, B. Fornberg, Painlevé IV with both parameters zero: a numerical](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref30) study, Stud. Appl. Math. 130 (2013) 108–133.
- <span id="page-12-10"></span>[31] S. [Olver, A general framework for solving Riemann-Hilbert problems](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref31) numerically, Numer. Math. 122 (2012) 305–340.
- [32] A.S[. Fokas, A.R. Its, A.A. Kapaev, V.Y. Novokshenov, Painlevé Transcendents:](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref32) The Riemann-Hilbert Approach, AMS, Providence, 2006.

- [33] M. [Jimbo, T. Miwa, Monodromy preserving deformation of linear ordinary](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref33) differential euqations with rational coefficients. III, Physica D 4 (1981) 26–46.
- [34] A.V[. Kitaev, Self-similar solutions of the modified nonlinear Schrödinger](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref34) equation, Theoret. Math. Phys. 64 (1985) 878–894.
- [35] A.E[. Milne, P.A. Clarkson, A.P. Bassom, Application of the isomonodromy](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref35) deformation method to the fourth Painlevé equation, Inverse Problems 12 (1997) 421–439.
- <span id="page-12-11"></span>[36] A.A. Kapaev, Global asymptotics of the fourth Painlevé transcendent, PDMI Preprint-6/1996 (1996).
- <span id="page-12-12"></span>[37] A.A. Kapaev, Scaling limits in the fourth Painlevé transcendent, PDMI Preprint-15/1996, 1996.
- <span id="page-12-13"></span>[38] S. Olver, RHPackage, alpha version 0.4 (2011) [http://www.comlab.ox.ac.uk/](http://www.comlab.ox.ac.uk/people/Sheehan.Olver/projects/RHPackage.html) [people/Sheehan.Olver/projects/RHPackage.html.](http://www.comlab.ox.ac.uk/people/Sheehan.Olver/projects/RHPackage.html)
- <span id="page-12-14"></span>[39] V.I. [Gromak, I. Laine, S. Shimomura, Painlevé Differential Equations in the](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref39) Complex Plane, Walter de Gruyter, Berlin, 2002.
- [40] A.S[. Fokas, U. Mugan, M.J. Ablowitz, A method of linearization for Painlevé](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref40) equations: Painlevé IV V, Physica D 30 (1988) 247–283.
- [41] N.[A. Lukashevich, Theory of the fourth Painlevé equation, Differ. Equ. 3 \(1967\)](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref41) 395–399.
- [42] Y. [Murata, Rational solutions of the second and fourth Painlevé equations,](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref42) Funkcial. Ekvac. 28 (1985) 1–32.
- <span id="page-12-15"></span>[43] J.A. Reeger, A computational study of the fourth Painlevé equation and a discussion of adams predictor–corrector methods, Ph.D. Thesis, University of Colorado, Boulder, Colorado, 2013.
- <span id="page-12-16"></span>[44] P.A[. Clarkson, Special polynomials associated with rational solutions of the](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref44) defocusing nonlinear Schrödinger equation and the fourth Painlevé equation, European J. Appl. Math. 17 (2006) 293–322.
- <span id="page-12-17"></span>[45] K. [Okamoto, Studies on the Painlevé equations III. Second and fourth Painlevé](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref45) equations, PII and PIV, Math. Ann. 275 (1986) 221–255.
- <span id="page-12-20"></span>[46] C.[M. Bender, S.A. Orszag, Advanced Methods for Scientists and Engineers,](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref46) McGraw-Hill, New York, 1978.
- <span id="page-12-21"></span>[47] P.A[. Clarkson, The fourth Painlevé equation and associated special polynomi](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref47)als, J. Math. Phys. 44 (2003).
- <span id="page-12-22"></span>[48] G.F[. Corliss, Integrating ODE's in the complex plane-Pole vaulting, Math. Comp.](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref48) 35 (1980) 1181–1189.
- [49] A.P[. Bassom, P.A. Clarkson, A.C. Hicks, Numerical studies of the fourth Painlevé](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref49) equation, IMA J. Appl. Math. 50 (1993) 167–193.
- [50] A.A[. Abramov, L.F. Yukhno, A method for the numerical solution of the Painlevé](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref50) equations, Comp. Math. Math. Phys. 53 (2013) 702–726.
- <span id="page-12-23"></span>[51] V.Y[. Novokshenov, Padé approximations for Painlevé I and II transcendants,](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref51) Theoret. Math. Phys. 159 (2009) 853–862.
- <span id="page-12-24"></span>[52] L. Peltonen, Numerical Solution of ODEs with Poles, Master's thesis, Worcester College, Walton Street, Oxford, Oxfordshire, OX1 2HB, 01865 278300, University of Oxford, 2011.
- <span id="page-12-25"></span>[53] D. Bermúdez, D.J. Fernández, Solution hierarchies for the Painlevé IV equation, in: Proceedings of the XXXI Workshop on Geometric Methods in Physics, Bialowieza, Poland.
- <span id="page-12-26"></span>[54] A.R[. Its, A.A. Kapaev, Connection formulae for the fourth Painlevé tran](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref54)scendent: Clarkson-McLeod solution, J. Phys. A: Math. Gen. 31 (1998) 4073–4113.
- <span id="page-12-27"></span>[55] A.P[. Bassom, P.A. Clarkson, A.C. Hicks, J.B. McLeod, Integral equations and exact](http://refhub.elsevier.com/S0167-2789(14)00083-9/sbref55) solutions of the fourth Painlevé equation, Proc. R. Soc. A 437 (1992) 1–24.