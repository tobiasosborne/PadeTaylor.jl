![](_page_0_Picture_1.jpeg)

Contents lists available at [ScienceDirect](http://www.elsevier.com/locate/physd)

### Physica D

journal homepage: [www.elsevier.com/locate/physd](http://www.elsevier.com/locate/physd)

![](_page_0_Picture_5.jpeg)

## A computational overview of the solution space of the imaginary Painlevé II equation

![](_page_0_Picture_7.jpeg)

Bengt Fornberg [a](#page-0-0) , J.A.C. Weideman[b,](#page-0-1)[∗](#page-0-2)

- <span id="page-0-0"></span><sup>a</sup> *Department of Applied Mathematics, University of Colorado, Boulder, CO 80309, USA*
- <span id="page-0-1"></span><sup>b</sup> *Department of Mathematical Sciences, Stellenbosch University, Stellenbosch 7600, South Africa*

#### h i g h l i g h t s

- The full 3-parameter solution space to the Imaginary Painlevé II equation is surveyed.
- Numerous illustrations show the equation's pole field dynamics.
- Full agreement with known asymptotic results is demonstrated.
- Novel pole dynamics that can encourage new analytical investigations are presented.

#### g r a p h i c a l a b s t r a c t

![](_page_0_Figure_17.jpeg)

### a r t i c l e i n f o

*Article history:* Received 17 February 2015 Received in revised form 23 May 2015 Accepted 23 July 2015 Available online 30 July 2015 Communicated by P.D. Miller

*Keywords:* Imaginary Painlevé II Modified Painlevé II Ablowitz–Segur solutions Tronquée solutions Pole field solver

#### a b s t r a c t

The six Painlevé equations were first formulated about a century ago. Since the 1970s, it has become increasingly recognized that they play a fundamental role in a wide range of physical applications. A recently developed numerical pole field solver (Fornberg and Weideman, 2011) now allows their complete solutions spaces to be surveyed across the complex plane. Following such surveys of the *P<sup>I</sup>* , *PII* and *PIV* equations, we consider here the case of the imaginary *PII* equation (the standard *PII* equation, with a change of sign for its nonlinear term). Solutions to this equation share many features with other classes of Painlevé transcendents, including a rich variety of pole field configurations, with connection formulas linking asymptotic behaviors in different directions of the complex plane.

© 2015 Elsevier B.V. All rights reserved.

### **1. Introduction**

Among the six Painlevé equations, the second one has a particularly large number of applications, as surveyed in [\[1](#page-10-0)[,2\]](#page-10-1). It is defined by

$$y'' = 2y^3 + zy + \beta, \tag{1}$$

*E-mail addresses:* [fornberg@colorado.edu](mailto:fornberg@colorado.edu) (B. Fornberg), [weideman@sun.ac.za](mailto:weideman@sun.ac.za) (J.A.C. Weideman).

<span id="page-0-2"></span><sup>∗</sup> Corresponding author.

where the primes denote differentiation with respect to the independent variable, z, and  $\beta$  is a constant. Attention has also been given in the literature to the special case when y(z) is purely imaginary along the real z-axis [1,3-5], with an application to nonlinear optics discussed in [6]. With the change of variables

$$y(z) = i u(z)$$
 and  $\beta = i \alpha$ , (2)

one arrives at the modified or imaginary  $P_{II}$  equation (Im  $P_{II}$  for short)

$$u'' = -2u^3 + zu + \alpha, \quad z = x + iy,$$
 (3)

with  $\alpha$  real and u(z) real valued along the real axis. (These conditions on  $\alpha$  and u(z) will be assumed throughout the rest of the paper.) Since the two equations only differ by a trivial variable change,  $\operatorname{Im} P_{II}$  will share many features of the regular  $P_{II}$  equation (which we denote by  $\operatorname{Re} P_{II}$  when assuming real solutions for z real). These include having first order poles as the only possible singularities. As noted in Section 2.1,  $\operatorname{Im} P_{II}$  differs in other ways, such as an absence of Bäcklund-type transformations and, with that, of non-trivial closed form solutions. Along the real axis, solutions are singularity free, but will always be oscillatory far out in at least one of the two directions.

Like for the regular  $P_{II}$  case, solutions are usually dominated by dense pole fields across the full complex plane. However, tronquée solutions (featuring pole free sectors extending to infinity) are also for Im  $P_{II}$  of particular interest, and will be highlighted in this study. The primary goal however is to provide the first full survey of the complete three parameter solution space to (3), taking as parameters  $\{\alpha, u(0), u'(0)\}$ .

The present study is made practical by the pole field solver, first described in [7] in connection with  $P_I$ , and since used to explore the solution spaces of Re  $P_{II}$  [8] and Re  $P_{IV}$  [9,10]. It is based on a numerical ODE initial value solver [11] which, by using Padé (or continued fraction) approximations in each step, becomes 'immune' to the presence of poles even in the immediate vicinity of its steps. By bringing up the order of accuracy to the 30-60 range, the computation becomes very efficient. In MATLAB, using a standard notebook computer, execution times of around a second are typical for each of the pole field plots shown later (all based on  $161 \times 161$  grid calculations, as seen in Fig. 3). Another key numerical component is the use of a quasi-randomized branching network of integration paths, to effectively cover full regions of the complex plane (rather than following single paths). It should also be noted that, when using this 'pole-friendly' ODE stepping process, Painlevé type initial value problems are numerically well conditioned within pole fields, while they are (mathematically) far more sensitive to infinitesimal perturbations in smooth solution regions. If very large smooth regions are present, it may be preferable to keep the integration paths mostly within the pole fields, and then treat smooth regions as boundary value problems [7].

It is straightforward to calculate solutions to  $Im P_{II}$  along the real axis with any standard ODE solver, since such solutions are singularity free. The only previously published instances of numerical solutions in a complex vicinity of the origin for  $Im P_{II}$  appear to be two illustrations in [12] (also published elsewhere by the same author), and two in [13]. All of these were obtained by using a single very high order Padé expansion in the special case

#### 2. Overview of the main features of solutions to $\operatorname{Im} P_{II}$

#### <span id="page-1-0"></span>2.1. General properties

Since (3) is invariant under the variable changes  $u(z) \to -u(z)$ ,  $\alpha \to -\alpha$ , the discussion here is limited to  $\alpha \geq 0$ . Analogously to the Re  $P_{II}$  equation having as its only possible singularities first

order poles with residues +1 or -1, it follows from the variable change (2) (or from direct substitution of a Laurent expansion into (3)) that all poles for  $\operatorname{Im} P_{II}$  also are first order, but instead with residues +i or -i. An immediate consequence is that no poles can be located on the real axis, since this would conflict with u(z) being real for z real.

<span id="page-1-2"></span>This lack of poles on the real axis influences how we characterize the solution space in this study. While pole counts along the positive and negative real axes played a key role in the Re  $P_{II}$  survey [8], asymptotic solution properties for  $x \to \pm \infty$  will achieve the same goal here. (By the notation  $x \to \pm \infty$  it is understood that y = 0, i.e., the limit is considered along the real z-axis.)

<span id="page-1-1"></span>Recursive use of the Bäcklund transformation for  $P_{II}$  [2, Sect. 32.7]

$$y(\beta + 1; z) = -y(\beta; z) - \frac{2\beta + 1}{2v'(\beta; z) + 2v(\beta; z)^2 + z}$$
(4)

<span id="page-1-3"></span>provides classes of closed form solutions to Re  $P_{II}$  when starting from the trivial solution y(0;z)=0 or from the one-parameter Airy solution family for  $y(\frac{1}{2};z)$ . For Im  $P_{II}$ , u(0;z)=0 is again a solution, but the counterpart to (4) will then only produce solutions violating u(z) being real for z real. In fact, no closed form solution other than u(0;z)=0 is known in the class of solutions considered here (u real along the real axis).

#### 2.2. Leading asymptotic behaviors for $x \to \pm \infty$

The only possibilities for smooth (non-oscillatory) asymptotic behavior of solutions to (3) as  $x \to \pm \infty$  can be found by setting u'' to zero, i.e. by considering

<span id="page-1-5"></span><span id="page-1-4"></span>
$$-2u^3 + xu + \alpha = 0. \tag{5}$$

Asymptotic balances then give

$$u \to 0$$
 as  $x \to -\infty$ , and  $u \to 0$ ,  $u \sim \pm \sqrt{\frac{1}{2}x}$ , as  $x \to +\infty$ . (6)

The corresponding full asymptotic expansions will be discussed in Section 3.1.1  $(x \to -\infty)$  and Section 3.1.2  $(x \to +\infty)$ . The various behaviors determined by (6) are represented schematically in Fig. 1, with the different cases indicated by the letters (A)–(D). For each of the cases (A), (B), (C), there exists a family of solutions that converge with oscillations of decreasing amplitude to the indicated asymptotic limit. For each of these three families, our numerical investigations showed that there further exists a unique set of initial conditions  $\{u(0), u'(0)\}$  that yields a solution with non-oscillatory convergence. (We are unaware of any theoretical results related to this observation.) The absence of oscillations will represent a pole free sector surrounding the real axis in the corresponding direction, known as a tronquée (or truncated) solution. The case (D) represents the family of solutions that are bounded across the entire real line. Known as the Ablowitz-Segur (AS) solutions [14], their convergence to the indicated limit is always non-oscillatory.

#### 2.3. An overview of some solution features

For a given  $\alpha$ , there remains two real parameters in the solution space, which we take as u(0) and u'(0). As a first step towards characterizing such a solution space, we note that, as  $x \to \pm \infty$ , no other possibilities exist than those sketched out in Fig. 1. For  $x \to -\infty$ , all solutions converge to 0. For  $x \to +\infty$ , convergence will occur to one of the three limits 0,  $\pm \sqrt{x/2}$ . The latter three limits are displayed in the diagrams shown in Fig. 2, with a shaded region

<span id="page-2-0"></span>![](_page_2_Figure_2.jpeg)

<span id="page-2-1"></span>**Fig. 1.** The diagrams show the curves  $-2u^3 + xu + \alpha = 0$  (cf. (5)) that represent the various asymptotic behaviors of (3) on the real axis, for the special value  $\alpha = 0$  and for a generic  $\alpha > 0$ , respectively.

![](_page_2_Figure_4.jpeg)

Fig. 2. Character of solutions u(x) for different choices of  $\{u(0), u'(0)\}$  in the cases  $\alpha=0, \alpha=\frac{1}{2}$  and  $\alpha=4$  (with  $\alpha=\frac{1}{2}$  and  $\alpha=4$  representative of intermediate and large values of  $\alpha$ , respectively; cf. additional illustrations in Fig. 12). The regions with convergence to  $\pm\sqrt{x/2}$  as  $x\to+\infty$  are shown as blank and shaded, respectively. Numerical values of u(0) and u'(0) for the marked points are given in Table 1. The points represent: (a) Unique case of non-oscillatory convergence to 0 for  $x\to-\infty$ , (b, c) Unique cases of non-oscillatory convergence to  $\pm\sqrt{x/2}$ , respectively, for  $x\to+\infty$ , (d) Example of solution along the Ablowitz–Segur curve, (e, f) Examples of solutions with oscillatory convergence to  $\pm\sqrt{x/2}$ , respectively, for  $x\to+\infty$ . In the  $\alpha=4$  case, the tip of the shaded region has retreated about  $1\frac{3}{4}$  turns from its position for  $\alpha=0$ .

<span id="page-2-2"></span>**Table 1**Numerical values of the points marked (a)–(f) in Fig. 2.

|     | $\alpha = 0$ |             | $\alpha = \frac{1}{2}$ |             | $\alpha = 4$ |            |
|-----|--------------|-------------|------------------------|-------------|--------------|------------|
|     | u(0)         | u'(0)       | <u>u(0)</u>            | u'(0)       | <u>u(0)</u>  | u'(0)      |
| (a) | 0            | 0           | 0.54353821             | 0.29696044  | 1.25994238   | 0.13258479 |
| (b) | 0.39507520   | 0.40522552  | 0.64216297             | 0.26816749  | 1.25994182   | 0.13258178 |
| (c) | -0.39507520  | -0.40522552 | 0.13395730             | -0.86583605 | -0.19549371  | 3.41118861 |
| (d) |              |             | 0                      | 1.56004311  |              |            |
| (e) |              |             | -0.6                   | 0           |              |            |
| (f) |              |             | -1.4                   | 0           |              |            |

representing the limit  $-\sqrt{x/2}$  and a blank region representing the limit  $+\sqrt{x/2}$ . The curve separating these regions represents the limit 0, constituting the AS family of solutions.

Some special points are marked in Fig. 2. The unique initial condition  $\{u(0), u'(0)\}$  that generates a non-oscillatory solution converging to 0 as  $x \to -\infty$  is marked by (a). In the case  $\alpha = 0$ , this solution is the trivial u = 0. The points (b) and (c) correspond to the unique cases of non-oscillatory convergence to the limits  $\pm \sqrt{x/2}$  as  $x \to +\infty$ , respectively. The  $\alpha = 0$  solution corresponding to point (b) is shown in Fig. 3. One observes a wide pole free region for  $|\arg(z)| \le 2\pi/3$ , and a dense pole field to the left of this. Note that point (c) corresponds to a sign reversal of the solution

represented by point (b). This symmetry is lost as soon as  $\alpha$  is increased from 0, however; cf. Fig. 2.

For the  $\alpha=\frac{1}{2}$  case, Fig. 4 shows the solutions corresponding to all the points (a)–(f). The cases (e) and (f) represent the 'typical' solution behavior—featuring pole fields extending in all directions of the complex plane. By contrast, each one of the cases (a)–(d) represents a solution with large pole free sectors.

Examining Figs. 2 and 4 and the data in Table 1, one notices that the points marked (a) and (b) approach each other rapidly in the  $\{u(0), u'(0)\}$ -plane when  $\alpha$  increases. (An empirical fit suggests the distance between the two points decreases exponentially in  $\alpha$ .) Yet, they represent completely different pole configurations in the

<span id="page-3-0"></span>![](_page_3_Figure_2.jpeg)

**Fig. 3.** The solution of  $\operatorname{Im} P_{II}$  that exhibits non-oscillatory convergence to  $+\sqrt{x/2}$  as  $x \to +\infty$  when  $\alpha = 0$ . Initial conditions correspond to the point marked (b) in the left diagram of Fig. 2. The real part of u is shown, with its graph along the real axis shown as the thicker curve. (Because the imaginary part is similar it is not displayed here.) In the first subplot of Fig. 8 two alternate views of the same solution are displayed, namely its pole locations in the complex plane and its profile on the real axis.

<span id="page-3-1"></span>![](_page_3_Figure_4.jpeg)

**Fig. 4.** Solutions of  $\operatorname{Im} P_{II}$  for the case  $\alpha = \frac{1}{2}$ . The label on each subplot corresponds to the cases (a)–(f) as marked and described in Fig. 2. The diagrams show the pole locations at the top, with the dark (blue) and light (yellow) dots representing residues +i and -i, respectively (a convention followed in all figures shown below). Immediately below each pole diagram the solution profile on the real axis is displayed. (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

<span id="page-4-1"></span>![](_page_4_Figure_2.jpeg)

**Fig. 5.** Pole fields and solutions along the real axis for the case  $\alpha = \frac{1}{2}$ , with initial conditions interpolated between cases (a) and (b) of Figs. 2 and 4. The label on each pole field refers to the value of the interpolation parameter t in (7). As the initial conditions vary between cases (a) and (b), the pole field switches completely from the right half-plane to the left, with some complicated interaction patterns in between. Comparing the pole fields of the first and last subplots, we note that they are not simply left-right reflections of each other. The residue patterns are entirely different, in a way that allows the solution on the real axis to remain almost the same.

sense that the poles of cases (a) and (b) are restricted to sectors in the right and left half-planes, respectively. This is therefore quite an unstable situation as far as pole locations are concerned.

To examine what happens when the data  $\{u(0),u'(0)\}$  vary between the cases (a) and (b), we computed for  $\alpha=\frac{1}{2}$  a set of pole configurations with initial conditions linearly interpolated between the two cases. That is,

$$u(0) = u_{a}(0) + t(u_{b}(0) - u_{a}(0)),$$
  

$$u'(0) = u'_{a}(0) + t(u'_{b}(0) - u'_{a}(0)),$$
(7)

where  $\{u_a(0), u'_a(0)\}$  and  $\{u_b(0), u'_b(0)\}$  represent the  $\{u(0), u'(0)\}$  values of cases (a) and (b) given in Table 1. The results are shown in Fig. 5, for six different values of t, with t=0 representing case (a) as already shown in Fig. 4 and t=1 representing case (b).

The trends just noted become even more pronounced for still larger  $\alpha$ , as will be considered later in Section 3.2. The case  $\alpha=4$  is shown in the right diagram of Fig. 2. It shows the points (a) and (b) as indistinguishable on the scale of that figure, as also confirmed by the data in Table 1.

#### 2.4. Oscillations along the real axis

In Fig. 2 the blank and shaded regions represent oscillatory convergence to the asymptotic  $u \sim \pm \sqrt{x/2}$ , respectively, as  $x \to +\infty$ . To leading order [15]

$$u(x) \mp \sqrt{x/2} \sim x^{-1/4} d_+ \widetilde{u}(x), \quad x \to +\infty,$$
 (8)

<span id="page-4-0"></span>where  $\widetilde{u}(x)$  is oscillatory, with unit amplitude. The number  $d_+$ , which is a function of  $\alpha$  and the initial conditions, has a non-zero value everywhere in the  $\{u(0), u'(0)\}$  plane apart from at the points denoted by (b) and (c) above. The oscillatory convergence on the negative real axis satisfies a similar expression

<span id="page-4-2"></span>
$$u(x) \sim (-x)^{-1/4} d_{-}\widetilde{u}(x), \quad x \to -\infty,$$
 (9)

with  $d_-$  non-zero apart from at the point marked (a). Fig. 6 shows these amplitudes  $d_-$  in the  $\{u(0), u'(0)\}$ -plane for a few values of  $\alpha$ . The quantity  $d_-$  will play a key role in the discussion of the AS solutions in Section 3.2.

#### 3. Tronquée solutions

As illustrated already in Fig. 2, tronquée solutions fall into two classes: Isolated cases, in the figure represented by the dots marked

<span id="page-5-2"></span>![](_page_5_Figure_2.jpeg)

<span id="page-5-3"></span>**Fig. 6.** Amplitude  $d_-$  of the asymptotic oscillation as  $x \to -\infty$  (cf. (9)), displayed over the  $\{u(0), u'(0)\}$ -plane. The white dots indicate where the amplitude vanishes, corresponding to the non-oscillatory unique cases labeled (a) in Fig. 2.

![](_page_5_Figure_4.jpeg)

Fig. 7. Tronquée solutions associated with  $u \to 0$ ,  $x \to -\infty$ . Initial conditions correspond to the points marked (a) in Fig. 2. (Note that the middle diagram also features in Fig. 4.) The  $\alpha=0$  member of this family is the u=0 solution, which is pole free. As  $\alpha$  is increased from 0 a group of poles enters from  $+\infty$ , with a widening of the pole free gap around the real axis. Correspondingly, the solution u(x) on the indicated real interval become less oscillatory with increasing  $\alpha$ .

(a), (b), (c), and the one-parameter AS family, represented by the boundary curve between the blank and shaded regions. We next describe these two classes in more detail.

#### 3.1. Isolated cases

#### <span id="page-5-0"></span>3.1.1. Non-oscillatory convergence to 0 for $x \to -\infty$

The solutions that correspond to points marked (a) in Fig. 2 have the full asymptotic behavior

$$u(x) \sim -\frac{\alpha}{x} \sum_{n=0}^{\infty} \frac{a_n}{x^{3n}}, \quad x \to -\infty,$$

with coefficients  $a_0 = 1$  and

$$a_1 = 2(1 + \alpha^2),$$
  $a_2 = 4(3\alpha^2 + 10)(1 + \alpha^2),$   
 $a_3 = 8(1 + \alpha^2)(12\alpha^4 + 117\alpha^2 + 280),$  etc. (10)

To compute the corresponding pole fields, this expression was used to compute accurate solutions at some value of x=-L, at which point the pole field solver was initialized. (A suitable value of L was determined by experimentation, with values in the range 7–10 usually satisfactory.) The computed value of  $\{u(0), u'(0)\}$  is given in the top row of Table 1 and shown as the point (a) in Fig. 2  $(\alpha = \frac{1}{2})$ .

The corresponding solutions for a range of  $\alpha$ -values can be seen in Fig. 7, featuring wide pole free sectors surrounding  $\mathbb{R}^-$ . For increasing  $\alpha$ , the poles of this family of solutions enter from the right, with a corresponding widening of the pole free gap around  $\mathbb{R}^+$ 

#### <span id="page-5-1"></span>3.1.2. Non-oscillatory convergence to $\pm \sqrt{x/2}$ for $x \to +\infty$

The solutions that correspond to the unique points marked (b,c) in Fig. 2 have the full asymptotic behavior

$$u(x) \sim \pm \sqrt{\frac{x}{2}} \sum_{n=0}^{\infty} \frac{(\pm 1)^n b_n}{x^{3n/2}}, \quad x \to +\infty.$$

The plus and minus sign choices correspond to the cases (b) and (c), respectively. The coefficients are  $b_0 = 1$  and

<span id="page-5-4"></span>
$$b_1 = \frac{\alpha}{\sqrt{2}}, \quad b_2 = \frac{1 - 6\alpha^2}{8}, \quad b_3 = \frac{\alpha(16\alpha^2 - 11)}{8\sqrt{2}}, \text{ etc.}$$

The corresponding pole fields were computed by the same strategy used above, namely the pole field solver was initialized at a point x = L, with  $\{u(L), u'(L)\}$  values provided by the asymptotic expansion. The results are shown in Figs. 8 (plus sign) and 9 (minus sign)

With either sign choice there are wide pole free sectors that surround  $\mathbb{R}^+$ . The plus sign generates a relatively large pole free

<span id="page-6-0"></span>![](_page_6_Figure_2.jpeg)

**Fig. 8.** Tronquée solutions associated with  $u \sim \sqrt{x/2}$ ,  $x \to +\infty$ . Initial conditions correspond to the points marked (b) in Fig. 2. (Note that the middle diagram also features in Fig. 4.) As  $\alpha$  is increased from 0 the group of poles in the left half-plane move slowly to  $-\infty$ , with a widening of the pole free gap around the real axis. Correspondingly, the solution u(x) on the indicated real interval becomes less oscillatory with increasing  $\alpha$ .

<span id="page-6-2"></span>![](_page_6_Figure_4.jpeg)

**Fig. 9.** Tronquée solutions associated with  $u \sim -\sqrt{x/2}$ ,  $x \to +\infty$ . Initial conditions correspond to the points marked (c) in Fig. 2. (Note that the middle diagram also features in Fig. 4.) As  $\alpha$  is increased the pole dynamics is in a sense the opposite of that in Fig. 8: the group of poles in the left half-plane move slowly to  $+\infty$ , with a narrowing of the pole free gap around the real axis. Correspondingly, the solution u(x) on the indicated real interval becomes more oscillatory with increasing  $\alpha$ .

region around  $\mathbb{R}^-$  as well, particularly for large  $\alpha$ , which results in relatively smooth solutions across  $\mathbb{R}$ . The minus sign, by contrast, produces solutions with poles just off the negative axis, which results in oscillations with large amplitude on  $\mathbb{R}^-$ .

Recalling that  $\operatorname{Im} P_{II}$  is invariant when changing the signs of both  $\alpha$  and u(z), we can view the pole field sequence in Fig. 9 as a direct continuation of the one in Fig. 8. Both would then represent  $u \sim \sqrt{x/2}$ , but with  $\alpha$  continuing from positive to negative values.

Along the real axis, the last case shown in Figs. 7 and 8 look very similar, and superficially like the everywhere non-oscillatory Hastings–McLeod solutions to Re  $P_{II}$ . The similarity between the two shown curves is not surprising, given how the points marked (a) and (b) in Fig. 2 approach each other for increasing  $\alpha$ . The difference to the Hastings–McLeod case is profound, however. The present Im  $P_{II}$  solutions eventually become oscillatory sufficiently far out in one of the two directions. The contrast in the pole fields of the solutions (a) and (b) reflects this strikingly. (This provides a good example of the insight that can be gained by considering

pole fields even in cases when the primary motivation comes from solution properties along the real axis. For a related example, see [8, Fig. 10].)

#### <span id="page-6-1"></span>3.2. Ablowitz-Segur (AS) solutions

All solutions of (3) are bounded on  $\mathbb{R}^-$ . By contrast, all solutions are unbounded on  $\mathbb{R}^+$ , with the exception of those that satisfy  $u \to 0$  as  $x \to +\infty$ . These AS solutions can therefore be characterized by the fact that they are bounded across all of  $\mathbb{R}$ . In the space of initial conditions, these solutions are represented by the curves that separate the blank and shaded regions in Fig. 2.

#### 3.2.1. Parametrization of the AS curves

The AS solutions are parameterized by a single parameter, say k, and can be expressed in the form

<span id="page-6-3"></span>
$$u(z) = B(\alpha; z) + e(\alpha, k; z), \tag{11}$$

<span id="page-7-2"></span>![](_page_7_Figure_2.jpeg)

**Fig. 10.** A smaller section of the contour plot of the asymptotic oscillation  $d_-$  shown in Fig. 6, with superimposed on it the Ablowitz–Segur curve from Fig. 2. The white dot marked (a) is where  $d_-$  vanishes, and the black dot represents the solution that corresponds to k=0 in (13). The diagram illustrates that the k=0 point is the unique point on the Ablowitz–Segur curve (i.e.,  $u\to 0$  as  $x\to +\infty$ ) that minimizes the amplitude  $d_-$  of the oscillation on  $\mathbb{R}^-$ .

where

$$B(\alpha;z) \sim -\frac{\alpha}{z} \sum_{n=0}^{\infty} \frac{a_n}{z^{3n}}, \quad z \to +\infty,$$
 (12)

with the coefficients  $a_n$  defined in (10). Further,

$$e(\alpha, k; z) \sim \frac{k}{2\sqrt{\pi}} \frac{e^{-(2/3)z^{3/2}}}{z^{1/4}} \sum_{n=0}^{\infty} \frac{c_n}{z^{3n/2}} + \left(\frac{k}{2\sqrt{\pi}}\right)^2 \frac{e^{-(4/3)z^{3/2}}}{z^{4/4}} \sum_{n=0}^{\infty} \frac{d_n}{z^{3n/2}} + O\left(k^3 \frac{e^{-(6/3)z^{3/2}}}{z^{7/4}}\right),$$
(13)

with coefficients again readily obtained recursively by substituting into (3) and equating coefficients (or by changing  $\alpha \rightarrow i\alpha$  in the corresponding formulas for Re  $P_{II}$ ). The rest of this section

is focused on showing how the AS solutions change when the parameter k ranges over  $-\infty < k < \infty$ .

Like in the Re  $P_{II}$  case, direct evaluation of (12) poses a challenge for  $\alpha>0$  and z>0, since the sum is so rapidly divergent that neither optimal truncation, nor convergence acceleration appear able to overcome this. Quite apart from this growth in its coefficients, the quantity  $e(\alpha,k;z)$  is asymptotically smaller than every term in (12), indicating that evaluation of this expansion alone might be insufficient for associating k-values to different AS solutions. Note however the further discussion on this issue in the Appendix.

The sums in (13) are also divergent, but less so, and are possible to sum sufficiently accurately for  $z \ge L$ . As first noted in [8], connection formulas [15] then allow  $B(\alpha;z)$  to be obtained via (11) rather than via (12). In the Re  $P_{II}$  case, there is for each  $\alpha$  a unique Hastings–McLeod solution that is readily computed by a standard ODE boundary value solver, and which in [15, Sect. 4] is shown to correspond to  $k = \cos \pi \alpha$ . With then both u(z) and  $e(\alpha, k; z)$  in (11) available at some z = L, the two quantities  $\left. \left\{ B(\alpha; z), \frac{d}{dz} B(\alpha; z) \right\} \right|_{z=L}$  can be calculated. With these as initial conditions, the pole field solver will generate the Re  $P_{II}$  solution  $B(\alpha; z)$ . These solutions were referred to as k = 0 solutions in [8], and again we will do so here.

<span id="page-7-1"></span>In the  $\operatorname{Im} P_{II}$  case, there is no counterpart to the Hastings–McLeod solution, and in fact no specific non-zero k value that features any similar readily distinguishable solution feature. An examination of the above quoted connection formula, however, additionally reveals that k=0 corresponds to a unique local minimum of the oscillation amplitude  $d_-$ , as  $z\to-\infty$ . To demonstrate this in a  $\{u(0),u'(0)\}$ -display, a smaller section of the middle case in Fig. 6 is reproduced in Fig. 10, with the AS curve superimposed. The minimum point along this curve is the point at which the AS curve is parallel to a contour line of  $d_-$ , or equivalently, where it is orthogonal to the gradient of  $d_-$ .

<span id="page-7-0"></span>With the k=0 point thus identified, we can compute corresponding  $\{u(0), u'(0)\}$ -values by a univariate optimization algorithm that minimizes  $d_-$  along the AS curve. Using this strategy, the  $B(\alpha; z)$  function becomes again available. Fig. 11 shows these k=0 solutions for some  $\alpha$ -values. Note the relatively wide pole free regions in the neighborhood of  $\mathbb{R}^-$ , which accounts for the oscillations of minimum amplitude as  $x \to -\infty$ .

With the  $B(\alpha; z)$  values available, one can now also identify  $k \neq 0$  solutions via (13). We summarize this in the  $\{u(0), u'(0)\}$ -planes shown in Fig. 12 by placing markers for different k-values

<span id="page-7-3"></span>![](_page_7_Figure_15.jpeg)

**Fig. 11.** Ablowitz–Segur solutions corresponding to the special k=0 case in the asymptotic formula (13), for various values of  $\alpha$ .

<span id="page-8-0"></span>![](_page_8_Figure_2.jpeg)

**Fig. 12.** Diagrams that represent asymptotic behavior, similar to the ones shown in Fig. 2. That is, the blank and shaded regions represent convergence to  $\pm \sqrt{x/2}$ , respectively, as  $x \to +\infty$ . The black dots correspond to values of k in the asymptotic formula (13), with the bigger white circle representing the special k=0 case. Solutions that lie in the segment between the black crosses in the middle diagram will be shown in Fig. 15.

<span id="page-8-2"></span>![](_page_8_Figure_4.jpeg)

Fig. 13. Examples of Ablowitz–Segur solutions in the case  $\alpha = \frac{1}{2}$ , for various choices of k in the asymptotic expression (13). The choice k = 0 minimizes the amplitude of the oscillations on  $\mathbb{R}^-$ , and also pushes the onset of the oscillatory motion as far to the left as possible. Corresponding pole diagrams are shown in Fig. 14.

<span id="page-8-1"></span>![](_page_8_Figure_6.jpeg)

**Fig. 14.** Pole diagrams of some of the Ablowitz–Segur solutions shown in Fig. 13. The choice k = 0 in (13) results in a pole configuration with the widest possible pole free neighborhood of the real axis.

along the AS curves. Some of these  $k \neq 0$  solutions are shown in comparison with the k = 0 solution along the real axis in Fig. 13, and as pole fields in the complex plane in Fig. 14. Again the amplitude minimizing property on  $\mathbb{R}^-$  of the k = 0 solution is evident.

In summary, for  $\operatorname{Im} P_{II}$  the AS family of solutions is present for all  $\alpha$ , and extends for  $|k| < \infty$ . This can be contrasted to the Re  $P_{II}$  case, where it is limited to  $|\alpha| \leq \frac{1}{2}$ , and then further to  $|k| \leq \cos \pi \alpha$ .

To conclude this section we illustrate in Fig. 15 the pole dynamics as initial conditions are varied between two AS solutions. We chose the two sets of initial conditions represented by the black crosses on the positive u(0) axis in the middle subplot of

Fig. 12, and display the solutions as u(0) is varied (in nonuniform increments) between the two points, keeping u'(0) = 0. The pole dynamics is described in the figure caption.

#### 4. Conclusions

In the absence of any explicit solutions, previous analytical studies of  $\operatorname{Im} P_{II}$  centered on asymptotic results in the  $|z| \to \infty$  limit, notably in the form of connection formulas. These formulas link behaviors in different directions far out in the complex plane. Numerical studies, on the other hand, centered on solutions along the real axis. These are relatively straightforward to compute, because of the absence of poles on the real axis.

<span id="page-9-1"></span>![](_page_9_Figure_2.jpeg)

**Fig. 15.** Pole dynamics as initial conditions vary between two Ablowitz–Segur solutions, in the case  $\alpha = \frac{1}{2}$ . In each subplot u'(0) = 0 and u(0) varies (in nonuniform increments) in the interval between the two crosses in the middle diagram of Fig. 12. Note how a pole field enters from  $+\infty$ , meets up with the field to the left, and then leaves behind one curve of poles as it returns to  $+\infty$ .

With the recently developed pole field solver [7] as our primary tool, we can now also illustrate the pole dynamics that arise as the three parameters in the solution space of  $\operatorname{Im} P_{II}$  are varied. While many of these features are reminiscent of those of the  $\operatorname{Re} P_{II}$  equation, the  $\operatorname{Im} P_{II}$  case is sufficiently distinct to justify separate attention. Notable differences to the  $\operatorname{Re} P_{II}$  case include: (i) all solutions of  $\operatorname{Im} P_{II}$  are pole free along the real axis, (ii) poles in the complex plane are of first order, with residue either +i or -i, (iii) the solutions do not change qualitative features when  $\alpha$  takes certain values (in contrast to the special cases arising for  $\operatorname{Re} P_{II}$  when  $\alpha = \frac{1}{2}, \frac{3}{2}, \frac{5}{2}, \ldots$ ), (iv) there are no counterparts to the Hastings–McLeod family of solutions, and (v) the Ablowitz–Segur family is present for all values of  $\alpha$  (rather than only for  $|\alpha| < \frac{1}{2}$ ).

One particular finding involves the solutions described by (a)  $u \to 0$  as  $x \to -\infty$  and (b)  $u \sim +\sqrt{x/2}$  as  $x \to +\infty$ . As  $\alpha$  increases, these two solutions rapidly approach each other

along the real axis, in spite of having their pole fields in very different locations in the complex plane and with different residue patterns. These pole fields, which are restricted to wedge-like regions around the real axis, are in turn fundamentally different from those of the Hastings–McLeod solutions to  $Re\ P_{II}$ , which are restricted to wedge-like regions around the imaginary axis (cf. [8, Fig. 10]). Despite these differences in the pole properties, all of these cases have quite similar solution features along the real axis

# <span id="page-9-0"></span>Appendix. Additional observations on the $B(\alpha; z)$ asymptotic expansion

The discussion above relating to (11)–(13) raises a general question about how much information is contained in a rapidly divergent expansion, such as (12).

<span id="page-10-16"></span>![](_page_10_Figure_2.jpeg)

**Fig. A.16.** Imaginary part of the function s(z) defined by (A.1). Note the similarity to the pole fields in Fig. 11 (which are displayed over the smaller domain  $[-10, 10] \times [-10, 10]$ ).

In the case of  $\alpha=0$ , the coefficients in (12), as given by (10), simplify to  $a_n=\frac{\Gamma(3n)}{3^{n-1}\Gamma(n)}$ . For  $\alpha\geq 0$ , the leading order approximation  $a_n\sim\frac{\sinh\pi\alpha}{\alpha}\left(\frac{3}{2}\right)^{2n+\frac{1}{2}}\Gamma(2n+\frac{1}{2})$  follows from [16, p. 428] (with no further correction terms available). The  $\alpha$ -dependence has here been reduced to a scalar factor. Briefly omitting this factor, one might thus consider the expansion

$$-\frac{1}{z} \sum_{n=0}^{\infty} \left(\frac{3}{2}\right)^{2n+\frac{1}{2}} \Gamma\left(2n+\frac{1}{2}\right) / z^{3n}.$$

By the Borel–Ritt theorem [17, p. 28], there will exist functions that obey this to all orders. The choice

$$s(z) = -\frac{\pi}{2z^{1/4}} \left\{ e^{\frac{2}{3}z^{3/2}} \operatorname{Erfc}\left(\sqrt{\frac{2}{3}}z^{3/4}\right) + e^{-\frac{2}{3}z^{3/2}} \operatorname{Erfi}\left(\sqrt{\frac{2}{3}}z^{3/4}\right) \right\}$$
(A.1)

is particularly simple algebraically, featuring a branch cut along  $\mathbb{R}^-$ . Here,  $\operatorname{Erfc}(z)=1-\operatorname{Erf}(z)$  and  $\operatorname{Erfi}(z)=\operatorname{Erf}(iz)/i$  are both entire functions, real along  $\mathbb{R}^+$ . Fig. A.16 displays  $\operatorname{Im} s(z)$  in the vicinity of the origin of the complex z-plane. (The figure for  $\operatorname{Re} s(z)$  is similar and not displayed here.) This function s(z) shares a number of features with all the k=0 solutions (cf. Fig. 11), such as being smooth within a sector with boundaries at  $\operatorname{arg}(z)=\pm\frac{\pi}{3}$ , a pattern of ridges that is reminiscent of the distribution of poles of residues  $\pm i$ , etc. When including the factor  $(\sin\pi\alpha)/\alpha$ , it provides in the smooth sector good approximations to  $B(\alpha;z)$  for low values of  $\alpha$ . However, the lack of any obvious means for error estimation makes it unclear how/if (A.1) can be utilized computationally.

#### References

- <span id="page-10-0"></span>[1] P.A. Clarkson, Painlevé equations—nonlinear special functions, in: F. Marcellàn, W. van Assche (Eds.), Orthogonal Polynomials and Special Functions, in: Lecture Notes in Math., vol. 1883, Springer, Berlin, 2006, pp. 331–411.
- <span id="page-10-1"></span>[2] P.A. Clarkson, Painlevé transcendents, in: F.W.J. Olver, D.W. Lozier, R.F. Boisvert, C.W. Clark (Eds.), NIST Handbook of Mathematical Functions, US Dept. Commerce, Washington, DC, 2010, pp. 723–740.
- <span id="page-10-2"></span>[3] A.S. Abdullayev, Pure imaginary solutions of the second Painlevé equation, Int. J. Math. Math. Sci. 53–56 (2004) 2989–3009.
- [4] A.R. Its, A.S. Fokas, A.A. Kapaev, On the asymptotic analysis of the Painlevé equations via the isomonodromy method, Nonlinearity 7 (5) (1994) 1291–1325
- 1291–1325.
   A.R. Its, V.Y. Novokshenov, The Isomonodromic Deformation Method in the Theory of Painlevé Equations, in: Lecture Notes in Mathematics, vol. 1191, Springer-Verlag, Berlin, 1986.
- <span id="page-10-3"></span>[6] J.A. Giannini, R.I. Joseph, The role of the second Painlevé transcendent in nonlinear optics, Phys. Lett. A 141 (8–9) (1989) 417–419.
- <span id="page-10-4"></span>[7] B. Fornberg, J.A.C. Weideman, A numerical methodology for the Painlevé equations, J. Comput. Phys. 230 (2011) 5957–5973.
- <span id="page-10-13"></span><span id="page-10-5"></span>[8] B. Fornberg, J.A.C. Weideman, A computational exploration of the second Painlevé equation. Found. Comput. Math. 14 (2014) 985–1016.
- <span id="page-10-6"></span>[9] J.A. Reeger, B. Fornberg, Painlevé IV with both parameters zero: A numerical study, Stud. Appl. Math. 130 (2013) 108–133.
- <span id="page-10-7"></span>[10] J.A. Reeger, B. Fornberg, Painlevé IV: A numerical study of the fundamental domain and beyond, Physica D 280–281 (2014) 1–13.
- <span id="page-10-8"></span>[11] I.M. Willers, A new integration algorithm for ordinary differential equations based on continued fraction approximations, Comm. ACM. 17 (1974) 504–508.
- <span id="page-10-9"></span>[12] V.Y. Novokshenov, Tronquée solutions of the Painlevé II equation, Theoret. Math. Phys. 172 (2012) 1136–1146.
- <span id="page-10-10"></span>[13] M. Bertola, On the location of poles for the Ablowitz–Segur family of solutions of the second Painlevé equation, Nonlinearity 25 (2012) 1179–1185.
- <span id="page-10-11"></span>[14] M.J. Ablowitz, H. Segur, Asymptotic solutions of the Korteweg-deVries equation, Stud. Appl. Math. 57 (1) (1976–1977) 13–44.
- <span id="page-10-12"></span>[15] B.M. McCoy, S. Tang, Connection formulae for Painlevé functions, Physica D 18 (1986) 190–196.
- <span id="page-10-14"></span>[16] A.S. Fókas, A.R. Its, A.A. Kapaev, V.Y. Novokshenov, Painlevé Transcendents: The Riemann-Hilbert Approach, American Mathematical Society, Providence, RI, 2006.
- <span id="page-10-15"></span>[17] A. Fruchard, R. Schäfke, Composite Asymptotic Expansions, in: Lecture Notes in Mathematics, vol. 2066, Springer, Heidelberg, 2013.