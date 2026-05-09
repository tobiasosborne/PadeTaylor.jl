# **A Computational Exploration of the Second Painlevé Equation**

**Bengt Fornberg · J.A.C. Weideman**

Received: 3 September 2012 / Revised: 26 March 2013 / Accepted: 9 April 2013 /

Published online: 16 May 2013

© SFoCM 2013

**Abstract** The pole field solver developed recently by the authors (J. Comput. Phys. 230:5957–5973, [2011\)](#page-30-0) is used to survey the space of solutions of the second Painlevé equation that are real on the real axis. This includes well-known solutions such as the Hastings–McLeod and Ablowitz–Segur type of solutions, as well as some novel solutions. The speed and robustness of this pole field solver enable the exploration of pole dynamics in the complex plane as the parameter and initial conditions of the differential equation are varied. Theoretical connection formulas are also verified numerically.

**Keywords** Painlevé transcendents · PII equation · Connection formulas · Pole field solver

**Mathematics Subject Classification (2010)** 33E17 · 34M55 · 34M40 · 33F05

Communicated by Elizabeth Mansfield.

B. Fornberg (-)

Department of Applied Mathematics, University of Colorado, 526 UCB, Boulder, CO 80309, USA e-mail: [fornberg@colorado.edu](mailto:fornberg@colorado.edu)

B. Fornberg

National Institute for Theoretical Physics (NITheP), Western Cape, South Africa

J.A.C. Weideman

Department of Mathematical Sciences, University of Stellenbosch, Private Bag X1, Matieland 7602, South Africa

e-mail: [weideman@sun.ac.za](mailto:weideman@sun.ac.za)

![](_page_0_Picture_19.jpeg)

![](_page_0_Picture_20.jpeg)

## **1 Introduction**

The six Painlevé equations were introduced a little over a century ago [\[28,](#page-31-0) [29\]](#page-31-1). They define transcendental functions that have become firmly established in mathematical physics especially over the last few decades. The first two of these are defined by

<span id="page-1-0"></span>
$$P_{I}: \frac{d^{2}u}{dz^{2}} = 6u^{2} + z, \tag{1}$$

and

$$P_{\text{II}}$$
:  $\frac{d^2 u}{dz^2} = 2u^3 + zu + \alpha,$  (2)

where *u(z)* is a meromorphic function of *z* and *α* is a constant. Where it is necessary to make explicit the dependence of the solution of PII on *α*, we shall write *u* = *u(α*; *z)*.

For definitions of the remaining four equations in the Painlevé family see [[8,](#page-30-1) [10](#page-30-2), [13,](#page-30-3) [16\]](#page-31-2), where lists of various applications can also be found. These include nonlinear wave motion (where PII arises as a reduction of the Korteweg–de Vries equation for water waves [[1,](#page-30-4) [25](#page-31-3), [31](#page-31-4)]), combinatorics, random matrices and statistical physics (where PII appears in the well-known Tracy–Widom distribution [[33\]](#page-31-5)), and electrostatic theory [\[21](#page-31-6), [22\]](#page-31-7). Short summaries of the history of the Painlevé equations are given in the introductions of [\[13](#page-30-3)] and [\[16](#page-31-2)].

The pole distribution of the Painlevé transcendents is an active research area. Not only do the poles carry physical significance in some applications, but the location of the poles in the complex plane also provides information about the solution characteristics (oscillations, decay) on the real line. A lot is known about the distribution of poles as |*z*|→∞, based on various asymptotic approaches to the problem [[4,](#page-30-5) [13,](#page-30-3) [23\]](#page-31-8). Much less is known for finite |*z*|, where numerical computation seems to be the only recourse.

Until recently, the presence of vast pole fields in the complex plane was perceived to be a numerical challenge. The methodology presented in [\[14](#page-30-0)] overcomes this issue. It exploits in two key ways the fact that the solutions are meromorphic:

- (i) the absence of branch-points allows for a flexible and effective path selection strategy even through large and dense pole fields, and
- (ii) a "pole-friendly" Padé-based ODE stepping scheme [[35\]](#page-31-9) along these paths maintains high accuracy both when integrating near to as well as right up to poles of low order.

When this new technique was applied to PI, numerical solutions reported in [\[14](#page-30-0)] yielded errors on the order of 10−<sup>10</sup> after integrating over distances as long as 104 through dense pole fields with average spacing on the order of 1.

Having established an effective numerical procedure, we are in a position to apply it to other Painlevé equations to look for solution features and for pole field patterns that have not been observed before. Here we focus on PII, for the case *α* real and with *u* restricted to be real-valued on the real axis. In particular, Hastings–McLeod solutions (smooth and nonoscillatory along the real axis), Ablowitz–Segur solutions (oscillatory and bounded), and tronquée-type solutions (featuring pole free sectors in

![](_page_1_Picture_15.jpeg)

the complex plane) are of special interest. A similar investigation for  $P_{IV}$  (with both its parameters zero) is reported in [30].

The body of literature for numerical results of  $P_{II}$  is not as extensive as it is for analytical aspects, but it is growing. Solutions on the real axis are presented in [9, 10, 21, 22, 27], and some pole fields in the complex plane, in the special case  $\alpha=0$ , are shown in [4, 26]. The key difference between the present paper and earlier numerical studies is that we investigate here how these pole fields evolve with  $\alpha$ . For such investigations the methodology of [14] offers speed and flexibility, which facilitate experimentation and create the possibility of numerical animations. A selection of these will be made available on the web page [34].

The outline of the paper is as follows: After a brief overview of some general relations and explicit solutions in Sect. 2, we focus in Sect. 3 on pole counting diagrams and the use of our pole field solver as the main tools for surveying the three-parameter  $P_{\rm II}$  solution space (the parameter  $\alpha$  and two initial conditions for the ODE). During this process, we come across all previously identified solution types, and also find some generalizations of these. These generalizations include what we call secondary Hastings–McLeod solutions, studied in Sect. 3, and in the section following that also quasi-Hastings–McLeod and quasi-Ablowitz–Segur solutions. In Sect. 4 we also develop effective approaches for calculating what we have named "k=0 solutions" and for numerically verifying connection formulas for arbitrary  $\alpha$ . We conclude with some additional solution illustrations in Sect. 5 and final remarks in Sect. 6.

## <span id="page-2-0"></span>2 General Relations and Explicit Solutions

Although the solutions to  $P_{\rm II}$  cannot, in general, be expressed in terms of classical special functions, a few special cases are known. They are briefly summarized in this section, in order to place the new computed solutions in the context of existing knowledge. In addition, several series expansions are available. We start this section, however, by listing a few transformations that enable the construction of new solutions from known ones.

#### 2.1 Some Useful Relations

It is readily seen that solutions for  $\pm \alpha$  are connected via the symmetry transformation

<span id="page-2-2"></span><span id="page-2-1"></span>
$$u(\alpha; z) = -u(-\alpha; z). \tag{3}$$

This makes it sufficient to consider only  $\alpha \ge 0$ , which we assume throughout this paper.

The Bäcklund transformation is

$$u(\alpha + 1; z) = -u(\alpha; z) - \frac{2\alpha + 1}{2u'(\alpha; z) + 2u(\alpha; z)^2 + z},$$
(4)

which defines a useful recurrence relation for generating solutions; see [10, Sect. 32.7], [15], [16, Sect. 19]. The prime denotes differentiation with respect

![](_page_2_Picture_14.jpeg)

to z. (A similar relation is available for  $u(\alpha - 1; z)$  in terms of  $u(\alpha; z)$ ; see [10, Sect. 32.7].) Solutions corresponding to  $\alpha = 0$  and  $\frac{1}{2}$  are connected via

<span id="page-3-0"></span>
$$u\left(\frac{1}{2};z\right) = 2^{-1/3} \frac{u'(0;-2^{-1/3}z)}{u(0;-2^{-1/3}z)},\tag{5}$$

a relation that will be discussed in Sect. 5.2. In principle it is therefore only necessary to consider  $0 \le \alpha < \frac{1}{2}$ , as solutions for arbitrary  $\alpha$  can then be reconstructed from (3) and (4) (or (5)). We shall not limit our investigations to this fundamental interval in parameter space, however, as it is not straightforward to see how pole locations and other solution features of  $u(\alpha; z)$  and  $u(\alpha + 1; z)$  are connected by (4).

Near a pole, say at  $z = z_0$ , the Laurent expansion of  $P_{II}$  solutions follows immediately from substitution into (2), namely

$$u(z) = \frac{c_{-1}}{z - z_0} + c_0 + c_1(z - z_0) + c_2(z - z_0)^2 + c_3(z - z_0)^3 + c_4(z - z_0)^4 + O((z - z_0)^5),$$

where

$$c_{-1} = \pm 1,$$
  $c_0 = 0,$   $c_1 = \mp \frac{z_0}{6},$   $c_2 = \frac{\mp 1 - \alpha}{4},$   $c_3 = \gamma,$   $c_4 = \frac{z_0}{72}(\pm 1 + 3\alpha);$ 

<span id="page-3-2"></span>see [2, Sect. 3.6] and [16, Sect. 17]. The only free parameters (in the infinite expansion) are the pole location  $z_0$ , the sign choice in  $c_{-1}$  and the coefficient  $\gamma$  (first appearing in  $c_3$ ). All poles are therefore of first order, with residues of either +1 or -1. This is in contrast to the situation for  $P_I$  where all poles are of second order, with strength of 1 and residue 0.

In the figures below, poles with residues +1 or -1 are denoted by dark (blue) and light (yellow) circles, respectively. Where zeros are shown, a smaller (red) square is used.

#### 2.2 Rational Solutions

There exists one rational solution for each  $\alpha = n$ , n = 1, 2, 3, ..., which will be denoted by  $u_n(z)$  [8, 10, 16]. (Because of (3), there is a corresponding rational solution for each  $\alpha = -n$ .) One way to generate them is to start with  $u_0 = 0$  and then apply the recursion (4). This yields the sequence

<span id="page-3-1"></span>
$$u_1(z) = -\frac{1}{z},$$
  $u_2(z) = \frac{4 - 2z^3}{4z + z^4},$   $u_3(z) = \frac{3z^2(160 + 8z^3 + z^6)}{320 - 24z^6 - z^9},$  (6)

etc., with some further members listed explicitly in [8] and [16, Sect. 20]. Figure 1 shows pole and zero locations of these solutions; similar plots were displayed previously in [8, 11]. These rational solutions are the only solutions of  $P_{II}$  to have finitely many poles in the complex plane [16, Sect. 20]. The vast pole-free (tronquée) regions in Fig. 1 are therefore rather special and are in fact completely destroyed by perturbations in the data (for example, changes in  $\alpha$  or changes in the values of u(0), u'(0)).

![](_page_4_Figure_2.jpeg)

<span id="page-4-3"></span><span id="page-4-0"></span>Fig. 1 Poles and zeros of the first six rational solutions of  $P_{II}$ , in the complex plane z = x + iy. The *dark* (*blue*) and *light* (*yellow*) *circles* represent poles with residues +1 or -1, respectively. Zeros are represented by the *smaller* (*red*) *squares* (Color figure online)

Note in Fig. 1 how the zeros of the rational solutions interlace the poles. The same thing is observed in Fig. 2 of the next section. To avoid clutter, we shall therefore cease to plot zeros along with the poles from there on.

## 2.3 Airy Solutions

There is a one-parameter family of Airy solutions when  $\alpha = n + \frac{1}{2}$ , which will be denoted by  $u_{n+\frac{1}{2}}(z)$  [8, 10, 16]. To define them, let  $\phi$  be the solution to  $\phi'' = -\frac{1}{2}z\phi$ , i.e.,

<span id="page-4-2"></span><span id="page-4-1"></span>
$$\phi(z) = c_1 \operatorname{Ai}\left(\frac{-z}{2^{1/3}}\right) + c_2 \operatorname{Bi}\left(\frac{-z}{2^{1/3}}\right),$$

with  $c_1$  and  $c_2$  arbitrary constants. Now define  $\Phi(z) = \phi'(z)/\phi(z)$  and combine  $c_1$  and  $c_2$ , by setting

$$c_1 = \cos\frac{\theta}{2}, \qquad c_2 = \sin\frac{\theta}{2}, \quad 0 \le \theta \le 2\pi,$$
 (7)

making  $\Phi(z)$   $2\pi$ -periodic in the parameter  $\theta$ .

The Airy solution corresponding to  $\alpha = \frac{1}{2}$  is then given by

$$u_{\frac{1}{2}}(z) = -\Phi,\tag{8}$$

![](_page_4_Picture_13.jpeg)

<span id="page-5-0"></span>![](_page_5_Figure_2.jpeg)

Fig. 2 Poles and zeros of the first three Airy solutions of  $P_{II}$  in the case  $\theta = 0$ . For typical variation with  $\theta$ , see Fig. 3 (Color figure online)

and the others follow from (4) together with the relation  $\Phi'(z) = -\frac{1}{2}z - \Phi(z)^2$ . The next two are

<span id="page-5-1"></span>
$$u_{\frac{3}{2}}(z) = \frac{2\Phi^3 + z\Phi - 1}{2\Phi^2 + z} \tag{9}$$

and

$$u_{\frac{5}{2}}(z) = \frac{4z\Phi^4 + 6\Phi^3 + 4z^2\Phi^2 + 3z\Phi + z^3 - 1}{(4\Phi^3 + 2z\Phi - 1)(2\Phi^2 + z)},$$
 (10)

with some further members listed in [8] and [16, Sect. 21]. Because of the free parameter  $\theta$ , the Airy solutions define a much bigger family of solutions than the rational solutions.

Figure 2 shows locations of poles and zeros of the solutions (8)–(10), corresponding to the case  $\theta=0$ . These follow single lines aligned with the real axis when  $\alpha=\frac{1}{2}$ , triple lines when  $\alpha=\frac{3}{2}$ , etc. Figure 3 shows similar pole fields but with  $\theta$  ranging over its  $2\pi$ -cycle, in the case  $\alpha=\frac{5}{2}$ . Notice how five curves of poles enter from  $-\infty$ , meet up with another group of poles aligned along the positive real axis, and then carry with them five poles from this group back to  $-\infty$ . For the general Airy solution with  $\alpha=n+\frac{1}{2}$ , a total of 2n+1 poles are transferred in this manner.

The pole fields displayed in Figs. 1–3 were computed from explicit solution formulas. In the sequel, such formulas are unavailable and all pole fields shown were computed with the pole field solver of [14]. The first of these is displayed in Fig. 4, which shows pole dynamics as  $\alpha$  is varied across  $[\frac{5}{2}, \frac{7}{2}]$ , for the special case u(0) = u'(0) = 0. This sequence starts and ends with Airy solutions, with a rational solution halfway between. During the transition the central group of nine poles remains relatively intact as curves of poles alternately enter from infinity, and recede back to infinity.

The class of solutions of  $P_{II}$  that satisfy u(0) = u'(0) = 0 for  $\alpha \neq 0$  was analyzed in [23], where asymptotic formulas are given in the limits  $x \to \pm \infty$  (here, and below, we use x in place of z to indicate that the variable is real). Also noted in [23] is the fact that u(3n; z) give rational solutions for  $n = 1, 2, 3, \ldots$ , while both  $u(3n - \frac{1}{2}; z)$  and  $u(3n + \frac{1}{2}; z)$  give Airy solutions. The first, middle and last subplots of Fig. 4

![](_page_5_Picture_12.jpeg)

![](_page_6_Figure_2.jpeg)

<span id="page-6-1"></span>Fig. 3 Poles of the Airy family of solutions of  $P_{II}$  ( $\alpha=\frac{5}{2}$ ). Note in the  $\theta=\frac{1}{3}\pi$  and  $\frac{5}{3}\pi$  cases the exact symmetries, and also in these cases the resemblance of the central pole groups to the rational solutions for  $\alpha=2$  and  $\alpha=3$ , respectively (Other cases, such as  $\theta=\frac{2}{3}\pi$  and  $\frac{4}{3}\pi$ , are similar to the  $\theta=\pi$  case in terms of lacking the 3-fold symmetry. For  $\theta=2\pi$ , the solution has returned to the  $\theta=0$  case) (Color figure online)

represent the n=1 case. Similar results hold when  $\alpha$  is varied and the origin is a pole of either residue  $\pm 1$ .

Figures 3 and 4 also illustrate a type of pole behavior that is often seen near tronquée solutions of  $P_{II}$ . As the parameter ( $\alpha$  or  $\theta$ , for example) is varied through its critical value, curves of poles move out to infinity (leaving behind the tronquée solution), and then move back in. At infinity, there can be a vertical shift in alignment of the poles, which affects solution features on either side of the tronquée case. In Fig. 4, there is a tronquée solution right in the middle ( $\alpha$  = 3). Note the change in alignment of the poles in the curves in the right half-plane as they move out to  $+\infty$  and back in as  $\alpha$  is varied through the critical value 3 (cf. the middle row of subplots in Fig. 4). Similar behavior is seen in Fig. 3, where the first subplot ( $\theta$  = 0) is a tronquée case. Note the alignment shift in the poles in the curves in the left half-plane as  $\theta$  is varied through the critical value 0 (cf. the sixth subplot, which is the same as  $\theta$  =  $-10^{-5}$ , followed by the first subplot and then the second).

#### <span id="page-6-0"></span>3 Pole Counting Diagrams

The rational and the Airy solutions discussed in Sect. 2 are the only known explicit (or classical) solutions of  $P_{\rm II}$ . There are furthermore solutions associated with the

![](_page_6_Picture_8.jpeg)

![](_page_7_Figure_2.jpeg)

<span id="page-7-0"></span>**Fig. 4** Poles of  $P_{II}$  in the case u(0) = u'(0) = 0 as  $\alpha$  is varied (in nonuniform increments) across  $\left[\frac{5}{2}, \frac{7}{2}\right]$ . The *first* and *last subplots* show Airy solutions, namely  $u_{\frac{5}{2}}(z)$   $(\theta = \frac{5}{3}\pi)$  and  $u_{\frac{7}{2}}(z)$   $(\theta = \frac{1}{3}\pi)$ . The *middle subplot* shows the rational solution  $u_3(z)$  (Color figure online)

names of Ablowitz, Boutroux, Hastings, McLeod, Segur, and others. No expressions for these solutions in terms of known special functions are available, however.

In the remainder of this paper we shall explore these and other solutions for  $\alpha \geq 0$ . Following [14, 30], we base the discussion on pole counting diagrams. These diagrams display the number of poles on the positive and negative real axes (denoted by  $\mathbb{R}^+$  and  $\mathbb{R}^-$ , respectively) as the initial data (u(0), u'(0)) are varied. This type of display covers for each  $\alpha$  all possible solutions and it is therefore particularly effective in identifying and relating solutions of different types, such as those referred to in the paragraph above.

A sequence of pole counting diagrams is shown in Fig. 5. Initial values that generate a finite number of poles along  $\mathbb{R}^+$  define curves in the (u(0), u'(0))-plane. These are labeled  $n^+$ , where n is the pole count. When similarly counting poles along  $\mathbb{R}^-$ ,

![](_page_7_Picture_7.jpeg)

![](_page_8_Figure_2.jpeg)

<span id="page-8-0"></span>**Fig. 5** Pole counting diagrams in the (u(0), u'(0))-plane. *Curves* labeled  $n^+$  denote initial conditions that generate solutions with n poles on  $\mathbb{R}^+$ . *Curves* and *regions* labeled  $n^-$  represent n poles on  $\mathbb{R}^-$ 

![](_page_8_Picture_4.jpeg)

both curves (zero width) and regions (finite width) are possible. As we shall see in Sect. 5.1, the count  $n^-$  within a region need not agree with the count along its edges. All unmarked (white) regions feature an infinity of poles along both  $\mathbb{R}^+$  and  $\mathbb{R}^-$ .

In each diagram in Fig. 5, there is a single curve marked  $0^+$  that runs diagonally between top left and bottom right, roughly. These  $0^+$  curves represent solutions that are pole free on  $\mathbb{R}^+$ . On either side of the  $0^+$  curves there are two curves running top to bottom, labeled  $1^+$ . These denote solutions with a single pole on  $\mathbb{R}^+$ . The adjacent pair of curves farther out denote solutions with two poles on  $\mathbb{R}^+$ , etc. (although we suppressed their labels to avoid clutter). The pole count increases by 1 from the inside out.

With regard to  $\mathbb{R}^-$ , there is a central curve marked  $0^-$  that runs diagonally between top right and bottom left, with an accompanying  $0^-$  region that remains visible in the first five subplots in Fig. 5 (but which does not vanish entirely until  $\alpha = \frac{3}{2}$ ). These curves/regions represent solutions that are pole free on  $\mathbb{R}^-$ . Within regions marked  $0^-$ ,  $1^-$ , etc., the indicated number of poles is generally followed by oscillations in u(x) as  $x \to -\infty$  (examples can be seen in Figs. 12, 13, 18 and 19); the rational solutions are exceptional in that the oscillation amplitudes vanish.

In the sections below, we shall see that many of the special solutions (and further generalizations of these) can be identified with intersections of curves or of curves with regions in these diagrams. For example, the intersection of the  $0^+$  curves with either the  $0^-$  curves or the  $0^-$  regions, represent solutions that are pole free on the entire real axis. Before looking at this in more detail, we describe how the diagrams in Fig. 5 were created.

<span id="page-9-0"></span>Using the pole field solver of [14], we integrated  $P_{\Pi}$  as an initial-value problem, first along the interval [0, L] using a large number of initial values (u(0), u'(0)). For each set of initial values we recorded the number of poles on [0, L]. To obtain these counts, a strategy based on the residue theorem was used. This was checked with a strategy based on computing zeros of Padé denominators. The value of  $L \gg 1$  was adjusted until we reached confidence that the pole count thus obtained is accurate for  $\mathbb{R}^+$ . The process was then repeated for  $\mathbb{R}^-$ .

## 3.1 Re-visiting the Rational and Airy Solutions

We begin the survey of the  $P_{II}$  solution space by noting where the explicit solutions of Sect. 2 fit into the pole counting diagrams. For example, the rational solution  $u_3(z)$  defined in (6) is generated by initial data (u(0), u'(0)) = (0, 0). Looking at the  $\alpha = 3$  diagram in Fig. 5, we note that the origin is located on a  $1^+$  curve and inside a  $2^-$  region. The one pole on  $\mathbb{R}^+$  and two poles on  $\mathbb{R}^-$  can be seen in the third subplot of Fig. 1.

The rational solutions  $u_1(z)$  and  $u_2(z)$ , by contrast, have poles at the origin and thus (u(0), u'(0)) are both infinite. To represent such solutions, the pole counting diagrams have to be extended to infinity. Consider, for example, the rational solution  $u_1(z) = -1/z$ , with  $u_1'(z) = 1/z^2$ . If this solution is perturbed so that this pole crosses the origin in the direction of  $\mathbb{R}^+$ , the initial data (u(0), u'(0)) switch from  $(-\infty, +\infty)$  to  $(+\infty, +\infty)$  and the pole count on  $\mathbb{R}^+$  increases from 0 to 1. This means that the  $0^+$  and  $1^+$  curves are connected at infinity, as indicated by the dash-dot line

![](_page_9_Picture_10.jpeg)

![](_page_10_Figure_2.jpeg)

![](_page_10_Figure_3.jpeg)

<span id="page-10-0"></span>**Fig. 6** Extended pole counting diagrams that show (in *dash-dot line* type) the location at infinity of the  $u_1(z)$  and  $u_2(z)$  rational solutions defined in (6). For clarity, only pole counting details for  $\mathbb{R}^+$  are included here

segment in the first subplot of Fig. 6. Likewise, when the pole at z = 0 of  $u_2(z) = (4-2z^3)/(4z+z^4) \sim 1/z$ ,  $z \to 0$ , crosses the origin from left to right, the initial data switch from  $(+\infty, -\infty)$  to  $(-\infty, -\infty)$  while the pole count on  $\mathbb{R}^+$  increases from 0 to 1. The corresponding connection is shown in the second subplot of Fig. 6.

In general, if a pole of residue  $\pm 1$  crosses the origin in this manner, then  $u' \sim \mp u^2$  shows that the connecting curve will be at  $u'(0) = \mp \infty$  in the diagram. The relationship  $u' \sim \mp u^2$  also confirms the parabolic nature of the curves where  $|u'(0)| \gg 1$ . In Sect. 5.1, similar arguments will be applied to pole counts on  $\mathbb{R}^-$ .

Next, we consider the Airy solutions. Because they are one-parameter families of solutions, they define not single points but curves in the pole counting diagrams when the parameter  $\theta$  defined in (7) is varied across its  $2\pi$ -cycle. Figure 7 shows two such Airy curves, superimposed on the corresponding diagrams from Fig. 5. Observe how the Airy solution curves trace out edges that exist in the pole counting diagrams. In the case  $\alpha = \frac{1}{2}$ , this curve is the parabola  $u'(0) = u(0)^2$ . Four additional Airy curves, without pole counting diagrams, are shown in Fig. 8.

<span id="page-10-1"></span>There is one detail in Fig. 7 that deserves pointing out. In the  $\alpha=\frac{1}{2}$  diagram (first subplot) there is a small but visible gap between the intersection point of the  $0^+$  and  $0^-$  curves and the Airy initial conditions corresponding to  $\theta=0$ . This gap seems to have disappeared in the  $\alpha=\frac{3}{2}$  diagram (second subplot), but it still exists. This can be seen in the magnified version of the diagram shown in the third subplot of Fig. 7. (The gap gets smaller still for  $\alpha=\frac{5}{2},\frac{7}{2}$ , etc.) We postpone an examination of solution features associated with the third subplot of Fig. 7 until Sect. 5.3.

### 3.2 Hastings–McLeod Solutions

Figure 9 shows four pole counting diagrams, three of which were displayed on a larger domain in Fig. 5. Only the  $0^+$  curves and  $0^-$  curves/regions are now shown, because our present interest is in solutions that are pole free on the entire real axis.

![](_page_10_Picture_11.jpeg)

![](_page_11_Figure_2.jpeg)

<span id="page-11-0"></span>Fig. 7 The first two subplots show (in thicker line type) curves of Airy initial conditions, superimposed on the  $\alpha$  half-integer pole counting diagrams. The third subplot shows the detail in a small neighborhood of the  $\theta = 0$  point in the second subplot. Solutions in this region will be examined in more detail in Sect. 5.3

<span id="page-11-1"></span>![](_page_11_Figure_4.jpeg)

Fig. 8 Curves of initial conditions of Airy solutions as the parameter  $\theta$  in (7) is varied over its  $2\pi$ -cycle

When  $\alpha = 0$ , the  $0^+$  curve is seen to intersect the  $0^-$  region along a curve segment that connects the two points

$$(u(0), u'(0)) \approx (\pm 0.3670615515480784, \mp 0.2953721054475501).$$

The corresponding two solutions (which, by (3), only differ in sign), are therefore pole free and nonoscillatory on the real axis, i.e., this corresponds to the well-known Hastings–McLeod solution [17]. When  $\alpha$  increases from 0 the symmetry is broken and the two dots in Fig. 9 move closer together, coalescing when  $\alpha = \frac{1}{2}$ . This implies that if  $0 < \alpha < \frac{1}{2}$ , then there are two (nonsymmetric) solutions that are pole free and nonoscillatory on the entire real axis. Only one of these has been described in any detail in the literature, as a generalized Hastings–McLeod solution [3, 7, 13].

![](_page_11_Picture_9.jpeg)

<span id="page-12-0"></span>![](_page_12_Figure_2.jpeg)

**Fig. 9** Pole counting diagrams on  $[-2, 2] \times [-2, 2]$  showing the location of Hastings–McLeod initial conditions as *dots* 

To distinguish between these two solutions, consider for  $\alpha > 0$  the asymptotic boundary conditions [7, 10, 13, 17, 18]

<span id="page-12-2"></span><span id="page-12-1"></span>
$$u(x) \sim \begin{cases} \pm \sqrt{-\frac{1}{2}x}, & x \to -\infty; \\ -\alpha/x, & x \to +\infty. \end{cases}$$
 (11)

These are obtained as asymptotic balances between the right-hand side terms in (2) when solutions are smooth and the second derivative becomes negligible, giving

$$2u^3 + xu + \alpha = 0. (12)$$

(In the case  $\alpha = 0$ , this argument needs refinement; see (13) and the discussion that follows it.)

The plus and minus signs in (11) are identified with the lower and upper dots in Fig. 9, respectively. Corresponding solutions in the case  $\alpha=0.495$  are shown in the middle subplot of the first column of Fig. 10. The lower solution (minus sign in (11)), is pole free, nonoscillatory, and monotone, and is the solution described in [3, 7, 13]. The upper solution (plus sign in (11)) is likewise pole free and nonoscillatory but not monotone. These features of the upper solution have not been noted in the literature as far as we know, although its asymptotic properties are recorded in [13, Theorem 11.7]. We refer to these two solutions as the primary and secondary Hastings–McLeod solutions, respectively.

As  $\alpha$  increases from 0, the primary solution changes slowly as can be seen in Fig. 10. By contrast, the secondary solution develops a steep gradient on the negative real axis, which moves out to  $-\infty$  as  $\alpha \to \frac{1}{2}$ . When  $\alpha = \frac{1}{2}$ ,  $u \sim -\sqrt{-\frac{1}{2}x}$  as  $x \to -\infty$ , and the two solutions coalesce pointwise.

The pole fields shown in the second and third columns of Fig. 10 clarify the process. They were obtained using the values of (u(0), u'(0)) for which approximations are given in the left column of the figure. (Similar pole fields for the case  $\alpha = 0$  have previously been displayed in [26].)

In the second column, we note that the poles for the primary Hastings–McLeod solutions are located in two wedge-like regions, well separated from the real axis. As  $\alpha$  is increased, these wedges move slowly to the right and away from the real axis, thereby preserving the smoothness of the primary solution. (We look at this in more detail in Sect. 4.3.) By contrast, the third column shows that, as  $\alpha$  increases,

![](_page_12_Picture_13.jpeg)

![](_page_13_Figure_2.jpeg)

<span id="page-13-0"></span>**Fig. 10** Hastings—McLeod solutions and their corresponding pole fields. The *thinner curves* in the *first column* are the branches of the cubic equation (12) that defines the asymptotic boundary conditions (11) (One of these is completely obscured behind the *thick solution curves*, but is better visible in Figs. 12, 13 and 18 below) (Color figure online)

the leftmost curve of poles dislodges from the two pole wedges. A conjugate pair of these poles approaches the real axis as  $\alpha \to \frac{1}{2}$ , which causes the steepening of the gradient of the secondary solution already noted.

### 3.3 Ablowitz–Segur Solutions for $\alpha = 0$

This class of solutions was first described in [1]. Like the Hastings–McLeod solutions the Ablowitz–Segur solutions are pole free on the entire real axis, but in addition they are bounded on all of  $\mathbb{R}$  [2, Sect. 3.7].

In the pole counting diagrams, initial conditions corresponding to the Ablowitz–Segur solutions lie on  $0^+$  curves that connect the two Hastings–McLeod points. This can be seen in Fig. 11, which displays two pole counting diagrams, once again showing only  $0^+$  and  $0^-$  detail. We examine the case  $\alpha>0$  in Sect. 4, and focus on  $\alpha=0$  for now.

![](_page_13_Picture_8.jpeg)

![](_page_14_Figure_2.jpeg)

<span id="page-14-1"></span>**Fig. 11** Pole counting diagrams showing the location of Ablowitz–Segur initial conditions as the *thicker curve* sections that connect the Hastings–McLeod points. The *shaded region* is the  $0^-$  region in Fig. 5 and the curve is the  $0^+$  curve. The values of k that correspond to the Hastings–McLeod points are also shown; cf. (21)

The Ablowitz–Segur solutions decay to zero as  $x \to +\infty$ . To leading order Eq. (2), with  $\alpha = 0$ , then reduces to the linear Airy equation, i.e., it will hold that

<span id="page-14-4"></span><span id="page-14-2"></span><span id="page-14-0"></span>
$$u(x) \sim k \operatorname{Ai}(x), \quad x \to +\infty,$$
 (13)

where k is a free parameter.

The question addressed and solved in [1, 3, 12, 13, 17, 24, 25, 31, 32] and elsewhere, is to find the asymptotic behavior as  $x \to -\infty$  of the solution that satisfies (13) as  $x \to +\infty$ . The answer depends on k and one needs to distinguish between |k| < 1, |k| = 1 and |k| > 1.

In the case |k| < 1, the asymptotic behavior on the negative real axis is

$$u(x) \sim d(-x)^{1/4} \sin\left(\frac{2}{3}(-x)^{3/2} - \frac{3}{4}d^2\log(-x) - \theta_0\right), \quad x \to -\infty,$$
 (14)

where d and  $\theta_0$  are constants. The connection formulas

<span id="page-14-3"></span>
$$d^2 = -\pi^{-1}\log(1 - k^2) \tag{15}$$

and

$$\theta_0 = \frac{3}{2}d^2\log 2 + \arg\Gamma\left(1 - \frac{1}{2}id^2\right) + \frac{1}{4}\pi(1 - 2\operatorname{sgn}k) \tag{16}$$

hold when the two parts of the solution, (13) and (14), are connected smoothly.

As  $k \to \pm 1$  it was proven in [17] that the oscillatory behavior of (14) turns into the square root behavior  $u \sim \pm \sqrt{-\frac{1}{2}x}$ . (This defines the Hastings–McLeod solution, which means the second boundary condition in (11) should be modified to  $u(x) \sim \pm \mathrm{Ai}(x)$ ,  $x \to +\infty$  in the case  $\alpha = 0$ .) When |k| > 1, poles appear on the negative real axis. (The physical significance of the value |k| = 1, in the water wave application, is that it distinguishes between nonlinearity dominating dispersion and vice versa [25, 31].)

![](_page_14_Picture_16.jpeg)

![](_page_15_Figure_2.jpeg)

<span id="page-15-0"></span>**Fig. 12** Solutions of PII (*α* = 0) corresponding to two values of *k* in the asymptotic boundary condition *u(x)* ∼ *k*Ai*(x)*, *x* → +∞; cf. ([17\)](#page-16-1). The solution in the *first subplot* is an Ablowitz–Segur solution. (The in-between case, *k* = 1, is one of the Hastings–McLeod solutions shown in the first subplot of Fig. [10.](#page-13-0)) The *figures on the right* are the pole fields of the solutions on the *left* (Color figure online)

Figure [12](#page-15-0) shows solutions for values of *k* near the critical value 1, i.e., near the secondary Hastings–McLeod solution. (By ([3\)](#page-2-1), there exist similar solutions near the primary Hastings–McLeod solution, *k* = −1.) When *k* is slightly less than the critical value we observe an Ablowitz–Segur solution, with oscillatory tail described by ([14\)](#page-14-2). (In fact, when graphs of ([13\)](#page-14-0) and ([14](#page-14-2)) are superimposed on the graph shown here, no discrepancy can be seen in the regions *x >* 0 and *x <* −3, respectively. Similar agreements were previously shown in [\[25](#page-31-3), [31](#page-31-4)].) When *k* slightly exceeds the critical value, the solution features an infinity of poles on R−.

Similar solutions (along the real axis) were displayed previously in [\[8](#page-30-1), [10](#page-30-2), [21,](#page-31-6) [22](#page-31-7)]; in fact, the values *k* = 1 ± 0*.*001 were chosen so that the solution curves in Fig. [12](#page-15-0) match the curves in [\[10,](#page-30-2) Fig. 32.3.6]. The corresponding pole fields are displayed in the right column of Fig. [12](#page-15-0).

The pole fields consist of the wedge-like sectors already observed in Fig. [10](#page-13-0), supplemented by a separate field farther into the left half-plane. When *k* is slightly less than 1 some of these poles are located just off the real axis, which causes the oscillatory features of the solution seen in the first subplot of Fig. [12](#page-15-0) and described asymptotically by [\(14](#page-14-2)). As *k* increases this field moves to the left, and at the criti-

![](_page_15_Picture_7.jpeg)

<span id="page-16-0"></span>cal k=1 it reaches  $z=-\infty$ , leaving only the two wedges. As k increases further, the poles enter again from  $z=-\infty$ , but this time vertically aligned such that there are poles along the real axis. This is the familiar behavior near a tronquée solution, already noted at the end of Sect. 2.3.

#### 4 Generalizations to $\alpha > 0$

In addition to the primary and secondary Hastings–McLeod solutions defined in Sect. 3.2, we describe below a class of solutions for the case  $\alpha > \frac{1}{2}$ , which has the same asymptotic properties as the secondary Hastings–McLeod solution, but with a finite number of poles on  $\mathbb{R}$ . These will be referred to as quasi-Hastings–McLeod solutions. To keep the discussion uncluttered, we introduce here the abbreviations pHM, sHM, and qHM for primary, secondary, and quasi-Hastings–McLeod solutions, respectively. Similarly, we abbreviate the regular and quasi-Ablowitz–Segur solutions by AS and qAS, respectively.

### 4.1 Refinements of the Boundary Condition for $x \to +\infty$

Here we derive higher order expansions for the boundary conditions (11) and (13) that will be necessary for our numerical work. We first note that any solution  $u(\alpha; x)$  that is smooth for  $x \to +\infty$  can be written as

<span id="page-16-2"></span><span id="page-16-1"></span>
$$u(\alpha; x) = B(\alpha; x) + e(\alpha, k; x), \tag{17}$$

where  $B(\alpha; x)$  satisfies an asymptotic series of the form

$$B(\alpha; x) \sim -\frac{\alpha}{x} \sum_{n=0}^{\infty} \frac{b_n}{x^{3n}},\tag{18}$$

and  $e(\alpha, k; x) \sim k \operatorname{Ai}(x)$  is independent of  $\alpha$  to leading order and contains only exponentially small terms in x [5, 6]. (The expansions (17)–(18) hold, in fact, in the sector  $\operatorname{arg} x \in (-\frac{1}{3}\pi, \frac{1}{3}\pi)$ ). This is a result due to Boutroux, who also considered similar expansions in larger sectors but which are not real on the real axis, and therefore outside the scope of the present paper.)

In (18),  $b_0 = 1$  and a recurrence formula for the other  $b_n$  is given in [8, 15]. It transpires that each  $b_n$  is a polynomial of degree 2n in even powers of  $\alpha$ .

Besides this dependence on  $\alpha$ , the function  $B(\alpha; x)$  contains no free parameters, while  $e(\alpha, k; x)$  contains the additional parameter k. In the case  $\alpha = 0$ , the function  $B(\alpha; x)$  vanishes and (17) matches (13) in the limit  $x \to +\infty$ . When  $\alpha = 1, 2, 3, \ldots$ , the series (18) converges and reproduces the rational solutions described in Sect. 2.2 [9]. In this case (17) matches the type of boundary condition considered in [2, 9].

The function  $e(\alpha, k; x)$  can likewise be expanded in the limit  $x \to +\infty$  as

$$e(\alpha, k; x) \sim \frac{k}{2\sqrt{\pi}} \frac{e^{-(2/3)x^{3/2}}}{x^{1/4}} \sum_{n=0}^{\infty} \frac{c_n}{x^{3n/2}} + \left(\frac{k}{2\sqrt{\pi}}\right)^2 \frac{e^{-(4/3)x^{3/2}}}{x^{4/4}} \sum_{n=0}^{\infty} \frac{d_n}{x^{3n/2}} + O\left(k^3 \frac{e^{-(6/3)x^{3/2}}}{x^{7/4}}\right), \quad (19)$$

<span id="page-16-3"></span>E°L''  with coefficients

$$c_0 = 1,$$
  $c_1 = \frac{1}{48} (-5 + 96\alpha^2),$   
 $c_2 = \frac{1}{4608} (385 - 7872\alpha^2 + 9216\alpha^4),$  etc.,

and

$$d_0 = 0$$
,  $d_1 = -2\alpha$ ,  $d_2 = \frac{1}{12} (77\alpha - 96\alpha^3)$ , etc.

These (as well as the  $b_n$ ) are readily generated by substituting (17) (expanded according to (18) and (19)) into (2) and equating coefficients.

### 4.2 Connection Formulas

Connection formulas analogous to (14)–(16) for arbitrary  $\alpha$  were derived in [20, 24]. We do not reproduce these formulas here, except to note that (15) should be modified to

<span id="page-17-0"></span>
$$d^{2} = -\pi^{-1} \log(\cos^{2}(\pi \alpha) - k^{2}). \tag{20}$$

(For the meaning of d we refer to [24].) The critical values of k that define pHM and sHM solutions are therefore

<span id="page-17-1"></span>
$$k_p = -\cos\pi\alpha, \qquad k_s = \cos\pi\alpha.$$
 (21)

The main challenge in verifying this theory numerically is to compute the k=0 function  $B(\alpha;x)$ . When  $\alpha$  is not an integer, the series (18) diverges violently and neither truncation at the optimal point nor sequence acceleration can be used to compute  $B(\alpha;x)$  accurately. This is not surprising since (18) alone, even with all the  $b_n$  known, does not define  $B(\alpha;x)$  uniquely. The expansion (19), however, allows this ambiguity to be bypassed, also leading to an effective computational procedure to obtain a unique  $B(\alpha;x)$  (for any  $\alpha>0$ ), which we will call the k=0 solution.

Based on the two critical choices for k defined in (21), it follows from (17) that

$$u_p(\alpha; x) = B(\alpha; x) + e(\alpha, k_p; x),$$
  $u_s(\alpha; x) = B(\alpha; x) + e(\alpha, k_s; x),$  (22)

where  $u_p(\alpha; x)$  and  $u_s(\alpha; x)$  are the pHM and sHM solutions. Because  $u_p(\alpha; x)$  is unique and pole free along the entire real axis for every  $\alpha$ , it is computed without difficulty as the solution to an ODE boundary-value problem. The expansion (19) provides accurate values of  $e(\alpha, k_p; x)$  and  $e(\alpha, k_s; x)$  and their x-derivatives at sufficiently large x = X. Using, for example,  $c_0, c_1, \ldots, c_{12}$  and  $d_0, d_1, \ldots, d_8$ , near machine accuracy is achieved with X = 8. Because (22) implies

<span id="page-17-2"></span>
$$u_s(\alpha; x) = u_p(\alpha; x) + e(\alpha, k_s; x) - e(\alpha, k_p; x)$$
(23)

the sHM solution can now be obtained by initiating the pole field solver with the now known values for the right-hand side of (23), i.e., for  $u_s(\alpha; X)$  and  $u'_s(\alpha; X)$ . A similar procedure is used for computing  $B(\alpha; x)$ , based on the first equation in (22). (We remark that in the case  $0 \le \alpha < \frac{1}{2}$  an alternative procedure for computing  $u_s(\alpha; x)$  is to solve it as a boundary-value problem, for it too is then a pole free solution albeit

![](_page_17_Picture_18.jpeg)

![](_page_18_Figure_2.jpeg)

<span id="page-18-0"></span>**Fig. 13** Solutions of  $P_{\text{II}}$  ( $\alpha=\frac{1}{3}$ ) corresponding to four values of k in the asymptotic boundary condition  $u(\frac{1}{3};x)\sim B(\frac{1}{3};x)+k\mathrm{Ai}(x),\ x\to +\infty;$  cf. (17). Values  $|k|<\frac{1}{2}$  correspond to the smooth, oscillatory solutions (Ablowitz–Segur) and values  $|k|>\frac{1}{2}$  correspond to the solutions with poles on  $\mathbb{R}^-$ 

not as smooth as  $u_p(\alpha; x)$ ; cf. Fig. 10. Being able to compute  $u_s(\alpha; x)$  in two distinct ways offers a useful numerical check.)

Before we have a closer look at the solutions  $u_s(\alpha; x)$  and  $B(\alpha; x)$ , let us return briefly to the pole counting diagram for  $\alpha = \frac{1}{3}$  shown in the second subplot of Fig. 11. The objective here is to verify the validity of the critical k formula (21), which turns into  $k_p = -\frac{1}{2}$  and  $k_s = \frac{1}{2}$  in the  $\alpha = \frac{1}{3}$  case. The expansion (19) was used to compute  $e(\frac{1}{3}, k; x)$  for four values near  $k = \pm \frac{1}{2}$  and the corresponding solutions were computed as described in the previous paragraph.

Figure 13 shows two AS solutions corresponding to  $k = \pm 0.4999$ , both oscillatory on  $\mathbb{R}^-$ . It also shows two solutions, corresponding to  $k = \pm 0.5001$ , which exhibit a string of poles on  $\mathbb{R}^-$ .

<span id="page-18-1"></span>The HM solutions for  $k=\pm\frac{1}{2}$  are not shown in Fig. 13, but they continue to follow the shown asymptotic curves  $\pm\sqrt{-x/2}$  as  $x\to-\infty$ . The corresponding pole fields are not shown either, but they are similar to the pole fields shown in Fig. 12. When |k| is just less than the critical value  $\frac{1}{2}$ , the poles are located just off the real axis, leading to the oscillations in u(x) on  $\mathbb{R}^-$ . At the critical  $|k|=\frac{1}{2}$  these poles have cleared to  $-\infty$ , to reappear precisely on the real axis when |k| just exceeds  $\frac{1}{2}$ .

## 4.3 Quasi-Hastings–McLeod and Quasi-Ablowitz–Segur Solutions

Consider again the pole counting diagrams displayed in Fig. 11 and the first three subplots of Fig. 9. These are typical diagrams for the case  $0 \le \alpha < \frac{1}{2}$ , with the two HM points located on two opposite edges of the  $0^-$  region, and the connecting curve representing the AS solutions. As  $\alpha$  increases from 0 to  $\frac{1}{2}$ , the secondary point moves closer to the primary one. They are both located on the  $0^+$  curve, at opposite edges of the  $0^-$  region. The two points coalesce at  $\alpha = \frac{1}{2}$ , at which time the  $0^-$  region has shrunk to zero width in that neighborhood.

When  $\alpha$  is increased beyond  $\frac{1}{2}$ , the sHM point passes through the primary one, while continuing along the  $0^+$  curve but now entering a  $1^-$  region. (Visualize the

![](_page_18_Picture_11.jpeg)

![](_page_19_Figure_2.jpeg)

<span id="page-19-0"></span>Fig. 14 The curves represent initial conditions that correspond to the quasi-Ablowitz–Segur solutions. The endpoints of the curves represent the primary Hastings–McLeod solution, marked  $k_p$ , and the secondary, quasi-Hastings–McLeod solution, marked  $k_s$ ; cf. (21). The curve labels  $n^-$  and  $n^+$  denote, respectively, n poles on  $\mathbb{R}^-$  and n poles on  $\mathbb{R}^+$ , same as in Fig. 5

morphing of the third subplot into the fourth in Fig. 5.) The first subplot in Fig. 14 shows a typical case. Here  $\alpha=\frac{2}{3}$ , with  $k_p=\frac{1}{2}$  and  $k_s=-\frac{1}{2}$ , and the corresponding HM points are located on opposite edges of the  $1^-$  region. Where the pHM point is located, the  $1^-$  region is bounded by a  $0^-$  curve, so the pHM solution remains pole free. The sHM solution has picked up a solitary pole on  $\mathbb{R}^-$ , however. Aside from the existence of this pole, the solution has the same asymptotic properties as the regular sHM solution (as defined by (11), with the plus sign choice). To emphasize the existence of at least one pole on  $\mathbb{R}$ , we call these qHM solutions. Before showing graphs of such solutions, let us consider in Fig. 14 what happens as  $\alpha$  is increased further.

At approximately  $\alpha=0.77$  (not shown) the qHM point approaches  $(u(0),u'(0))=(-\infty,+\infty)$  in the pole counting diagram. Similar to the situation described in Fig. 6, this means a pole of residue -1 will be passing through the origin and will enter the (u(0),u'(0))-plane from the direction  $(+\infty,+\infty)$ , this time on a  $1^+$  curve. The second subplot in Fig. 14 shows a typical case, corresponding to  $\alpha=1$ . The third and fourth subplots in Fig. 14 follow similar patterns, but more poles appear on  $\mathbb{R}$ . These diagrams suggest that if  $\frac{1}{2} < \alpha < \frac{3}{2}$  then the qHM solution has one pole on  $\mathbb{R}$ , if  $\frac{3}{2} < \alpha < \frac{5}{2}$  then it has two poles, etc. At the half-integer  $\alpha$  values the pHM and qHM solutions coincide, so no poles. Our pole field solver enables us to find the locations and residues of these poles, and the results are displayed in Fig. 15.

A few qHM solutions and their pole fields are shown in Fig. 16, which is an extension of Fig. 10. As these figures illustrate, increasing  $\alpha$  causes the pHM pole field wedges to move smoothly to the right and away from the real axis. The pole fields for

![](_page_19_Picture_7.jpeg)

![](_page_20_Figure_2.jpeg)

<span id="page-20-0"></span>Fig. 15 Pole locations on  $\mathbb{R}$  for the quasi-Hastings–McLeod solutions for various  $\alpha$ . The *dash-dot* and *solid curves* represent residues -1 and +1, respectively. The *horizontal lines* represent half-integer values of  $\alpha$ 

the sHM solution agree with these when  $\alpha=\frac{1}{2},\frac{3}{2},\frac{5}{2}$ , etc., but Figs. 14–16 indicate a fundamentally different process at intermediate  $\alpha$  values. An example of this is seen in Fig. 17, where  $\alpha$  is increased from  $\frac{5}{2}$  to  $\frac{7}{2}$ . A band made up of five curves of poles enters from the left, and joins the pole wedges. As  $\alpha$  approaches  $\frac{7}{2}$ , seven curves of poles dislodge from the wedges on the left sides, and exit together towards minus infinity (to be compared to the dislodging of the single curve of poles seen in one of the subplots of Fig. 10). During this process, the HM wedges have actually moved slightly towards the real axis and to the left, but the end result gives the opposite impression, because its two leftmost rows of poles got separated off and transported out to minus infinity. Throughout this process, the pole configuration along the real axis remains reminiscent of the one for the rational solution at  $\alpha=3$ , sliding in from the left and then returning out to the left again.

As noted by [19], the characteristics of the qHM solutions can be partially explained by the Bäcklund transformation (4). Consider the second boundary condition in (11): if  $u(\alpha; x) \sim -\alpha/x$ ,  $x \to +\infty$ , is substituted into the right-hand side of (4) one gets  $u(\alpha+1;x) \sim -(\alpha+1)/x$ , i.e., this boundary condition is preserved by the transformation. The first boundary condition in (11), extended by one term, reads  $u(\alpha;x) \sim \pm \sqrt{-x/2} + \alpha/(2x) + \cdots$ ,  $x \to -\infty$ . Then (4) produces  $u(\alpha+1;x) \sim \mp \sqrt{-x/2}$ , i.e., this boundary condition is preserved but for a sign change. That is, the upper solution  $u(\alpha;x)$  maps to the lower solution  $u(\alpha+1;x)$  and vice versa. Next, check the sign of the denominator in (4), say  $v = 2u' + 2u^2 + x$ , in the limits  $x \to \pm \infty$ . When the plus sign is considered in (11), v is negative in both x limits but when the negative sign is considered v has opposite signs. This sign change

![](_page_20_Picture_6.jpeg)

![](_page_21_Figure_2.jpeg)

<span id="page-21-0"></span>Fig. 16 Same as Fig. 10, but now including the quasi-Hastings–McLeod solutions. These solutions have the same asymptotics as the regular Hastings–McLeod solutions as  $x \to \pm \infty$ , but they have one or more poles on  $\mathbb R$  (Color figure online)

in the denominator of (4) makes it plausible that the secondary HM solution picks up an additional pole on the real axis when the value of  $\alpha$  is increased by one. A proof that one and only one pole is added in this manner will require a deeper analysis, however.

The curves in Fig. 14 represent values of k between  $k_p$  and  $k_s$ , and are associated with solutions similar to the regular AS solutions (oscillatory as  $x \to -\infty$ , smoothly decaying as  $x \to +\infty$ ), except for the existence of poles in between. We refer to these as qAS solutions (previously noted in [2] in the special case of integer  $\alpha$ ). They coincide with the pHM solution when  $\alpha$  is half-integer. Otherwise, like the qHM solutions they have  $[\alpha + \frac{1}{2}]$  poles along  $\mathbb{R}$  (where  $[\cdot]$  denotes the integer part).

In the  $\alpha=1$  case in Fig. 14, the k=0 point is located at  $(\pm\infty, +\infty)$  in the diagram, as indicated by the dash-dot line segment. This is the rational solution  $u_1(z)$ , defined in (6) and displayed in a pole counting diagram in Fig. 6. Similarly, in the  $\alpha=2$  subplot, the rational solution  $u_2(z)$  is located at  $(\pm\infty, -\infty)$ , again shown as a dash-dot line segment. In the  $\alpha=3$  subplot, the point k=0 corresponds to the rational solution  $u_3(z)$ .

![](_page_21_Picture_7.jpeg)

![](_page_22_Figure_2.jpeg)

<span id="page-22-0"></span>Fig. 17 Pole dynamics of secondary, quasi-Hastings–McLeod solutions as  $\alpha$  is varied across  $[\frac{5}{2}, \frac{7}{2}]$  (Color figure online)

Figure 18 shows some of the qAS solutions in the special case of  $\alpha = \frac{2}{3}$  (cf. the first subplot of Fig. 14). In this case, the two HM solutions are obtained for  $k = \pm \frac{1}{2}$ , with the qAS solutions arising for the in-between values of k. Subplots 3–7 of Fig. 18 show qAS solutions, all decaying on  $\mathbb{R}^+$ , oscillating on  $\mathbb{R}^-$ , and featuring a single pole on  $\mathbb{R}^-$ . The central subplot shows the k=0 solution. Subplots 2 and 8 correspond to the sHM and pHM cases, respectively, with associated pole fields similar to those shown in the top row of Fig. 16 (for  $\alpha=1$  rather than for  $\alpha=\frac{2}{3}$ ). The main feature of the pole dynamics between these two HM cases is that a pole field enters from the left, and then exits again to the left, bringing with it the pole that was originally located on the real axis.

#### 4.4 k = 0 Solutions

The k=0 solutions include as special cases the rational solutions when  $\alpha$  is an integer and the HM solutions when  $\alpha$  is a half-integer. A sequence of these solutions with  $\alpha$  varying across [0,2] is displayed in Fig. 19, the rightmost column of which shows the rational and HM solutions. The  $\alpha=0$  case corresponds to the trivial u=0 solution (not shown). Because these k=0 solutions are special cases of qAS solutions, they show the familiar pattern of oscillations on  $\mathbb{R}^-$ , decay on  $\mathbb{R}^+$ , and  $[\alpha+\frac{1}{2}]$  poles along  $\mathbb{R}$ . (The exceptions are the  $\alpha$  integer cases, where the oscillations are not present, and the  $\alpha$  half-integer cases, where neither oscillations nor poles are present.)

We continue our display of k = 0 solutions in Fig. 20, but now showing pole fields instead as  $\alpha$  varies across [2, 3]. All these solutions are pole free in large sectors of

![](_page_22_Picture_8.jpeg)

![](_page_23_Figure_2.jpeg)

<span id="page-23-2"></span>**Fig. 18** Solutions of  $P_{II}$  ( $\alpha = \frac{2}{3}$ ) corresponding to nine values of k in the asymptotic boundary condition  $u(x) \sim B(x) + k \operatorname{Ai}(x)$ ,  $x \to +\infty$ ; cf. (17). Values  $|k| < \frac{1}{2}$  correspond to quasi-Ablowitz–Segur solutions,  $k = \frac{1}{2}$  corresponds to the regular Hastings–McLeod solution, and  $k = -\frac{1}{2}$  to a quasi-Hastings–McLeod solution. Values  $|k| > \frac{1}{2}$  feature an infinity of poles on  $\mathbb{R}^-$  (All solutions are displayed on  $[-8, 8] \times [-10, 10]$ )

<span id="page-23-1"></span><span id="page-23-0"></span>the right half-plane, and therefore of tronquée-type. Observe also how the group of four poles associated with the rational solution  $u_2(z)$  move along the negative real axis as  $\alpha$  is increased, and disappears to  $-\infty$  as  $\alpha \to \frac{5}{2}$ . As  $\alpha$  is increased further, the group of nine poles associated with the rational solution  $u_3(z)$  enters from the left, and settles in a symmetric configuration around the origin as  $\alpha \to 3$ .

#### 5 Additional Solution Illustrations

## 5.1 Pole Counting Diagrams: Domain Edges and Poles Crossing the Origin

When we follow a family of solutions for which a pole passes the origin, we noted in Sect. 3.1 (cf. Fig. 6) how the corresponding initial conditions at z=0 would be affected. These were seen to follow trajectories in the (u(0), u'(0))-plane, which connect at either  $u'(0) = +\infty$  or  $u'(0) = -\infty$  according to if the pole passing the origin has residue -1 or +1. If finite, both the counts  $n^+$  and  $n^-$  change by one. When counting poles on  $\mathbb{R}^+$ , the regions with finite count  $n^+$  always form curves (of zero width) in the pole counting diagrams.

![](_page_23_Picture_8.jpeg)

![](_page_24_Figure_2.jpeg)

<span id="page-24-1"></span>Fig. 19 Solutions (17) with k = 0 as  $\alpha$  is varied (in nonuniform increments) across [0, 2]. In the rightmost column these solutions reduce to the rational solutions in the case  $\alpha$  integer, and to the Hastings–McLeod solutions in the case  $\alpha$  half-integer (All solutions are displayed on  $[-10, 6] \times [-4, 4]$ )

Turning to similar counts on  $\mathbb{R}^-$ , Fig. 5 showed that finite counts  $n^-$  can occur also within entire regions, not only along curves. Additionally, the counts along the edges of these regions need not agree with counts for within the regions. As an example, the left subplot of Fig. 21 shows all these different counts for  $\mathbb{R}^-$  and the associated connection lines in the case of  $\alpha = 3$ . The right subplot of that figure shows corresponding pole counts on  $\mathbb{R}^+$ .

The  $\alpha$  half-integer cases reveal still an additional feature in that some of the  $n^-$  regions join each other side-to-side, whereas others lose their width and turn into  $n^-$  curves. Figure 22 gives the different  $n^-$  counts and connection lines in such a case  $(\alpha = \frac{3}{7}; \text{ cf. also the last two subplots in Fig. 7}).$ 

<span id="page-24-0"></span>In all cases, there is (at least) one  $0^-$  curve between bottom left and top right in the pole counting diagrams (either as a domain boundary or, in the half-integer cases, partly a domain boundary) and also a single  $0^+$  curve between bottom right and top left. As discussed in more detail already, this ensures the existence of Hastings–McLeod solution(s) for all  $\alpha$ .

## 5.2 Relating u(0; z) and $u(\frac{1}{2}; z)$ Solutions

Here we consider the connection between  $\alpha = 0$  and  $\alpha = \frac{1}{2}$  solutions as described by (5). It follows from that representation that both poles and zeros of u(0; z) can generate poles of  $u(\frac{1}{2}; z)$ .

![](_page_24_Picture_9.jpeg)

![](_page_25_Figure_2.jpeg)

<span id="page-25-0"></span>Fig. 20 Pole dynamics for the k = 0 solutions as  $\alpha$  is varied across [2, 3] (Color figure online)

For example, if  $z \rightarrow z_0$  with  $z_0$  the location of a pole, then

$$u(0;z) \sim \frac{\pm 1}{z - z_0} \implies \frac{u'(0;z)}{u(0;z)} \sim \frac{-1}{z - z_0} \implies u\left(\frac{1}{2};z\right) \sim \frac{1}{z - (-2^{1/3}z_0)}.$$

On the other hand, if  $z_0$  is a zero of single multiplicity (as all zeros of  $P_{II}$  with  $\alpha = 0$  are) then for some nonzero constant c

$$u(0;z) \sim c(z-z_0) \implies \frac{u'(0;z)}{u(0;z)} \sim \frac{1}{z-z_0} \implies u\left(\frac{1}{2};z\right) \sim \frac{-1}{z-(-2^{1/3}z_0)}.$$

Therefore, a pole or a zero of u(0; z) at  $z = z_0$  generates a pole at  $z = -2^{1/3}z_0$  of  $u(\frac{1}{2}; z)$ . When generated by a pole (of either residue) of u(0; z), the pole of  $u(\frac{1}{2}; z)$  has a residue +1, while if generated by a zero of u(0; z) the pole of  $u(\frac{1}{2}; z)$  has residue -1.

![](_page_26_Figure_2.jpeg)

<span id="page-26-0"></span>Fig. 21 Pole counting diagrams in the case  $\alpha = 3$ , which show (in the *dash-dot* connections) how pole counts increase/decrease when a pole crosses the origin. The *left subplot* refers to  $\mathbb{R}^-$ , and illustrates how pole counts on edges of the  $n^-$  domains are related to pole counts in the interior. The *right subplot* refers to  $\mathbb{R}^+$  (For clarity, the aspect ratio of the diagrams here is not the same as in Fig. 5)

To see where these solutions fit into the pole counting diagrams, define for simplicity  $(u(0;0),u'(0;0))=(v_0,v_1)$  and  $(u(\frac{1}{2};0),u'(\frac{1}{2};0))=(w_0,w_1)$ . Then it follows from (5) that

<span id="page-26-2"></span><span id="page-26-1"></span>
$$w_0 = (2^{-1/3})\frac{v_1}{v_0}, \qquad w_1 = (2^{-2/3})\frac{v_1^2 - 2v_0^4}{v_0^2},$$
 (24)

with inverse formula

$$v_0 = \pm (2^{-1/6})\sqrt{w_0^2 - w_1}, \qquad v_1 = \pm (2^{1/6})w_0\sqrt{w_0^2 - w_1}.$$
 (25)

One such relationship between a specific pair of u(0; x) and  $u(\frac{1}{2}; x)$  solutions is shown in Fig. 23. To plot these solutions, we chose initial conditions for the  $u(\frac{1}{2}; x)$  solution where the  $3^+$  curve in the pole counting diagram intersects the horizontal axis. This point has coordinates approximately  $(w_0, w_1) = (1.670027, 0)$ , as indicated by the point marked (a) in the second subplot in Fig. 24. Using the mapping (25), the initial conditions for the u(0; x) solution are then computed to be approximately  $(v_0, v_1) = (1.487825, 3.130535)$ , or the negatives of these values. These points, marked (a) in the first pole counting diagram of Fig. 24, are located on  $1^-$  edges of the  $1^-$  region. These results are consistent with the fact that the u(0; x)

![](_page_26_Picture_9.jpeg)

![](_page_27_Figure_2.jpeg)

<span id="page-27-0"></span>**Fig. 22** Equivalent illustration to Fig. 21 but for  $\alpha = \frac{3}{2}$ 

![](_page_27_Figure_4.jpeg)

<span id="page-27-1"></span>Fig. 23 The two zeros and one pole of u(0; x) on  $\mathbb{R}^-$  (*left subplot*) generate two poles of residue +1 and one pole of residue -1 of  $u(\frac{1}{2}; x)$  on  $\mathbb{R}^+$  (*right subplot*)

and  $u(\frac{1}{2}; x)$  solutions shown in Fig. 23 have one pole on  $\mathbb{R}^-$  and three poles on  $\mathbb{R}^+$ , respectively.

Some additional points of interest, including Hastings–McLeod points, are also marked in Fig. 24.

Another insight that can be gained from (5) relates to the fact that, for  $\alpha$  half-integer, the different  $n^-$  regions inside the Airy boundaries (cf. Figs. 7 and 8) feature

![](_page_28_Figure_2.jpeg)

<span id="page-28-1"></span>**Fig. 24** The points in the *left* and *right subplots* are connected via the mappings (24) and (25). The *points* marked (a) define the initial conditions for the two solutions shown in Fig. 23. The *points* marked (b) correspond to Hastings–McLeod solutions. The *points* marked (c) show a  $1^+$  solution of u(0; x) that generates a  $1^-$  solution of  $u(\frac{1}{2}; x)$ 

<span id="page-28-0"></span>no in-between gaps, and outside these they have collapsed to curves (of zero width). We focus on the  $\alpha=\frac{1}{2}$  case, seen in Figs. 5 and 24 (right subplot). According to (24) and (25), the exterior of the Airy parabola in the  $\alpha=\frac{1}{2}$  diagram (i.e.,  $w_1< w_0^2$ ) maps to the entire  $(v_0,v_1)$  plane in the  $\alpha=0$  diagram. The collapsed  $n^-$  regions are consistent with the fact that the  $\alpha=0$  diagram has no  $n^+$  regions (i.e., there are  $n^+$  curves only). The interior of the Airy parabola in the  $\alpha=\frac{1}{2}$  diagram, on the other hand, will map to the entire  $(v_0,v_1)$  plane of the imaginary (or modified)  $P_{II}$  equation,  $u''=-2u^3+zu$ , which is obtained from (2) by  $u\to iu$ . The absence of gaps is consistent with Theorem 3.1 in [9], which asserts that real solutions of this equation lacks poles on  $\mathbb R$  and has at most a finite number of zeros along  $\mathbb R^+$ .

## 5.3 Re-visiting the $u(\frac{3}{2}; z)$ Solution

We return to the pole counting diagrams for the case  $\alpha = \frac{3}{2}$  shown in Fig. 7, and the fact that the HM and Airy ( $\theta = 0$ ) points are located in close proximity of each other. To examine the solution features in this region, consider the points marked (a)–(f) in the third subplot of that figure. Pole fields corresponding to these six points are displayed in Fig. 25, showing the rich diversity of pole dynamics that can arise

![](_page_28_Picture_7.jpeg)

![](_page_29_Figure_2.jpeg)

<span id="page-29-0"></span>**Fig. 25** Pole fields of PII *(α* <sup>=</sup> <sup>3</sup> <sup>2</sup> *)* generated by the initial conditions labeled (**a**)–(**f**) in the third subplot of Fig. [7](#page-11-0) (Color figure online)

already from a tiny region in parameter space. Although these six pole fields superficially might look similar, there are several significant differences. In particular, we can note:

- Between subplots (a), (b), and between (d), (e), the leftmost, the top right and bottom right pole fields exit simultaneously and then return with changed alignment, featuring at one instant the Airy-type solution.
- Between subplots (a), (d), and again (b), (e) and (c), (f), both the far leftmost pole field and the central triple band exit and then return from minus infinity. However, only the triple band undergoes a change in vertical alignment (similar to what was seen in subplots 4–6 of Fig. [20](#page-25-0) when varying *α* for *k* = 0 solutions).
- Changes in the pole field to the right occur between subplots (b), (c) and also between (e), (f) (with no other significant changes).

We have described in Sect. [2.3](#page-4-3) how the poles move when we follow the Airy curve through its period, including for when *θ* passes zero. When following the curve marked 0<sup>+</sup> through the point marked HM in the third subplot of Fig. [7](#page-11-0), the upper and lower HM pole wedges remain in place, the rightmost pole field remains at plus infinity, and the leftmost pole field as well as the triple pole band exits to and then returns from minus infinity.

![](_page_29_Picture_9.jpeg)

#### <span id="page-30-7"></span>6 Conclusions

With the aid of the pole field solver of [14] and pole counting diagrams, we have surveyed the solution space of the  $P_{\rm II}$  equation ( $\alpha \ge 0$ ; u(z) real for z real). Previously described solution types have been revisited and extended throughout the complex plane.

New solutions (or at least solutions that have not been identified in any detail in the literature) have been described. These include what we call the secondary Hastings—McLeod solutions, which are nonoscillatory and pole free on the entire real axis when  $0 \le \alpha < \frac{1}{2}$ , but not monotone like the primary Hastings—McLeod solutions. The solutions corresponding to k = 0 in (17) have likewise not been computed before.

In addition to the discussions and illustrations in the present work, animations of solutions (as the parameter and ODE initial conditions are varied) can be found on the web site [34].

<span id="page-30-9"></span><span id="page-30-4"></span>Acknowledgements Financial support for this work was provided by NSF grant DMS-0914647 (first author) and the National Research Foundation in South Africa (second author), as well as the National Institute for Theoretical Physics (NITheP) in South Africa. Communications with Andrew Bassom, Peter Clarkson, Percy Deift, Alexander Its, Andrei Kapaev, Jonah Reeger and Harvey Segur are also acknowledged. The workshop "Numerical solution of the Painlevé Equations", held in May 2010 at the International Center for the Mathematical Sciences (ICMS), in Edinburgh, stimulated [14] and the present study.

#### <span id="page-30-15"></span><span id="page-30-14"></span><span id="page-30-11"></span><span id="page-30-5"></span>References

- M.J. Ablowitz, H. Segur, Exact linearization of a Painlevé transcendent, *Phys. Rev. Lett.* 38, 1103– 1106 (1977).
- <span id="page-30-12"></span>2. M.J. Ablowitz, H. Segur, Solitons and the inverse scattering transform (SIAM, Philadelphia, 1981).
- <span id="page-30-1"></span> A.P. Bassom, P.A. Clarkson, C.K. Law, J.B. McLeod, Application of uniform asymptotics to the second Painlevé transcendent, Arch. Ration. Mech. Anal. 143, 241–271 (1998).
- <span id="page-30-6"></span> M. Bertola, On the location of poles for the Ablowitz–Segur family of solutions of the second Painlevé equation, Nonlinearity 25, 1179–1185 (2012).
- <span id="page-30-2"></span> P. Boutroux, Remarques sur les singularités transcendantes des fonctions de deux variables, Bull. Soc. Math. Fr. 39, 296–304 (1911).
- <span id="page-30-10"></span> P. Boutroux, Recherches sur les transcendantes de M. Painlevé et l'étude asymptotique des équations différentielles du second ordre (suite), Ann. Sci. École Norm. Super. (3) 31, 99–159 (1914).
- <span id="page-30-13"></span> T. Claeys, A.B.J. Kuijlaars, M. Vanlessen, Multi-critical unitary random matrix ensembles and the general Painlevé II equation, *Ann. Math.* (2) 168, 601–641 (2008).
- <span id="page-30-3"></span>8. P.A. Clarkson, Painlevé equations—nonlinear special functions, in *Orthogonal polynomials and special functions*. Lecture notes in mathematics, vol. 1883 (Springer, Berlin, 2006), pp. 331–411.
- P.A. Clarkson, Asymptotics of the second Painlevé equation, in Special functions and orthogonal polynomials. Contemporary mathematics, vol. 471 (Am. Math. Soc., Providence, 2008), pp. 69–83.
- <span id="page-30-8"></span><span id="page-30-0"></span>P.A. Clarkson, Painlevé transcendents, in NIST handbook of mathematical functions (U.S. Dept. Commerce, Washington, 2010), pp. 723–740.
- P.A. Clarkson, E.L. Mansfield, The second Painlevé equation, its hierarchy and associated special polynomials, *Nonlinearity* 16, R1–R26 (2003).
- P.A. Clarkson, J.B. McLeod, A connection formula for the second Painlevé transcendent, Arch. Ration. Mech. Anal. 103, 97–138 (1988).
- A.S. Fokas, A.R. Its, A.A. Kapaev, V.Y. Novokshenov, Painlevé transcendents: the Riemann–Hilbert approach (Am. Math. Soc., Providence, 2006).
- B. Fornberg, J.A.C. Weideman, A numerical methodology for the Painlevé equations, J. Comput. Phys. 230, 5957–5973 (2011).
- 15. V.I. Gromak, Solutions of the second Painlevé equation, Differ. Urayn. 18, 753-763, 914-915 (1982).

![](_page_30_Picture_23.jpeg)

- <span id="page-31-19"></span><span id="page-31-18"></span><span id="page-31-15"></span><span id="page-31-14"></span><span id="page-31-6"></span><span id="page-31-2"></span>16. V.I. Gromak, I. Laine, S. Shimomura, *Painlevé differential equations in the complex plane* (de Gruyter, Berlin, 2002).
- <span id="page-31-8"></span><span id="page-31-7"></span>17. S.P. Hastings, J.B. McLeod, A boundary value problem associated with the second Painlevé transcendent and the Korteweg–de Vries equation, *Arch. Ration. Mech. Anal.* **73**, 31–51 (1980).
- 18. A.R. Its, A.A. Kapaev, Quasi-linear Stokes phenomenon for the second Painlevé transcendent, *Nonlinearity* **16**, 363–386 (2003).
- <span id="page-31-16"></span><span id="page-31-3"></span>19. A.A. Kapaev, Private communication.
- <span id="page-31-12"></span>20. A.A. Kapaev, Global asymptotics of the second Painlevé transcendent, *Phys. Lett. A* **167**, 356–362 (1992).
- <span id="page-31-11"></span>21. A.V. Kashevarov, The second Painlevé equation in electrostatic probe theory: numerical solutions, *Comput. Math. Math. Phys.* **38**, 950–958 (1998).
- <span id="page-31-0"></span>22. A.V. Kashevarov, The second Painlevé equation in the electrostatic probe theory: numerical solutions for the partial absorption of charged particles by the surface, *Tech. Phys.* **49**, 3–9 (2004).
- <span id="page-31-1"></span>23. A.V. Kitaev, Symmetric solutions for the first and the second Painlevé equation, *J. Math. Sci.* **73**, 494–499 (1995).
- 24. B.M. McCoy, S. Tang, Connection formulae for Painlevé functions, *Physica D* **18**, 190–196 (1986).
- <span id="page-31-10"></span><span id="page-31-4"></span>25. J.W. Miles, On the second Painlevé transcendent, *Proc. R. Soc. Lond. Ser. A* **361**, 277–291 (1978).
- 26. V.Y. Novokshenov, Padé approximations for Painlevé I and II transcendents, *Theor. Math. Phys.* **159**, 853–862 (2009).
- <span id="page-31-17"></span>27. S. Olver, Numerical solution of Riemann–Hilbert problems: Painlevé II, *Found. Comput. Math.* **11**, 153–179 (2011).
- <span id="page-31-5"></span>28. P. Painlevé, Mémoire sur les équations différentielles dont l'intégrale générale est uniforme, *Bull. Soc. Math. Fr.* **28**, 201–261 (1900).
- <span id="page-31-13"></span><span id="page-31-9"></span>29. P. Painlevé, Sur les équations différentielles du second ordre et d'ordre supérieur dont l'intégrale générale est uniforme, *Acta Math.* **25**, 1–85 (1902).
- 30. J.A. Reeger, B. Fornberg, Painlevé IV with both parameters zero: a numerical study, *Stud. Appl. Math.* **130**, 108–133 (2013).
- 31. R.R. Rosales, The similarity solution for the Korteweg–de Vries equation and the related Painlevé transcendent, *Proc. R. Soc. Lond. Ser. A* **361**, 265–275 (1978).
- 32. H. Segur, M.J. Ablowitz, Asymptotic solutions of nonlinear evolution equations and a Painlevé transcendent, *Physica D* **3**, 165–184 (1981).
- 33. C.A. Tracy, H. Widom, Painlevé functions in statistical physics, *Publ. Res. Inst. Math. Sci.* **47**, 361– 374 (2011).
- 34. J.A.C. Weideman, [dip.sun.ac.za/~weideman/PAINLEVE/](http://dip.sun.ac.za/~weideman/PAINLEVE/).
- 35. I.M. Willers, A new integration algorithm for ordinary differential equations based on continued fraction approximations, *Commun. ACM* **17**, 504–508 (1974).

![](_page_31_Picture_22.jpeg)