# Heun Functions and some of their applications in Physics <sup>∗</sup>

M. Horta¸csu †

Mimar Sinan Fine Arts University, Department of Physics, Istanbul, Turkey

February 25, 2019

#### Abstract

Most of the theoretical physics known today is described by using a small number of differential equations. For linear systems, different forms of the hypergeometric or the confluent hypergeometric equations often suffice to describe the system studied. These equations have power series solutions with simple relations between consecutive coefficients and/ or can be represented in terms of simple integral transforms. If the problem is nonlinear, one often uses one form of the Painlev´e equations. There are important examples, however, where one has to use higher order equations. Heun equation is one of these examples, which recently is often encountered in problems in general relativity and astrophysics. Its special and confluent forms take names as Mathieu, Lam´e and Coulomb spheroidal equations. For these equations whenever a power series solution is written, instead of a two way recursion relation between the coefficients in the series, we find one between three or four different ones. An integral transform solution using simpler functions also is not obtainable.The use of this equation in physics and mathematical literature exploded in the later years, more than doubling the number of papers with these solutions in the last decade, compared to time period since this equation was introduced in 1889 up to 2008. We use SCI data to conclude this statement, which is not precise, but in the correct ballpark. Here this equation will be introduced and examples for its use, especially in general relativity literature will be given.

<sup>∗</sup>This paper is a revised and many times updated version of the conference talk by the same author, published in Proceedings of the 13th Regional Conference on Mathematical Physics, Antalya, Turkey,27-31 October 2010, edited by U˘gur Camcı and Ibrahim Semiz, pp. 23-39, World Scientific (2013).

<sup>†</sup>E-mail: hortacsu@itu.edu.tr

## 1 Introduction

Most of the theoretical physics known today is described by using a small number of differential equations. If we study only linear systems, different forms of the hypergeometric or the confluent hypergeometric equations often suffice to describe the system studied. These equations have power series solutions with simple relations between consecutive coefficients and/ or can be represented in terms of simple integral transforms. If the problem is described in terms of nonlinear differential equations, then one often uses one form of the Painlev´e equations.

There are important examples, however, where one has to use higher order equations. Such an equation was proposed by Karl Heun in 1889 [\[1\]](#page-13-0). This equation and its confluent forms becomes indispensable in general relativity if one studies exact solutions of wave equations in the background of certain metrics. A well known example is the Kerr metric [\[2\]](#page-13-1). Although it is possible to solve the wave equations in the background of some metrics in terms of hypergeometric functions or its confluent forms, this is not possible for the much studied Kerr metric. If we also study even the trivially extended forms of some metrics by adding a flat dimension to the existing metric, we may have to solve the Heun equation to obtain the exact solution.

Here we will find introduce the Heun equation and its confluent forms and mention some of the properties of the Heun equation. Then we will give some examples in physics, mainly in gravitational physics, where one can find many recent papers. This part is meant to be a survey of the work done in the field of General Relativity and Quantum Gravity concentrating in the last decades. In another section we will give an example where the Heun equation emerges from a trivial extension of a wave equation in the background of the Eguchi-Hanson instanton metric [\[3\]](#page-14-0). We will end with some concluding remarks.

## 2 Heun equation

Let us review some well known facts about second order differential equations. Differential equations are classified according to their singularity structure [\[4,](#page-14-1) [5\]](#page-14-2). If a differential equation has no singularities over the full complex plane, it can only be a constant. Singularities are classified as regular singular and irregular singular points. If the coefficient of the first derivative has at most single poles, and the coefficient of the term without a derivative has at most double poles when the coefficient of the second derivative is unity, this second order differential equation has regular singularities, which gives us one regular solution while expanding around this singular point. In general the second solution has a pole or a branch point singularity. If the poles of these coefficients are higher, we have irregular singularities and the general solution has an essential singularity [\[6\]](#page-14-3).

As stated in Morse and Feshbach [\[4\]](#page-14-1) an example of a second order differential equation with one regular singular point is

$$\frac{d^2w}{dz^2} = 0. (1)$$

This equation has one solution which is constant. The second solution blows up at infinity. The differential equation

$$\frac{d^2w}{dz^2} + k^2w = 0. (2)$$

has one irregular singularity at infinity which gives an essential singularity at this point. The equation

$$z\frac{d^2w}{dz^2} + (1+a)\frac{dw}{dz} = 0. (3)$$

has two regular singular points, at zero and at infinity.

In physics an often used equation is the hypergeometric equation

$$z(1-z)\frac{d^2w}{dz^2} + [c - (1+a+b)z]\frac{dw}{dz} - abw = 0.$$
(4)

This equation has three regular singular points, at zero, one and infinity. Jacobi, Legendre, Gegenbauer, Tchebycheff equations are special forms of this equation. When the singular points at z=1 and z equals infinity are "coalesced" at infinity, we get the confluent hypergeometric equation

$$z\frac{d^2w}{dz^2} + (c - z)\frac{dw}{dz} - aw = 0. (5)$$

with an essential singularity at infinity and a regular singularity at zero. Bessel, Laguerre, Hermite equations can be reduced to this form.

An important property of all these equations is that they allow infinite series solutions about one of their regular singular points where a recursion relation can be found between two consecutive coefficients. This fact allows one to have an idea about the general properties of the solution, as the asymptotic behaviour at distant points, the radius of convergence of the series, etc.

A new equation was introduced in 1889 by Karl M. W. L. Heun [1]. This is an equation with four regular singular points at zero, one, an arbitrary point f between zero and one and infinity. This equation is discussed in the book edited by Ronveaux [7]. Most of the general information we give below is taken from this book. As discussed there, any equation with four regular singular points can be transformed to the equation given below:

$$\frac{d^2w}{dz^2} + \left[\frac{c}{z} + \frac{d}{z-1} + \frac{e}{z-f}\right] \frac{dw}{dz} - \frac{abz - q}{z(z-1)(z-f)} w = 0.$$
 (6)

There is a relation between the constants given as a + b + 1 = c + d + e. This relation is not related to the regularity of the singularity at infinity. It just gives the exponents of the term multiplying the series solution around infinity in terms of u = 1/z as a, b.

If we try to obtain a solution in terms of a power series, one can not get a recursion relation between two consecutive coefficients. We have a relation at least between three coefficients.

It is known that [8] any second order differential equation with n regular singular points has a family of  $2^{n-1}n!$  local solutions, which splits into 2n sets of

$$2^{n-2}(n-1)!$$

equivalent expressions, each set defining one of the two Frobenius solutions in the neighborhood of a singular point. The n! factor comes from permuting the n singular points and the  $2^{n-1}$  factor from negating exponent differences. Maier [8] gave the list of 192 local solutions for the Heun equation.

The set of transformations that can be applied to the a Fuchian equation with n singular points to generate alternative expressions for this equation has order  $2^{n-1}n!$  and acts on the parameter space of the equation. This group of transformations is isomorphic to the Coxeter group  $D_n$ . These transformations generate  $2^{n-2}(n-1)!$  solutions. For the Heun case n=4, and this group is isomorphic to  $D_4$ , a group of order 192. These transformations will be the combination of Mobius transformations and transformations which multiply the desired solution by powers.

It turns out that the Mobius group PGL(2,C), which takes x to  $\frac{(Ax+B)}{(Cx+D)}$ , for nonvanishing AD-BC, can be used where x takes values from the different singular points. For Heun equation with four regular singular points, this transformation takes each singular point to five other points, which have zeroes at the same value. These points are given below:

$$x, x/(x-1), x/f, x/(x-f), (1-f)x/(x-f), (f-1)x/f(x-1),$$

$$1-x, (x-1)/x, (x-1)/(x-f), (x-1)/(f-1), d(x-1)/(x-f), f(x-1)/(f-1)x,$$

$$1/x, 1/(1-x), f/x, f/(f-x), (f-1)/(x-1), (1-f)/(x-f),$$

$$(x-f)/x, (f-x)/a, (x-f)/(x-1), (f-x)/(f-1), (x-f)/f(x-1), (f-x)/(f-1)x.$$

Any one of these transformations maps three of the four points, 0,1,f, infinity, into 0,1, infinity, but generally changes the value of f, which takes one of the six possible values:  $f_1 = f$ ,  $f_2 = 1 - f$ ,  $f_3 = 1/f$ ,  $f_4 = 1/(1-f)$ ,  $f_5 = f/(f-1)$ ,  $f_6 = (f-1)/f$ . Each value is taken four times.

Just recall the Heun equation:

$$\frac{d^2w}{dx^2} + \left[\frac{c}{x} + \frac{d}{x-1} + \frac{e}{x-f}\right]\frac{dw}{dx} - \frac{abx - q}{x(x-1)(x-f)}w = 0,\tag{7}$$

written in terms of the real variable x. One writes the solution to the Heun equation in the form:

$$y(x) = x^{r}(x-1)^{s}(1-x/f)^{t}u(x).$$

This changes the form of the differential equation. For (i) r = 0 or 1 - c, (ii) s = 0 or 1 - d, (iii) t = 0 or 1 - e, however, the resulting equation has the Heun form. The values given above are the exponents at the singularities [9], [10].

Of course, the parameters of the equations change. For each such combination, say for r = 0, there are four possible values s and t can take, namely both equal to zero; s = 1 - d, t = 0; s = 0, t = 1 - d; s = 1 - d, t = 1 - e. Thus we get three more solutions for each solution. Another factor of six comes from the six different possible values f can take. In total for expansions around a single regular singular point, we have twenty four equivalent solutions, obtained by simply transforming the original equation.

The presence of two different indices for expansion around each singular point doubles the number of equivalent solutions, resulting in 48 solutions for expansions around each singular point. Four singular points multiplies this number by four giving the total of 192 local solutions.

It turns out that for infinite set of values of the parameter q, there are solutions which are analytic at 0 and at 1. These are called *Heun functions*, whereas those which are analytic only at one point are called *local Heun functions* [11].

For integer values of one of a, c-a, d-a, e-a, and for special finite values of q, solutions analytic at three singularities exist, the so called *Heun polynomials*. A special case is for a = -n, n = 0, 1, 2 and  $q_{n,m}, m = 0, 1, ..., n$ , where  $q_{n,m}$  are eigenvalues of a tridiagonal matrix, we get the solution as a polynomial of degree n, which is analytic at three singular points, 0,1 and f. [12].

"No example has been given of a solution of Heun's equation expressed in the form of a definite integral or contour integral involving only functions which are, in some sense, simpler" [13]. This statement does not exclude the possibility of having an infinite series of integrals with "simpler" integrands.

One can obtain different confluent forms of this equation. When we "coalesce" two regular singular points, we get the confluent Heun equation: The standard form of the confluent form equation is given as [14].

$$\frac{d^2w}{dz^2} + \left(\alpha + \frac{\gamma+1}{z-1} + \frac{\beta+1}{z}\right)\frac{dw}{dz} + \left(\frac{\nu}{z-1} + \frac{\mu}{z}\right)w = 0.$$
 (8)

with solution

$$\begin{split} HeunC(\alpha,\beta,\gamma,\delta,\eta,z). \\ \delta &= \mu + \nu - \alpha(\frac{\beta + \gamma + 2}{2}), \\ \eta &= \frac{\alpha(\beta + 1)}{2} - \mu - (\frac{\beta + \gamma + \beta\gamma}{2}). \end{split}$$

Another version of this equation can be written as

$$\frac{d}{dz}((z^2-1)\frac{dw}{dz}) + [-p^2(z^2-1) + 2p\beta z - \lambda - \frac{m^2 + s^2 + 2msz}{(z^2-1)}]w = 0.$$
(9)

Special forms of this equation are obtained in problems with two Coulombic centers,

$$\frac{d}{dz}((z^2 - 1)\frac{dw}{dz}) + [-p^2(z^2 - 1) + 2p\beta z - \lambda - \frac{m^2}{(z^2 - 1)}]w = 0,$$
(10)

whose special form, when b = 0, is the spheroidal equation,

$$\frac{d}{dz}((z^2 - 1)\frac{dw}{dz}) + [-p^2(z^2 - 1) - \lambda - \frac{m^2}{(z^2 - 1)}]w = 0.$$
(11)

Another form is the algebraic form of the Mathieu equation:

$$\frac{d}{dz}((z^2-1)\frac{dw}{dz}) + [-p^2(z^2-1) - \lambda - \frac{1}{4(z^2-1)}]w = 0.$$
(12)

If we coalesce two regular singular points pairwise, we obtain the double confluent form:

$$D^{2}w + (\alpha_{1}z + \frac{\alpha_{-1}}{z})Dw + \left[ (B_{1} + \frac{\alpha_{1}}{2})z + (B_{0} + \frac{\alpha_{1}\alpha_{-1}}{2}) + (B_{-1} - \frac{\alpha_{-1}}{2})\frac{1}{z} \right]w = 0.$$
(13)

Here  $D = z \frac{d}{dz}$ . We can reduce the new equation to the Mathieu equation, an equation with two irregular singularities at zero and at infinity if we reduce this equation to the form:

$$D^{2}y + (Bz^{2} + B_{0} + Bz^{-2})y = 0. (14)$$

Another form is the biconfluent form, where three regular singularities are coalesced. The result is an equation with a regular singularity at zero and an irregular singularity at infinity of higher order:

$$z^{2}\frac{d^{2}w}{dz^{2}} + z\frac{dw}{dz}w + (A_{0} + A_{1}z + A_{2}z^{2} + A_{3}z^{3} - z^{4})w = 0.$$
(15)

The anharmonic equation in three dimensions can be reduced to this equation:

$$\frac{d^2w}{dz^2} + (E - \frac{\nu}{r^2} - \mu r^2 - \lambda r^4 - \eta r^6)w = 0.$$
 (16)

In the triconfluent case, all regular singular points are "coalesced" at infinity which gives the equation below:

$$\frac{d^2w}{dz^2} + (A_0 + A_1z + A_2z^2 - \frac{9}{4}z^4)w = 0. (17)$$

These different forms are used in different problems in physics.

## 3 Some Examples of the Heun equation in Physical Applications

In SCI we found about one hundred thirty papers when Heun functions were searched in the summer of 2010. Now, at the end of April 2018, the number exceeded 330. The number of published articles in SCI more than doubled in the last eight years. More than three fourths of these papers were published in the last ten years. The rest of the papers were published between 1990 and 2005, except a single paper in 1986 [15]. These numbers may differ depending on the institution where one uses the SCI, since different universities in Turkey start their search from different dates. We think we are still in the correct ball park. This shows that although the Heun equation was found in 1889, it was largely neglected in the physics literature until recently. Earlier papers on this topic are mostly articles in mathematics journals. If one looks for books on this topic, published before the year 2000, one finds out the list of books is not very long. There is a book edited by A.Ronveaux, which is a collection of papers presented in the "Centennial Workshop on Heun's Equations: Theory and Application. Sept.3-8 1989, Schloss Ringberg". It was published by the Oxford University Press in 1995 by the title Heun's Differential Equations [7]. There are two books on functions which are special cases of the Heun Equation: Mathieusche Funktionen und Sphaeroidfunktionen mit anwendungen auf physikalische und technische Probleme by Joseph Meixner and Friedrich Wilhelm Schaefke, published by Springer Verlag in 1954 [16] and a Dover reprint of a book first published in 1946, Theory and Applications of Mathieu Functions by N.W. McLachlan in 1963 [17]. Classical mathematical physics books, such as Morse and Feshbach [4], Whittaker and Watson [18] or the Batemann Manuscript [19], have sections or chapters on the special forms of the Heun equation like Mathieu, Lamé or spheroidal functions. Some papers on different mathematical properties of these functions can be found in references [20]- [26].

A reason why more physicists are interested in the Heun equation recently may be, perhaps, a demonstration of the fact that we do not have simple problems in theoretical physics anymore. Mathematical physicists have to tackle more difficult problems, either with more difficult metrics or in higher dimensions. Both of these extensions may necessitate the use of the Heun functions among the solutions. We can give the Eguchi-Hanson case as an example. The wave equation for the scalar particle in the background of the Eguchi-Hanson metric [3] in four dimensions has hypergeometric functions as solutions [27], whereas the Nutku helicoid [28, 29] metric, the next higher one, gives us Mathieu functions [30], a member of the Heun function set, if the method of separation of variables is used to get a solution. We also find that the scalar particle, in the background of the Eguchi-Hanson metric, trivially extended to five dimensions gives Heun type solutions. [31].

Note that the problem does not need to be very complicated to work with these equations. We encounter Mathieu functions if we consider two dimensional problems with elliptical shapes [32]. Let us use  $x = \frac{1}{2}a\cosh\mu\cos\theta$ ,  $y = \frac{1}{2}a\sinh\mu\sin\theta$ , where a is the distance from the origin to the focal point. Then the Helmholtz equation can be written as

$$\partial_{\mu\mu}\psi + \partial_{\theta\theta}\psi + \frac{1}{4}a^2k^2[\cosh^2\mu - \cos^2\theta]\psi = 0$$
(18)

which separates into two equations

$$\frac{d^2H}{d\theta^2} + (b - h^2\cos^2\theta)H = 0,\tag{19}$$

$$-\frac{d^2M}{d\mu^2} + (b - h^2 \cosh^2 \mu)M = 0.$$
 (20)

The solutions to these two equations can be represented as Mathieu and modified Mathieu functions.

If we combine different inverse powers of r, starting from first up to the fourth, or if we combine the quadratic potentials with inverse even powers of two, four and six, we see that the solution of the Schrodinger equation involves Heun functions [33]. Solution to symmetric double Morse potentials also needs these functions, like  $V(x) = B^2/4\sinh 2x - (s+1/2)B\cosh x$  where s = (0, 1/2, 1, ...) [33]. Similar problems are treated in references [34], [35] and [36]

o In atomic physics further problems such as separated double wells, Stark effect, hydrogen molecule ion use these functions. Physics problems which end up with these equations are given in the book by S.Y. Slavyanov and S. Lay [37]. Here we see that even the Stark effect, hydrogen atom in the presence of an external electric field, gives rise to this equation. As described in page 166 of Slavyanov's book, cited above (original reference is Epstein [38], also treated by S.Yu Slavyanov [39]), when all the relevant constants, namely Planck constant over  $2\pi$ , electron mass and electron charge are set to unity, the Schrodinger equation for the hydrogen atom in a constant electric field of magnitude F in the z direction is given by

$$\left(\Delta + 2\left[E - \left(Fz - \frac{1}{r}\right)\right]\right)\Psi = 0. \tag{21}$$

Here  $\Delta$  is the laplacian operator. Using parabolic coordinates, where the cartesian ones are given in terms of the new coordinates by  $x=\sqrt{\xi\eta}cos\phi, y=\sqrt{\xi\eta}sin\phi, z=\frac{\xi-\eta}{2}$  and writing the wave function in the product form

$$\Psi = \frac{1}{\sqrt{\xi \eta}} V(\xi) U(\eta) \exp(im\phi), \tag{22}$$

we get two separated equations:

$$\frac{d^2V}{d\xi^2} + (\frac{E}{2} + \frac{\beta_1}{\xi} - \frac{F}{4}\xi + \frac{1 - m^2}{4\xi^2})V(\xi) = 0,$$
(23)

$$\frac{d^2U}{d\eta^2} + (\frac{E}{2} + \frac{\beta_2}{\eta} + \frac{F}{4}\eta + \frac{1 - m^2}{4\eta^2})U(\eta) = 0.$$
 (24)

Here  $\beta_1$  and  $\beta_2$  are separation constants that must add to one. We note that these equations are of the biconfluent Heunform.

The hydrogen molecule also is treated in reference [40]. When the hydrogen-molecule ion is studied in the Born-Oppenheimer approximation, where the ratio of the electron mass to the proton mass is very small, one gets two singly

confluent Heun equations if the prolate spheroidal coordinates  $\xi = \frac{r_1 + r_2}{2c}$ ,  $\eta = \frac{r_1 + r_2}{2c}$  are used. Here c is the distance between the two centers. Assuming

$$\psi == \sqrt{\xi \eta} V(\xi) U(\eta) \exp(im\phi), \tag{25}$$

we get two confluent Heun equations:

$$\frac{d}{d\xi} \left( (1 - \xi^2) \frac{dV}{d\xi} \right) + \left( \lambda^2 \xi^2 - \kappa \xi - \frac{m^2}{1 - \xi^2} + \mu \right) V = 0, \tag{26}$$

$$\frac{d}{d\eta}\left((1-\eta^2)\frac{dU}{d\eta}\right) + \left(\lambda^2\eta^2 - \frac{m^2}{1-\eta^2} + \mu\right)U = 0. \tag{27}$$

Some additional physics papers with Heun type solutions include:

Three relatively recent papers which treat atoms in magnetic fields:

- o Exact low-lying states of two interacting equally charged particles in a magnetic field are studied by Truong and Bazzali [41]
- o The energy spectrum of a charged particle on a sphere under a magnetic field and Coulomb force are studied by Ralko and Truong [42]
- o B.S. Kandemir presented an analytical analysis of the two-dimensional Schrodinger equation for two interacting electrons subjected to a homogeneous magnetic field and confined by a two-dimensional external parabolic potential. Here a biconfluent Heun (BHE) equation is used [43]
- o Arda and Sever, in one instance with Aydoğdu studied Schrodinger equation with different potentials and in two cases found Heun and confluent Heun solutions [44,45].

In two papers Hammann et al [46,47] solved the one-dimensional Schrodinger equation for position-dependent masses, and obtained Heun solutions. The importance of these papers is the derivation and use of a relations between Heun functions which are functions of z and 1-z, which can be used for obtaining the reflection and transition amplitudes for scattering problems for waves described in terms of Heun functions.

o Recently Ishkhanyan showed that the solution of the Schrodinger equation for the  $V_0/\sqrt{x}$  can be given as a derivative of a triconfluent Heun function [48]. In another paper, solution for the same potential is given [49] as a linear combination of two confluent hypergeometric functions. For another potential which is an inverse square root near the origin and vanishes exponentially at infinity, solution is given in terms of linear combination of Gauss hypergeometric functions [50]. These potentials belong to the Heun class.

Downing showed that the solution to the one dimensional Schrodinger equation with a hyperbolic double well potential is obtained by a transformation of the confluent Heun equation [51].

Hartmann and Portnoi calculated the bound modes of two-dimensional massless Dirac fermions confined within a hyperbolic secant potential [52]

Portnoi et al. continued studying the two-dimensional Dirac particles in two papers, first confined in nonuniform magnetic fields, and second in Poschl-Teller waveguide [53,54] in terms of confluent Heun functions.

- o In a relatively recent work P. Dorey, [55] showed that equations in finite lattice systems also reduce to Heun equations.
- o Dislocation movement in crystalline materials, quantum diffusion of kinks along dislocations are some solid state applications of this equation. The book by S.Y. Slavyanov and S. Lay [37] is a general reference on problems solved before 2000.

We also cite a recent mathematical application by A.M. Ishkhanyan et al. where "total fifteen potentials for which the stationary Klein-Gordon equation is solvable in terms of the confluent Heun functions are presented.. Only nine of

the potentials are independent due to the transposition symmetry of regular singular points of the equation. Four of these equations can be reduced to the hypergeometric form. The remaining five independent Heun potentials are four-parametric and have solutions only in terms of irreducible confluent Heun functions [56]. Prof Ishkhanyan expands the Heun solution in terms of hypergeometric functions and shows that the sum has only finite number of terms in his cases. Prof. Ishkhanyan wrote additional papers after this one using the same method for other potentials. We will not comment on them, however, since from this point on, we will confine ourselves only to papers on general relativity and cosmology.

Among the papers in general relativity, we also will not be able to comment on all the works of some experts like Prof Fiziev on this field, who wrote scores of papers on Heun equations. We will give only the earlier papers and leave the reader to investigate the later ones in the ArXiv.

o In general relativity, in a relatively early work, Teukolsky studied the perturbations of the Kerr metric [57]. If we take

$$\Psi = \exp(-i\omega t) \exp(im\phi) S(\theta) R(r),$$

for the scalar particle we get two equations.

$$\frac{d}{dr} \left( \Delta \frac{dR}{d\theta} \right) + \left( [(r^2 + a^2)^2 \omega^2 - 4aMr\omega m + a^2 m^2] \Delta^{-1} - A - a^2 \omega^2 \right) R = 0, \tag{28}$$

$$\frac{1}{\sin\theta} \left( \frac{d}{d\theta} \sin\theta \frac{dS}{d\theta} \right) + \left( a^2 \omega^2 \cos^2\theta - \frac{m^2}{\sin^2\theta} + A \right) S = 0. \tag{29}$$

Here A is the separation constant,  $\Delta = r^2 - 2Mr + a^2$ .

Teukolsky just stated these equations [57]. Later these equations were found to be two coupled singly confluent Heun equations [58].

o Quasi-normal modes of rotational gravitational singularities were also studied by solving these equations by E.W. Leaver [59].

In recent applications in general relativity, Heun type equations become indispensable when one studies phenomena in higher dimensions, or in different geometries. We must note that even the simplest black hole metric, the Schwarzschild, has solutions in the Heun form [60] [61].

Some other references for general relativity applications are:

o D. Batic, H. Schmid, M. Winklmeier where the Dirac equation in the Kerr-Newman metric and static perturbations of the non-extremal Reisner-Nordstrom solution are studied [62]. D. Batic and H. Schmid also studied the Dirac equation for the Kerr-Newman metric and looked for its propagator [63]. They found that the equation satisfied is a form of a general Heun equation described in Reference [62]. In later work Batic, with collaborators continued studying Heun equations and their generalizations [64]. In his most recent paper Batic, with collaborators studied *Semi commuting and commuting operators for the Heun family* [65].

Prof. P.P. Fiziev studied problems whose solutions are Heun equations extensively.

o In a paper published in gr-qc/0603003, he studied the exact solutions of the Regge-Wheeler equation in the Schwarschild black hole interior [60].

o He presented a novel derivation of the Teukolsky-Starobinsky identities, based on properties of the confluent Heun functions [66]. These functions define analytically all exact solutions to the Teukolsky master equation, as well as to the Regge-Wheeler and Zerilli ones.

o In a talk given at 29th Spanish Relativity Meeting (ERE 2006), he depicted in more detail the exact solutions of Regge-Wheeler equation, which described the axial perturbations of Schwarzschild metric in linear approximation, in the Schwarzschild black hole interior and on Kruskal-Szekeres manifold in terms of the confluent Heun functions [67].

- o All classes of exact solutions to the Teukolsky master equation were described in terms of confluent Heun functions in Reference [\[68\]](#page-16-3) [\[69\]](#page-16-4).
- o In reference [\[70\]](#page-16-5) he reveals important properties of the confluent Heun's functions by deriving a set of novel relations for confluent Heun's functions and their derivatives of arbitrary order. Specific new sub classes of confluent Heun's functions are introduced and studied. A new alternative derivation of confluent Heun's polynomials is presented.
- o In another paper [\[71\]](#page-16-6) he, with a collaborator, noted that weak gravitational, electromagnetic, neutrino and scalar fields, considered as perturbations on Kerr background satisfied Teukolsky Master Equation. The two non-trivial equations were obtained after separating the variables, one equation only with the polar angle and another using only the radial variable. These were solved by transforming each one into the form of a confluent Heun equation.
- o Fiziev is an expert in this topic. Two further articles by him and his collaborator are Solving systems of transcendental equations involving the Heun functions, [\[72\]](#page-16-7) and Application of the confluent Heun functions for finding the quasinormal modes of non rotating black holes [\[73\]](#page-16-8).

We also cite one of the last papers of Fiziev on the mathematical properties of this subject which can have applications in physics. In [\[74\]](#page-16-9), the author "introduces and studies a novel type of solutions to the general Heun equation". His approach is based " on the symmetric form of the Heun differential equation yielded by development of the Papperitz-Klein symmetric form of the Fuchsian equations with an arbitrary number of regular singular points greater than 4. The symmetry group of these equations turns to be a proper extension of the Mobius group". He also introduces and studies" new series solutions and derives solutions for the four singular point case which treats simultaneously and on an equal footing all singular points."

Among other papers on this subject one may cite the following papers:

- o R.Manvelyan, H.J.W. Muller Kirsten, J.Q. Liang and Y. Zhang, calculated the absorption rate of a scalar by a D3 brane in ten dimensions in terms of modified Mathieu functions, and obtained the S-matrix in reference [\[75\]](#page-16-10).
- o T.Oota and Y.Yasui studied the scalar laplacian on a wide class of five dimensional toric Sasaki-Einstein manifolds, ending in two Heun's differential equations in reference [\[76\]](#page-16-11).
- o S.Musiri and G. Siopsis found out that the wave equation, obtained in calculating the asymptotic form of the quasi-normal frequencies for large AdS black holes in five dimensions, reduces to a Heun equation, in reference [\[77\]](#page-16-12).
- o A. Al-Badawi and I. Sakalli studied the Dirac equation in the rotating Bertotti-Robinson spacetime [\[78\]](#page-16-13) ending up with a Heun type equation.
- I first encountered this type of equation when we tried to solve the scalar wave in the background of the Nutku helicoid instanton [\[30\]](#page-14-22). In this case, for a scalar particle in this background metric, one gets the Mathieu equation which is a special case of the Heun equation. In the same paper, the solutions in four dimensions involve the product of two exponentials and two Heun functions. These solutions can be summed to give the Green's function for this problem in a closed form. We could not succeed to obtain a closed form solution for the Greens function when the similar problem is studied in five dimensions [\[79,](#page-16-14) [80\]](#page-16-15).
- o The helicoid instanton is a double-centered solution. As remarked above, for the simpler instanton solution of Eguchi-Hanson [\[3\]](#page-14-0) hypergeometric solutions are sufficient [\[27\]](#page-14-19). Here one must remark that another paper using the Eguchi-Hanson metric ends up with the confluent Heun equation [\[81\]](#page-16-16). These two papers show that sometimes judicious choice of the coordinate system and separation ansatz matters.
- o Sucu and Unal also obtained closed solutions for the spinor particle written in t ¨ he background of the Nutku helicoid instanton [\[27\]](#page-14-19), whereas using the separation of variables method gives us an infinite series of product of two Mathieu functions [\[80\]](#page-16-15)
- o One can show that the solutions of Sucu and Unal can be expanded in terms of Mathieu functions if one attempts ¨ to use the separation of variables method. as described by L.Chaos-Cador and E. Ley-Koo [\[82\]](#page-16-17).

o Tolga Birkandan and I also found an extension of the Heun equation with five singular points [\[31\]](#page-14-23), and calculated the solution of a scalar field in the background of the Eguchi-Hanson equation trivially extended to five dimensions [\[31\]](#page-14-23). Then the solution for the radial component turned out to be given in terms of the confluent Heun equation.

o Mirjam Cvetiˇc and Finn Larsen studied grey body factors and event horizons for rotating black holes with two rotation parameters and five charges in five dimensions. When the Klein-Gordon equation for a scalar particle in this background is written, one gets a confluent Heun equation. In the asymptotic region this equation turns into the hypergeometric form [\[83\]](#page-16-18). When they studied the similar problem for the rotating black hole with four U(1) charges, they again obtained a confluent Heun equation for the radial component of the Klein-Gordon equation, which they reduce to the hypergeometric form by making approximations [\[84\]](#page-16-19). These two papers are partly repeated in [\[85\]](#page-16-20) . Same equations were obtained which were reduced to approximate forms which gave solutions in the hypergeometric form.

M. Cvetiˇc encounters this function in several of her publications and reduces them to the hypergeometric form by giving physical arguments to drop certain terms in the equation. The hypergeometric solution points to the presence of conformal symmetry in the reduced model [\[86,](#page-16-21)[87\]](#page-16-22). The method is going to the extreme and near extreme (Kerr/CFT correspondence) limits, going to the boundary and in some cases using a "subtracted metric" using a warp factor which preserves all the near horizon properties of the black hole such as the entropy and the thermodynamic potentials, and if necessary dropping certain terms which are negligible in these limits [\[88–](#page-16-23)[90\]](#page-16-24).

"In general, conformal symmetry does not exist in the non-extremal cases. The solutions often turn out to be of the Heun form. In the extremal case two horizons overlap. In the near extremal case they are very close to each other. In these two cases and in the near horizon limit, we find conformal symmetry, resulting in solutions which are hypergeometric functions, or one of its confluent forms. If we want conformal symmetry without going to the extremal or the near horizon limit, we have to change the 'warp factor'. When the warp factor is changed, the rest of the metric preserves its initial form. The thermodynamic potentials and entropy do not change. You have to drop some terms resulting in solutions in the hypergeometric form. This is equivalent to putting the black hole into a conic box. If you go to the asymptotic or to the scaling limit, this is seen clearly. In these limits the Einstein equations are not satisfied unless the energy-momentum tensor, on the right side of the Einstein equations are also changed, to account for putting the system into the conic box". [\[91\]](#page-16-25)

Cvetiˇc also studied black holes in supergravity with Birkandan. Heun solutions also exist for the Wu Black Hole which is the most general solution of maximally supersymmetric gauged supergravity in D=5 [\[92\]](#page-16-26). Here they did not study the limiting cases. For the massless Klein-Gordon equation in the background of the most general black hole in four dimensions and N=2 gauge supersymmetry with U(1)<sup>2</sup> gauge symmetry (Chow-Compere solution [\[93\]](#page-16-27)), the angular equation gives Heun type solutions. The radial equation has five regular singularities, which reduce to hypergeometric functions in the near horizon extremal limit [\[94\]](#page-16-28).

o We should also mention two papers by H.R. Christiansen and M.S. Cunha with Heun type solutions. These are: Confluent Heun functions in gauge theories on thick braneworlds [\[95\]](#page-16-29), and Kalb-Ramond excitations in a thick-brane scnario with dilaton [\[96\]](#page-17-0). In the first paper, the propagation modes of gauge fields in an infinite Randall-Sundrum scenario are investigated. Here a sine-Gordon soliton represents the thick four dimensional braneworld while an exponentially coupled scalar field acts for the dilaton. For the gauge field motion a differential equation is found which can be transformed into a confluent Heun equation. In the second paper a similar scenario is used. Here a bulk Kalb- Ramond field is coupled to a dilaton, in a warped space-time in the presence of a brane field in five dimensions. Full spectrum and eigenstates are studied. In the general case, the solution to the field equations are given in terms of the confluent Heun function, which reduces to the confluent hypergeometric function for special values of the parameters.

Other relevant references I could find, are listed as references [\[97\]](#page-17-1) - [\[104\]](#page-17-2).

o The more recent papers on this subject include The quantum treatment of the 5D-warped Friedman-Robertson-Walker universe in Schrodinger Picture [\[105\]](#page-17-3). Here the time-evolving Schrodinger version of the Wheeler-De Witt equation, written for the five dimensional warped k=0-FRW Universe is studied. For small values of the cosmological scale factor, a, the wave function of the Universe is expressed in terms of the Heun Double Confluent functions, whereas for large values of this parameter the solution becomes the Hermite associated functions. Two papers by the same authors using Heun type functions are Fermions in magnestar's crust in terms of Heun double confluent functions [\[106\]](#page-17-4), and The approximative analytic study of fermions in magnetar's crust; ultra-relativistic plane waves, Heun and Mathieu solutions and beyond [\[107\]](#page-17-5).

- o In Fermi surfaces and analytic Green's functions from conformal gravity [\[109\]](#page-17-6), T2-symmetric charged AdS black holes are constructed in conformal gravity. The most general solution up to an overall conformal factor contains three nontrivial parameters: the mass, electric charge and a quantity that can be identified as the massive spin-2 hair. The Dirac equation for the charged massless spinor in this background can be solved in terms of the general Heun's function for generic frequency ω and wave number k. This allows us to obtain the analytic Green's function G(ω, k) for both extremal and non-extremal black holes. For some special choice of black hole parameters, the Green's function reduces to simpler hypergeometric or confluent hypergeometric functions.
- o Two of the authors of the paper quoted above had calculated the Greens's functions in terms of the Heun function in an earlier paper, Exact Green's functions from conformal gravity [\[108\]](#page-17-7).
- o Another paper is: Quantized black hole and Heun function by D. Momeni, K. Yerzhanov and R. Myrzakulov [\[110\]](#page-17-8) where a black hole is quantized using the Bohr method. The solution turns to be of the Heun type .
- o In On an approach to constructing static ball models in general relativity by A.M. Baranov, some solutions of the Einstein equation were described by Heun functions [\[111\]](#page-17-9).
- o In an unpublished paper, On analytic solutions of wave equations in regular coordinate systems on Schwarzschild background Dennis Philipp and Volker Perlick claim that they find " The wave equation for the propagation of (massless) scalar, electromagnetic and gravitational waves on fixed Schwarzschild background spacetime, which is described by the general time-dependent Regge-Wheeler equation, canbe transformed to usual Schwarzschild, Eddington-Finkelstein, Painleve Gullstrand and Kruskal-Szekeres coordinates. In the first three cases, but not in the last one, it is possible to separate a harmonic time-dependence. Then the resulting radial equations belong to the class of confluent Heun equations" [\[112\]](#page-17-10).

Among additional papers we can also cite the article of Bezerra et al, Exact solutions of the Klein-Gordon equation in the Kerr-Newman background and Hawking Radiation where both the radial and angular solutions are given in terms of confluent Heun functions [\[113\]](#page-17-11). In the particular case corresponding to an extreme Kerr-Newman black hole, the solution is given by the double confluent Heun functions [\[114\]](#page-17-12). Biconfluent Heun functions were obtained for the exact solution of the Schrodinger equation for a particle (galaxy) moving in a Newtonian universe with a cosmological constant [\[115\]](#page-17-13).

Other papers on general relativity written in 2015 also include New results for electromagnetic quasinormal and quasibound modes of Kerr black holes, by D.Staicova and P.Fiziev [\[116\]](#page-17-14), where the authors solve Teukolsky equations with confluent Heun solutions numerically. In Heun functions describing fermions evolving in paralel and magnetic fields, by C. Dariescu and M.A. Dariescu, [\[117\]](#page-17-15), the solutions are in terms of double confluent Heun functions. Same authors also published Quantum analysis of k=-1 Robert-Walker Universe, where they solved the Wheeler-DeWitt equation [\[118\]](#page-17-16). The solutions turned out to be Heun functions. M.C. E. Cedeno and J.C.N. Araujo show that for Master equation solutions in the linear regime of characteristic formulation of general relativity, the solution is in terms of confluent Heun's functions for radiative case in the Schwarzschild's background [\[119\]](#page-17-17). In Massless Dirac particles in the vacuum C-metric, D.Bini et al show that the Dirac equation, written in the background of the C- metric can be reduced to a radial and an angular equation, both of which can be solved in terms of general Heun functions [\[120\]](#page-17-18). Vieira et al [\[121\]](#page-17-19) show that for Charged massive scalar fields are considered in the gravitational and electromagnetic field produced by a dyonic black hole with a cosmic string along its axis of symmetry " exact solutions of both angular and radial parts of the covariant Klein-Gordon equation in this background can be obtained, and are given in terms of the confluent Heun functions". In [\[122\]](#page-17-20) , Kofron separates test fields equations on the non-rotating C- metric background. He finds that the resulting equations are of the Heun or confluent Heun form for the general case. These equations, however, can be reduced to hypergeometric functions in the static, axial symmetric and the extremal case where the inner and outer horizons coalesce. In another paper [\[123\]](#page-17-21), the same author studies the similar phenomena on the background of the rotating C- metric. For the general case, the radial equation has five regular singularities. In the extremal, static and axial symmetric cases, one obtains a polynomial solution.

Some other papers published in 2016 in the field of general relativity where solutions to field equations in the background

of different metrics are as follows:

Valtancoli [\[124\]](#page-17-22) found Heun solutions for the radial part of the Klein-Gordon equation when the scalar field is conformally coupled to a charged BTZ black hole.

Vieira and Bezerra [\[125\]](#page-17-23) study "resonant frequencies, Hawking radiation and scattering of scalar waves...", and find confluent Heun solutions. They also study [\[126\]](#page-17-24) the class of solutions of the Wheeler-DeWitt equation in the Friedmann-Robertson-Walker universe. In still another paper [\[127\]](#page-17-25), these authors find confluent Heun solutions for the massless Klein-Gordon equation in the background metric of the three dimensional rotating and four dimensional canonical acoustic black holes.

Sakalli [\[128\]](#page-17-26) finds analytical solutions in rotating linear dilaton black holes.

Kraniotis [\[129\]](#page-18-0) studies the Klein-Gordon equation in the background metric of the Kerr-Newman (anti) de Sitter black hole. He first reduces the radial and angular equations to the Heun form, writes the solution in terms of local Heun and confluent Heun functions. In my opinion this paper should be also praised for the introduction of the "false singular point" concept, which reduces the solution to hypergeometric functions for certain values of the physical parameters in the equation.

Since we updated this paper in February 2017, we find close to thirty new publications if one searches the word Heun Functions in the index Web of Knowledge ln April 2018. Many of these papers are on the mathematical aspects of the equation and solving Schrodinger equations for different new potentials in terms of Heun or linear combinations of Heun functions. There are also solutions in terms of Heun functions for equations used in different branches of physics. Here we will attempt to review only the papers for applications in physics related to general relativity.

In [\[130\]](#page-18-1), Arda et al solve the energy relations obtained with the help of the quantization rule for the Klein-Gordon equation with a linear plus an inverse-linear potential in terms of bi-confluent Heun equations. Vieira wrote two papers [\[131,](#page-18-2) [132\]](#page-18-3) where he first studied Resonant frequencies of a hydrodynamic vortex. The radial equation has solutions in terms of double confluent Heun functions. In the second paper, analytic solutions for sound perturbations in the presence of a rotating acoustic black hole which is an analogue of the conical Kerr metric were studied. In the massless case, the radial equation has Heun type solutions. Vieira also wrote another paper with co-authors [\[133\]](#page-18-4), where massive scalar fields are considered in the gravitational field produced by a Schwarzschild black hole with a global monopole in f(R) gravity. The exact solution of the radial part of the Klein-Gordon equation in this background is given in terms of the general Heun functions. The properties of the general Heun functions are applied to study the Hawking radiation and the resonant frequencies of scalar particles.

Ciprian Dariescu wrote two papers with collaborators [\[134,](#page-18-5) [135\]](#page-18-6). In the the first paper, using a perturbative method, Klein-Gordon equation for a charged massive field in the background of a magnetar is solved both in the interior solution and outside the star. Equations can be seperated with general and confluent Heun function solutions. With special conditions on parameters, polynomial solutions can be found and first order transition amplitudes are computed [\[134\]](#page-18-5). In the second paper, for the spatially open Friedmann-Robertson-Walker (FRW) Universe with stiff matter and radiation as non-interacting matter sources, the scale function coming from the integration of the Friedmann equation is expressed in terms of elliptic integrals. For a negative cosmological constant, the allowed ranges for the models parameters are identified. Within the quantum analysis, the Wheeler DeWitt (WDW) equation turns into a modified Morse equation whose solutions are Mathieu and Heun functions. [\[135\]](#page-18-6).

Sobhani et al [\[136\]](#page-18-7) wrote a paper where the thermodynamicl properties of the anharmonic oscillator cosmic string framework are studied. The Schrodinger equation is written in the cosmic string framework and anharmonic oscillations are investigated. The wave function and energy spectrum are derived using confluent Heun functions.

Birkandan was also active in this period. He wrote four papers. In the first paper, with Bouaziz, he showed that the deformed Schrodinger equation for a singular inverse square potential in coordinate space with a minimal length is solved in terms of Heun functions [\[137\]](#page-18-8). In his second paper with a collaborator, confluent Heun solutions to the radial equations of two Halilsoy-Badawi metrics are found. For the first metric, the radial part of the massless Dirac equation and for the second case, the radial part of the massless Klein-Gordon equation are studied [\[138\]](#page-18-9), both with Heun type solutions. In the third paper, he and his collaborator showed that Heun-type exact solutions emerged for both the radial and the angular equations for the case of a scalar particle coupled to the zero mass limit of both the Kerr and Kerr-(anti)de-Sitter space times. Since any type D metric has Heun-type solutions, it is interesting that this property is retained when the black hole has a zero mass limit. This work further refuted the claims that mass of the black hole, going to zero limit of the Kerr metric was both locally and globally the same as the Minkowski metric [139]. We comment on the fourth paper in the Conclusion section.

A comprehensive bibliography can be found at the bibliography section of http://tcpa.uni-sofia.bg/heun/home.html, compiled by Profs. Plamen Fiziev and Denitsa Staicova.

Just to give an example of how the Heun function is emerges in a simple problem, in the next section, our work in [31] for the scalar particle in the background metric of the extended Eguchi-Hanson solution will be sketched..

#### 4 Scalar field in the background of the extended Eguchi-Hanson solution

To go to five dimensions, we can add a time component to the Eguchi-Hanson [3] metric so that we have

$$ds^{2} = -dt^{2} + \frac{1}{1 - \frac{a^{4}}{4}}dr^{2} + r^{2}(\sigma_{x}^{2} + \sigma_{y}^{2}) + r^{2}(1 - \frac{a^{4}}{r^{4}})\sigma_{z}^{2}$$
(30)

where

$$\sigma_x = \frac{1}{2}(-\cos\xi d\theta - \sin\theta\sin\xi d\phi) \tag{31}$$

$$\sigma_y = \frac{1}{2} (\sin \xi d\theta - \sin \theta \cos \xi d\phi) \tag{32}$$

$$\sigma_z = \frac{1}{2}(-d\xi - \cos\theta d\phi). \tag{33}$$

This is a vacuum solution. If we take

$$\Phi = e^{ikt}e^{in\phi}e^{i(m+\frac{1}{2})\xi}\varphi(r,\theta),\tag{34}$$

we find the scalar equation as

$$\varphi(r,\theta) = \left(\frac{r^4 - a^4}{r^2}\partial_{rr} + \frac{3r^4 + a^4}{r^3}\partial_r + k^2r^2 + \frac{4a^4m^2}{a^4 - r^4} + 4\cot\theta\partial_\theta + \frac{8mn\cos\theta - 4(m^2 + n^2)}{\sin^2\theta}\right)\varphi(r,\theta).$$
(35)

If we take  $\varphi(r,\theta) = f(r)g(\theta)$ , the solution of the radial part is expressed in terms of confluent Heun  $(H_C)$  functions.

$$f(r) = \left(-a^4 + r^4\right)^{\frac{1}{2}m} H_C\left(0, m, m, \frac{1}{2}k^2a^2, \frac{1}{2}m^2 - \frac{1}{4}\lambda - \frac{1}{4}k^2a^2, \frac{a^2 + r^2}{2a^2}\right)$$

$$+ \left(a^2 + r^2\right)^{-\frac{1}{2}m} \left(r^2 - a^2\right)^{\frac{1}{2}m} H_C\left(0, -m, m, \frac{1}{2}k^2a^2, \frac{1}{2}m^2 - \frac{1}{4}\lambda - \frac{1}{4}k^2a^2, \frac{a^2 + r^2}{2a^2}\right)$$

$$(36)$$

If the variable transformation  $r = a\sqrt{\cosh x}$  is made, one solution can be expressed as

$$f(x) = \left(\sinh(x)\right)^m H_C\left(0, m, m, \frac{1}{2}k^2a^2, \frac{1}{2}m^2 - \frac{1}{4}\lambda - \frac{1}{4}k^2a^2, \frac{1}{2}\cosh^2(x/2)\right). \tag{37}$$

We tried to express the equation for the radial part in terms of  $u = \frac{a^2 + r^2}{2a^2}$  to see the singularity structure more clearly. Then the radial differential operator reads

$$4\frac{d^2}{du^2} + 4\left(\frac{1}{u-1} + \frac{1}{u}\right)\frac{d}{du} + k^2a^2\left(\frac{1}{u-1} + \frac{1}{u}\right) + \frac{m^2}{u^2(1-u)^2}.$$
 (38)

This operator has two regular singularities at zero and one, and an irregular singularity at infinity, the singularity structure of the confluent Heun equation. This is different from the hypergeometric equation, which has regular singularities at zero, one and infinity.

The solution of the angular equation which is regular at  $\theta = \pi$  for m greater than n is given below in terms of hypergeometric functions.

$$g(\theta) = \sin(\theta)^m \cot(\theta/2)^n \times_2 F_1(([m + \frac{1}{2}\sqrt{\lambda + 1} + \frac{1}{2}, m - \frac{1}{2}\sqrt{\lambda + 1} + \frac{1}{2}],) [1 + n + m], \frac{1}{2}\cos^2(\theta)/2).$$

#### 5 Conclusion

In this paper, first the Heun function is introduced, then some its uses in physics, especially in the field of general relativity and gravitation are demonstrated. We have to note that most of the physicists that bluntly state their solution is in terms of Heun functions are mainly from the third world. We see physicists from Bulgaria, Romania, Brazil, Armenia, even Turkey in this group. There are mathematicians from the western world, though, who are experts in this field. Batic, a mathematician, although he now works in U.A.E. may be considered from the western world. Ronveaux from Belgium, and many other mathematicians are from the western world.

They are not really many exceptions to this observation. Cvetič and Larsen demonstrate what the physicists from the western world do. They try to express their solutions in terms of hypergeometric functions, by going to the asymptotics, to the extremal or to the near extremal limit, or putting the solution into a conic box, by changing the energy momentum term if necessary, but keeping the thermodynamic potentials same. A long endeavor was necessary to label *Teukolsky Master Equations* as belonging to the Heun class [58]. Only recently the equation given by 't Hooft [140] was shown to belong to the Heun class if it were not modified [141]. When modified the solution is the manageable hypergeometric function. We agree that this impression may be wrong, but it is just an observation.

The first version of this paper was submitted to the 13th Regional Conference on Mathematical Physics, which was held in Antalya, Turkey on 27-31 October 2010 and printed in [142].

## 6 Acknowledgement

I am grateful to Prof.s Cemsinan Deliduman and Kayhan Ülker for providing me a shelter at Mimar Sinan Fine Arts University during my days in retirement. I indebted to Tolga Birkandan for collaboration and technical assistance. I am grateful to Prof. Dr. André Ronveaux for informing me of a slight error in my reference 9. I thank Science Academy, Turkey for support.

#### <span id="page-13-0"></span>References

- <span id="page-13-1"></span>[1] K. Heun, Math. Annalen **33**(1889) 161.
- [2] R.P. Kerr, Phys. Rev. Lett. 11, 237 (1963)238.

- <span id="page-14-1"></span><span id="page-14-0"></span>[3] T. Eguchi, A.J. Hanson, Physics Letters 74B (1978) 249.
- <span id="page-14-2"></span>[4] Philip M. Morse and Herman Feshbach Methods of Theoretical Physics McGraw Hill Company, New York (1953).
- <span id="page-14-3"></span>[5] E.L. Ince, Ordinary Differential Equations (Dover Publications) (1926, 1956).
- [6] Philip M. Morse and Herman Feshbach Methods of Theoretical Physics McGraw Hill Company, New York 1953, pp532.
- <span id="page-14-5"></span><span id="page-14-4"></span>[7] A. Ronveaux (ed.), Heun's Differential Equations (Oxford University Press) (1995).
- <span id="page-14-6"></span>[8] R.S. Maier, J. Diff. Equations 213, (2005)171-203 , [math.CA/0203264.](http://arxiv.org/abs/math/0203264)
- <span id="page-14-7"></span>[9] F.M.Arscott, in Heun's Differential Equations, ed. by A. Ronveaux, Oxford University Press(1995).
- <span id="page-14-8"></span>[10] A. Ronveaux, Applied Math. and Comp., 141(2005) 177.
- <span id="page-14-9"></span>[11] F.M.Arscott, in Heun's Differential Equations, ed. by A. Ronveaux, Oxford University Press(1995). pp ,11-12, 39-41.
- <span id="page-14-10"></span>[12] F.M.Arscott, in Heun's Differential Equations, ed. by A. Ronveaux, Oxford University Press(1995). pp 12, 41-44.
- <span id="page-14-11"></span>[13] F.M. Arscott: in Heun's Differential Equation, ed.by A. Ronveaux, Oxford University Press (1995), p. 65.
- <span id="page-14-12"></span>[14] P.P. Fiziev: Class. Quant. Grav. 27 (2010) 135001.
- <span id="page-14-13"></span>[15] G.Valent,Siam J. Math. Analysis 17 (1986) 688 .
- [16] J. Meixner and F. W. Schaefke , Mathieusche Funktionen und Sphaeroidfunktionen mit anwendungen auf physikalische und technische Probleme , Springer Verlag (1954).
- <span id="page-14-15"></span><span id="page-14-14"></span>[17] N.W. McLachlan, Theory and Applications of Mathieu Functions, Dover reprint from 1946 edition (1963).
- <span id="page-14-16"></span>[18] A Course of Modern Analysis, E.T.Whittaker and G.N.Watson, Cambridge University Press,(1963).
- [19] Bateman Manuscript, Higher Transcendental Functions, Vol. III, A.Erdelyi, W.Magnus,F.Oberhettinger and F.G.Tricomi, Mc.Graw Hill (1955).
- <span id="page-14-17"></span>[20] R. Schafke, D. Schmidt, SIAM J. Math. Analysis 11 (1980) 848.
- [21] R.S. Maier, R.S. Maier, J. Diff. Equations 213(2005) 171, [math.CA/0203264,](http://arxiv.org/abs/math/0203264) [math.CA/0408317.](http://arxiv.org/abs/math/0408317)
- [22] N. Gurappa, P.K. Panigrahi, Journal of Physics A: Math. General 37(2004) L605-L608 .
- [23] K.Kuiken, SIAM J. Math. Analysis 10 (1979) 655.
- [24] B.D.B. Figueiredo, J. Math. Phys., 46 , 113503 (2005), [math-phys/0509013.](http://arxiv.org/abs/math-phys/0509013)
- <span id="page-14-18"></span>[25] S.Bellucci, V.Yeghikyan,J. Math. Phys., 54 (2013) 082103, [arXiv:1302.0798.](http://arxiv.org/abs/1302.0798)
- [26] P.Aydiner and T.Birkandan, Physical problems admitting Heun-to-hypergeometric reduction, Proceedings of the International Conference DAYS on DIFFRACTION 2015, pp.27-33 [,arXiv:1401.0449](http://arxiv.org/abs/1401.0449) [math-ph](2015).
- <span id="page-14-20"></span><span id="page-14-19"></span>[27] Y. Sucu, N. Unal, Class. Quant. Grav. ¨ 21 (2004) 1443.
- <span id="page-14-21"></span>[28] Y. Nutku, Phys. Rev. Lett. 77 (1996) 4702.
- <span id="page-14-22"></span>[29] D. Lorenz-Petzold, J. Math. Phys. 24(1983) 2632 .
- <span id="page-14-23"></span>[30] A.N. Aliev, M. Horta¸csu, J. Kalaycı , Y. Nutku, Class. Quant. Grav. 16 (1999) 631.
- <span id="page-14-24"></span>[31] T. Birkandan, M. Horta¸csu, J. Math. Phys. 49 (2008) 00054101 .
- [32] Philip M. Morse and Herman Feshbach Methods of Theoretical Physics McGraw Hill Company, New York 1953, p. 1407.

- <span id="page-15-1"></span><span id="page-15-0"></span>[33] B.D.B. Figueiredo, [math-ph 10402071](http://arxiv.org/abs/math-ph/1040207).
- <span id="page-15-2"></span>[34] B.D.B.. Figueiredo, J. Math. Phys. 48(2007) 013503 .
- <span id="page-15-3"></span>[35] L.J. El-Jaick, B.D.B. Figueiredo, J. Math. Phys. 49 (2008) 083508 .
- <span id="page-15-4"></span>[36] L.J. El-Jaick, B.D.B. Figueiredo,J. Math. Phys. 50 (2008) 123511 .
- <span id="page-15-5"></span>[37] S.Y. Slavyanov, W. Lay, Special Functions, A Unified Theory Based on Singularities(Oxford University Press) (2000).
- <span id="page-15-6"></span>[38] P.S. Epstein, Phys. Rev. 28 (1926) 695.
- [39] S.Y. Slavyanov, Asymptotic Solutions of the One-dimensional Schrodinger Equation (Leningrad University Press) (in Russian) (1991), Translation into English: S.Y. Slavyanov, Asymptotic Solutions of the One-dimensional Schrodinger Equation (Amer. Math. Soc. Trans. of Math. Monographs) 151 (1996).
- <span id="page-15-8"></span><span id="page-15-7"></span>[40] A.H. Wilson, Proc. Royal Soc. London, A118 (1928) 617, also reference [\[37\]](#page-15-4), p.167.
- <span id="page-15-9"></span>[41] T.T. Truong, D. Bazzali, Phys. Lett. A 269 (2000) 186.
- <span id="page-15-10"></span>[42] A. Ralko A, T.T. Truong, J. Phys. A 35 (2002) 9573.
- <span id="page-15-11"></span>[43] B.S. Kandemir, J. Math. Phys. 46 (2005) 032110.
- <span id="page-15-12"></span>[44] A.Arda, O.Aydo˘gdu and R. Sever, J. of Physics A, 43 (2010) 425204 .
- <span id="page-15-13"></span>[45] A.Arda and R.Sever, Comm. Theor. Physics 64 (2015) 269.
- <span id="page-15-14"></span>[46] L.Dekar, Lyazid Chetouani, T.F. Hammann, Jour. Mathematical Physics 39 (1998) 2551.
- <span id="page-15-15"></span>[47] L.Dekar, Lyazid Chetouani, T.F. Hammann, Phys.Rev A 59 (1999) 107.
- <span id="page-15-16"></span>[48] A.M. Ishkhanyan, Eur. Phys. Lett. 112(2015) 10006.
- <span id="page-15-17"></span>[49] A.M. Ishkhanyan, Physics Letters A, 380(2016) 3786.
- <span id="page-15-18"></span>[50] A.M. Ishkhanyan, Europ.Phys. J. Plus 133 (2018) Article Number: 83.
- <span id="page-15-19"></span>[51] C.A. Downing, J. Math. Phys, 54(201389 072101.
- <span id="page-15-20"></span>[52] R.R.Hartmann, M.E. Portnoi, Phys. Rev. A89 (2014) 012101.
- <span id="page-15-21"></span>[53] C.A. Downing, M.E. Portnoi,Phys. Rev B, 94 (2016) Article Number: 165404.
- <span id="page-15-22"></span>[54] R.R.Hartmann, M.E. Portnoi,Scientific Reports, 7(2017) Article Number:11599
- <span id="page-15-23"></span>[55] P.Dorey, J. Suzuki and R.Tateo,J. Phys. A, 37(2004) 2047.
- <span id="page-15-24"></span>[56] A.S. Tarloyan, T.A.Ishkhanyan and A.M.Ishkhanyan, AM, Annalen der Physik 528 (2016)264.
- <span id="page-15-25"></span>[57] S.A. Teukolsky, Phys. Rev. Lett. 29 (1972) 1114.
- <span id="page-15-26"></span>[58] D. Batic, J. Schmid, J. Math. Phys. 48 (2007) 042502.
- <span id="page-15-27"></span>[59] E.W. Leaver, Proc. Royal Soc. London A 402 (1985) 285 .
- <span id="page-15-28"></span>[60] P. Fiziev, Class. Quant. Grav. 23 (2006) 2447.
- <span id="page-15-29"></span>[61] D. Philipp , V. Perlick (2015) [arXiv:1503.08101.](http://arxiv.org/abs/1503.08101)
- <span id="page-15-30"></span>[62] D. Batic, H. Schmid, M.. Winklmeier, [gr-qc/0607017,](http://arxiv.org/abs/gr-qc/0607017) J. Phys. A, 39(2006) 12559 .
- <span id="page-15-31"></span>[63] D. Batic, H. Schmid, Prog. Theor. Phys. , 116 (2006)5 17.
- [64] D. Batic, M. Sandoval, e-Print: [arXiv:0805.4399](http://arxiv.org/abs/0805.4399) [gr-qc].

- <span id="page-16-1"></span><span id="page-16-0"></span>[65] D.Batic, D. Mills, M. Nowakowski, Theoret. Mathem. Physics, 195 (2018) 6.
- <span id="page-16-2"></span>[66] P.P. Fiziev, Phys. Rev. D 80 (2009) 124001 .
- <span id="page-16-3"></span>[67] P.P. Fiziev,J.Phys.Conf.Ser.66(2007)012016,. e-Print: [gr-qc/0702014.](http://arxiv.org/abs/gr-qc/0702014)
- <span id="page-16-4"></span>[68] P.P. Fiziev, Class. Quant. Grav. 27 (2010) 135001.
- <span id="page-16-5"></span>[69] D.R.Staicova and P.P. Fiziev , Astrophys. Space Sci. 332 (2011) 385 .
- <span id="page-16-6"></span>[70] P.P. Fiziev, J. Phys. A,43 (2010) 035203.
- <span id="page-16-7"></span>[71] R.S. Borissov, P.P. Fiziev, e-Print: [arXiv:0903.3617](http://arxiv.org/abs/0903.3617) [gr-qc].
- [72] P. Fiziev, and D. Staicova, Towards New Paradigms : Proceeding of the Spanish Relativity Meeting 2011, ed. by I.B. Jimenez,J.S.R. Cembranos, A. Dobado, et. Al, Book Series: AIP Conference Proceedings, Vol 1458, pp. 395-398, (2012).
- <span id="page-16-9"></span><span id="page-16-8"></span>[73] P. Fiziev and D. Staicova, Phys. Rev. D, 84(2011) 127502.
- <span id="page-16-10"></span>[74] P. Fiziev,Intern. Frontier Science Lett. 7 11 (2016).
- [75] R.Manvelyan, H.J.W. Muller Kirsten, J.Q. Liang and Y. Zhang, Nucl. Phys. B579 (2000) 177 , e-Print: [hep-th/0001179.](http://arxiv.org/abs/hep-th/0001179)
- <span id="page-16-12"></span><span id="page-16-11"></span>[76] T. Oota, Y. Yasui, Nucl. Phys. B 742 (2006) 275 .
- <span id="page-16-13"></span>[77] S Musiri, G. Siopsis, Phys. Lett. B 576 (2003) 309, e-Print[:hep-th/0308196.](http://arxiv.org/abs/hep-th/0308196)
- <span id="page-16-14"></span>[78] A.Al-Badawi, I. Sakalli, J. Math. Phys. 49 (2008) 052501 .
- <span id="page-16-15"></span>[79] T.Birkandan, M. Horta¸csu, J. Math. Phys. 48 (2007) 092301.
- <span id="page-16-16"></span>[80] T. Birkandan, M. Horta¸csu, J. Phys. A, 40(2007) 1105.
- <span id="page-16-17"></span>[81] A.Malmendier, J. of Math. Phys., 44(2003) 4308.
- <span id="page-16-18"></span>[82] L. Chaos-Cador, E. Ley-Koo, Rev. Mex. Fis. 48(2002) 67.
- <span id="page-16-19"></span>[83] M. Cvetiˇc and F. Larsen, Phys. Rev. D 56(1997) 4994 . e-Print[:hep-th/9705192.](http://arxiv.org/abs/hep-th/9705192)
- <span id="page-16-20"></span>[84] M. Cvetiˇc and F. Larsen, Nucl. Phys. B 506(1997) 107. e-Print: [hep-th/9706071](http://arxiv.org/abs/hep-th/9706071)
- <span id="page-16-21"></span>[85] M.Cvetiˇc and F. Larsen, J.High Energy Phys. 0909 (2009) 088. e-Print: [arXiv:0908.1136](http://arxiv.org/abs/0908.1136) [hep-th]
- <span id="page-16-22"></span>[86] T. Birkandan and M.Cvetiˇc, Phys. Rev. D, 84(2011) 044018.
- <span id="page-16-23"></span>[87] M.Cvetiˇc and F. Larsen, J. High Energy Phys., 1202 (2012) 122.
- [88] M.Cvetiˇc and F. Larsen, J. High Energy Phys. , 1209 (2012) 076.
- <span id="page-16-24"></span>[89] M. Cvetiˇc and F. Larsen, J. High Energy Phys. ,1411 (2014) 033.
- <span id="page-16-25"></span>[90] M.Cvetiˇc and G. Gibbons , J. High Energy Phys. 1207 (2012) 014.
- <span id="page-16-26"></span>[91] Tolga Birkandan (Private communication).
- <span id="page-16-27"></span>[92] T. Birkandan, M. Cvetiˇc J.High Energy Phys. 1409 (2014) 121.
- <span id="page-16-28"></span>[93] D. D. K. Chow and G. Compere, Class. Quant. Grav. 31 (2014) 022001.
- <span id="page-16-29"></span>[94] T.Birkandan and M.Cvetiˇc, Class.Quant.Grav. 32 (2015) 085007.
- [95] M.S. Cunha and H.R.Christiansen, Phys. Rev. D, 84(2011) 085002 , [arXiv:1109.3486.](http://arxiv.org/abs/1109.3486)

- <span id="page-17-1"></span><span id="page-17-0"></span>[96] M.S.Cunha and H.R. Christiansen, European Phys. J. C 72 (2012)1942, [arXiv:1203.2172.](http://arxiv.org/abs/1203.2172)
- [97] G. Siopsis, [hep-th/0407157,](http://arxiv.org/abs/hep-th/0407157) Nucl. Phys. B715 (2005) 483.
- [98] L.Anguelova, P. Langfelder, J. High Energy Phys., 0303 (2003) 057.
- [99] S.R. Lau, Class. Quant. Grav., 21 (2004) 4147, also J.Computational Physics, 199 (2004) 376.
- [100] H.Suzuki, E.Takasugi and H.Umetsu, Prog. Theor. Phys., 100(1998) 491.
- [101] A.Zecca, Nuovo Cimento B, 125(2010) 191.
- [102] A.Ensico and N. Kamran, Comm. Math. Physics 290 (2009) 105.
- <span id="page-17-2"></span>[103] G.Esposito and R. Roychowdhury, Gen. Rel. Grav. 42 (2010)1221.
- <span id="page-17-3"></span>[104] S.Yoshida, N. Uchikata and T.Fumatase, Phys Rev. D 81(2010) 044005.
- <span id="page-17-4"></span>[105] C.Dariescu, M.A. Dariescu and C.Cretu, Int. J. Theor. Phys.,52 (2013) 1345.
- <span id="page-17-5"></span>[106] C.Dariescu and M.A. Dariescu, Modern Phys. Lett., 27 (2012) 1250184.
- <span id="page-17-7"></span>[107] C.Dariescu and M.A. Dariescu, Astrophysics and Space Science, 341 (2012) 429.
- <span id="page-17-6"></span>[108] H. Lu and, Zhao-Long Wang, Phys. Lett. B,718 (2013) 1536.
- <span id="page-17-8"></span>[109] Jun Li, Hai-Shan Liu, H. Lu ,Zhao-Long Wang, J. High Energy Phys., 1302 (2013) 109.
- <span id="page-17-9"></span>[110] D. Momeni, K. Yerzhanov and R. Myrzakulov, Canadian J. Phys., 90 Issue: 9 (2012) 877.
- <span id="page-17-10"></span>[111] A.M. Baranov, Gravitation and Cosmology, 18 (2012) 201.
- <span id="page-17-11"></span>[112] Dennis Philipp and Volker Perlick, [arXiv:1503.0810](http://arxiv.org/abs/1503.0810) (2015).
- <span id="page-17-12"></span>[113] H.S. Vieira, V.B. Bezerra, C.R. Muniz , Annals Phys.(N.Y.) 350 (2014) 14.
- <span id="page-17-13"></span>[114] V.B.Bezerra, H.S.Vieira , A.Costa, Class.Quant.Grav. 31 (2014) 045003.
- <span id="page-17-14"></span>[115] H.S.Vieira and V.B. Bezerra, J.Math.Phys. 56 (2015) 092501.
- <span id="page-17-15"></span>[116] D.Staicova and P.P.Fiziev, Astrophysics and Space Science,358,issue 1, (2015) 10.
- <span id="page-17-16"></span>[117] C. Dariescu and M.A. Dariescu, Chinese Physics Lett.32 (2015) 071101.
- <span id="page-17-17"></span>[118] C. Dariescu and M.A. Dariescu, Found. Phys,45 (2015) 1495.
- <span id="page-17-18"></span>[119] M.C.E.Cedeno and J.C.N. de Araujo, Phys. Rev. D 92 (2014)124015.
- <span id="page-17-19"></span>[120] D. Bini, E. Bittencour, A. Geralico, Class. Quant. Gravity 32 (2015) 215010.
- <span id="page-17-20"></span>[121] H.S.Vieira, V.B.Bezerra, G.V.Silva, Annals of Physics 362 (2015)576.
- <span id="page-17-21"></span>[122] D. Kofron, Phys. Rev. D92 (2015)124064, ArXiv: 1603.01452.
- <span id="page-17-22"></span>[123] D. Kofron, Phys. Rev. D93 (2016) 104012, ArXiv: 1604.05638.
- <span id="page-17-23"></span>[124] P.Valtancoli, Annals of Physics, 369, (2016) 161.
- <span id="page-17-24"></span>[125] H.S.Vieira and V.B. Bezerra, General Relativity and Gravitation 48, (2016) Issue 7, Article 88
- <span id="page-17-25"></span>[126] H.S.Vieira and V.B. Bezerra, Physical Review D, 94 (2016) 023511.
- <span id="page-17-26"></span>[127] H.S.Vieira and V.B. Bezerra, Annals of Physics, 373 (2016) 28.
- [128] I. Sakallı, Physical Review D, 94 (2016) 084040.

- <span id="page-18-1"></span><span id="page-18-0"></span>[129] G.V. Kraniotis, Class. and Quant. Gravity, 33 (2016) 225011.
- <span id="page-18-2"></span>[130] A. Arda, C. Tezcan, R. Sever, PRAMANA-Journal of Physics, 88(2017) Article number 39.
- <span id="page-18-3"></span>[131] H.S.Vieira,Int. J. Modern Phys. D 26 (2017) Article Number 1750035.
- <span id="page-18-4"></span>[132] H.S.Vieira, Chinese Physics C 41 (2017) Article Number 043105.
- <span id="page-18-5"></span>[133] H.S.Vieira, J.P.Morais Gra¸ca, V.B. Bezerra, Chinese Physics C 41 (2017) Article Number 095102.
- <span id="page-18-6"></span>[134] C. Dariescu, M.A. Dariescu, C.Stelea, Gen. Rel. Grav. 49 (2017) no.12 153.
- <span id="page-18-7"></span>[135] M.A. Dariescu, C. Dariescu, Int. J. Theor. Phys. 57 (2018) no.3 652.
- <span id="page-18-8"></span>[136] H.Sobhani, H.Hassanabdi, W.S.Chung, Eur. Phys. J C 78 (2018), no.2, 106.
- <span id="page-18-9"></span>[137] D.Bouaziz, T. Birkandan, Ann. of Phys. (N.Y.) 387 2017 62.
- <span id="page-18-10"></span>[138] T.Birkandan, M.Horta¸csu, Euro. Phys. Lett.119 no.2, 20002.
- <span id="page-18-11"></span>[139] T.Birkandan, M.Horta¸csu, Gen. Rel. Grav. 50 (2018) no.3, 28.
- <span id="page-18-12"></span>[140] G. 't Hooft, Phys. Rev. 14 (1976) 3432.
- <span id="page-18-13"></span>[141] T. Birkandan, M. Horta¸csu, Reports Math. Phys. 79 (2017) 81, [arXiv:1605.07848](http://arxiv.org/abs/1605.07848) [hep-th].
- [142] M. Horta¸csu, Proceedings of the 13th Regional Conference on Mathematical Physics, Antalya, Turkey,27-31 October 2010, ed. by U˘gur Camcı and ˙Ibrahim Semiz, pp. 23-39, World Scientific (2013).