# A software pa kage for the numeri al integration of ODE by means of high-order Taylor methods

# Angel Jorba (1) and Maorong Zou(2)

- (1) Departament de Matemati
  a Apli
  ada i Analisi, Universitat de Bar
  elona, Gran Via 585, 08007 Bar
  elona, Spain. E-mail: angelmaia.ub.es
- (2) Department of Mathemati
  s, The University of Texas at Austin,

To the memory of Wil liam F. S
helter

### Abstra t

This paper revisits the Taylor method for the numeri
al integration of initial value problems of Ordinary Dierential Equations (ODEs). The main issue is to present a 
omputer program that, given a set of ODEs, produ
es the 
orresponding Taylor numeri
al integrator. The step size 
ontrol adaptively sele
ts both order and step size to a
hieve a pres
ribed error, and trying to minimize the global number of operations. The pa
kage provides support for several extended pre
ision arithmeti
s, in
luding user-dened types.

The paper dis
usses the performan
e of the resulting integrator in some examples. As it 
an sele
t the order of the approximation used, it has a very good behaviour for high a

ura
y 
omputations. In fa
t, if we are interested in a very a

urate 
omputation in extended pre
ision arithmeti
, it be
omes the best 
hoi
e by far. The main drawba
k is that the Taylor method is an expli
it method, so it

sti systems.

| 1          |                              | Introdu tion                               | 3  |
|------------|------------------------------|--------------------------------------------|----|
| 2          |                              | A short summary on automati  dierentiation | 5  |
|            | 2.1                          | An example: The Van der Pol equation .<br> | 7  |
| 3          | Degree and step size  ontrol |                                            | 8  |
|            | 3.1                          | On the optimal step size sele tion<br><br> | 8  |
|            | 3.2                          | Dangerous step sizes<br><br>               | 9  |
|            | 3.3                          | A more pra ti al approa h<br><br>          | 10 |
|            |                              | 3.3.1<br>First step size  ontrol .<br><br> | 11 |
|            |                              | 3.3.2<br>Se ond step size  ontrol<br><br>  | 12 |
|            |                              | 3.3.3<br>User dened step size  ontrol<br> | 12 |
|            | 3.4                          | High a  ura y  omputations<br><br>         | 12 |
| 4          |                              | Software implementation                    | 13 |
|            | 4.1                          | Extended arithmeti<br><br>                 | 13 |
|            | 4.2                          | Step size  ontrol<br><br>                  | 14 |
| 5          |                              | Examples                                   | 15 |
|            | 5.1                          | The three-body problem<br><br>             | 15 |
|            | 5.2                          | Extended pre ision  al ulations<br><br>    | 18 |
|            | 5.3                          | Some  omparisons .<br><br>                 | 18 |
| Referen es |                              | 19                                         |    |

Let us 
onsider the following initial value problem: nd a smooth fun
tion <sup>x</sup> : [a; b℄ ! <sup>R</sup><sup>m</sup>

$$\begin{cases} x'(t) &= f(t, x(t)), \\ x(a) &= x_0, \end{cases}$$
 (1)

where <sup>f</sup> : [a; b℄ - <sup>R</sup> - <sup>R</sup><sup>m</sup> ! <sup>R</sup><sup>m</sup> and <sup>m</sup> 1. There is a 
lassi
al result of the theory of ODE that ensures the existen
e and uniqueness of a fun
tion x(t) satisfying (1). However, the ee
tive 
omputation of su
h a fun
tion x(t) is a mu
h more diÆ
ult question.

The sear
h of good numeri
al methods for (1) is one of the 
lassi
al problems in numeri
al analysis. The usual pro
edures for solving numeri
ally this problem are based in approximating the values x(t), for a suitable mesh of values of t. For the moment, and to simplify the presentation, we will use an equispa
ed mesh: if <sup>N</sup> <sup>2</sup> N, then we dene

$$h = \frac{b-a}{N}$$
,  $t_n = a + nh$ ,  $0 \le n \le N$ .

So, the problem is to nd approximations to the values x(tn). In what follows, we will denote these approximated values by xn.

In this paper we will revisit one of the oldest numeri
al pro
edures for this problem: the Taylor method. For simpli
ity, we will assume that the fun
tion <sup>f</sup> is analyti
 for (t; x) <sup>2</sup> [a; b℄ -. The idea of the method is very simple: given the initial 
ondition x(tn) = xn, the value x(tn+1) is approximated from the Taylor series of x(t) at <sup>t</sup> = tn. The algorithm is then,

$$x_0 = x(a),$$
  
 $x_{n+1} = x_n + x'(t_n)h + \frac{x''(t_n)}{2!}h^2 + \dots + \frac{x^{(p)}(t_n)}{p!}h^p, \quad n = 0, \dots, N-1.$  (2)

We refer to [But87℄ for a dis
ussion of the basi
 properties of the method. One of the key points for a pra
ti
al implementation is the ee
tive 
omputation of the values of the derivatives x(j) (tn). A rst pro
edure to obtain them is to dierentiate the equation in (1) w.r.t. t, at the point <sup>t</sup> = tn. Hen
e,

$$x'(t_n) = f(t_n, x(t_n)), \quad x''(t_n) = f_t(t_n, x(t_n)) + f_x(t_n, x(t_n))x'(t_n),$$

and so on. Hen
e, the rst step to apply this method is, for a given <sup>f</sup> , to 
ompute these derivatives up to a suitable order. Then, for ea
h step of the integration (see (2)), we have to evaluate these expressions to obtain the 
oeÆ
ients of the power series of x(t) at <sup>t</sup> = tn. Usually, these expressions will be very 
umbersome, so it will take a signi
ant amount of time to evaluate them numeri
ally. This, jointly with the initial eort to 
ompute those derivatives, is the main drawba
k of the Taylor method. As an example to illustrate these ideas, let us 
onsider the Van der Pol equation:

$$\begin{cases}
 x' &= y \\
 y' &= (1 - x^2)y - x
 \end{cases}
 .
 (3)$$

The second and third order derivatives of x and y with respect to time are

$$x'' = (1 - x^{2})y - x$$

$$y'' = x^{3} - x - 2xy^{2} + (x^{4} - 2x^{2})y$$

$$x''' = x^{3} - x - 2xy^{2} + (x^{4} - 2x^{2})y$$

$$y''' = 2x^{3} - x^{5} + (-1 + 5x^{2} + 3x^{4} - x^{6})y + (-8x + 4x^{3})y^{2} - 2y^{3}$$
(4)

Hence a third order Taylor method for the initial value problem (3) is

$$\begin{pmatrix} x_{n+1} \\ y_{n+1} \end{pmatrix} = \begin{pmatrix} x_n \\ y_n \end{pmatrix} + \begin{pmatrix} y_n \\ (1 - x_n^2)y_n - x_n \end{pmatrix} h$$

$$+ \frac{1}{2!} \begin{pmatrix} (1 - x_n^2)y_n - x_n \\ x_n^3 - x_n - 2x_n y_n^2 + (x_n^4 - 2x_n^2)y_n \end{pmatrix} h^2$$

$$+ \frac{1}{3!} \begin{pmatrix} x_n^3 - x_n - 2x_n y_n^2 + (x_n^4 - 2x_n^2)y_n \\ 2x_n^3 - x_n^5 + (-1 + 5x_n^2 + 3x_n^4 - x_n^6)y_n + (-8x_n + 4x_n^3)y_n^2 - 2y_n^3 \end{pmatrix} h^3$$

As one can see from these equations, expressions for higher order derivatives are quite complicated, and the complexity increases dramatically as the order increases.

This difficulty can be overcomed by using the so-called automatic differentiation ([BKSF59], [Wen64], [Moo66], [Ral81], [GC91], [BCCG92], [BBCG96], [Gri00]). This is a procedure that allows for a fast computation of the derivatives of a given function, up to arbitrarily high orders. As far as we know, these ideas were first used in Celestial Mechanics problems ([Ste56], [Ste57]; see also [Bro71]).

A first inconvenience of this method is that the function f has to belong to a special class; fortunately, this class is large enough to contain the functions that appear in most applications. A second difficulty is that the program to compute these derivatives has to be specifically coded for each function f.

The goal of this work is to present a software that, given a function f (belonging to a suitable class), generates the code that computes the jet of derivatives needed by the Taylor method up to any order.

For an efficient numerical integration, we need some knowledge of the order  $p \in \mathbb{N}$  up to which the derivatives have to be computed, and an estimate of the step size h, in order to have a truncation error below a given threshold value  $\varepsilon$ . We note that, as we have to select the value of two parameters (p and h), we can ask for a second condition besides the size of the truncation error. Here we have chosen to minimize the number of operations needed to advance the independent variable t in one unit. We have also coded the algorithms to do these tasks so that the output of the program is, in fact, a numerical integrator for the initial value problem (1).

We have tested this Taylor integrator against some well-known integration methods. The results show that Taylor method is a very competitive method to integrate with the standard double precision arithmetic of the computer. However, the main motivation for writing this software is to address the need of highly accurate computations in some problems of Dynamical Systems and Mechanics (see, for instance, [MS99], [SV01] and [Sim]). The standard methods of the numerical integration of ODEs become extremely slow for more accurate computations that require extended precision arithmetic. Here we

show that, if we need extended pre
ision for our 
al
ulations, Taylor method be
omes the best 
hoi
e by far. The reason for this advantage is, basi
ally, that Taylor method does not need to redu
e the step size to in
rease a

ura
y; it 
an simply in
rease the order (see Se
tion 3.4). As we will see, this allows to greatly redu
e the number of arithmeti operations during the numeri
al integration.

As any expli
it s
heme, the Taylor method is not suitable for sti equations be
ause, in this 
ase, the errors grow too fast. However, there are modi
ations of the Taylor method to deal with these situations ([Bar80℄, [JZ85℄, [KC92℄ and [CGH<sup>+</sup> 97℄). These

In the paper we present the main details of our implementation. We have tried to produ
e an eÆ
ient pa
kage, in the sense that the produ
ed Taylor integrator be as fast as possible. Moreover, we have also in
luded support for quadruple pre
ision and multiple pre
ision arithmeti
. We have done several test to 
ompare the eÆ
ien
y and a

ura
y of the generated Taylor routine against other numeri
al integrators.

There are several papers that fo
us on 
omputer implementations of the Taylor method in dierent 
ontexts; see, for instan
e, [BWZ70℄, [CC82℄, [CC94℄. A good survey is [NJC99℄ (see also [Cor95℄).

This pa
kage has been released under the GNU Publi
 Li
ense, so anybody with Internet a

ess is free to get it and to redistribute it. To obtain a 
opy, visit the URLs http://www.ma.utexas.edu/~mzou/taylor/ or http://www.maia.ub.es/~angel/. We note that the a
tual version of the pa
kage is written to run under the GNU/Linux operating system. We do not expe
t major problems to run it under any version of Unix, but we do not plan to write ports for other operating systems.

The paper has been split as follows: Se
tion 2 
ontains a survey about automati dierentiation, Se
tion 3 is devoted to the sele
tion of step size and trun
ation degree, Se
tion 4 gives some details about the software and Se
tion 5 provides some examples of the use of the pa
kage.

# 2 A short summary on automati dierentiation

Automati
 dierentiation is a re
ursive pro
edure to 
ompute the value of the derivatives of 
ertain fun
tions at a given point (see [Moo66, Ral81℄). The 
onsidered fun
tions are those that 
an be obtained by sum, produ
t, quotient, and 
omposition of elementary fun
tions.1 We will see this in detail.

To simplify the dis
ussion let us introdu
e the following notation: if <sup>a</sup> : <sup>t</sup> <sup>2</sup> <sup>I</sup> <sup>R</sup> 7! <sup>R</sup> denotes a smooth fun
tion, we 
all its normalized j-th derivative to the value

$$a^{[j]}(t) = \frac{1}{j!}a^{(j)}(t). \tag{5}$$

where <sup>a</sup>(j) (t0) denotes the <sup>j</sup> derivative w.r.t. t. In what follows, we will fo
us on the omputation of the values a[j℄ (t).

Elementary fun
tions in
lude polynomials, trigonometri
 fun
tions, real powers, exponentials and logarithms.

Assume now that a(t) = <sup>F</sup> (b(t); 
(t)) and that we know the values <sup>b</sup> [j℄ (t) and [j℄ (t), <sup>j</sup> = 0; : : : ; n, for a given t. The next proposition gives the n-th derivative of <sup>a</sup> at <sup>t</sup> for

Proposition 2.1 If the fun
tions <sup>b</sup> and are of 
lass <sup>C</sup><sup>n</sup> , and <sup>2</sup> <sup>R</sup> n f0g, we have

1. If 
$$a(t) = b(t) \pm c(t)$$
, then  $a^{[n]}(t) = b^{[n]}(t) \pm c^{[n]}(t)$ .

2. If 
$$a(t) = b(t)c(t)$$
, then  $a^{[n]}(t) = \sum_{i=0}^{n} b^{[n-i]}(t)c^{[i]}(t)$ .

3. If 
$$a(t) = \frac{b(t)}{c(t)}$$
, then  $a^{[n]}(t) = \frac{1}{c^{[0]}(t)} \left[ b^{[n]}(t) - \sum_{i=1}^{n} c^{[i]}(t) a^{[n-i]}(t) \right]$ 

4. If 
$$a(t) = b(t)^{\alpha}$$
, then  $a^{[n]}(t) = \frac{1}{nb^{[0]}(t)} \sum_{i=0}^{n-1} (n\alpha - i(\alpha+1)) b^{[n-i]}(t) a^{[i]}(t)$ .

5. If 
$$a(t) = e^{b(t)}$$
, then  $a^{[n]}(t) = \frac{1}{n} \sum_{i=0}^{n-1} (n-i) a^{[i]}(t) b^{[n-i]}(t)$ .

6. If 
$$a(t) = \ln b(t)$$
, then  $a^{[n]}(t) = \frac{1}{b^{[0]}(t)} \left[ b^{[n]}(t) - \frac{1}{n} \sum_{i=1}^{n-1} (n-i)b^{[i]}(t)a^{[n-i]}(t) \right]$ 

7. If a(t) = 
os (t) and b(t) = sin (t), then

$$a^{[n]}(t) = -\frac{1}{n} \sum_{i=1}^{n} ib^{[n-i]}(t)c^{[i]}(t)$$

$$b^{[n]}(t) = \frac{1}{n} \sum_{i=1}^{n} i a^{[n-i]}(t) c^{[i]}(t)$$

Proof: These proofs 
an be found in the literature, so we only give some hints about

- 2. It follows from Leibniz's formula:

$$a^{[n]}(t) = \frac{1}{n!}a^{(n)}(t) = \frac{1}{n!}\sum_{i=0}^{n} \binom{n}{i}b^{(n-i)}(t)c^{(i)}(t) = \sum_{i=0}^{n}b^{[n-i]}(t)c^{[i]}(t).$$

- 3. Apply item 2 to a(t)
  (t) = b(t).
- 4. Take logarithms and derivatives to obtain a0 (t)b(t) = a(t)b (t). Use item 2 and (5).

(t) = a(t)b

(t). Use item 2 and (5).

(t). Use item 2 and (5).

- 5. Take logarithms and derivatives to obtain a0
- 6. Take derivatives to obtain a0

(t)b(t) = <sup>b</sup>

7. Take derivatives to obtain a0 (t) = b(t) (t) and <sup>b</sup> (t) = a(t) (t). Use item 2 and (5).

Remark 2.1 It is possible to derive similar formulas for other fun
tions, like inverse trigonometri
 fun
tions.

Corollary 2.1 The number of arithmeti
 operations to 
ompute the derivatives up to order <sup>n</sup> of a ve
tor eld is of O(n2 ).

Proof: The number of operations to obtain the j-th derivative on
e we know the pre-<sup>P</sup> O(j) (see the previous formulas). Hen
e, the total number of operations is j=1 O(j) = O(n2 ).

Although these methods only allow for the derivation of a redu
ed set of fun
tions, we note that the fun
tions that 
an be 
oded in a 
omputer belong to this set.

## 2.1 An example: The Van der Pol equation

These rules 
an be applied re
ursively so that we 
an obtain re
ursive formulas for the derivatives of a fun
tion des
ribed by 
ombination of these basi
 fun
tions. As an example, we 
an apply them to the Van der Pol equation (3). To this end we de
ompose the righthand side of these equations in a sequen
e of simple operations:

$$\begin{array}{rcl}
 u_1 & = & x \\
 u_2 & = & y \\
 u_3 & = & u_1 u_1 \\
 u_4 & = & 1 - u_3 \\
 u_5 & = & u_4 u_2 \\
 u_6 & = & u_5 - u_1 \\
 x' & = & u_2 \\
 y' & = & u_6
 \end{array}$$

$$(6)$$

Then, we 
an apply the formulas given in Proposition 2.1 (items 1 and 2) to ea
h of the equations in (6) to derive re
ursive formulas for <sup>u</sup> , <sup>j</sup> = 1; : : : ; 6. Then, from the last two equations of (6) we have

$$x^{[n+1]} = \frac{1}{n+1}u_2^{[n]},$$
  
$$y^{[n+1]} = \frac{1}{n+1}u_6^{[n]}.$$

<sup>n</sup> + 1

The fa
tor <sup>1</sup> n+1 omes from the denition given in equation (5). Then, we 
an apply re
ursively these formulas for <sup>n</sup> = 0; 1; : : :, up to a suitable degree p, to obtain the jet of normalized derivatives for the solution of the ODE. Note that is not ne
essary to sele
t the value <sup>p</sup> in advan
e. Compare these formulas with the result of using formal derivatives on the right hand side of (3), as it has been done in the Introdu
tion.

The main task of the software we present is to read read the system of ODEs, to de ompose it in a sequen
e of basi
 operations, and to apply the formulas in Proposition 2.1 to this de
omposition. The result is an ANSI C routine that, given an initial 
ondition

and a nal degree, 
omputes normalized derivatives of the solution up to that degree.

# 3 Degree and step size control

Let us now focus on a single step of the integration procedure for equation (1). In what follows, to simplify the discussion we will assume that the step size is always positive (a similar discussion also holds for negative step sizes).

Assume that, for a given time  $t_n$ , the solution is at the point  $x_n$ , and that we want to compute the position of the trajectory for a new time  $t_{n+1} = t_n + h_n$  within a given error  $\varepsilon$ . We will not assume that  $t_{n+1}$  is fixed in advance, so we have to determine not only the degree of the Taylor expansion to be used but also the value  $h_n$ .

So, let us denote by  $\{x_n^{[j]}(t_n)\}_j$  the jet of normalized derivatives of the solution  $x_n$  of (1) that satisfies  $x_n(t_n) = x_n$ . Then, if  $h = t - t_n$  is small enough, we have

$$x_n(t) = \sum_{j=0}^{\infty} x_n^{[j]}(t_n) h^j.$$

The idea is to select a sufficiently small value  $h_n$  and a sufficiently large value p, such that the values

$$t_{n+1} \equiv t_n + h_n, \qquad x_{n+1} \equiv \sum_{j=0}^{p} x_n^{[j]}(t_n) h_n^j,$$

satisfy

$$||x_n(t_{n+1}) - x_{n+1}|| \le \varepsilon.$$

## 3.1 On the optimal step size selection

We are interested in values p and h such that a) the error is below a given value  $\varepsilon$ ; and b) the total number of operations of the numerical integration is as small as possible. To determine such values, we need some assumptions on the analyticity properties of the solution x(t). The following result can be found in [Sim01].

**Proposition 3.1** Assume that the function  $h \mapsto x(t_n + h)$  is analytic on a ball of radius  $\rho_n$ , and that there exists positive constants  $M_1$  and  $M_2$  such that

$$\frac{M_1}{\rho_n^j} \le |x_n^{[j]}| \le \frac{M_2}{\rho_n^j}, \qquad \forall n \in \mathbb{N}.$$
 (7)

Then, if the required accuracy  $\varepsilon$  tends to 0, the optimal value of h that minimizes the number of operations tends to

$$h_n = \frac{\rho_n}{e^2}$$
.

**Remark 3.1** Note that the optimal step size does not depend on the level of accuracy. The optimal order p is, in fact, the order that guarantees the required precision once the step size has been selected. If  $\varepsilon$  tends to 0, this order p behaves like

$$p_n = -\frac{1}{2} \ln \varepsilon.$$

Remark 3.2 The upper bound in (7) is the usual Cau
hy bound on the Taylor 
oeÆ
ients of an analyti
 fun
tion. Note that this upper bound is dire
tly related to the 
omputational eort, be
ause when the normalized derivatives are larger, we need more terms in the series to a
hieve a given a

ura
y. On the other hand, the lower bound in (7) is needed for the optimality of <sup>h</sup> w.r.t. the (total) number of operations: if the normalized derivatives 
ould be very smal l, then the step size obtained using the upper bound would be far from optimal.

Proof: Although the proof 
an be found in [Sim01℄, we will in
lude it here for the

The error introdu
ed when 
utting a Taylor series to a given degree <sup>p</sup> is of the order of the rst negle
ted term,

$$E \approx M \left(\frac{h}{\rho}\right)^{p+1}$$

$$h \approx \rho \left(\frac{\varepsilon}{M}\right)^{\frac{1}{p+1}},\tag{8}$$

On the other hand, the 
omputational eort to obtain the jet of normalized derivatives up to order <sup>p</sup> is O(p2 ) (p + 1)<sup>2</sup> (see Corollary 2.1). So, the (instantaneous) number of oating point operations per unit of time is given, in order of magnitude, by

<sup>M</sup>

$$\phi(p) = \frac{c(p+1)^2}{\rho\left(\frac{\varepsilon}{M}\right)^{\frac{1}{p+1}}},$$

and, solving <sup>0</sup> (p) = 0, we obtain

$$p = -\frac{1}{2} \ln \left( \frac{\varepsilon}{M} \right) - 1$$

So, we will use this order with the largest step size that keeps the lo
al error below ": inserting this value of <sup>p</sup> in (8) we have

$$h = \frac{\rho}{e^2}$$

Finally, note that the optimal order <sup>p</sup> behaves like <sup>1</sup> ln " when " goes to zero.

There are strategies to use step sizes that are larger than the radius of 
onvergen
e of the series, but we will not dis
uss them here. See [CC82℄ for more details.

# 3.2 Dangerous step sizes

An important point to take into a

ount is that, after the step size and the order have been sele
ted, we have to sum the Taylor series to obtain the nal result. Although adding a series seems a trivial numeri
al operation, there are 
ases in whi
h this is a very diÆ
ult problem. For instan
e, 
onsider the initial value problem

$$\dot{x} = -x, \qquad x(0) = 1,\tag{9}$$

and assume we are interested in 
omputing x(10) by numeri
al integration of the ODE (as this problem 
an be easily solved by hand, we already know that the solution is given by x(t) = exp(t), so x(10) = exp(10) 0:0000454). Using (9), the Taylor series at x(0) is very simple to obtain:

$$x(h) = 1 + \sum_{n>1} (-1)^n \frac{h^n}{n!}$$

Due to the entire 
hara
ter of this fun
tion, the optimal step size is <sup>h</sup> = 10, and the degree is sele
ted to have a trun
ation error smaller than a given pre
ision. However, using su
h a large <sup>h</sup> introdu
es strong 
an
ellations in the sum of the series, that lead to the loss of many signi
ant digits when the sum is 
arried out with 
oating point arithmeti
.

The solution for this problem is to use smaller values of h, to avoid these 
an
ellations. The key point is to use a step size su
h that the series is de
reasing in modulus. For instan
e, in this 
ase we should start with <sup>h</sup> = 1. For this reason, on
e the step size has to redu
e <sup>h</sup> when ne
essary. This will be dis
ussed again in Se
tion 3.3.2.

## 3.3 A more pra ti al approa h

The main drawba
k of the previous approa
h is that it requires information that we 
annot obtain easily, like the radius of 
onvergen
e of the Taylor series. Hen
e, in this se
tion we will dis
uss some modi
ations to simplify the numeri
al implementation.

relative error. This is a question we 
annot answer in a general way, so we have provided is a

omplished by asking the user for both toleran
es, and the 
ode looks for an step size that tries to satisfy both of them.

When 
ontrolling the relative error, there is a dangerous situation that we should avoid. Note that the usual denition of relative error "r,

$$\varepsilon_r = \frac{\varepsilon_a}{|x_n|},$$

assumes that the solution xn is dierent from zero. Hen
e, if some of the points xn are zero (or very small), we should avoid the use of the relative error for these points, sin
e this 
ould result in an undened (or ridi
ulously small) step size. Therefore, we have used the following estimate for the relative error,

$$\varepsilon_r = \frac{\varepsilon_a}{\max\{|x_n|, |\dot{x}_n|\}},\tag{10}$$

where <sup>j</sup>:j denotes the sup norm of a ve
tor, and \_xn <sup>f</sup> (tn; xn). There is still the possibility that the origin be an equilibrium point and that our solution is approa
hing it. We do not deal with this 
ase, sin
e there is no a unique answer for all the possible situations error 
ontrol that he wants (see also Se
tion 4.2).

Finally let us note that, as all the step size controls presented are local, the values used to control the absolute and relative errors used in one part of the phase space can be meaningless in a different region. As before, the user is responsible for setting these values properly, and to change them when needed.

### 3.3.1 First step size control

The simplest step size control we have implemented is based on a very simple estimation of the radius of convergence of the Taylor series from the size of the last two terms (a more sophisticated approach can be found in [CC82]). Then, the order and step size are estimated using the formulas of Proposition 3.1. To explain the details, let us define  $(t_n, x_n)$  as the actual point of the orbit, and let us call  $\varepsilon_a$  and  $\varepsilon_r$  the absolute and relative required accuracies. Moreover, we define

$$\varepsilon_n = \min \left\{ \frac{\varepsilon_a}{\max\{|x_n|, |\dot{x}_n|\}}, \ \varepsilon_r \right\}, \qquad p_n = -\frac{1}{2} \ln \varepsilon_n + 1.$$
(11)

Then, if we denote

$$ar{
ho}_j = \left(|x_n^{[j]}|\right)^{-rac{1}{j}}, \quad \hat{
ho}_j = \left(rac{|x_n^{[j]}|}{\max\{|x_n|,|\dot{x}_n|\}}
ight)^{-rac{1}{j}}, \quad j=1,\ldots,p_n,$$

we estimate the radius of convergence of the series as

$$\rho_n = \min\{\bar{\rho}_{p_n-1}, \, \bar{\rho}_{p_n} \, \hat{\rho}_{p_n-1}, \, \hat{\rho}_{p_n}\},\,$$

so the computed step size is

$$h_n = \pm \frac{\rho_n}{e^2}.$$

The sign "-" is for a backwards integration.

As these values are approximations to the asymptotics given in Proposition 3.1, we are not sure if they guarantee the required precision. For this reason, we want to note that, for  $j = p_n - 1$  and  $j = p_n$ ,

$$\begin{aligned} |x_n^{[j]}h_n^j| & \leq & \frac{|x_n^{[j]}|\rho_n^j}{e^{2j}} \leq \frac{1}{e^{2j}}, \\ \frac{|x_n^{[j]}h_n^j|}{\max\{|x_n|,|\dot{x}_n|\}} & \leq & \frac{|x_n^{[j]}|\rho_n^j}{|x_n|e^{2j}} \leq \frac{1}{e^{2j}}, \end{aligned}$$

and, if  $j = p_n - 1$ ,

$$\frac{1}{e^{2j}} = \frac{1}{e^{2(p_n - 1)}} \approx \frac{1}{e^{-\ln \varepsilon_n}} = \varepsilon_n.$$

This means that the terms of order  $p_n - 1$  in the Taylor series have a contribution of order  $\varepsilon_n$  (either in absolute value or relative to the size of the initial condition) while the terms of order  $p_n$  (the last terms to be considered) have a contribution of order  $\varepsilon_n/e^2$ . Hence, this shows that the proposed strategy is similar to the more straightforward method of looking for an  $h_n$  such that the last terms in the series are of the order of the error wanted.

### 3.3.2 Se ond step size ontrol

This is based in 
orre
ting the previous method to avoid too large step sizes, as it has been explained in Se
tion 3.2. There are several possibilities to implement this approa
h. First, we want to note that it 
an be very dangerous to ask for a stri
tly de
reasing series, sin
e this 
ould lead to a very drasti
 (and unne
essary) step redu
tions. A typi
al example of this situation o

urs when the solution x(t) has all the odd (or even) Taylor terms equal to zero. Hen
e, to avoid these situations we have used a weaker 
riterion: let h0 be the step size 
ontrol obtained using the methods of Se
tion 3.3.1, and let <sup>h</sup> be the largest value su
h that <sup>h</sup> h0 and that

$$|x^{[0]}| + h_0|x^{[1]}| \ge |x^{[j]}|h^j, \qquad j = 2, \dots, p.$$
 (12)

Note that, in many 
ases, it is enough to take <sup>h</sup> <sup>=</sup> h0 to meet this 
ondition. On the other hand, when dealing with entire fun
tions, this avoids sele
ting too large step sizes that 
ould lead to a dramati
 loss of pre
ision (we re
all the example in Se
tion 3.2).

### 3.3.3 User dened step size ontrol

All the implemented step sizes are based on rough estimations of the error and, for this reason, they 
annot be taken as \rigorous".

One of the main uses of Taylor methods is the soalled veried integration, or integration with guaranteed error bounds. This means that it is possible to produ
e rigorous bounds on the trun
ation error so that, working with interval arithmeti
, the tra je
tory is omputed with true error bounds (see, for instan
e, [CC82℄, [CC94℄, [Cor95℄ and [Hoe01℄). This is not the approa
h we have taken here. The proposed step size 
ontrols are more in the lines of the 
lassi
al numeri
al analysis, in whi
h the error is estimated by means of asymptoti
 formulas ([But87℄). The main reason for that is to 
ompete, both in speed and a

ura
y, with the 
lassi
al numeri
al integrators that one 
an nd in numeri
al libraries.

However, our implementation 
an also be used for dierent purposes. In our opinion, the main point of this pa
kage is that, given a ve
tor eld, produ
es an eÆ
ient C routine that 
omputes the jet of normalized derivatives at a given initial 
ondition (this is usually the most tedious part of building a Taylor series integrator). So, it is not diÆ
ult to use this routine with an interval arithmeti
, plus a suitable step size 
ontrol to integrate a given ODE with guaranteed error bounds. We also oer to the user the possibility of using his own error 
ontrol algorithm. This 
an be very useful in some spe
i

ases where some spe
ial properties of the solution are known (for instan
e, when the solution is known to be a polynomial). See the do
umentation that 
omes with the software for more details

# 3.4 High a ura y omputations

An important property of high order Taylor integrators is their suitability for high a

ura
y 
omputations.

Assume that we are solving an IVP like (1) and that, at a given step, we are using an step size h << 1 and an order <sup>p</sup> to obtain a lo
al error " << 1. The number of operations

needed to 
ompute all the derivatives is of O(p2 ) (see Corollary 2.1). As the number of operations to sum the power series is only O(p), the total operation 
ount for a single step of the Taylor method is of O(p2 ). Hen
e, if we want to in
rease the a

ura
y to, say, " we 
an simply in
rease the order of the Taylor method to 2p so the number of operations is in
reased by a fa
tor 4. Note that, if we want to a
hieve the same level of a

ura
y not by in
reasing the order but by redu
ing the step size h, we have to use an step size of <sup>h</sup><sup>2</sup> . This means that we will have to use h1 steps (of size <sup>h</sup><sup>2</sup> ) to 
ompute the orbit after <sup>h</sup> units of time so the total number of operations is, roughly speaking, in
reased by a fa
tor of h1 . For instan
e, assume that the trun
ation error is exa
tly <sup>h</sup><sup>p</sup> . If " = 1016 and <sup>p</sup> = 8, then the step size has to be <sup>h</sup> = 0:01. Note that, if <sup>p</sup> is xed, to a
hieve an a

ura
y of 1032 we have to use <sup>h</sup> = 104 , that for
es to use 100 times more steps (hen
e, 100 times more operations) than for the " = 1016 ase. Changing the value of <sup>p</sup> from 8 to 16 allows to keep the same step size <sup>h</sup> = 0:01 while the 
omputational eort required to obtain the derivatives is only in
reased by a fa
tor 4. If the required pre
ision were higher, these dieren
es would be even more dramati
.

Hen
e, it requires mu
h less work to in
rease the order rather than to redu
e the step size (this observation was already impli
it in Proposition 3.1, where it was shown that the optimal step size was independent from the level of a

ura
y required).

If we are using a numeri
al integrator of xed order (for instan
e, a Runge-Kutta or Adams-Bashford method), the only possibility to in
rease a

ura
y is to redu
e the step size. Hen
e, xed order methods are strongly penalized for high a

ura
ies, 
ompared with varying order methods. For this reason, if the required a

ura
y is high enough, Taylor methods {with varying order{ are the best option by far. This will be illustrated

# 4 Software implementation

The basi
 operations of the taylor translator are rst to parse the dierential equations to redu
e them to a sequen
e of binary operations and 
alls to mathemati
al fun
tions, and then to apply the rules of automati
 dierentiation as it has been explained in Se
tion 2.

A point we think is worth noting is the way we deal with the basi
 arithmeti
. When the Taylor translator writes the 
ode for the jet of derivatives and/or the step size 
ontrol, it de
lares all the real variables as MY FLOAT (we will explain this type later on), and substitutes ea
h mathemati
al operation by a suitable ma
ro; the name of these ma
ros is independent from the arithmeti
. In other words, if we use the taylor program to only generate the jet of derivatives and/or the step size 
ontrol, we do not need to spe
ify the kind of arithmeti
 we want, sin
e the output does not depend on that.

The denition of the type MY FLOAT and the body of the ma
ros is 
ontained in the header le. This le is produ
ed by the 
ag -header plus a 
ag spe
ifying the arithmeti wanted. For instan
e, the ma
ro for the produ
t of two gmp numbers is

```
/* multipli
ation r=a*b */
#define MultiplyMyFloatA(r,a,b) mpf_mul(r,(a), (b))
```

Here, mpf mul is a gmp fun
tion that multiplies the two numbers a and b and stores the result in r. Then, the C prepro
essor will substitute the ma
ros by the 
orresponding 
alls to the arithmeti
 library. Hen
e, to use an arithmeti
 dierent than the ones provided here we only have to modify the header le. For more details, see the manual that 
omes

The pa
kage in
ludes support for several extended pre
ision arithmeti
s, namely

- doubledouble This is a C++ library that denes a extended 
  oat type, in whi
  h ea
  h number is stored as the sum of two double numbers. The a

  ura
  y is then of nearly 30 de
  imal digits. The standard way of using this library is by means of overloading. See http://www.btexa
  t.
  om/people/briggsk2/doubledouble.html
- dd real, qd real This is also a C++ library, similar to doubledouble, that denes the types dd real (2 doubles) and qd real (4 doubles), providing a

  ura
  ies of nearly 32 and 64 de
  imal digits, respe
  tively. For more details, visit the URL http://www.ners
  .gov/~dhbailey/mpdist/mpdist.html
- GNU Multiple Pre
  ision Library (gmp) This is the standard GNU library for extended pre
  ision. This library allows to dene arbitrary long integer, rational and real types, and to operate on them by means of fun
  tion 
  alls (more details on the library 
  an be found in http://www.swox.
  om/gmp/). Unfortunately, this library does not provide trans
  endental fun
  tions so, in prin
  iple, we are restri
  ted to ve
   tor elds that 
  an be written with the basi
   arithmeti
   fun
  tions plus square root (those are the only fun
  tions provided for 
  oating point types). For more sophisti
  ated ve
  tor elds, the user has to 
  ode the 
  orresponding fun
  tions.

None of these 
oating point libraries is in
luded in our pa
kage.

# 4.2 Step size ontrol

We have in
luded two implementations of the previously mentioned step size methods. Moreover, the user 
an also provide his own step size 
ontrol.

The rst step 
ontrol we provide is the dire
t 
oding of the method explained in Se
tion 3.3.1. The Taylor series is generated up to the order pn dened in (11), and the step size is 
omputed as explained in Se
tion 3.3.1. Finally, the series is added using the Horner method. The se
ond step size 
ontrol is based on a modi
ation of the previous ontrol, as explained in Se
tion 3.3.2. The main dieren
e is a loop that runs over all the oeÆ
ients of the Taylor series and 
he
ks for the 
ondition (12). When it is not satised, <sup>h</sup> is redu
ed a

ordingly.

We have also in
luded a user dened step size 
ontrol. In this 
ase, the user must 
ode its own step 
ontrol. More details 
an be found in the do
umentation that 
omes with

```
/* ODE spe
ifi
ation: rtbp */
mu=0.01;
umu=1-mu;
r2=x1*x1+x2*x2+x3*x3;
rpe2=r2-2*mu*x1+mu*mu;
rpe3i=rpe2^(-3./2);
rpm2=r2+2*(1-mu)*x1+(1-mu)*(1-mu);
rpm3i=rpm2^(-3./2);
diff(x1, t)= x4+x2;
diff(x2, t)= x5-x1;
diff(x3, t)= x6;
diff(x4, t)= x5-(x1-mu)*(umu*rpe3i)-(x1+umu)*(mu*rpm3i);
diff(x5, t)=-x4-x2*(umu*rpe3i+mu*rpm3i);
diff(x6, t)=-x3*(umu*rpe3i+mu*rpm3i);
```

Figure 1: Input le for the restri
ted three-body problem.

# 5 Examples

In this se
tion we have in
luded a few examples to show how Taylor method works. We have assumed that the taylor translator is already installed in your 
omputer. For details about how to do this, 
he
k the do
umentation that 
omes with the pa
kage.

# 5.1 The three-body problem

This is, probably, the best known problem in Celestial Me
hani
s. Here we will use a simplied version of it, the soalled restri
ted three-body problem (for details about this problem see, for instan
e, [Sze67℄). The problem boils down to des
ribe the dynami
s of the Hamiltonian system

$$H_{RTBP} = \frac{1}{2}(p_x^2 + p_y^2 + p_z^2) + yp_x - xp_y - \frac{1-\mu}{r_{PS}} - \frac{\mu}{r_{PJ}}$$

being a mass parameter, <sup>r</sup> P S = (<sup>x</sup> ) + <sup>y</sup><sup>2</sup> + <sup>z</sup> P J = (<sup>x</sup> + 1)2 <sup>+</sup> <sup>y</sup><sup>2</sup> <sup>+</sup> <sup>z</sup> input le for the 
orresponding ve
tor eld is shown in Figure 1.

There are several kind of instru
tion in this le. First of all, anything between /\* \*/ is ignored by the program, so we 
an use them to put 
omments in the le. Next, we have some lines to dene numeri
al 
onstants, plus some operations with the variables of the system. The variables of the equation are labeled x1, x2 and so on, and the independent variable is labeled t. Finally, the last 6 lines are the denition of the dierential equations.

Although we think that the notation used is 
lear and self-expli
ative, we want to make some 
omments about it. First, the taylor translator does not make any kind of optimization on the input des
ription of the ve
tor eld, with the ex
eption of 
ommon expression eliminations. If one of your main 
on
erns is the eÆ
ien
y of the 
ode generated by taylor, you should apply other kind of optimizations \by hand" in your input le (for instan
e, to simplify algebrai
 expressions to minimize the number of operations).

A se
ond point we want to 
omment on is the use of the exponent \3:0=2" in the expressions. There are several ways of introdu
ing su
h an exponent. If we use the expression \1:5", the program will use the exp and ln fun
tions to dene it (this is true for any real exponent). If we use \3:0=2", then we 
an use the 
ag \-sqrt" of the translator to for
e the program to use the square root fun
tion instead of the exp and ln fun
tions. Without this 
ag, the value \3:0=2" is treated as \1:5".

The input le supports more features than the ones showed here (like the use of extern variables to re
eive parameters from the user's programs); for details 
he
k the manual that 
omes with the pa
kage.

To produ
e a numeri
al integrator for this ve
tor eld, assume that we have the 
ode of Figure 1 in the le rtbp.in. Then, you 
an type

```
taylor -name rtbp -o taylor_rtbp.
 -step -jet -sqrt rtbp.in
taylor -name rtbp -o taylor.h -header
```

(we have assumed that the taylor binary is in a dire
tory 
ontained in your path; otherwise you should spe
ify its lo
ation). The rst line outputs the le taylor rtbp. with the 
ode for the step size 
ontrol and the jet of derivatives. The se
ond line produ
es the header le; it is needed for the le taylor rtbp.
, and the user may also want to in
lude it in the 
alling routine, sin
e it 
ontains the prototype for the 
all to the integrator. There are more options to 
ontrol the output of taylor, look at the do
umentation for

Fortran 77 users 
an use a single instru
tion:

```
taylor -name rtbp -o taylor_rtbp.
 -step -jet -f77 -header -sqrt rtbp.in
```

This will put everything in the le taylor rtbp.
, so you 
an simply 
ompile and link it with your (Fortran) 
alling routine. For details about how to 
all the Taylor integration

Then, you 
an 
all the routine taylor\_step\_rtbp (with suitable parameters), to perform a numeri
al integration of the previous ve
tor eld. As it is usual in one step expli
it methods, ea
h 
all advan
es the independent variable in some amount that depends on the level of a

ura
y required. The details about this 
all (parameters, et
.) 
an be found

For instan
e, let us show an example of numeri
al integration with this routine by sele
ting the initial 
ondition x1=-0.45, x2=0.80, x3=0.00, x4=-0.80, x5=-0.45 and x6=0.58. We will perform rst a numeri
al integration with the standard double pre
ision of the 
omputer, for 1 unit of time. As the Hamiltonian fun
tion is a rst integral of the system, we will 
he
k the preservation of its value as an indi
ator of the a

ura
y of the integration. So, we have 
oded a small main program that uses this initial 
ondition to all the Taylor integrator till the time has advan
ed in one unit, asking for a relative error less than 1017. The results are shown in Figure 2. The three 
olumns that appear in the middle of the table are the following: the rst 
olumn is a 
ounter that is in
reased

```
numeri
al integration starts...
```

Figure 2: Integration, for one unit of time, of the Restri
ted Three-Body Problem. See

![](_page_16_Figure_4.jpeg)

Figure 3: Long term behaviour of the energy. Horizontal axis: time. Verti
al axis: variation of the value of the Hamiltonian, in multiples of the epsilon of the ma
hine.

ea
h time the main program 
alls the taylor integrator; the se
ond 
olumn is the value of the time variable, the third 
olumn is the order used in the Taylor method, and the fourth value is the dieren
e between the a
tual value of <sup>H</sup> and the 
orresponding value for the initial 
ondition, in multiples of the epsilon of the ma
hine. the rst step of the integrator, the energy has 
hanged in only one bit. In the next three alls, the energy 
hanges in an amount of two times the ma
hine epsilon, whi
h gives an idea of the extremely good behaviour of the rounding errors. The program stops exa
tly

An interesting experiment is to perform a longer integration, to see the long term behaviour of the error. This is shown in Figure 3 for 106 units of time. Although the total amount of deviation is very small, we note that it does not behave as a random walk, sin
e the plot shows a 
lear drift on the variation of the Hamiltonian. This is a weakness that should be taken into a

ount in very long integrations.

The epsilon of the ma
hine is the smallest number we 
an add to 1 su
h that the result is dierent from 1. It is 
ommonly used as a measure of the level of a

ura
y of a given arithmeti
.

```
numeri
al integration starts...
```

Figure 4: Numeri
al integration of the ve
tor eld in Figure 1, using gmp with a 512 bits mantissa and asking for a relative error of 10150. See the text for more details.

## 5.2 Extended pre ision al ulations

We will use the same example as in the previous se
tion. The main dieren
e when generating the 
ode for the taylor integrator with extended pre
ision is that we have to tell to the taylor translator that we want to use an extended arithmeti
. For instan
e, to generate 
ode using gmp, we 
an do

```
taylor -name rtbp -o taylor_rtbp.
 -step -jet -sqrt rtbp.in
taylor -name rtbp -o taylor.h -gmp -header
```

(note that we only need to spe
ify the kind of arithmeti
 when generating the header le). As a rst test, we 
an repeat the run shown in Figure 2 using a 512 bits mantissa and asking for a relative error of 10150. The results are displayed in Figure 4 where, in the fourth 
olumn, we have repla
ed the ma
hine epsilon by the dieren
e in the value of <sup>H</sup> with respe
t to its initial value.

The output for a relative error of 10300 (and using a mantissa of 1024 bits) is also shown in Figure 5. We note that the step sizes are very similar { the dieren
es 
ome from dierent estimates of the radius of 
onvergen
e of the Taylor series { and the most relevant 
hange is the order of the method.

## 5.3 Some omparisons

Let us start by 
omparing the implementation of the Taylor method presented here with the well-known Runge-Kutta-Fehlberg of orders 7 and 8 (from now on, rk78), using the standard double pre
ision of the 
omputer. First, let us repeat the 
al
ulation showed in Figure 3. We have asked the rk78 for an a

ura
y of 1016. Note that the step sizes ontrols of rk78 and Taylor method are totally dierent; from our numeri
al experiments, it seems that the best results (i.e., the best a

ura
y) for the rk78 in this problem are obtained by asking for a lo
al error 1016. The results, both for the rk78 and the Taylor method, 
an be summarized as follows: the error in the value of <sup>H</sup> after 10<sup>6</sup> time, in number of multiples of the ma
hine epsilon, is 13412 for the rk78, and 217 for the Taylor method. The time taken for the rk78 was of 9m and 34s, while the Taylor

numeri
al integration starts...

Figure 5: As Figure 4, but using a 1024 bits mantissa and asking for a relative error of 10300

method needed 4m and 4s (these runs were done in a GNU/Linux environment, running on a Pentium II pro
essor at 500 MHz). If we ask the rk78 for an a

ura
y of 1014 then the time taken goes down to 4m and 58s, but the nal error is 649368 times the epsilon of the ma
hine (that is, 1:44 -1010).

Next we des
ribe a se
ond ben
hmark using the standard quadruple pre
ision (the type long double in C) of a HP 9000/712 
omputer, with a 100 MHz PA-RISC 1.1 pro
essor. In this 
ase, the long double type is 16 bytes long, holding nearly 33 de
imal digits { we note that this is not the usual situation for many C 
ompilers, in whi
h the type long double is only 10 bytes long, holding nearly 19 de
imal digits. We have used the ve
tor eld of the Restri
ted Three-Body Problem, with the same initial 
ondition and mass parameter as before. The integration time has been restri
ted to 10 units, to avoid long testing times. We have asked for a lo
al error of 1032 for the rk78, and of the 1033 for the Taylor method. The total 
pu time taken for the rk78 is of 3m and 48s, while the Taylor method only needed 4.1 se
onds.

# A knowledgements

This work has been supported by the Comision Conjunta Hispano Norteameri
ana de Coopera
ion Cient
a y Te
nologi
a. A.J. has also been supported by the Spanish CICYT grant BFM2000{0623, the Catalan CIRIT grant 2000SGR{00027 and DURSI.

[Bar80℄ D. Barton. On Taylor series and sti equations. ACM Trans. Math. Software, 6(3):280{294, 1980.

- [BBCG96℄ M. Berz, C. Bis
  hof, G.F. Corliss, and A. Griewank, editors. Computational Dierentiation: Te
  hniques, Appli
  ations, and Tools. SIAM, Philadelphia,
- [BCCG92℄ C.H. Bis
  hof, A. Carle, G.F. Corliss, and A. Griewank. ADIFOR: Automati dierentiation in a sour
  e translation environment. In Paul S. Wang, editor, Pro
  eedings of the International Symposium on Symboli
   and Algebrai
   Computation, pages 294{302, New York, 1992. ACM Press.
- [BKSF59℄ L.M. Beda, L.N. Korolev, N.V. Sukkikh, and T.S. Frolova. Programs for automati
   dierentiation for the ma
  hine BESM. Te
  hni
  al Report, Institute for Pre
  ise Me
  hani
  s and Computation Te
  hniques, A
  ademy of S
  ien
  e, Mos
  ow, USSR, 1959. (In Russian).
- [Bro71℄ R. Brou
  ke. Solution of the N-Body Problem with re
  urrent power series. Celestial Me
  h., 4(1):110{115, 1971.
- [But87℄ J.C. But
  her. The Numeri
  al Analysis of Ordinary Dierential Equations. Wiley, 1987.
- [BWZ70℄ D. Barton, I.M. Willers, and R.V.M. Zahar. The automati
   solution of ordinary dierential equations by the method of Taylor series. Computer J., 14(3):243{248, 1970.
- [CC82℄ G.F. Corliss and Y.F. Chang. Solving ordinary dierential equations using Taylor series. ACM Trans. Math. Software, 8(2):114{144, 1982.
- [CC94℄ Y.F. Chang and G.F. Corliss. ATOMFT: Solving ODEs and DAEs using Taylor series. Computers and Mathemati
  s with Appli
  ations, 28:209{233,
- [CGH<sup>+</sup> 97℄ G.F. Corliss, A. Griewank, P. Henneberger, G. Kirlinger, F.A. Potra, and H.J. Stetter. High-order sti ODE solvers via automati
   dierentiation and rational predi
  tion. In Numeri
  al analysis and its appli
  ations (Rousse, 1996), pages 114{125. Springer, Berlin, 1997.
- [Cor95℄ G.F. Corliss. Guaranteed error bounds for ordinary dierential equations. In M. Ainsworth, J. Levesley, W. A. Light, and M. Marletta, editors, Theory of Numeri
  s in Ordinary and Partial Dierential Equations, pages 1{75. Oxford University Press, Oxford, 1995. Le
  ture notes for a sequen
  e of ve le
  tures at the VI-th SERC Numeri
  al Analysis Summer S
  hool, Lei
  ester University, 25 - 29 July, 1994.
- [GC91℄ A. Griewank and G.F. Corliss, editors. Automati
   Dierentiation of Algorithms: Theory, Implementation, and Appli
  ation. SIAM, Philadelphia,
- [Gri00℄ A. Griewank. Evaluating Derivatives. SIAM, Philadelphia, Penn., 2000.

[Hoe01℄ J. Hoefkens. Rigorous numeri
al analysis with high order Taylor methods. PhD thesis, Mi
higan State University, 2001.

- [JZ85℄ F. Jalbert and R.V.M. Zahar. A highly pre
  ise Taylor series method for sti ODEs. In Pro
  eedings of the fourteenth Manitoba 
  onferen
  e on numeri
  al mathemati
  s and 
  omputing (Winnipeg, Man., 1984), volume 46, pages 347{
- [KC92℄ G. Kirlinger and G.F. Corliss. On impli
  it Taylor series methods for sti ODEs. In Computer arithmeti
   and en
  losure methods (Oldenburg, 1991), pages 371{379. North-Holland, Amsterdam, 1992.
- [Moo66℄ R.E. Moore. Interval Analysis. Prenti
  e-Hall, Englewood Clis, N.J., 1966.
- [MS99℄ R. Martnez and C. Simo. Simultaneous binary 
  ollisions in the planar fourbody problem. Nonlinearity, 12(4):903{930, 1999.
- [NJC99℄ N.S. Nedialkov, K.R. Ja
  kson, and G.F. Corliss. Validated solutions of initial value problems for ordinary dierential equations. Appl. Math. Comput., 105(1):21{68, 1999.
- [Ral81℄ L.B. Rall. Automati
   Dierentiation: Te
  hniques and Appli
  ations, volume 120 of Le
  ture Notes in Computer S
  ien
  e. Springer Verlag, Berlin, 1981.
- [Sim℄ C. Simo. Dynami
  al properties of the gure eight solution of the three{body problem. To appear in A. Chen
  iner, R. Cushman, C. Robinson and Z. Xia, editors, Pro
  eedings of the Chi
  ago Conferen
  e dedi
  ated to Don Saari.
- [Sim01℄ C. Simo. Global dynami
  s and fast indi
  ators. In H.W. Broer, B. Krauskopf, and G. Vegter, editors, Global analysis of dynami
  al systems, pages 373{389, Bristol, 2001. IOP Publishing.
- [Ste56℄ J.F. Steensen. On the restri
  ted problem of three bodies. Danske Vid. Selsk. Mat.-Fys. Medd., 30(18):17, 1956.
- [Ste57℄ J.F. Steensen. On the problem of three bodies in the plane. Mat.-Fys. Medd. Danske Vid. Selsk., 31(3):18, 1957.
- [SV01℄ C. Simo and C. Valls. A formal approximation of the splitting of separatri
  es in the 
  lassi
  al Arnold's example of diusion with two equal parameters. To appear in Nonlinearity, 2001.
- [Sze67℄ V. Szebehely. Theory of Orbits. A
  ademi
   Press, 1967.
- [Wen64℄ R. E. Wengert. A simple automati
   derivative evaluation program. Comm. ACM, 7(8):463{464, 1964.