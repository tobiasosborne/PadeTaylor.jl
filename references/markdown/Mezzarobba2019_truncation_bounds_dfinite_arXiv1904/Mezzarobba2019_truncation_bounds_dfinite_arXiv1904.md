# MULTIPLIERS AND INTEGRATION OPERATORS BETWEEN CONFORMALLY INVARIANT SPACES

# DANIEL GIRELA AND NOEL MERCHAN´

Abstract. In this paper we are concerned with two classes of conformally invariant spaces of analytic functions in the unit disc D, the Besov spaces B<sup>p</sup> (1 ≤ p < ∞) and the Q<sup>s</sup> spaces (0 < s < ∞). Our main objective is to characterize for a given pair (X, Y ) of spaces in these classes, the space of pointwise multipliers M(X, Y ), as well as to study the related questions of obtaining characterizations of those g analytic in D such that the Volterra operator T<sup>g</sup> or the companion operator I<sup>g</sup> with symbol g is a bounded operator from X into Y .

# 1. Introduction

Let D = {z ∈ C : |z| < 1} denote the open unit disc of the complex plane C and let Hol(D) be the space of all analytic functions in D endowed with the topology of uniform convergence on compact subsets.

If 0 < r < 1 and f ∈ Hol(D), we set

$$M_p(r,f) = \left(\frac{1}{2\pi} \int_0^{2\pi} |f(re^{it})|^p dt\right)^{1/p}, \quad 0 
$$M_{\infty}(r,f) = \sup_{|z|=r} |f(z)|.$$$$

If 0 < p ≤ ∞ the Hardy space H<sup>p</sup> consists of those f ∈ Hol(D) such that

$$||f||_{H^p} \stackrel{\text{def}}{=} \sup_{0 < r < 1} M_p(r, f) < \infty.$$

We mention [\[18\]](#page-21-0) for the theory of H<sup>p</sup> -spaces.

If 0 < p < ∞ and α > −1, the weighted Bergman space A<sup>p</sup> α consists of those f ∈ Hol(D) such that

$$||f||_{A^p_\alpha} \stackrel{\text{def}}{=} \left( (\alpha + 1) \int_{\mathbb{D}} (1 - |z|)^{\alpha} |f(z)|^p dA(z) \right)^{1/p} < \infty.$$

The unweighted Bergman space A p 0 is simply denoted by A<sup>p</sup> . Here, dA(z) = π dx dy denotes the normalized Lebesgue area measure in D. We refer to [\[19\]](#page-21-1), [\[36\]](#page-21-2) and [\[58\]](#page-22-0) for the theory of these spaces.

<sup>2010</sup> Mathematics Subject Classification. 30H25, 47B38.

Key words and phrases. M¨obius invariant spaces and Besov spaces and Q<sup>s</sup> spaces and multipliers and integration operators and Carleson measures.

This research is supported in part by a grant from "El Ministerio de Ciencia, Innovaci´on y Universidades", Spain (PGC2018-096166-B-I00) and by grants from la Junta de Andaluc´ıa (FQM-210 and UMA18-FEDERJA-002).

We let  $\operatorname{Aut}(\mathbb{D})$  denote the set of all disc automorphisms, that is, of all one-toone analytic maps  $\varphi$  from  $\mathbb{D}$  onto itself. It is well known that  $\operatorname{Aut}(\mathbb{D})$  coincides with the set of all Möbius transformations from  $\mathbb{D}$  onto itself:

$$\operatorname{Aut}(\mathbb{D}) = \{ \lambda \varphi_a : |\lambda| = 1, \ a \in \mathbb{D} \},\,$$

where 
$$\varphi_a(z) = (a-z)/(1-\overline{a}z) \ (z \in \mathbb{D}).$$

A linear space X of analytic functions in  $\mathbb D$  is said to be *conformally invariant* or  $M\ddot{o}bius\ invariant$  if whenever  $f\in X$ , then also  $f\circ\varphi\in X$  for any  $\varphi\in \operatorname{Aut}(\mathbb D)$  and, moreover, X is equipped with a semi-norm  $\rho$  for which there exists a positive constant C such that

$$\rho(f \circ \varphi) \leq C\rho(f)$$
, whenever  $f \in X$  and  $\varphi \in \operatorname{Aut}(\mathbb{D})$ .

The articles [8] and [44] are fundamental references for the theory of Möbius invariant spaces which have attracted much attention in recent years (see, e.g., [3, 16, 17, 30, 47, 57, 58]).

The Bloch space  $\mathcal{B}$  consists of all analytic functions f in  $\mathbb{D}$  such that

$$\rho_{\mathcal{B}}(f) \stackrel{\text{def}}{=} \sup_{z \in \mathbb{D}} (1 - |z|^2) |f'(z)| < \infty.$$

The Schwarz-Pick lemma easily implies that  $\rho_{\mathcal{B}}$  is a conformally invariant seminorm, thus  $\mathcal{B}$  is a conformally invariant space. It is also a Banach space with the norm  $\|\cdot\|_{\mathcal{B}}$  defined by  $\|f\|_{\mathcal{B}} = |f(0)| + \rho_{\mathcal{B}}(f)$ . The little Bloch space  $\mathcal{B}_0$  is the set of those  $f \in \mathcal{B}$  such that  $\lim_{|z| \to 1} (1 - |z|^2) |f'(z)| = 0$ . Alternatively,  $\mathcal{B}_0$  is the closure of the polynomials in the Bloch norm. A classical reference for the theory of Bloch functions is [7]. Rubel and Timoney [44] proved that  $\mathcal{B}$  is the largest "reasonable" Möbius invariant space. More precisely, they proved the following result.

<span id="page-1-0"></span>**Theorem A.** Let X be a Möbius invariant linear space of analytic functions in  $\mathbb{D}$  and let  $\rho$  be a Möbius invariant seminorm on X. If there exists a non-zero decent linear functional L on X which is continuous with respect to  $\rho$ , then  $X \subset \mathcal{B}$  and there exists a constant A > 0 such that  $\rho_{\mathcal{B}}(f) \leq A\rho(f)$ , for all  $f \in X$ .

Here, a linear functional L on X is said to be decent if it extends continuously to  $\mathcal{H}ol(\mathbb{D})$ .

The space BMOA consists of those functions f in  $H^1$  whose boundary values have bounded mean oscillation on the unit circle  $\partial \mathbb{D}$  as defined by F. John and L. Nirenberg. There are many characterizations of BMOA functions. Let us mention the following:

If  $f \in \mathcal{H}ol(\mathbb{D})$ , then  $f \in BMOA$  if and only if  $||f||_{BMOA} \stackrel{def}{=} |f(0)| + \rho_*(f) < \infty$ , where

$$\rho_*(f) = \sup_{a \in \mathbb{D}} \|f \circ \varphi_a - f(a)\|_{H^2}.$$

It is well known that  $H^{\infty} \subset BMOA \subset \mathcal{B}$  and that BMOA equipped with the seminorm  $\rho_*$  is a Möbius invariant space. The space VMOA consists of those  $f \in BMOA$  such that  $\lim_{|a| \to 1} ||f \circ \varphi_a - f(a)||_{H^2} = 0$ , it is the closure of the

polynomials in the BMOA-norm. We mention [28] as a general reference for the space BMOA.

Other important Möbius invariant spaces are the Besov spaces and the  $Q_s$  spaces.

For  $1 , the analytic Besov space <math>B^p$  is defined as the set of all functions f analytic in  $\mathbb D$  such that  $f' \in A^p_{p-2}$ . All  $B^p$  spaces  $(1 are conformally invariant with respect to the semi-norm <math>\rho_{B^p}$  defined by

$$\rho_{{\scriptscriptstyle B}^p}(f) \stackrel{\text{def}}{=} \|f'\|_{A^p_{p-2}}$$

(see [8, p. 112] or [16, p. 46]) and Banach spaces with the norm  $\|\cdot\|_{B^p}$  defined by  $\|f\|_{B^p} = |f(0)| + \rho_{B^p}(f)$ . An important and well-studied case is the classical Dirichlet space  $B^2$  (often denoted by  $\mathcal{D}$ ) of analytic functions whose image has a finite area, counting multiplicities.

The space  $B^1$  requires a special definition: it is the space of all analytic functions f in  $\mathbb{D}$  for which  $f'' \in A^1$ . Although the semi-norm  $\rho$  defined by  $\rho(f) = ||f''||_{A^1}$  is not conformally invariant, the space itself is. An alternative definition of  $B^1$  with a conformally invariant semi-norm is given in [8], where it is also proved that  $B^1$  is contained in any Möbius invariant space. A lot of information on Besov spaces can be found in [8, 16, 17, 37, 56, 58]. Let us recall that

$$VMOA \subsetneq \mathcal{B}_0, \quad BMOA \subsetneq \mathcal{B},$$
  
 $B^1 \subsetneq B^p \subsetneq B^q \subsetneq VMOA \subsetneq BMOA, \quad 1$ 

If  $0 \le s < \infty$ , we say that  $f \in Q_s$  if f is analytic in  $\mathbb{D}$  and

$$\sup_{a\in\mathbb{D}} \int_{\mathbb{D}} |f'(z)|^2 g(z,a)^s \, dA(z) < \infty,$$

where  $g(z, a) = \log(|1 - \overline{a}z|/|a - z|)$  is the Green function of  $\mathbb{D}$ . These spaces were introduced by Aulaskari and Lappan [12] while looking for characterizations of Bloch functions (see [50] for the case s = 2). For s > 1,  $Q_s$  is the Bloch space,  $Q_1 = BMOA$ , and

$$\mathcal{D} \subsetneq Q_{s_1} \subsetneq Q_{s_2} \subsetneq BMOA, \qquad 0 < s_1 < s_2 < 1.$$

It is well known [14, 46] that for every s with  $0 \le s < \infty$ , a function  $f \in \mathcal{H}ol(\mathbb{D})$  belongs to  $Q_s$  if and only if

$$\rho_{Q_s}(f) \stackrel{\text{def}}{=} \left( \sup_{a \in \mathbb{D}} \int_{\mathbb{D}} |f'(z)|^2 (1 - |\varphi_a(z)|^2)^s dA(z) \right)^{1/2} < \infty.$$

All  $Q_s$  spaces  $(0 \le s < \infty)$  are conformally invariant with respect to the seminorm  $\rho_{Q_s}$ . They are also Banach spaces with the norm  $\|\cdot\|_{Q_s}$  defined by  $\|f\|_{Q_s} = |f(0)| + \rho_{Q_s}(f)$ . We mention [52, 53] as excellent references for the theory of  $Q_s$ -spaces.

Let us recall the following two facts which were first observed in [10].

<span id="page-2-0"></span>If 
$$0 , then  $B^p \subset Q_s$  for all  $s > 0$ . (1.1)$$

<span id="page-3-1"></span>If 
$$2 , then  $B^p \subset Q_s$  if and only if  $1 - \frac{2}{p} < s$ . (1.2)$$

For g analytic in D, the Volterra operator T<sup>g</sup> is defined as follows:

$$T_g(f)(z) \stackrel{\text{def}}{=} \int_0^z g'(\xi) f(\xi) d\xi, \ f \in \mathcal{H}ol(\mathbb{D}), \ z \in \mathbb{D}.$$

We define also the companion operator I<sup>g</sup> by

$$I_g(f)(z) \stackrel{\text{def}}{=} \int_0^z g(\xi) f'(\xi) d\xi, \ f \in \mathcal{H}ol(\mathbb{D}), \ z \in \mathbb{D}.$$

The integration operators T<sup>g</sup> and I<sup>g</sup> have been studied in a good number of papers. Let us just mention here that Pommerenke [\[43\]](#page-22-9) proved that T<sup>g</sup> is bounded on H<sup>2</sup> if and only if g ∈ BMOA and that Aleman and Siskakis [\[4\]](#page-20-6) characterized those g ∈ Hol(D) for which T<sup>g</sup> is bounded on H<sup>p</sup> (p ≥ 1), while Aleman and Cima characterized in [\[1\]](#page-20-7) those g ∈ Hol(D) for which T<sup>g</sup> maps H<sup>p</sup> into H<sup>q</sup> . Aleman and Siskakis [\[5\]](#page-20-8) studied the operators I<sup>g</sup> and T<sup>g</sup> acting on Bergman spaces.

For g ∈ Hol(D), the multiplication operator M<sup>g</sup> is defined by

$$M_g(f)(z) \stackrel{\text{def}}{=} g(z)f(z), \quad f \in \mathcal{H}ol(\mathbb{D}), \ z \in \mathbb{D}.$$

If X and Y are two Banach spaces of analytic function in D continuously embedded in Hol(D) and g ∈ Hol(D) then g is said to be a multiplier from X to Y if Mg(X) ⊂ Y . The space of all multipliers from X to Y will be denoted by M(X, Y ) and M(X) will stand for M(X, X). Using the closed graph theorem we see that for the three operators Tg, Ig, Mg, we have that if one of them maps X into Y then it is continuous from X to Y . We remark also that

$$T_g(f) + I_g(f) = M_g(f) - f(0)g(0).$$
 (1.3)

Thus if two of the operators Tg, Ig, M<sup>g</sup> are bounded from X to Y so is the third one.

It is well known that if X is nontrivial then M(X) ⊂ H<sup>∞</sup> (see, e. g., [\[2,](#page-20-9) Lemma 1. 1] or [\[48,](#page-22-10) Lemma 1. 10]), but M(X, Y ) need not be included in H<sup>∞</sup> if Y 6⊂ X. However, when dealing with M¨obius invariant spaces we have the following result.

<span id="page-3-0"></span>Proposition 1.1. *Let* X *and* Y *be two M¨obius invariant spaces of analytic functions in* D *equipped with the seminorms* ρ<sup>X</sup> *and* ρ<sup>Y</sup> *, respectively. Suppose that there exists a non-trivial decent linear functional* L *on* Y *which is continuous with respect to* ρ<sup>Y</sup> *. Let* g ∈ Hol(D)*. Then the following statements hold.*

- (i) *If* M<sup>g</sup> *is continuous from* (X, ρ<sup>X</sup> ) *into* (Y, ρ<sup>Y</sup> )*, then* g ∈ H<sup>∞</sup>*.*
- (ii) *If* I<sup>g</sup> *is continuous from* (X, ρ<sup>X</sup> ) *into* (Y, ρ<sup>Y</sup> )*, then* g ∈ H<sup>∞</sup>*.*

Before embarking into the proof of Proposition [1.1,](#page-3-0) let us mention that, as usual, throughout the paper we shall be using the convention that C = C(p, α, q, β, . . .) will denote a positive constant which depends only upon the displayed parameters p, α, q, β . . . (which sometimes will be omitted) but not necessarily the same at different occurrences. Moreover, for two real-valued functions E1, E<sup>2</sup> we write

 $E_1 \lesssim E_2$ , or  $E_1 \gtrsim E_2$ , if there exists a positive constant C independent of the arguments such that  $E_1 \leq CE_2$ , respectively  $E_1 \geq CE_2$ . If we have  $E_1 \lesssim E_2$  and  $E_1 \gtrsim E_2$  simultaneously then we say that  $E_1$  and  $E_2$  are equivalent and we write  $E_1 \approx E_2$ . Also, if 1 , <math>p' will stand for its conjugate exponent, that is,  $\frac{1}{p} + \frac{1}{p'} = 1$ .

Proof of Proposition 1.1. Since X is conformally invariant,  $\operatorname{Aut}(\mathbb{D}) \subset X$  [8, p. 114] and

<span id="page-4-0"></span>
$$\rho_{\scriptscriptstyle X}(\varphi_a) \simeq 1, \quad a \in \mathbb{D}.$$
(1.4)

Suppose that  $M_g$  is continuous from  $(X, \rho_X)$  into  $(Y, \rho_Y)$ . Using this, Theorem A, and (1.4) we obtain

$$\rho_{\mathcal{B}}(g\,\varphi_a) \lesssim \rho_{\mathcal{Y}}(g\,\varphi_a) \lesssim \rho_{\mathcal{X}}(\varphi_a) \lesssim 1, \quad a \in \mathbb{D}.$$

This implies that

$$(1 - |a|^2) |g'(a)\varphi_a(a) + g(a)\varphi'_a(a)| \lesssim 1, \quad a \in \mathbb{D}.$$

Since  $\varphi(a) = 0$  and  $\varphi'_a(a) = -(1 - |a|^2)^{-1}$ , it follows that

$$|g(a)| \lesssim 1, \quad a \in \mathbb{D},$$

that is,  $g \in H^{\infty}$ .

Similarly, if we assume that  $I_g$  is continuous from  $(X, \rho_X)$  into  $(Y, \rho_Y)$ , we obtain

$$\rho_{\mathcal{B}}\left(I_g(\varphi_a)\right) \lesssim 1, \quad a \in \mathbb{D}.$$

This implies that

$$(1-|a|^2) |(I_g(\varphi_a))'(a)| = (1-|a|^2)|\varphi_a'(a)||g(a)| = |g(a)| \lesssim 1, \quad a \in \mathbb{D}.$$

For notational convenience, set

$$\mathcal{BQ} = \{Q_s : 0 \le s < \infty\} \cup \{B^p : 1 \le p < \infty\}.$$

The main purpose of this paper is characterizing, for a given pair of spaces  $X, Y \in \mathcal{BQ}$ , the functions  $g \in \mathcal{H}ol(\mathbb{D})$  such that the operators  $M_g$ ,  $T_g$  and/or  $I_g$  map X into Y. When X and Y are Besov spaces this question has been extensively studied (see, e.g. [9, 26, 32, 45, 49, 59]). Thus we shall concentrate ourselves to study these operators when acting between a certain Besov space  $B^p$  and a certain  $Q_s$  space and when acting between  $Q_{s_1}$  and  $Q_{s_2}$  for a certain pair of positive numbers  $s_1, s_2$ .

# 2. Multipliers and integration operators from Besov spaces into $Q_s$ -spaces

For  $\alpha > 0$ , the  $\alpha$ -logarithmic Bloch space  $\mathcal{B}_{\log,\alpha}$  is the Banach space of those functions  $f \in \mathcal{H}ol(\mathbb{D})$  which satisfy

$$||f||_{\log,\alpha} \stackrel{\text{def}}{=} |f(0)| + \sup_{z \in \mathbb{D}} (1 - |z|^2) \left( \log \frac{2}{1 - |z|^2} \right)^{\alpha} |f'(z)| < \infty.$$
 (2.1)

For simplicity, the space  $\mathcal{B}_{\log,1}$  will be denoted by  $\mathcal{B}_{\log}$ .

It is clear that  $B_{\log,\alpha} \subset \mathcal{B}_0$ , for all  $\alpha > 0$ . Using the characterization of VMOAin terms of Carleson measures [28, p. 102], it follows easily that

$$B_{\log,\alpha} \subset VMOA$$
, for all  $\alpha > 1/2$ .

In particular,  $\mathcal{B}_{log} \subset VMOA$ .

Brown and Shields [15] showed that  $M(\mathcal{B}) = \mathcal{B}_{\log} \cap H^{\infty}$ . The spaces  $M(B^p, \mathcal{B})$  $(1 \le p < \infty)$  were characterized in [25]. Namely, Theorem 1 of [25] asserts that  $M(B^1, \mathcal{B}) = H^{\infty}$  and

<span id="page-5-0"></span>
$$M(B^p, \mathcal{B}) = H^{\infty} \cap \mathcal{B}_{\log, 1/p'}, \quad 1$$

where p' is the exponent conjugate to p, that is,  $\frac{1}{p} + \frac{1}{p'} = 1$ .

In this section we extend these results. In particular, we shall obtain for any pair (p,s) with  $2 and <math>0 < s < \infty$  a complete characterization of the space of multipliers  $M(B^p, Q_s)$ .

Let us start with the case  $s \ge 1$  which is the simplest one.

**Theorem 2.1.** Let  $g \in \mathcal{H}ol(\mathbb{D})$ . Then:

- (i)  $I_g$  maps  $B^1$  into  $\mathcal{B}$  if and only if  $g \in H^{\infty}$ .
- (ii)  $M_g$  maps  $B^1$  into  $\mathcal{B}$  if and only if  $g \in H^{\infty}$ . (iii)  $T_g$  maps  $B^1$  into  $\mathcal{B}$  if and only if  $g \in \mathcal{B}$ .

*Proof.* If  $I_q(B^1) \subset \mathcal{B}$  then, using Proposition 1.1, it follows that  $g \in H^{\infty}$ .

To prove the converse it suffices to recall that  $B^1 \subset \mathcal{B}$ . Indeed, suppose that  $g \in H^{\infty}$  and take  $f \in B^1$ . Then

$$(1-|z|^2)\left|\left(I_g(f)\right)'(z)\right| = (1-|z|^2)|f'(z)||g(z)| \le ||f||_{\mathcal{B}}||g||_{H^{\infty}}.$$

Thus  $I_q(f) \in \mathcal{B}$ .

Hence (i) is proved. Now, (ii) is contained in [25, Theorem 1].

It remains to prove (iii). If  $T_q(B^1) \subset \mathcal{B}$  then  $T_q(1) = g - g(0) \in \mathcal{B}$  and, hence  $g \in \mathcal{B}$ . Conversely, if  $g \in \mathcal{B}$  and  $f \in \mathcal{B}^1$  then, using the fact that  $B^1 \subset H^{\infty}$ , we obtain

$$(1 - |z|^2) |(T_a(f))'(z)| = (1 - |z|^2) |g'(z)| |f(z)| \le ||g||_{\mathcal{B}} ||f||_{H^{\infty}}.$$

Thus  $T_q(f) \in \mathcal{B}$ . Hence (iii) is also proved.  $\square$ 

<span id="page-5-1"></span>**Theorem 2.2.** Suppose that  $1 , <math>\frac{1}{p} + \frac{1}{p'} = 1$ , and let  $g \in \mathcal{H}ol(\mathbb{D})$ . Then:

- (i)  $I_q$  maps  $B^p$  into  $\mathcal{B}$  if and only if  $g \in H^{\infty}$ .
- (ii)  $M_g$  maps  $B^p$  into  $\mathcal{B}$  if and only if  $g \in H^{\infty} \cap \mathcal{B}_{\log,1/p'}$ .
- (iii)  $T_q$  maps  $B^p$  into  $\mathcal{B}$  if and only if  $g \in \mathcal{B}_{\log 1/p'}$ .

*Proof.* If  $I_q$  maps  $B^p$  into  $\mathcal{B}$  then Proposition 1.1 implies that  $g \in H^{\infty}$ . Conversely, using that  $B^p \subset \mathcal{B}$ , we see that if  $q \in H^{\infty}$  and  $f \in B^p$  then

$$(1 - |z|^2) \left| (I_q(f))'(z) \right| = (1 - |z|^2) |f'(z)| |g(z)| \le ||f||_{\mathcal{B}} ||g||_{H^{\infty}}.$$

Hence,  $I_q(f) \in \mathcal{B}$ . Thus (i) is proved and (ii) reduces to (2.2).

Finally, (iii) follows from the following more precise result.

<span id="page-6-0"></span>**Theorem 2.3.** Suppose that  $1 , <math>\frac{1}{p} + \frac{1}{p'} = 1$ , and let  $g \in \mathcal{H}ol(\mathbb{D})$ . Then the following conditions are equivalent.

- (a)  $T_q$  maps  $B^p$  into  $\mathcal{B}$ .
- (b)  $g \in \mathcal{B}_{\log,1/p'}$ .
- (c)  $T_a$  maps  $B^p$  into  $\mathcal{B}_0$ .

Proof of Theorem 2.3. (a)  $\Rightarrow$  (b) Suppose (a). By the closed graph theorem  $T_g$  is a bounded operator from  $B^p$  into  $\mathcal{B}$ , hence

<span id="page-6-1"></span>
$$(1-|z|^2)|g'(z)f(z)| \lesssim ||f||_{B^p}, \quad z \in \mathbb{D}, \ f \in B^p.$$
 (2.3)

For  $a \in \mathbb{D}$  with  $a \neq 0$ , set

<span id="page-6-3"></span>
$$f_a(z) = \left(\log \frac{1}{1 - |a|^2}\right)^{-1/p} \log \frac{1}{1 - \overline{a}z}, \quad z \in \mathbb{D}.$$
 (2.4)

It is readily seen that  $f_a \in B^p$  for all a and that  $||f_a||_{B^p} \approx 1$ . Using this and taking  $f = f_a$  and z = a in (2.3), we obtain

$$(1-|a|^2)|g'(a)|\left(\frac{1}{1-|a|^2}\right)^{1/p'} \lesssim 1,$$

that is  $g \in \mathcal{B}_{\log,1/p'}$ .

(b)  $\Rightarrow$  (c) Suppose (b) and take  $f \in B^p$ . It is well known that

$$|f(z)| = o\left(\left(\log \frac{1}{1 - |z|^2}\right)^{1/p'}\right), \text{ as } |z| \to 1,$$

(see, e.g., [37, 56]). This and (b) immediately yield that  $T_q(f) \in \mathcal{B}_0$ .

The implication (c)  $\Rightarrow$  (a) is trivial. Hence the proof of Theorem 2.3 is finished and, consequently, Theorem 2.2 is also proved.  $\Box$ 

Let us turn now to the case  $0 < s \le 1$ . We shall consider first the Volterra operators  $T_g$ . For  $0 < s < \infty$  and  $\alpha > 0$  we set

$$Q_{s,\log,\alpha} = \left\{ f \in \mathcal{H}ol(\mathbb{D}) : \sup_{a \in \mathbb{D}} \left( \log \frac{2}{1 - |a|} \right)^{2\alpha} \int_{\mathbb{D}} |f'(z)|^2 (1 - |\varphi_a(z)|^2)^s dA(z) < \infty \right\}.$$

We have the following results.

<span id="page-6-2"></span>**Theorem 2.4.** Suppose that  $0 < s \le 1$  and let  $g \in \mathcal{H}ol(\mathbb{D})$ . Then:

- (i)  $T_g$  maps  $B^1$  into  $Q_s$  if and only if  $g \in Q_s$ .
- (ii) If  $1 , <math>0 < s \le 1$ , and  $T_g$  maps  $B^p$  into  $Q_s$ , then  $g \in Q_{s,\log,1/p'}$ .
- (iii) If  $1 , then <math>T_g$  maps  $B^p$  into  $Q_1 = BMOA$  if and only if  $g \in Q_{1,\log,1/p'}$ .
- (iv) If 2 , <math>0 < s < 1, and  $1 \frac{2}{p} < s$  then  $T_g$  maps  $B^p$  into  $Q_s$  if and only if  $g \in Q_{s,\log,1/p'}$ .

Before we get into the proofs of these results we shall introduce some notation and recall some results which will be needed in our work.

If  $I \subset \partial \mathbb{D}$  is an interval, |I| will denote the length of I. The Carleson square S(I) is defined as  $S(I) = \{re^{it}: e^{it} \in I, 1 - \frac{|I|}{2\pi} \le r < 1\}$ . Also, for  $a \in \mathbb{D}$ , the Carleson box S(a) is defined by

$$S(a) = \left\{ z \in \mathbb{D} : 1 - |z| \le 1 - |a|, \left| \frac{\arg(a\bar{z})}{2\pi} \right| \le \frac{1 - |a|}{2} \right\}.$$

If s>0 and  $\mu$  is a positive Borel measure on  $\mathbb{D}$ , we shall say that  $\mu$  is an s-Carleson measure if there exists a positive constant C such that

$$\mu(S(I)) < C|I|^s$$
, for any interval  $I \subset \partial \mathbb{D}$ ,

or, equivalently, if there exists C > 0 such that

$$\mu(S(a)) \le C(1-|a|)^s$$
, for all  $a \in \mathbb{D}$ .

A 1-Carleson measure will be simply called a Carleson measure.

These concepts were generalized in [55] as follows: If  $\mu$  is a positive Borel measure in  $\mathbb{D}$ ,  $0 \le \alpha < \infty$ , and  $0 < s < \infty$ , we say that  $\mu$  is an  $\alpha$ -logarithmic s-Carleson measure if there exists a positive constant C such that

$$\frac{\mu\left(S(I)\right)\left(\log\frac{2\pi}{|I|}\right)^{\alpha}}{|I|^{s}} \leq C, \quad \text{for any interval } I \subset \partial \mathbb{D}$$

or, equivalently, if

$$\sup_{a \in \mathbb{D}} \frac{\mu\left(S(a)\right) \left(\log \frac{2}{1 - |a|^2}\right)^{\alpha}}{(1 - |a|^2)^s} < \infty.$$

Carleson measures and logarithmic Carleson measures are known to play a basic role in the study of the boundedness of a great number of operators between analytic function spaces. In particular we recall the Carleson embedding theorem for Hardy spaces which asserts that if  $0 and <math>\mu$  is a positive Borel measure on  $\mathbb{D}$  then  $\mu$  is a Carleson measure if and only if the Hardy space  $H^p$  is continuously embedded in  $L^p(d\mu)$  (see [18, Chapter 9]).

In the next theorem we collect a number of known results which will be needed in our work.

(i) If  $0 < s \le 1$  and  $f \in \mathcal{H}ol(\mathbb{D})$ , then  $f \in Q_s$  if and only if Theorem B. the Borel measure  $\mu$  on D defined by

$$d\mu(z) = (1 - |z|^2)^s |f'(z)|^2 dA(z)$$

is an s-Carleson measure.

(ii) If  $0 < \alpha < \infty$ ,  $0 < s < \infty$ , and  $\mu$  is a positive Borel measure on  $\mathbb{D}$  then  $\mu$  is an  $\alpha$ -logarithmic s-Carleson measure if and only if

$$\sup_{a\in\mathbb{D}} \left(\log\frac{2}{1-|a|^2}\right)^{\alpha} \int_{\mathbb{D}} \left(\frac{1-|a|^2}{|1-\overline{a}\,z|^2}\right)^s \, d\mu(z) \, < \, \infty.$$

- (iii) If  $1 then <math>B^p \subset Q_s$  for all s > 0. (iv) If  $2 and <math>1 \frac{2}{p} < s$ , then  $B^p \subset Q_s$ .

(v) For s > -1, we let  $\mathcal{D}_s$  be the space of those functions  $f \in \mathcal{H}ol(\mathbb{D})$  for which

$$||f||_{\mathcal{D}_s} \stackrel{def}{=} |f(0)| + \left( \int_{\mathbb{D}} (1 - |z|^2)^s |f'(z)|^2 dA(z) \right)^{1/2} < \infty.$$

Suppose that 0 < s < 1 and  $\alpha > 1$ , and let  $\mu$  be a positive Borel measure on  $\mathbb{D}$ . If  $\mu$  is an  $\alpha$ -logarithmic s Carleson measure, then  $\mu$  is a Carleson measures for  $\mathcal{D}_s$ , that is,  $\mathcal{D}_s$  is continuously embedded in  $L^2(d\mu)$ .

Let us mention that (i) is due to Aulaskari, Stegenga and Xiao [13], (ii) is due to Zhao [55], (iii) and (iv) were proved by Aulaskari and Csordas in [10], and (v) is due to Pau and Peláez [41, Lemma 1].

Using Theorem B (ii) and the fact that

$$1 - |\varphi(z)|^2 = \frac{(1 - |a|^2)(1 - |z|^2)}{|1 - \overline{a}z|^2},$$

we see that for a function  $f \in \mathcal{H}ol(\mathbb{D})$  we have that  $f \in Q_{s,\log,\alpha}$  if and only if the measure  $\mu$  defined by  $d\mu(z) = (1-|z|^2)^s |f'(z)|^2 dA(z)$  is a  $2\alpha$ -logarithmic s-Carleson measure.

Proof of Theorem 2.4 (i). Suppose that  $T_g$  maps  $B^1$  into  $Q_s$ . Since the constant functions belong to  $B^1$ , we have that  $T_g(1) = g - g(0) \in Q_s$  and, hence,  $g \in Q_s$ . To prove the converse, suppose that  $g \in Q_s$ . Then the measure  $\mu$  defined by

$$d\mu(z) = (1 - |z|^2)^s |g'(z)|^2 dA(z)$$

is an s-Carleson measure. Take now  $f \in B^1$ , then  $f \in H^{\infty}$  and, hence,

$$(1-|z|^2)^s \left| (T_g(f))'(z) \right|^2 = (1-|z|^2)^s |g'(z)|^2 |f(z)|^2 \le ||f||_{H^{\infty}}^2 (1-|z|^2)^s |g'(z)|^2.$$

Since  $\mu$  is an s-Carleson measure, it follows readily that the measure  $\nu$  given by  $d\nu(z) = (1-|z|^2)^s \left| (T_g(f))'(z) \right|^2 dA(z)$  is also an s-Carleson measure and, hence,  $T_g(f) \in Q_s$ .  $\square$ 

Proof of Theorem 2.4 (ii).

Suppose that  $0 < s \le 1$ ,  $1 , and that <math>T_g$  maps  $B^p$  into  $Q_s$ . For  $a \in \mathbb{D} \setminus \{0\}$ , set

$$f_a(z) = \left(\log \frac{1}{1 - |a|^2}\right)^{-1/p} \log \frac{1}{1 - \overline{a}z}, \quad z \in \mathbb{D},$$

as in (2.4). We have that  $||f_a||_{B^p} \approx 1$  and it is also clear that

$$|f_a(z)| \simeq \left(\log \frac{1}{1 - |a|^2}\right)^{1/p'}, \quad z \in S(a).$$

Using these facts, we obtain

$$\frac{\left(\log \frac{1}{1-|a|^2}\right)^{2/p'}}{(1-|a|^2)^s} \int_{S(a)} (1-|z|^2)^s |g'(z)|^2 dA(z)$$

$$\approx \frac{1}{(1-|a|^2)^s} \int_{S(a)} (1-|z|^2)^s |g'(z)f_a(z)|^2 dA(z)$$

$$= \frac{1}{(1-|a|^2)^s} \int_{S(a)} (1-|z|^2)^s |(T_g(f_a))'(z)|^2 dA(z).$$

The fact that  $T_g$  is a bounded operator from  $B^p$  into  $Q_s$ , implies that the measures  $(1-|z|^2)^s |(T_g(f_a))'(z)|^2 dA(z)$  are s-Carleson measures with constants controlled by  $||T_g||^2$ . Then it follows that the measure  $(1-|z|^2)^s |g'(z)|^2 dA(z)$  is a 2/p'-logarithmic s-Carleson measure and, hence,  $g \in Q_{s,\log,1/p'}$ .  $\square$ 

Proof of Theorem 2.4 (iii) and (iv). In view of (ii) we only have to prove that if  $g \in Q_{s,\log,1/p'}$  then  $T_g$  maps  $B^p$  into  $Q_s$ .

Hence, take  $g \in Q_{s,\log,1/p'}$  and set

$$K(g) = \sup_{a \in \mathbb{D}} \left( \log \frac{2}{1 - |a|} \right)^{2/p'} \int_{\mathbb{D}} |g'(z)|^2 (1 - |\varphi_a(z)|^2)^s dA(z),$$

and take  $f \in B^p$ . Set  $F = T_g(f)$ , we have to prove that  $F \in Q_s$  or, equivalently, that the measure  $\mu_F$  defined by

$$d\mu_{F}(z) = (1 - |z|^{2})^{s} |F'(z)|^{2} dA(z)$$

is an s-Carleson measure. Let  $a \in \mathbb{D}$ . Using the well known fact that

$$1 - |a|^2 \approx |1 - \overline{a}z|, \quad z \in S(a),$$

we obtain

$$\frac{1}{(1-|a|^2)^s} \int_{S(a)} |F'(z)|^2 (1-|z|^2)^s dA(z) \approx \int_{S(a)} |F'(z)|^2 \frac{(1-|z|^2)^s (1-|a|^2)^s}{|1-\overline{a}z|^{2s}} dA(z)$$

$$= \int_{S(a)} |f(z)|^2 |g'(z)|^2 (1-|\varphi_a(z)|^2)^s dA(z)$$

$$\leq 2 \int_{\mathbb{D}} |f(a)|^2 |g'(z)|^2 (1-|\varphi_a(z)|^2)^s dA(z)$$

$$+ 2 \int_{\mathbb{D}} |f(z)-f(a)|^2 |g'(z)|^2 (1-|\varphi_a(z)|^2)^s dA(z)$$

$$= 2T_1(a) + 2T_2(a). \tag{2.5}$$

Using the fact that

<span id="page-9-1"></span><span id="page-9-0"></span>
$$|f(a) - f(0)| \lesssim ||f||_{B^p} \left(\log \frac{2}{1 - |a|^2}\right)^{1/p'},$$
 (2.6)

we obtain

$$T_1(a) \lesssim \|f\|_{B^p}^2 \left(\log \frac{2}{1 - |a|^2}\right)^{2/p'} \int_{\mathbb{D}} |g'(z)|^2 (1 - |\varphi_a(z)|^2)^s dA(z) \lesssim K(g) \|f\|_{B^p}^2.$$
(2.7)

To estimate  $T_2(a)$  we shall treat separately the cases s = 1 and 0 < s < 1. Let us start with the case s = 1. Then

<span id="page-10-0"></span>
$$T_2(a) = \int_{\mathbb{D}} |f(z) - f(a)|^2 |g'(z)|^2 (1 - |\varphi_a(z)|^2) dA(z).$$

Making the change of variable  $w = \varphi(z)$  in the last integral, we obtain

$$T_2(a) = \int_{\mathbb{D}} |(f \circ \varphi_a)(w) - f(a)|^2 |(g \circ \varphi_a)'(w)|^2 (1 - |w|^2) dA(w).$$

Since  $Q_{1,\log,1/p'} \subset Q_1 = BMOA$ ,  $g \in BMOA$  and then it follows that, for all  $a \in \mathbb{D}$ ,  $g \circ \varphi_a \in BMOA$  and  $\rho_*(g \circ \varphi_a) = \rho_*(g)$ . This gives that all the measures  $(1 - |w|^2)|(g \circ \varphi_a)'(w)|^2 dA(w)$   $(a \in \mathbb{D})$  are Carleson measures with constants controlled by  $||g||^2_{BMOA}$ . Then, using the Carleson embedding theorem for  $H^2$  and the fact that  $B^p$  is continuously embedded in BMOA, it follows that

$$T_2(a) \lesssim \|g\|_{BMOA}^2 \|f \circ \varphi_a - f(a)\|_{H^2}^2 \lesssim \|g\|_{BMOA}^2 \|f\|_{BMOA}^2 \lesssim \|g\|_{BMOA}^2 \|f\|_{B^p}^2.$$

Putting together this, (2.5), and (2.7), we see that the measure  $\mu_F$  is a Carleson measure. This finishes the proof of part (iii).

To finish the proof of part (iv) we proceed to estimate  $T_2(a)$  assuming that 2 , <math>0 < s < 1, and  $1 - \frac{2}{p} < s$ . Notice that

$$T_2(a) = (1 - |a|^2)^s \int_{\mathbb{D}} \left| \frac{f(z) - f(a)}{(1 - \overline{a}z)^s} \right|^2 |g'(z)|^2 (1 - |z|^2)^s dA(z).$$

Since 0 < s < 1, 2/p' > 1, and the measure  $(1 - |z|^2)^s |g'(z)|^2 dA(z)$  is a 2/p'-logarithmic s-Carleson measure, using Theorem B (v), it follows that

$$T_2(a) \lesssim (1-|a|^2)^s \left(|f(a)-f(0)|^2 + \int_{\mathbb{D}} \left| \left(\frac{f(z)-f(a)}{(1-\overline{a}z)^s}\right)' \right|^2 (1-|z|^2)^s dA(z) \right).$$

The growth estimate (2.6) and simple computations yield

$$T_{2}(a) \lesssim \|f\|_{B^{p}}^{2} (1 - |a|^{2})^{s} \left(\log \frac{2}{1 - |a|^{2}}\right)^{2/p'} + \int_{\mathbb{D}} |f'(z)|^{2} (1 - |\varphi_{a}(z)|^{2})^{s} dA(z)$$

$$+ \int_{\mathbb{D}} \frac{|f(z) - f(a)|^{2}}{|1 - \overline{a} z|^{2}} (1 - |\varphi_{a}(z)|^{2})^{s} dA(z)$$

$$\lesssim \|f\|_{B^{p}}^{2} + \int_{\mathbb{D}} |f'(z)|^{2} (1 - |\varphi_{a}(z)|^{2})^{s} dA(z) + \int_{\mathbb{D}} \frac{|f(z) - f(a)|^{2}}{|1 - \overline{a} z|^{2}} (1 - |\varphi_{a}(z)|^{2})^{s} dA(z).$$

By Theorem B (iv), our assumptions on s and p imply that  $B^p$  is continuously embedded in  $Q_s$ . Hence,  $f \in Q_s$ . This implies that

$$\int_{\mathbb{D}} |f'(z)|^2 (1 - |\varphi_a(z)|^2)^s dA(z) \le ||f||_{Q_s}^2 \lesssim ||f||_{B^p}^2$$

and that

$$\int_{\mathbb{D}} \frac{|f(z) - f(a)|^2}{|1 - \overline{a}|^2} (1 - |\varphi_a(z)|^2)^s dA(z) \lesssim ||f||_{Q_s}^2 \lesssim ||f||_{B^p}^2,$$

by a result proved by Pau and Peláez in [41, pp. 551–552]. Consequently, we have proved that  $T_2(a) \lesssim ||f||_{B^p}^2$ . This, together with (2.5) and (2.7), shows that  $\mu_F$  is an s-Carleson measure as desired. Thus the proof is also finished in this case.

The case when 1 and <math>0 < s < 1 remains open. This is so because if we set  $\alpha = 2/p'$ , then  $\alpha \le 1$  and, hence,  $\alpha$  is not in the conditions of Theorem B (v). We can prove the following result.

<span id="page-11-2"></span>**Theorem 2.5.** Suppose that 1 and <math>0 < s < 1, and let  $g \in \mathcal{H}ol(\mathbb{D})$ . The following statements hold.

- (i) If  $T_g$  maps  $B^p$  into  $Q_s$  then  $g \in Q_{s,\log,1/p'}$ .
- (ii) If  $\alpha > 1/2$  and  $g \in Q_{s,\log,\alpha}$  then  $T_g$  maps  $B^p$  into  $Q_s$ .

*Proof.* (i) follows from part (ii) of Theorem 2.4.

Let us turn to prove (ii). Suppose that 0 < s < 1,  $\alpha > 1/2$ , and  $g \in Q_{s,\log,\alpha}$ . Set

$$K(g) = \sup_{a \in \mathbb{D}} \left( \log \frac{2}{1 - |a|} \right)^{2\alpha} \int_{\mathbb{D}} |g'(z)|^2 (1 - |\varphi_a(z)|^2)^s dA(z),$$

and take  $f \in B^p$ . Set  $F = T_g(f)$ , we have to prove the  $F \in Q_s$  or, equivalently, that the measure  $\mu_F$  defined by

$$d\mu_F(z) = (1 - |z|^2)^s |F'(z)|^2 dA(z)$$

is an s-Carleson measure. Now we argue as in the proof of Theorem 2.4 (iv). For  $a \in \mathbb{D}$ , we obtain

<span id="page-11-1"></span>
$$\frac{1}{(1-|a|^2)^s} \int_{S(a)} |F'(z)|^2 (1-|z|^2)^s dA(z) \lesssim 2T_1(a) + 2T_2(a), \qquad (2.8)$$

where  $T_1(a)$  and  $T_2(a)$  are defined as in the proof of Theorem 2.4. Using (2.6) and the fact that  $\frac{1}{p'} \leq \frac{1}{2} < \alpha$ , we obtain

<span id="page-11-0"></span>
$$|f(a) - f(0)| \lesssim ||f||_{B^p} \left(\log \frac{2}{1 - |a|^2}\right)^{\alpha}.$$

This yields

$$T_1(a) \lesssim \|f\|_{B^p}^2 \left(\log \frac{2}{1 - |a|^2}\right)^{2\alpha} \int_{\mathbb{D}} |g'(z)|^2 (1 - |\varphi_a(z)|^2)^s dA(z) \lesssim K(g) \|f\|_{B^p}^2.$$
(2.9)

To estimate  $T_2(a)$ , observe that the measure  $(1-|z|^2)^s|g'(z)|^2\,dA(z)$  is a  $2\alpha$ -logarithmic s-Carleson measure. Since  $2\alpha>1$ , using Lemma 1 of [41], this implies that the measure  $(1-|z|^2)^s|g'(z)|^2\,dA(z)$  is a Carleson measure for  $\mathcal{D}_s$ . Then, arguing as in the proof of Theorem 2.4 (iv), we obtain  $T_2(a)\lesssim \|f\|_{B^p}^2$ . This, together with (2.9) and (2.8), implies that the measure  $\mu_F$  is an s-Carleson measure.

Regarding the operators  $I_g$  and  $M_g$  we have the following results.

<span id="page-12-0"></span>**Theorem 2.6.** Let  $g \in \mathcal{H}ol(\mathbb{D})$ , then:

- (1) If  $1 and <math>0 < s \le 1$  then:
  - (1a)  $I_a$  maps  $B^p$  into  $Q_s$  if and only if  $g \in H^{\infty}$ .
  - (1b) If  $M_g$  maps  $B^p$  into  $Q_s$  then  $g \in Q_{s,\log 1/p'} \cap H^{\infty}$ .
- (1c) If  $g \in Q_{s,\log,\alpha} \cap H^{\infty}$  for some  $\alpha > 1/2$  then  $M_g$  maps  $B^p$  into  $Q_s$ . (2) If  $2 and <math>1 \frac{2}{p} < s \le 1$  then:
- - (2a)  $I_a$  maps  $B^p$  into  $Q_s$  if and only if  $g \in H^{\infty}$ .
  - (2b)  $M_g$  maps  $B^p$  into  $Q_s$  if and only if  $g \in Q_{s,\log,1/p'} \cap H^{\infty}$ .
- (3) If  $2 and <math>0 < s \le 1 \frac{2}{n}$  then:
  - (3a)  $I_q$  maps  $B^p$  into  $Q_s$  if and only if  $g \equiv 0$ .
  - (3b)  $M_a$  maps  $B^p$  into  $Q_s$  if and only if  $q \equiv 0$ .

Proof of Parts (1) and (2) of Theorem 2.6. Using Proposition 1.1 it follows that if either  $I_g$  or  $M_g$  maps  $B^p$  into  $Q_s$  for any pair (s,p) with  $0 < s < \infty$  and  $1 then <math>q \in H^{\infty}$ .

Suppose now that s and p are in the conditions of (1) or (2) and that  $g \in H^{\infty}$ . Take  $f \in B^p$ . We have to prove  $I_g(f) \in Q_s$  or, equivalently, that the measure

<span id="page-12-1"></span>
$$(1 - |z|^2)^s |f'(z)|^2 |g(z)|^2 dA(z)$$
 is an s-Carleson measure. (2.10)

Using (1.1) and (1.2), we see that  $B^p \subset Q_s$ . Hence  $f \in Q_s$  which is the same as saying that  $(1-|z|^2)^s |f'(z)|^2 dA(z)$  is an s-Carleson measure. This and the fact that  $g \in H^{\infty}$  trivially yield (2.10). Thus (1a) and (2a) are proved. Then (1b), (1c), and (2b) follow using Proposition 1.1, the fact that if two of the operators  $T_g$ ,  $I_g$ ,  $M_g$  map  $B^p$  into  $Q_s$  so does the third one, Theorem 2.4, and Theorem 2.5.

In order to prove Theorem 2.6 (3), for 2 we shall consider the function $F_p$  defined by

$$F_p(z) = \sum_{k=1}^{\infty} \frac{1}{k^{1/2} 2^{k/p}} z^{2^k}, \quad z \in \mathbb{D}.$$
 (2.11)

Using [10, Corollary 7] or [14, Theorem 6], we see that  $F_p \in B^p$  and  $F_p \notin Q_{1-\frac{2}{n}}$ . Hence

$$F_p \in B^p \setminus Q_s, \quad 0 < s \le 1 - \frac{2}{p}, \quad 2 < p < \infty.$$
 (2.12)

Let us estimate the integral means  $M_2(r, F'_n)$ . We have

$$zF'_p(z) = \sum_{k=1}^{\infty} 2^{k/p'} k^{-1/2} z^{2^k}, \quad z \in \mathbb{D}$$

and, hence,

$$M_2(r, F_p')^2 \gtrsim \sum_{k=1}^{\infty} 2^{2k/p'} k^{-1} r^{2^{k+1}}, \quad 0 < r < 1.$$

Set 
$$r_n = 1 - 2^{-n}$$
  $(n = 1, 2, ...)$ . Then

$$M_2(r_n, F_p')^2 \gtrsim \sum_{k=1}^{\infty} 2^{2k/p'} k^{-1} r_n^{2^{k+1}}$$
  
  $\gtrsim 2^{2n/p'} n^{-1} r_n^{2^{n+1}} \gtrsim 2^{2n/p'} n^{-1} \approx \frac{1}{(1-r_n)^{2/p'} \log \frac{2}{1-r_n}}, \quad n = 1, 2, \dots$ 

This readily yields

<span id="page-13-2"></span>
$$M_2(r, F_p')^2 \gtrsim \frac{1}{(1-r)^{2/p'} \log \frac{2}{1-r}}, \quad 0 < r < 1.$$
 (2.13)

Proof of part (3) of Theorem 2.6. Suppose that  $2 and <math>0 < s \le 1 - \frac{2}{p}$  and  $g \in \mathcal{H}ol(\mathbb{D})$  is not identically zero.

Suppose first that either  $I_g$  or  $M_g$  maps  $B^p$  into  $Q_s$ . We know that then  $g \in H^\infty$  and then, by Fatou's theorem and the Riesz uniqueness theorem, we know that g has a finite non-tangential limit  $g(e^{i\theta})$  for almost every  $\theta \in [0, 2\pi]$  and that  $g(e^{i\theta}) \neq 0$  for almost every  $\theta$ . Then it follows that there exist C > 0,  $r_0 \in (0, 1)$ , and a measurable set  $E \subset [0, 2\pi]$  whose Lebesgue measure |E| is positive such that

<span id="page-13-0"></span>
$$|g(re^{i\theta})| \ge C, \quad \theta \in E, \quad r_0 < r < 1.$$
 (2.14)

Since  $F_p$  is given by a power series with Hadamard gaps, Lemma 6. 5 in [60, Vol. 1, p. 203] implies that

<span id="page-13-1"></span>
$$\int_{E} |F_{p}'(re^{i\theta})|^{2} d\theta \approx M_{2}(r, F_{p}')^{2}, \quad 0 < r < 1.$$
(2.15)

Using the fact that  $s \leq 1 - \frac{2}{n}$ , (2.14), (2.15), and (2.13), we obtain

$$\int_{0}^{1} (1-r)^{s} M_{2}(r, F'_{p}g)^{2} dr \geq \int_{r_{0}}^{1} (1-r)^{1-\frac{2}{p}} M_{2}(r, F'_{p}g)^{2} dr$$

$$\gtrsim \int_{r_{0}}^{1} (1-r)^{1-\frac{2}{p}} \int_{E} |F'_{p}(re^{i\theta})|^{2} |g(re^{i\theta})|^{2} d\theta dr \gtrsim \int_{r_{0}}^{1} (1-r)^{1-\frac{2}{p}} \int_{E} |F'_{p}(re^{i\theta})|^{2} d\theta dr$$

$$\gtrsim \int_{r_{0}}^{1} (1-r)^{1-\frac{2}{p}} M_{2}(r, F'_{p})^{2} dr \gtrsim \int_{r_{0}}^{1} \frac{dr}{(1-r)\log\frac{2}{1-r}} = \infty. \tag{2.16}$$

If we assume that  $I_g$  maps  $B^p$  into  $Q_s$  then  $I_g(F_p) \in Q_s$  and then, using [11, Proposition 3. 1], it follows that

<span id="page-13-3"></span>
$$\int_0^1 (1-r)^s M_2(r, F_p' g)^2 dr < \infty.$$

This is in contradiction with (2.16).

Assume now that  $M_g$  maps  $B^p$  into  $Q_s$ . Since 1 and  $F_p$  belong to  $B^p$ , we have that g and  $F_p g$  belong to  $Q_s$  and then, by [11, Proposition 3.1],

<span id="page-13-4"></span>
$$\int_0^1 (1-r)^s M_2(r,g')^2 dr < \infty \tag{2.17}$$

and

$$\int_0^1 (1-r)^s M_2(r, (F_p g)')^2 dr < \infty.$$
 (2.18)

Notice that  $F_p \in H^{\infty}$  and then

$$M_2(r, F_p g') \lesssim M_2(r, g'), \quad 0 < r < 1.$$

This and (2.17) imply that

<span id="page-14-0"></span>
$$\int_0^1 (1-r)^s M_2(r, F_p' g)^2 dr < \infty.$$
 (2.19)

We have arrived to a contradiction because it is clear that (2.16) and (2.19) cannot be simultaneously true.  $\square$ 

In the other direction we have the following result.

**Theorem 2.7.** Suppose that  $0 < s < \infty$  and  $1 \le p < \infty$  and let  $g \in \mathcal{H}ol(\mathbb{D})$ . Then the following conditions are equivalent

- (i)  $M_g$  maps  $Q_s$  into  $B^p$ .
- (ii)  $g \equiv 0$ .

*Proof.* Suppose that  $g \not\equiv 0$ . Choose an increasing sequence  $\{r_n\}_{n=1}^{\infty} \subset (0,1)$  with  $\lim \{r_n\} = 1$  and a sequence  $\{\theta_n\}_{n=1}^{\infty} \subset [0,2\pi]$  such that

$$|g(r_n e^{i\theta_n})| = M_{\infty}(r_n, g), \quad n = 1, 2, \dots$$

For each n set

$$f_n(z) = \log \frac{1}{1 - e^{-i\theta_n} z}, \quad z \in \mathbb{D}.$$

Notice that  $M(r_1, g) > 0$  and that the sequence  $\{M(r_n, g)\}$  is increasing. Set

$$f_n(z) = \log \frac{1}{1 - e^{-i\theta_n z}}, \quad z \in \mathbb{D}, \quad n = 1, 2, \dots$$

We have that  $f_n \in Q_s$  for all n and

$$||f_n||_{Q_s} \asymp 1.$$

Assume that  $M_g$  maps  $Q_s$  into  $B^p$ . Then, by the closed graph theorem,  $M_g$  is bounded operator from  $Q_s$  into  $B^p$ . Hence the sequence  $\{g f_n\}_{n=1}^{\infty}$  is a bounded sequence on  $B^p$ , that is,

$$||g f_n||_{B^p} \lesssim 1.$$

Then it follows that

$$M(r_1, g) \log \frac{1}{1 - r_n} \le M(r_n, g) \log \frac{1}{1 - r_n} = |g(r_n e^{i\theta_n}) f_n(r_n e^{i\theta_n})|$$

$$\lesssim ||g f_n||_{B^p} \left(\log \frac{1}{1 - r_n}\right)^{1/p'} \lesssim \left(\log \frac{1}{1 - r_n}\right)^{1/p'}.$$

This is a contradiction.  $\square$ 

## 3. Multipliers and integration operators between $Q_s$ spaces

As we mentioned above the space of multipliers  $M(\mathcal{B}) = M(Q_s)$  (s > 1) was characterized by Brown and Shields in [15]. Ortega and Fàbrega [40] characterized the space  $M(BMOA) = M(Q_1)$ . Pau and Peláez [41] and Xiao [54] characterized the spaces  $M(Q_s)$  (0 < s < 1) closing a conjecture formulated in [51]. Indeed, Theorem 1 of [41] and Theorem 1.2 of [54] assert the following.

<span id="page-15-3"></span>**Theorem C.** Suppose that  $0 < s \le 1$  and let g be an analytic function in the unit disc  $\mathbb{D}$ . Then:

- (i)  $T_g$  maps  $Q_s$  into itself if and only if  $g \in Q_{s,\log,1}$ .
- (ii)  $I_g$  maps  $Q_s$  into itself if and only if  $g \in H^{\infty}$ .
- (ii)  $M_g$  maps  $Q_s$  into itself if and only if  $g \in Q_{s,\log,1} \cap H^{\infty}$ .

We shall prove the following results.

<span id="page-15-0"></span>**Theorem 3.1.** Suppose that  $0 < s_1 \le s_2 \le 1$  and let  $g \in \mathcal{H}ol(\mathbb{D})$ . Then:

- (i)  $T_g$  maps  $Q_{s_1}$  into  $Q_{s_2}$  if and only if  $g \in Q_{s_2,\log,1}$ .
- (ii)  $I_g$  maps  $Q_{s_1}$  into  $Q_{s_2}$  if and only if  $g \in H^{\infty}$ .
- (iii)  $M_g$  maps  $Q_{s_1}$  into  $Q_{s_2}$  if and only if  $g \in Q_{s_2,\log,1} \cap H^{\infty}$ .

<span id="page-15-4"></span>**Theorem 3.2.** Suppose that  $0 < s_1 < s_2 \le 1$  and let  $g \in \mathcal{H}ol(\mathbb{D})$ . Then the following conditions are equivalent:

- (i)  $I_g maps Q_{s_2} into Q_{s_1}$ .
- (ii)  $M_g$  maps  $Q_{s_2}$  into  $Q_{s_1}$ .
- (iii)  $q \equiv 0$ .

Proof of Theorem 3.1. For  $a \in \mathbb{D}$  we set

$$h_a(z) = \log \frac{2}{1 - \overline{a}z}, \quad z \in \mathbb{D}.$$

Then  $h_a \in Q_{s_1}$  for all  $a \in \mathbb{D}$  and

<span id="page-15-1"></span>
$$||h_a||_{Q_{s_1}} \asymp 1. \tag{3.1}$$

• If  $T_g$  maps  $Q_{s_1}$  into  $Q_{s_2}$  then  $T_g$  is a bounded operator from  $Q_{s_1}$  into  $Q_{s_2}$ . Using this and (3.1), it follows that for all  $a \in \mathbb{D}$  the measure  $(1 - |z|^2)^{s_2}|g'(z)|^2|h_a(z)|^2dA(z)$  is an  $s_2$ -Carleson measure and that

<span id="page-15-2"></span>
$$\int_{S(a)} (1 - |z|^2)^{s_2} |g'(z)|^2 |h_a(z)|^2 dA(z) \lesssim (1 - |a|^2)^{s_2}, \quad a \in \mathbb{D}.$$
 (3.2)

Since

$$|h_a(z)| \simeq \log \frac{2}{1 - |a|^2}, \quad z \in S(a),$$

(3.2) implies that

$$\left(\log \frac{2}{1-|a|^2}\right)^2 \int_{S(a)} (1-|z|^2)^{s_2} |g'(z)|^2 dA(z) \lesssim (1-|a|^2)^{s_2}.$$

This is the same as saying that the measure  $(1 - |z|^2)^{s_2} |g'(z)|^2 dA(z)$  is a 2-logarithmic  $s_2$ -Carleson measure or, equivalently, that  $g \in Q_{s_2,\log,1}$ .

If  $g \in Q_{s_2,\log,1}$  then, by Theorem C,  $T_g$  maps  $Q_{s_2}$  into itself. Since  $Q_{s_1} \subset Q_{s_2}$ , it follows trivially that  $T_g$  maps  $Q_{s_1}$  into  $Q_{s_2}$ . Hence (i) is proved

• Proposition 1.1 shows that if  $I_g$  maps  $Q_{s_1}$  into  $Q_{s_2}$  then  $g \in H^{\infty}$ .

Conversely, suppose that  $g \in H^{\infty}$ . In order to prove that  $I_g$  maps  $Q_{s_1}$  into  $Q_{s_2}$ , we have to prove that for any  $f \in Q_{s_1}$  the measure  $(1-|z|^2)^{s_2}|g(z)|^2|f'(z)|^2 dA(z)$  is an  $s_2$ -Carleson measure. So, take  $f \in Q_{s_1}$ . Then  $(1-|z|^2)^{s_1}|f'(z)|^2 dA(z)$  is an  $s_1$ -Carleson measure. Then it follows that

$$\int_{S(a)} (1 - |z|^2)^{s_2} |g(z)|^2 |f'(z)|^2 dA(z)$$

$$\leq ||g||_{H^{\infty}}^2 (1 - |a|^2)^{s_2 - s_1} \int_{S(a)} (1 - |z|^2)^{s_1} |f'(z)|^2 dA(z)$$

$$\lesssim (1 - |a|^2)^{s_2}.$$

This shows that  $(1-|z|^2)^{s_2}|g(z)|^2|f'(z)|^2 dA(z)$  is an  $s_2$ -Carleson measure as desired, finishing the proof of (ii).

• If  $M_g$  maps  $Q_{s_1}$  into  $Q_{s_2}$  then, Proposition 1.1,  $g \in H^{\infty}$ . Then (i) implies that  $I_g$  maps  $Q_{s_1}$  into  $Q_{s_2}$ . Since  $M_g(f) = I_g(f) + T_g(f) + f(0)g(0)$ , it follows that  $T_g$  maps  $Q_{s_1}$  into  $Q_{s_2}$ . Then (i) yields  $g \in Q_{s_2,\log,1}$ . Then we have that  $g \in Q_{s_2,\log,1} \cap H^{\infty}$ .

Conversely, if  $g \in Q_{s_2,\log,1} \cap H^{\infty}$  then (i) and (ii) immediately give that both  $T_g$  and  $I_g$  map  $Q_{s_1}$  into  $Q_{s_2}$  and then so does  $M_g$ .  $\square$ 

Some results from [11] will be used to prove Theorem 3.2. As we have already noticed if  $0 < s \le 1$  and  $f \in Q_s$  then  $\int_0^1 (1-r)^s M_2(r,f')^2 dr < \infty$ . Using ideas from [27], Aulaskari, Girela and Wulan [11, Theorem 3.1] proved that this result is sharp in a very strong sense.

<span id="page-16-0"></span>**Theorem D.** Suppose that  $0 < s \le 1$  and let  $\varphi$  be a positive increasing function defined in (0,1) such that

$$\int_0^1 (1-r)^s \, \varphi(r)^2 \, dr < \infty.$$

Then there exists a function  $f \in Q_s$  given by a power series with Hadamard gaps such that  $M_2(r, f') \ge \varphi(r)$  for all  $r \in (0, 1)$ .

Proof of Theorem 3.2. Suppose that  $g \not\equiv 0$  and that either  $I_g$  or  $M_g$  maps  $Q_{s_2}$  into  $Q_{s_1}$ . By Proposition 1.1,  $g \in H^{\infty}$  and then it follows that there exist C > 0,  $r_0 \in (0,1)$ , and a measurable set  $E \subset [0,2\pi]$  whose Lebesgue measure |E| is positive such that

$$|g(re^{i\theta})| \ge C, \quad \theta \in E, \quad r_0 < r < 1.$$

• Suppose that  $I_g$  maps  $Q_{s_2}$  into  $Q_{s_1}$ . Then we use Theorem D to pick a function  $F \in Q_{s_2}$  given by a power series with Hadamard gaps so that

<span id="page-16-1"></span>
$$M_2(r, F') \ge \frac{1}{(1-r)^{(1+s_1)/2}}, \quad 0 < r < 1.$$
 (3.3)

Since Ig(F) ∈ Q<sup>s</sup><sup>1</sup> ,

<span id="page-17-0"></span>
$$\int_0^1 (1-r)^{s_1} M_2(r, F'g)^2 dr < \infty.$$
 (3.4)

However, using Lemma 6. 5 in [\[60,](#page-22-16) Vol. 1, p. 203] and [\(3.3\)](#page-16-1), it follows that

$$\int_{0}^{1} (1-r)^{s_{1}} M_{2}(r, F'g)^{2} dr \gtrsim \int_{r_{0}}^{1} (1-r)^{s_{1}} \int_{E} |F'(re^{i\theta})|^{2} |g(re^{i\theta})|^{2} d\theta dr$$

$$\gtrsim \int_{r_{0}}^{1} (1-r)^{s_{1}} \int_{E} |F'(re^{i\theta})|^{2} d\theta dr$$

$$\asymp \int_{r_{0}}^{1} (1-r)^{s_{1}} M_{2}(r, F')^{2} dr$$

$$\gtrsim \int_{r_{0}}^{1} (1-r)^{-1} dr$$

$$= \infty.$$

This is in contradiction with [\(3.4\)](#page-17-0).

• Suppose now that M<sup>g</sup> maps Q<sup>s</sup><sup>2</sup> into Q<sup>s</sup><sup>1</sup> . Take ε > 0 with s<sup>2</sup> − s<sup>1</sup> − ε > 0 and use Theorem [D](#page-16-0) to pick a function H ∈ Q<sup>s</sup><sup>2</sup> given by a power series with Hadamard gaps so that

<span id="page-17-1"></span>
$$M_2(r, H') \ge \frac{1}{(1-r)^{(1+s_1+\varepsilon)/2}}, \quad 0 < r < 1.$$
 (3.5)

Since gH ∈ Q<sup>s</sup><sup>1</sup> we have that

<span id="page-17-2"></span>
$$\int_0^1 (1-r)^{s_1} M_2(r, g'H + gH')^2 dr < \infty.$$
 (3.6)

Using Lemma 6. 5 in [\[60,](#page-22-16) Vol. 1, p. 203] and [\(3.5\)](#page-17-1), we obtain as above that

$$\int_{0}^{1} (1-r)^{s_{1}+\varepsilon} M_{2}(r, H'g)^{2} dr \gtrsim \int_{r_{0}}^{1} (1-r)^{s_{1}+\varepsilon} \int_{E} |H'(re^{i\theta})|^{2} d\theta dr$$

$$\gtrsim \int_{r_{0}}^{1} (1-r)^{s_{1}+\varepsilon} M_{2}(r, H')^{2} dr$$

$$\gtrsim \int_{r_{0}}^{1} \frac{dr}{1-r}$$

$$= \infty. \tag{3.7}$$

Notice that g ∈ Q<sup>s</sup><sup>1</sup> . Using this and the fact that

<span id="page-17-4"></span><span id="page-17-3"></span>
$$|H(z)| \lesssim \log \frac{2}{1-|z|}, \quad z \in \mathbb{D},$$

it follows that

$$\int_{0}^{1} (1-r)^{s_{1}+\varepsilon} M_{2}(r, Hg')^{2} dr \lesssim \int_{0}^{1} (1-r)^{s_{1}+\varepsilon} \left(\log \frac{2}{1-r}\right)^{2} M_{2}(r, g')^{2} dr$$

$$\lesssim \int_{0}^{1} (1-r)^{s_{1}+\frac{\varepsilon}{2}} M_{2}(r, g') dr < \infty.$$
(3.8)

We have arrived to a contradiction because (3.6), (3.7), and (3.8) cannot hold simultaneously.  $\square$ 

Remark 3.3. The implication (ii)  $\Rightarrow$  (iii) in Theorem 3.2 was obtained by Pau and Peláez [42, Corollary 4] using the fact that there exists a function  $f \in Q_{s_2}$ ,  $f \not\equiv 0$ , whose sequence of zeros is not a  $Q_{s_1}$ -zero set.

This idea gives also the following:

$$M(\mathcal{B}, Q_s) = \{0\}, \quad 0 < s \le 1.$$

Indeed, it is well known that there exists a function  $f \in \mathcal{B}$ ,  $f \not\equiv 0$ , whose sequence of zeros does not satisfy the Blaschke condition [7, 31]. If  $g \not\equiv 0$  were a multiplier from  $\mathcal{B}$  into  $Q_s$  for some  $s \leq 1$  then the sequence of zeros of fg would satisfy the Blaschke condition. But this is not true because all the zeros of f are zeros of f.

### 4. Some further results

The inner-outer factorization of functions in the Hardy spaces plays an outstanding role in lots of questions. In many cases the outer factor  $O_f$  of f inherits properties of f. Working in this setting the following concepts arise as natural and quite interesting.

A subspace X of  $H^1$  is said to have the f-property (also called the property of division by inner functions) if  $h/I \in X$  whenever  $h \in X$  and I is an inner function with  $h/I \in H^1$ .

Given  $v \in L^{\infty}(\partial \mathbb{D})$ , the Toeplitz operator  $T_v$  associated with the symbol v is defined by

$$T_v f(z) = P(vf)(z) = \frac{1}{2\pi i} \int_{\partial \mathbb{D}} \frac{v(\xi)f(\xi)}{\xi - z} d\xi, \quad f \in H^1, \quad z \in \mathbb{D}.$$

Here, P is the Szegö projection.

A subspace X of  $H^1$  is said to have the K-property if  $T_{\overline{\psi}}(X) \subset X$  for any  $\psi \in H^{\infty}$ .

These notions were introduced by Havin in [34]. It was also pointed out in [34] that the K-property implies the f-property: indeed, if  $h \in H^1$ , I is inner and  $h/I \in H^1$  then  $h/I = T_{\overline{I}}h$ .

In addition to the Hardy spaces  $H^p$   $(1 many other spaces such as the Dirichlet space [34, 38], several spaces of Dirichlet type including all the Besov spaces <math>B^p$  (1 [20, 21, 22, 39], the spaces <math>BMOA and VMOA [35], and the  $Q_s$  spaces (0 < s < 1) [23] have the K-property. The Hardy space  $H^1$ ,  $H^\infty$  and  $VMOA \cap H^\infty$  are examples of spaces which have the f-property bur fail to have the K-property [35].

The first example of a subspace of  $H^1$  not possessing the f-property is due to Gurarii [33] who proved that the space of analytic functions in  $\mathbb{D}$  whose sequence of Taylor coefficients is in  $\ell^1$  does not have the f-property. Anderson [6] proved

that the space  $\mathcal{B}_0 \cap H^{\infty}$  does not have the f-property. Later on it was proved in [29] that if  $1 \leq p < \infty$  then  $H^p \cap \mathcal{B}$  fails to have the f-property also.

Since as we have already mentioned the Besov spaces  $B^p$  ( $1 ) and the <math>Q_s$  spaces ( $0 < s \le 1$ ) have the K-property (and, also, the f-property), it seems natural to investigate whether the spaces of multipliers and the spaces  $Q_{s,\log,\alpha}$  that we have considered in our work have also these properties. We shall prove the following results.

<span id="page-19-0"></span>**Theorem 4.1.** The spaces of multipliers  $M(B^p, Q_s)$   $(0 < s \le 1, 1 \le p < \infty)$ ,  $M(Q_{s_1}, Q_{s_2})$   $(0 < s_1, s_2 \le 1)$ , and  $M(B^p, B^q)$   $(1 \le p, q < \infty)$  have the f-property.

<span id="page-19-1"></span>**Theorem 4.2.** For  $\alpha > 0$  and 0 < s < 1 the space  $Q_{s,\log,\alpha}$  has the K-property.

Theorem 4.1 follows readily from the following result.

**Lemma 4.3.** Let X and Y be to Banach spaces of analytic functions which are continuously contained in  $H^1$ . Suppose that X contains the constants functions and that Y has the f-property. Then the space of multipliers M(X,Y) also has the f-property.

*Proof.* Since X contains the constants functions  $M(X,Y) \subset Y \subset H^1$ .

Suppose that  $F \in M(X,Y)$ , I is an inner function, and  $F/I \in H^1$ . Take  $f \in X$ . Then  $fF \in Y \subset H^1$  and then  $fF/I \in H^1$ . Since Y has the f-property, it follows that  $fF/I \in Y$ . Thus, we have proved that  $F/I \in M(X,Y)$ .  $\square$ 

Theorem 4.2 will follows from a characterization of the spaces  $Q_{s,\log,\alpha}$  in terms of pseudoanalytic continuation. We refer to Dyn'kin's paper [24] for similar descriptions of classical smoothness spaces, as well as for other important applications of the pseudoanalytic extension method.

Let,  $\mathbb{D}_{-}$  denotes the region  $\{z \in \mathbb{C} : |z| > 1\}$ , and write

$$z^* \stackrel{\text{def}}{=} 1/\overline{z}, \quad z \in \mathbb{C} \setminus \{0\}.$$

We shall need the Cauchy-Riemann operator

$$\overline{\partial} = \frac{\partial}{\partial \overline{z}} \stackrel{\text{def}}{=} \frac{1}{2} \left( \frac{\partial}{\partial x} + i \frac{\partial}{\partial y} \right), \quad z = x + iy.$$

The following result is an extension of [23, Theorem 1].

<span id="page-19-2"></span>**Theorem 4.4.** Suppose that 0 < s < 1,  $\alpha > 0$ , and  $f \in \bigcap_{0 < q < \infty} H^q$ . Then the following conditions are equivalent.

(i)  $f \in Q_{s,\log,\alpha}$ .

(ii) 
$$\sup_{|a|<1} \left( \log \frac{2}{1-|a|} \right)^{2\alpha} \int_{\mathbb{D}} |f'(z)|^2 \left( \frac{1}{|\varphi_a(z)|^2} - 1 \right)^s dA(z) < \infty.$$

(iii) There exists a function  $F \in C^1(\mathbb{D}_-)$  satisfying

$$\begin{split} F(z) &= O(1), \quad as \ z \to \infty, \\ \lim_{r \to 1^+} F(re^{i\theta}) &= f(e^{i\theta}), \quad a.e. \ and \ in \ L^q([-\pi, \pi]) \ for \ all \ q \in [1, \infty) \ , \\ \sup_{|a| < 1} \left( \log \frac{2}{1 - |a|} \right)^{2\alpha} \int_{\mathbb{D}_-} \left| \overline{\partial} F(z) \right|^2 \left( |\varphi_a(z)|^2 - 1 \right)^s \ dA(z) < \infty. \end{split}$$

Theorem 4.4 can be proved following the arguments used in the proof of [23, Theorem 1], we omit the details. Once Theorem 4.4 is established, noticing that if  $\alpha > 0$  and 0 < s < 1 then  $Q_{s,\log\alpha} \subset Q_s \subset BMOA$ , Theorem 4.2 can be proved following the steps in the proof of [23, Theorem 2]. Again, we omit the details.

**Acknowledgements.** We wish to thank the referees for reading carefully the paper and making a number of nice suggestions to improve it.

#### References

- <span id="page-20-7"></span>[1] A. Aleman and J. A. Cima, An integral operator on H<sup>p</sup> and Hardy's inequality, J. Anal. Math. **85** (2001), 157–176.
- <span id="page-20-9"></span>[2] A. Aleman, P. L. Duren, M. J. Martín and D. Vukotić, *Multiplicative isometries and isometric zero-divisors*, Canad. J. Math. **62** (2010), no. 5, 961-974.
- <span id="page-20-1"></span>[3] A. Aleman and A. Simbotin, *Estimates in Möbius invariant spaces of analytic functions*, Complex Var. Theory Appl. **49** (2004), no. 7-9, 487-510.
- <span id="page-20-6"></span>[4] A. Aleman and A. G. Siskakis, An integral operator on H<sup>p</sup>, Complex Variables Theory Appl. 28 (1995), no. 2, 149-158.
- <span id="page-20-8"></span>[5] A. Aleman and A. G. Siskakis, Integration operators on Bergman spaces, Indiana Univ. Math. J. 46 (1997), no. 2, 337–356.
- <span id="page-20-14"></span><span id="page-20-2"></span>[6] J. M. Anderson, On division by inner factors, Comment. Math. Helv. 54 (1979), no. 2, 309–317.
- [7] J. M. Anderson, J. Clunie and Ch. Pommerenke, On Bloch functions and normal functions, J. Reine Angew. Math. 270 (1974), 12–37.
- <span id="page-20-0"></span>[8] J. Arazy, S. D. Fisher and J. Peetre, Möbius invariant function spaces, J. Reine Angew. Math. 363 (1985), 110–145.
- <span id="page-20-10"></span>[9] N. Arcozzi, R. Rochberg and E. Sawyer, Carleson measures for analytic Besov spaces. Rev. Mat. Iberoamericana 18 (2002), no. 2, 443–510.
- <span id="page-20-5"></span>[10] R. Aulaskari and G. Csordas, Besov spaces and the  $Q_{q,0}$  classes, Acta Sci. Math. (Szeged) **60** (1995), no. 1–2, 31-48.
- <span id="page-20-13"></span>[11] R. Aulaskari, D. Girela and H. Wulan, Taylor coefficients and mean growth of the derivative of  $Q_p$  functions, J. Math. Anal. Appl. 258 (2001), no. 2, 415-428.
- <span id="page-20-3"></span>[12] R. Aulaskari and P. Lappan, Criteria for an analytic function to be Bloch and a harmonic or meromorphic function to be normal, Complex Analysis and its Applications (Harlow), Pitman Research Notes in Math, vol. 305, Longman Scientific and Technical, 1994, 136–146.
- <span id="page-20-12"></span>[13] R. Aulaskari, D. A. Stegenga and J. Xiao, Some subclasses of BMOA and their characterization in terms of Carleson measures, Rocky Mountain J. Math. 26 (1996), no. 2, 485-506.
- <span id="page-20-4"></span>[14] R. Aulaskari, J. Xiao and R. Zhao, On subspaces and subsets of BMOA and UBC, Analysis 15 (1995), no. 2, 101-121.
- <span id="page-20-11"></span>[15] L. Brown and A. L. Shields, Multipliers and cyclic vectors in the Bloch space, Michigan Math. J. 38 (1991), no. 1, 141–146.

- <span id="page-21-3"></span>[16] J. J. Donaire, D. Girela and D. Vukotić, On univalent functions in some Möbius invariant spaces, J. Reine Angew. Math. **553** (2002), 43–72.
- <span id="page-21-4"></span><span id="page-21-0"></span>[17] J. J. Donaire, D. Girela and D. Vukotić, On the growth and range of functions in Möbius invariant spaces, J. Anal. Math. 112 (2010), 237-260.
- [18] P. L. Duren, Theory of H<sup>p</sup> spaces, Academic Press, New York-London, 1970. Reprint: Dover, Mineola-New York, 2000.
- <span id="page-21-1"></span>[19] P. L. Duren and A. P. Schuster, Bergman Spaces, Math. Surveys and Monographs, Vol. 100, American Mathematical Society, Providence, Rhode Island, 2004.
- <span id="page-21-16"></span>[20] K. M. Dyakonov, Factorization of smooth analytic functions via Hilbert-Schmidt operators, (in Russian), Algebra i Analiz 8(1996), no. 4, 1–42. English translation in St. Petersburg Math. J. 8 (1997), no. 4, 543–569.
- <span id="page-21-17"></span>[21] K. M. Dyakonov, Equivalent norms on Lipschitz-type spaces of holomorphic functions, Acta. Math. 178 (1997), 143–167.
- <span id="page-21-18"></span>[22] K.M. Dyakonov, Holomorphic functions and quasiconformal mappings with smooth modulii, Adv. in Math. 187 (2004), 146–172.
- <span id="page-21-21"></span>[23] K. M. Dyakonov and D. Girela, On  $Q_p$  spaces and pseudoanalytic extension, Ann. Acad. Sci. Fenn. Math. **25** (2000), no. 2, 477-486.
- <span id="page-21-24"></span><span id="page-21-10"></span>[24] E. M. Dyn'kin, The pseudoanalytic extension, J. Anal. Math. 60 (1993), 45–70.
- [25] P. Galanopoulos, D. Girela and M. J. Martín, Besov spaces, multipliers and univalent functions, Complex Anal. Oper. Theory 7 (2013), no. 4, 1081-1116.
- <span id="page-21-8"></span>[26] P. Galanopoulos, D. Girela and J. A. Peláez, Multipliers and integration operators on Dirichlet spaces, Trans. Amer. Math. Soc. 363 (2011), no. 4, 1855-1886.
- <span id="page-21-12"></span>[27] D. Girela, Growth of the derivative of bounded analytic functions, Complex Variables Theory Appl. 20 (1992), no. 1–4, 221-227.
- <span id="page-21-6"></span>[28] D. Girela, Analytic functions of bounded mean oscillation, Complex Function Spaces, Mekrijärvi 1999 Editor: R. Aulaskari. Univ. Joensuu Dept. Math. Rep. Ser. 4, Univ. Joensuu, Joensuu (2001) pp. 61–170.
- <span id="page-21-23"></span>[29] D. Girela, C. González and J. A. Peláez, Multiplication and division by inner functions in the space of Bloch functions, Proc. Amer. Math. Soc. 134 (2006), no. 5, 1309-1314.
- <span id="page-21-5"></span>[30] D. Girela and N. Merchán, A generalized Hilbert operator acting on conformally invariant spaces, Banach J. Math. Anal. 12 (2018), no. 2, 374-398.
- <span id="page-21-13"></span>[31] D. Girela, M. Nowak and P. Waniurski, On the zeros of Bloch functions, Math. Proc. Cambridge Philos. Soc. 129 (2000), no. 1, 117-128.
- <span id="page-21-9"></span>[32] D. Girela and J. A. Peláez, Carleson measures, multipliers and integration operators for spaces of Dirichlet type, J. Funct. Anal. 241 (2006), no. 1, 334 358.
- <span id="page-21-22"></span>[33] V. P. Gurarii, On the factorization of absolutely convergent Taylor series and Fourier integrals, (in Russian) Zap. Naučn. Sem. Leningrad. Otdel. Mat. Inst. Steklov. (LOMI) 30 (1972), 15–32.
- <span id="page-21-14"></span>[34] V. P. Havin, On the factorization of analytic functions smooth up to the boundary (Russian), Zap. Naucn. Sem. Leningrad. Otdel. Mat. Inst. Steklov. (LOMI) 22 (1971) 202-205.
- <span id="page-21-20"></span>[35] H. Hedenmalm, On the f- and K-properties of certain function spaces, Contemporary Math. 91 (1989), 89–91.
- <span id="page-21-2"></span>[36] H. Hedenmalm, B. Korenblum and K. Zhu, Theory of Bergman Spaces, Graduate Texts in Mathematics, Vol. 199, Springer, New York, Berlin, etc. 2000.
- <span id="page-21-7"></span>[37] F. Holland and D. Walsh, Growth estimates for functions in the Besov spaces  $A_p$ , Proc. Roy. Irish Acad. Sect. A 88 (1988), 1–18.
- <span id="page-21-15"></span>[38] B. I. Korenblum, A certain extremal property of outer functions, (in Russian), Mat. Zametki 10 (1971), 53–56. English translation in Math. Notes 10 (1971), 456–458.
- <span id="page-21-19"></span>[39] B. I. Korenblum and V. M. Faĭvyševskiĭ, A certain class of compression operators that are connected with the divisibility of analytic functions, (in Russian), Ukrain. Mat. Z. 24 (1972), 692-695. English translation in Ukrainian Math. J. 24 (1973), 559–561.
- <span id="page-21-11"></span>[40] J. M. Ortega and J. Fàbrega, Pointwise multipliers and corona type decomposition in BMOA, Ann. Inst. Fourier (Grenoble) 46 (1996), no. 1, 111-137.

- <span id="page-22-19"></span><span id="page-22-15"></span>[41] J. Pau and J. A. Peláez, Multipliers of Möbius invariant  $Q_s$  spaces, Math. Z. **261** (2009), no. 3, 545-555.
- [42] J. Pau and J. A. Peláez, On the zeros of functions in Dirichlet-type spaces, Trans. Amer. Math. Soc. 363 (2011), no. 4, 1981-2002.
- <span id="page-22-9"></span>[43] Ch. Pommerenke, Schlichte Funktionen und analytische Funktionen von beschrakter mittlerer Oszillation, Comment. Math. Helv. **52** (1977), no. 4, 591-602.
- <span id="page-22-1"></span>[44] L. E. Rubel and R. M. Timoney, An extremal property of the Bloch space, Proc. Amer. Math. Soc. 75 (1979), no. 1, 45–49.
- <span id="page-22-11"></span><span id="page-22-6"></span>[45] D. Stegenga, Multipliers of the Dirichlet space. Illinois J. Math. 24 (1980), no. 1, 113–139.
- [46] K. Stroethoff, Besov-type characterisations for the Bloch space, Bull. Austral. Math. Soc. **39** (1989), no. 3, 405-420.
- <span id="page-22-10"></span><span id="page-22-2"></span>[47] R. M. Timoney, Natural function spaces, J. London Math. Soc. (2) 41 (1990), no. 1, 78-88.
- [48] S. A. Vinogradov, Multiplication and division in the space of analytic functions with area integrable derivative, and in some related spaces (in russian), Zap. Nauchn. Sem. S.-Peterburg. Otdel. Mat. Inst. Steklov. (POMI) 222 (1995), Issled. po Linein. Oper. i Teor. Funktsii 23, 45–77, 308. English translation in J. Math. Sci. (New York) 87, no. 5 (1997), 3806–3827.
- <span id="page-22-12"></span>[49] Z. Wu, Carleson measures and multipliers for Dirichlet spaces. J. Funct. Anal. 169 (1999), no. 1, 148–163.
- <span id="page-22-5"></span>[50] J. Xiao, Carleson measure, atomic decomposition and free interpolation from Bloch space, Ann. Acad. Sci. Fenn. Ser. A I Math. 19 (1994), 35–46.
- <span id="page-22-18"></span><span id="page-22-7"></span>[51] J. Xiao, The Q<sub>p</sub> corona theorem, Pacific J. Math. **194** (2000), no. 2, 491-509.
- <span id="page-22-8"></span>[52] J. Xiao, Holomorphic Q classes, Lecture Notes in Mathematics 1767, Springer-Verlag, 2001.
- <span id="page-22-17"></span>[53] J. Xiao, Geometric Q functions, Frontiers in Mathematics. Birkhäuser Verlag, 2006.
- <span id="page-22-14"></span>[54] J. Xiao, The Q<sub>p</sub> Carleson measure problem, Adv. Math. **217** (2008), no. 5, 2075-2088.
- <span id="page-22-4"></span>[55] R. Zhao, On logarithmic Carleson measures, Acta Sci. Math. (Szeged) 69 (2003), no. 3-4, 605-618.
- <span id="page-22-3"></span>[56] K. Zhu, Analytic Besov spaces, J. Math. Anal. Appl. 157 (1991), 318–336.
- [57] K. Zhu, A class of Möbius invariant function spaces, Illinois J. Math. 51 (2007), no. 3, 977-1002.
- <span id="page-22-0"></span>[58] K. Zhu, Operator Theory in Function Spaces, Marcel Dekker, New York, 1990. Reprint: Math. Surveys and Monographs, Vol. 138, American Mathematical Society, Providence, Rhode Island, 2007.
- <span id="page-22-13"></span>[59] N. Zorboska, Multiplication and Toeplitz operators on the analytic Besov spaces. In: More Progress in Analysis: Proceedings of the 5th International. Isaac Congress. Catania, Italy, 25–30 July 2005. Editors: H. G. W. Begehr and F. Nicolosi. World Scientific (2009), pp. 387–396.
- <span id="page-22-16"></span>[60] A. Zygmund, Trigonometric Series, Vol. I and Vol. II, Second edition, Camb. Univ. Press, Cambridge, 1959.

DEPARTAMENTO DE ANÁLISIS MATEMÁTICO, ESTADÍSTICA E INVESTIGACIÓN OPERATIVA, Y MATEMÁTICA APLICADA, UNIVERSIDAD DE MÁLAGA, 29071 MÁLAGA, SPAIN *E-mail address*: girela@uma.es

DEPARTAMENTO DE MATEMÁTICA APLICADA, UNIVERSIDAD DE MÁLAGA, 29071 MÁLAGA, SPAIN

E-mail address: noel@uma.es