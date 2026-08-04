from pathlib import Path
import json, hashlib, zipfile, shutil, textwrap, datetime
root=Path('/mnt/data/apn_checkpoint_work/APN_PROJECT_CHECKPOINT_2026-08-02')
chk=root/'checkpoint'; res=root/'results/current'; man=root/'manuscript'

priority = {
  "schema":"deep-priority-audit-v1",
  "date":"2026-08-02",
  "scope":"relative-trace rank-one and projective rank-two switchings of the Gold cube",
  "sources":[
    {"id":"BCL2009","type":"primary article/preprint","title":"Constructing new APN functions from known ones","authors":"L. Budaghyan, C. Carlet, G. Leander","year":2009,"locator":"Cryptology ePrint 2007/063; FFA 15(2), 150-159","audited":"full preprint","relevant":"constructs x^3+Tr(x^9); gives a general sufficient construction for adding a quadratic Boolean coordinate"},
    {"id":"EP2009","type":"primary article","title":"A new almost perfect nonlinear function which is not quadratic","authors":"Y. Edel, A. Pott","year":2009,"locator":"AMC 3(1), 59-81, Theorem 9 and dimension-8 catalogue","audited":"full article","relevant":"necessary and sufficient rank-one Boolean switching criterion; switching-equivalence catalogues"},
    {"id":"BFA2023","type":"primary extended abstract","title":"On quadratic APN functions F(x)+Tr(x)L(x)","authors":"H. Taniguchi","year":2023,"locator":"BFA 2023, pp. 123-127, Theorem 3","audited":"full extended abstract","relevant":"necessary and sufficient hyperplane criterion; computational counts"},
    {"id":"TPPA2025","type":"primary journal article","title":"Changing almost perfect nonlinear functions on affine subspaces of small codimensions","authors":"H. Taniguchi, A. Polujan, A. Pott, R. Arshad","year":2025,"locator":"EJC 32(4), #P4.61, Theorems 3.2, 3.5, 4.1 and Example 4.2","audited":"full article","relevant":"H-equivalence; all six-bit quadratic classes; codimension-two modifications; one overlapping 8-bit example"},
    {"id":"ARSHAD2018","type":"doctoral thesis","title":"Contributions to the theory of almost perfect nonlinear functions","authors":"R. Arshad","year":2018,"locator":"doi:10.25673/13406","audited":"metadata and attribution in TPPA2025 only; repository returned HTTP 429 for full text","relevant":"TPPA2025 states that parts of its work appeared in the thesis; exact theorem-level overlap remains unresolved"},
    {"id":"ZHENG2021","type":"primary preprint/article","title":"Constructing new APN functions through relative trace functions","authors":"L. Zheng, H. Kan, Y. Li, J. Peng, D. Tang","year":2021,"locator":"arXiv:2101.11535","audited":"abstract, introduction, constructions","relevant":"neighboring vector-valued relative-trace constructions; structurally different from the scalar product-of-traces family"}
  ],
  "claims":[
    {"id":"C1","claim":"centre-independent rank-r derivative-incidence/forbidden-set criterion","verdict":"qualified candidate novelty","confidence":"medium-high","boundary":"Rank-one Boolean switching is covered by Edel-Pott Theorem 9; the common trace-factor specialization is equivalent to BFA 2023 Theorem 3 / TPPA 2025 Theorem 3.2. Novelty may be claimed only for the arbitrary-kernel rank-r operator and precomputed forbidden-set formulation."},
    {"id":"C2","claim":"n=4 trace-zero coefficient characterization and EA reduction","verdict":"minor structural contribution","confidence":"medium","boundary":"No exact scalar statement located, but dimension four has one APN class and the conclusion creates no new class."},
    {"id":"C3","claim":"n=6 coefficients are exactly the six elements of order nine","verdict":"candidate novelty","confidence":"medium","boundary":"No explicit statement found in BCL2009, EP2009, BFA2023, or TPPA2025. Full Arshad thesis could not be audited; resulting global class is historical."},
    {"id":"C4","claim":"n=8 rank-one coefficient fibre is exactly GF(4)^*","verdict":"candidate novelty","confidence":"medium-high","boundary":"TPPA2025 Example 4.2 supplies one member of this fibre, not the complete coefficient classification. One associated class is historical EP No. 6."},
    {"id":"C5","claim":"no nonzero coefficient exists for every even n>=10","verdict":"strong candidate novelty","confidence":"high","boundary":"No dimension-wide nonexistence theorem for this exact scalar family located in audited sources."},
    {"id":"C6","claim":"complete P^1(GF(4)) rank-two displacement theorem and trace selector","verdict":"strong candidate novelty","confidence":"high","boundary":"TPPA2025 Theorem 4.1 gives a general codimension-two criterion, but not this normalized projective classification, selector, or exact bad-section counts."},
    {"id":"C7","claim":"involution-only obstruction: four bad sections for each of three involutions","verdict":"strong candidate novelty","confidence":"high","boundary":"No direct analogue located."},
    {"id":"C8","claim":"eight marked points partition into two EA/CCZ classes identified with EP Nos. 15 and 18","verdict":"new local organization of historical classes","confidence":"high","boundary":"The global classes are not new; explicit Frobenius EA cycles and their occurrence in the projective family are the contribution."},
    {"id":"C9","claim":"marked centre-relative switching geometry","verdict":"conceptual candidate novelty","confidence":"medium","boundary":"Must be formalized as an invariant before title-level use; EP2009 already defines switching equivalence."}
  ],
  "unresolved_risks":[
    "Full-text theorem-level audit of Arshad's 2018 thesis remains outstanding because the repository returned HTTP 429.",
    "EP2009 tables were audited for formulas and class IDs, but absence of an equivalent coefficient parametrization cannot be proved solely from keyword search; manuscript wording must remain 'we did not locate'.",
    "A final Crossref/MathSciNet/zbMATH citation search should be repeated immediately before submission."
  ],
  "safe_central_claim":"We completely classify a fixed scalar relative-trace switching subfamily of the Gold cube in all even dimensions, prove uniform nonexistence above dimension eight, and classify its exceptional dimension-eight rank-two extensions by a finite projective geometry. The reached global EA/CCZ classes are known; the contribution is the coefficient classification, dimension rigidity, projective selector, and centre-relative organization."
}
(res/'deep_priority_audit.json').write_text(json.dumps(priority,indent=2,ensure_ascii=False)+"\n")

audit_md=r'''# Stage 5W — deep priority audit

Date: 2026-08-02

## 1. Scope and method

The audit concerns the proposed paper on scalar relative-trace switchings of the Gold cube.  The question is not whether switching APN functions is known—it is—but whether the following exact statements have appeared previously:

1. the order-nine coefficient set in dimension six;
2. the \(\mathbb F_4^*\) coefficient fibre in dimension eight;
3. nonexistence for every even dimension at least ten;
4. the \(\mathbb P^1(\mathbb F_4)\) rank-two selector;
5. the partition of eight marked points into two EA/CCZ classes.

Primary sources were read at theorem level whenever full text was accessible.  A negative priority verdict is stated only when a direct overlap was located.  Otherwise the language is deliberately “no direct statement was located,” not “the result has never appeared.”

## 2. Source-by-source findings

### 2.1 Budaghyan–Carlet–Leander

The full ePrint version of *Constructing new APN functions from known ones* was audited.  Its principal switched-cube result is the APN family

\[
x^3+\operatorname{Tr}(x^9),
\]

and it gives a general sufficient construction for adding a quadratic Boolean coordinate to a known quadratic APN function.  It does not formulate the scalar family

\[
x^3+\theta\,\operatorname{Tr}(x)\operatorname{Tr}(\lambda x)
\]

as a dimension-classification problem, and no order-nine, \(\mathbb F_4^*\), all-even nonexistence, or projective-selector statement was located.

**Priority effect.**  Switched cubes and the idea of adding a Boolean coordinate are prior art.  Our paper must not present the family as the first trace switching of the Gold cube.

### 2.2 Edel–Pott

The full 2009 article was audited.  Proposition 3 and Theorem 9 give the general rank-one switching form

\[
F(x)+f(x)u
\]

and a necessary and sufficient APN condition in terms of four-tuples.  The paper also computes switching-equivalence classes in small dimensions; in dimension eight the EA switching class of \(x^3\) contains 17 CCZ-inequivalent functions.

This is the closest prior framework to our rank-one derivative-incidence lemma.  Therefore priority cannot be claimed for a general necessary-and-sufficient rank-one switching criterion.  What remains distinct is:

- the rank-\(r\), arbitrary-kernel matrix formulation;
- the forbidden-set preprocessing;
- the exact coefficient fibres for the chosen relative-trace kernels;
- the projective rank-two classification.

The global classes reached in dimension eight are already present in their catalogue.  In particular, our two rank-two projective classes correspond exactly to catalogue representatives 15 and 18.  This is historical identification, not a new-class claim.

No formula in the audited article states that the dimension-six coefficients are precisely the elements of order nine or that the dimension-eight scalar fibre is precisely \(\mathbb F_4^*\).  No \(\mathbb P^1(\mathbb F_4)\) selector was located.

### 2.3 Taniguchi, BFA 2023

The complete five-page extended abstract was audited.  Theorem 3 is the necessary and sufficient criterion for

\[
F(x)+\operatorname{Tr}(x)L(x)
\]

with \(F\) quadratic APN and \(L\) linear.  It also reports the counts 448 in dimension four, 4608 in dimension five, and approximately 40,000 in dimension six.

It contains neither the scalar coefficient classification nor the rank-two projective theorem.

### 2.4 Taniguchi–Polujan–Pott–Arshad 2025

The journal article was audited in full.

- Theorem 3.2 is the H-equivalence criterion for \(F+\operatorname{Tr}(x)L(x)\).
- Theorem 3.5 states that every quadratic six-bit APN EA-class has a representative \(x^3+\operatorname{Tr}(x)L(x)\).
- Theorem 4.1 gives a general necessary and sufficient condition for modifications by constants on the four cosets of a codimension-two subspace.
- Example 4.2 constructs
  \[
  x^3+\alpha^{85}\operatorname{Tr}_{1}^{8}(x)\operatorname{Tr}_{2}^{8}(x).
  \]

The last example is an exact overlap with one point of our rank-one dimension-eight geometry.  Internally it matches `CLASS-B`, already exactly CCZ-identified with Edel–Pott representative 6, \(x^9+\operatorname{Tr}(x^3)\).

Theorem 4.1 subsumes the *existence test* for arbitrary codimension-two constant modifications.  It does not supply our normalized projective parametrization, the diagonal synchronization condition, the trace selector, the involution-only obstruction, or the exact 8-to-2 class partition.

### 2.5 Arshad thesis

The thesis metadata and bibliographic record were located.  The 2025 article explicitly states that parts of its work appeared in the thesis.  However, the repository repeatedly returned HTTP 429 for the PDF, so a full theorem-level audit was not possible in this environment.

This is the main unresolved priority risk.  Claims concerning the hyperplane criterion and codimension-two modifications are already ceded to prior work regardless.  For the dimension-six order-nine set and the dimension-eight projective theorem, the safe wording is:

> We did not locate these explicit classifications in the accessible primary sources; the full Arshad thesis must still be checked before submission.

### 2.6 Relative-trace constructions

The Zheng–Kan–Li–Peng–Tang work studies vector-valued relative traces of quadratic functions and produces infinite APN families.  It is adjacent in terminology but structurally different from the scalar Boolean kernel

\[
Q_c(x)=\operatorname{Tr}_{K/\mathbb F_2}(cx)
       \operatorname{Tr}_{K/\mathbb F_2}(\lambda cx).
\]

It should be cited to delimit the terminology “relative trace,” not as a direct antecedent of the coefficient-fibre theorem.

## 3. Claim-level verdicts

| Claim | Verdict | Submission-safe wording |
|---|---|---|
| General rank-one switching criterion | Not novel | Cite Edel–Pott Theorem 9 and the H-equivalence literature. |
| Common-trace-factor criterion | Not novel | State exact equivalence with BFA 2023 Theorem 3 / 2025 Theorem 3.2. |
| Arbitrary-kernel rank-\(r\) forbidden-set formulation | Qualified candidate novelty | Present as a rank-update reformulation and extension, not as the first switching criterion. |
| \(n=4\) trace-zero characterization | Minor structural contribution | Useful base case; no new EA-class. |
| \(n=6\) order-nine characterization | Candidate novelty | “We determine explicitly…”; no new global class claim. |
| \(n=8\) coefficient fibre \(\mathbb F_4^*\) | Candidate novelty | Mention that one member is already Example 4.2 of the 2025 paper. |
| Every even \(n\ge10\): nonexistence | Strong candidate novelty | Headline dimension-rigidity theorem. |
| \(\mathbb P^1(\mathbb F_4)\) rank-two selector | Strong candidate novelty | Headline finite-geometric classification, while citing the general codimension-two theorem. |
| Involution-only obstruction | Strong candidate novelty | State as the mechanism behind the selector. |
| Eight points \(\to\) two EA/CCZ classes | New local organization | The classes themselves are historical representatives 15 and 18. |
| Marked centre-relative geometry | Conceptual candidate | Keep in discussion unless a formal invariant theorem is added. |

## 4. Final central claim

The following wording survives the audit:

> We completely classify a fixed scalar relative-trace switching subfamily of the Gold cube in every even dimension.  The family is nontrivial only in dimensions four, six and eight, and no nonzero coefficient is possible in any even dimension at least ten.  In dimension eight we classify the exceptional rank-two extensions by \(\mathbb P^1(\mathbb F_4)\), obtaining a trace selector and two centre-relative EA/CCZ orbits.  All reached global classes are known; the contribution is the coefficient classification, dimension rigidity, projective selector and marked local organization.

## 5. Mandatory pre-submission checks

1. Obtain and search the full Arshad thesis.
2. Repeat the exact-phrase and citation search in MathSciNet, zbMATH and Google Scholar shortly before submission.
3. Ask a specialist familiar with Edel–Pott switching tables to verify the priority boundary.
4. Avoid “first,” “new family,” and “new APN classes” unless a specific claim is independently confirmed.
'''
(chk/'STAGE5W_DEEP_PRIORITY_AUDIT_2026-08-02.md').write_text(audit_md)

matrix=r'''# Paper claims–evidence–bibliography matrix v3

Date: 2026-08-02

This version incorporates the full primary-source audit of BCL 2009, Edel–Pott 2009, BFA 2023 and the 2025 EJC paper.  The Arshad thesis remains a flagged full-text gap.

| ID | Manuscript claim | Evidence | Priority verdict | Exact literature boundary | Final treatment |
|---|---|---|---|---|---|
| C1 | Rank-\(r\) derivative-incidence criterion for \(F+U\circ q\) | Analytic | Qualified candidate novelty | EP09 Theorem 9 covers rank-one Boolean switching; BFA23 Theorem 3 and TPPA25 Theorem 3.2 cover the common trace factor. | Call it a rank-update reformulation and arbitrary-kernel extension. Do not claim the first switching criterion. |
| C2 | Precomputed forbidden sets reduce each candidate to \(2^r-1\) incidence tests | Exact algorithmic consequence | Candidate algorithmic contribution | No identical implementation located. | Proposition plus pseudocode and complexity. |
| C3 | \(n=4\): coefficients are the seven nonzero trace-zero elements and all outputs are EA-Gold | Analytic + exhaustive validation | Minor structural contribution | Small dimension has one APN class; no new class. | Short proposition establishing the base of the dimension ladder. |
| C4 | \(n=6\): coefficients are exactly the six elements of order nine | Compact 12-orbit computer-assisted theorem | Candidate novelty; medium confidence | No explicit statement located in audited papers. Thesis full text pending. Global class is historical Banff class 2. | Main exceptional-dimension theorem; include 12-row appendix certificate. |
| C5 | \(n=8\): rank-one fibre is exactly \(\mathbb F_4^*\) | Symbolic sufficiency + 33-orbit necessity certificate | Candidate novelty; medium-high confidence | TPPA25 Example 4.2 gives one member, not the full fibre. CLASS-B is EP representative 6. | Main theorem; explicitly cite the overlapping example. |
| C6 | Every even \(n\ge10\) has no nonzero coefficient | Artin–Schreier/Weil proof for \(n\ge14\); exact bridges \(n=10,12\) | Strong candidate novelty | No exact all-even theorem located. | Headline theorem. |
| C7 | Projective rank-two classification on \(\mathbb P^1(\mathbb F_4)\) | Analytic reduction + finite exact table | Strong candidate novelty | TPPA25 Theorem 4.1 is the general codimension-two criterion, not the normalized classification. | Headline finite-geometric theorem. |
| C8 | Only three semilinear involutions obstruct, four sections each | Analytic reduction + four-case table | Strong candidate novelty | No direct analogue located. | Key lemma in C7 proof. |
| C9 | Eight marked points form two EA/CCZ classes | Explicit EA cycles + distinguishing orthoderivative spectra | New local organization, historical global classes | Exactly EP representatives 15 and 18. | State both facts in the same theorem to avoid a new-class implication. |
| C10 | TPPA25 Example 4.2 matches CLASS-B / EP representative 6 | Exact stored CCZ witness and matching portable signature | Not novel; overlap anchor | Published example. | Relation-to-prior-work remark. |
| C11 | Centre-relative marked switching geometry refines global labels | Exact examples, no general invariant theorem yet | Conceptual candidate | EP09 already defines switching equivalence. | Discussion/future work unless formalized. |
| C12 | Current classes fail the complementary non-bent-subspace condition | Exact corrected exhaustive screen | Secondary negative result | Uses external permutation criterion. | Omit from main paper or place in supplement. |

## Claims that must not appear

- “We introduce switching of APN functions.”
- “We give the first criterion for \(F+\operatorname{Tr}(x)L(x)\).”
- “We construct new eight-bit APN classes.”
- “The eight marked points are eight classes.”
- “No earlier source contains the order-nine or projective results” before the Arshad thesis is checked.

## Safe abstract-level novelty sentence

> We determine the complete coefficient set of a fixed scalar relative-trace switching of the Gold cube in every even dimension, prove its uniform disappearance above dimension eight, and classify its exceptional dimension-eight rank-two extensions by a projective trace selector.
'''
(chk/'PAPER_CLAIMS_EVIDENCE_BIBLIOGRAPHY_MATRIX_V3_2026-08-02.md').write_text(matrix)

arch=r'''# Stage 5X — manuscript architecture

Date: 2026-08-02

## Working title

**Dimension Rigidity and Projective Geometry of Relative-Trace Switchings of the Gold Cube**

## 1. Fixed notation

Let \(n\) be even, \(K=\mathbb F_{2^n}\), \(E=\mathbb F_4\), and

\[
T=\operatorname{Tr}_{K/E},\qquad
\operatorname{tr}=\operatorname{Tr}_{K/\mathbb F_2}.
\]

Fix \(\lambda\in E\setminus\mathbb F_2\), so \(\lambda^2+\lambda+1=0\).  For \(c\in K^*\), define the Boolean quadratic form

\[
Q_c(x)=\operatorname{tr}(cx)\operatorname{tr}(\lambda cx).
\]

The normalized rank-one family is

\[
F_{n,\theta}(x)=x^3+\theta Q_1(x).
\]

For \(n=8\), \(K\) contains \(L=\mathbb F_{16}\).  For \(\rho\in L\setminus E\) and \(\eta,\zeta\in E^*\), define

\[
G_{\rho,\eta,\zeta}(x)
=x^3+\eta Q_1(x)+\zeta\rho^{-3}Q_\rho(x).
\]

The notation \([\rho]\) denotes a point of \(\mathbb P_E(L)\cong\mathbb P^1(E)\).  Literal Boolean forms are distinguished from their classes modulo linear functions whenever projective rescaling is used.

## 2. Theorem numbering and section plan

### Section 1. Introduction

- Prior switching constructions.
- Exact boundary with Edel–Pott Theorem 9, H-equivalence, and codimension-two modifications.
- Dimension-rigidity result and projective theorem.
- Explicit statement that no new global APN class is claimed.

### Section 2. Preliminaries and low-rank updates

**Lemma 2.1 (normalized polarization).**  State \(F(0)=q_i(0)=0\).

**Theorem 2.2 (rank-update criterion).**  Centre-independent derivative-incidence formulation.

**Corollary 2.3 (forbidden sets).**  Precomputation and \(2^r-1\) tests.

**Proposition 2.4 (relation to known switching criteria).**  Rank-one reduction to Edel–Pott and common-trace-factor equivalence to H-equivalence.

### Section 3. Relative-trace form and section products

**Lemma 3.1.**  Polar form of \(Q_c\).

**Lemma 3.2.**  Normalize \(c\) by input/output scaling.

**Theorem 3.3 (section-product criterion).**

\[
F_{n,\theta}\text{ is APN}
\iff
\theta\notin S_n,
\qquad
S_n=\{a^3+b^3:T(a)=1,T(b)=0\}.
\]

### Section 4. Dimension four

**Proposition 4.1.**  APN iff \(\operatorname{tr}(\theta)=0\), \(\theta\ne0\).

**Proposition 4.2.**  Explicit EA reduction to \(x^3\).

### Section 5. Dimension six

**Theorem 5.1.**  APN iff \(\theta^6+\theta^3+1=0\), equivalently \(\operatorname{ord}(\theta)=9\).

**Lemma 5.2.**  Kummer-coordinate reduction.

**Appendix reference.**  Twelve Frobenius-orbit witnesses.

### Section 6. Dimension eight, rank one

**Theorem 6.1.**  APN iff \(\theta\in E^*\).

**Lemma 6.2.**  Transverse two-space product formula.

**Lemma 6.3.**  Symbolic exclusion of products in \(E^*\) using the resultant/Bezout certificate.

**Appendix reference.**  Thirty-three Frobenius-orbit witnesses for necessity.

### Section 7. Uniform nonexistence above dimension eight

**Theorem 7.1.**  For every even \(n\ge10\) and \(\theta\ne0\), \(F_{n,\theta}\) is not APN.

**Lemma 7.2.**  Character expansion for \(h_n(\theta)\).

**Lemma 7.3.**  Geometry of the Fermat cubic and pole divisors.

**Lemma 7.4.**  Artin–Schreier character-sum bounds.

**Lemma 7.5.**  Exact bridge cases \(n=10,12\).

### Section 8. Projective rank-two geometry in dimension eight

**Theorem 8.1 (projective selector).**

\[
G_{\rho,\eta,\zeta}\text{ APN}
\iff
\eta=\zeta,
\quad
\operatorname{Tr}_{E/\mathbb F_2}(\delta([\rho])\eta)=1.
\]

**Lemma 8.2.**  Four non-base projective displacements.

**Lemma 8.3.**  Linear transitions are unobstructed.

**Lemma 8.4.**  Three semilinear involutions, four bad sections each.

### Section 9. EA/CCZ partition and historical identification

**Theorem 9.1.**  Eight marked points form exactly two EA/CCZ classes of four.

**Proposition 9.2.**  Explicit Frobenius EA cycle.

**Proposition 9.3.**  Identification with Edel–Pott representatives 15 and 18.

**Remark 9.4.**  The published codimension-two example belongs to a different rank-one class, Edel–Pott representative 6.

### Section 10. Computation and reproducibility

- Discovery versus proof.
- Certificates and hashes.
- What is computer-assisted.

### Section 11. Conclusion

- Dimension rigidity.
- Centre-relative geometry.
- Generalization to non-Gold quadratic centres deferred to a separate paper.

## 3. Proof-dependency summary

The logical spine is

\[
\text{Theorem 2.2}
\Longrightarrow
\text{Theorem 3.3}
\Longrightarrow
\begin{cases}
\text{Theorems 4.1, 5.1, 6.1},\\
\text{Theorem 7.1},\\
\text{Theorem 8.1}.
\end{cases}
\]

Theorem 9.1 depends on Theorem 8.1 only for completeness of the eight marked points; its equivalence proof uses explicit Frobenius EA maps and orthoderivative invariants independently.

## 4. Computer-assisted boundary

- Theorem 5.1 necessity: 12-row orbit certificate.
- Theorem 6.1 necessity: 33-row orbit certificate.
- Theorem 6.1 sufficiency: symbolic resultant/Bezout certificate, independently verifiable.
- Theorem 7.1: analytic for \(n\ge14\), exact finite enumeration for \(n=10,12\).
- Theorem 8.1: analytic reduction plus a four-case table over \(\mathbb F_{16}\).
- Theorem 9.1: identities checked on all 256 inputs, with formulas included.

## 5. Recommended paper length

- Main text: 24–30 pages.
- Appendix tables: 3–5 pages.
- Electronic supplement: scripts, JSON certificates, exact EA witnesses.

The Gröbner campaign chronology should occupy no more than one paragraph in the introduction and one short reproducibility subsection.
'''
(chk/'STAGE5X_MANUSCRIPT_ARCHITECTURE_2026-08-02.md').write_text(arch)

dep=r'''# Manuscript theorem dependency graph

```text
Definitions and normalization
        |
        v
Theorem 2.2: rank-r derivative-incidence criterion
        |
        +--> Corollary 2.3: forbidden sets
        |
        +--> Proposition 2.4: EP09 / H-equivalence specializations
        |
        v
Lemma 3.1 + Lemma 3.2
        |
        v
Theorem 3.3: section-product criterion
        |
        +--> Proposition 4.1 --> Proposition 4.2          (n=4)
        |
        +--> Lemma 5.2 + Appendix A --> Theorem 5.1      (n=6)
        |
        +--> Lemmas 6.2,6.3 + Appendix B --> Theorem 6.1 (n=8 rank 1)
        |
        +--> Lemmas 7.2--7.5 --> Theorem 7.1             (n>=10)
        |
        +--> Theorem 6.1 + Lemmas 8.2--8.4
                    |
                    v
              Theorem 8.1: projective rank-two selector
                    |
                    v
              Theorem 9.1: eight marked points -> two classes
                    |              \
                    |               +--> orthoderivative separation
                    +--> explicit Frobenius EA cycles
```

## External dependencies

- Edel–Pott Theorem 9: prior rank-one switching criterion.
- BFA 2023 Theorem 3 / TPPA 2025 Theorem 3.2: H-equivalence specialization.
- TPPA 2025 Theorem 4.1: prior general codimension-two modification criterion.
- Yoshiara: EA and CCZ equivalence coincide for quadratic APN functions.
- Standard Artin–Schreier/Weil bound on character sums over curves.

## Independence statements

- The dimension classification does not depend on the historical class identifications.
- The exact EA/CCZ partition does not depend on Sage or `sboxU`.
- The nonexistence theorem does not depend on the Gröbner campaign.
'''
(chk/'MANUSCRIPT_THEOREM_DEPENDENCY_GRAPH_2026-08-02.md').write_text(dep)

# Appendix tables
cert=json.load(open(res/'compact_section_orbit_certificates.json'))

def row_tex(row,n):
    rep=row.get('representative_hex', str(row.get('representative_code', row.get('representative','?'))))
    osz=row['orbit_size']; mult=row.get('representation_multiplicity','')
    w=row.get('section_witness',{})
    # n6 may use different keys
    parts=[]
    for key in ['a_hex','b_hex','p_hex','q_hex','r_hex','s_hex']:
        if key in w: parts.append(f"{key[:-4]}={w[key]}")
    wit=', '.join(parts) if parts else json.dumps(w,sort_keys=True)
    poly='+'.join('X' if e==1 else ('1' if e==0 else f'X^{{{e}}}') for e in row.get('minimal_polynomial_exponents',[]))
    return f"{rep} & {osz} & {mult} & ${poly}$ & \\texttt{{{wit}}} \\\\"

lines=[r'''% Auto-generated from results/current/compact_section_orbit_certificates.json
\section{Compact Frobenius-orbit certificates}
\label{app:orbit-certificates}

The section-product sets are invariant under the binary Frobenius map.  Thus one explicit section witness for every Frobenius orbit proves coverage of the entire orbit.

\subsection{Dimension six}
\begin{longtable}{llll}
\toprule
Orbit representative & Orbit size & Multiplicity & Section witness \\
\midrule
\endhead''']
for row in cert['n6']['rows']:
    rep=row.get('representative_hex', str(row.get('representative_code', row.get('representative','?'))))
    osz=row['orbit_size']; mult=row.get('representation_multiplicity','')
    w=row.get('section_witness',{})
    wit=', '.join(f"{k.replace('_hex','')}={v}" for k,v in w.items())
    lines.append(f"\\texttt{{{rep}}} & {osz} & {mult} & \\texttt{{{wit}}} \\\\")
lines += [r'''\bottomrule
\end{longtable}

The uncovered nonzero elements are exactly the six roots of $X^6+X^3+1$, equivalently the elements of multiplicative order nine.

\subsection{Dimension eight}
\begin{longtable}{lllll}
\toprule
Orbit representative & Orbit size & Multiplicity & Minimal polynomial support & Section witness \\
\midrule
\endhead''']
for row in cert['n8']['rows']:
    rep=row['representative_hex']; osz=row['orbit_size']; mult=row['representation_multiplicity']
    exps=row.get('minimal_polynomial_exponents',[])
    poly='\\{'+','.join(map(str,exps))+'\\}'
    w=row['section_witness']; wit=', '.join(f"{k.replace('_hex','')}={v}" for k,v in w.items())
    lines.append(f"\\texttt{{{rep}}} & {osz} & {mult} & ${poly}$ & \\texttt{{{wit}}} \\\\")
lines += [r'''\bottomrule
\end{longtable}

The only uncovered field elements are the four elements of $\mathbb F_4$; hence the nonzero APN coefficients are precisely $\mathbb F_4^*$.
''']
(man/'APPENDIX_ORBIT_CERTIFICATES_N6_N8.tex').write_text('\n'.join(lines)+'\n')

uniform=r'''\section{Uniform nonexistence above dimension eight}
\label{sec:uniform-nonexistence}

Let $q=2^n$, where $n$ is even, and let $K=\mathbb F_q$ and $E=\mathbb F_4$.  For $\theta\in K^*$ put
\[
h_n(\theta)=\#\{(a,b)\in K^2:T(a)=1,\ T(b)=0,\ a^3+b^3=\theta\}.
\]
By the section-product criterion, $F_{n,\theta}$ is APN if and only if $h_n(\theta)=0$.

\subsection{Character expansion}
Let $\chi_K$ and $\chi_E$ be the canonical additive characters.  Orthogonality gives
\[
h_n(\theta)=\frac1{16}\sum_{u,v\in E}\chi_E(u)
\sum_{a^3+b^3=\theta}\chi_K(ua+vb).
\]
The inner sum is taken over the affine part of the smooth projective Fermat cubic
\[
C_\theta:\quad X^3+Y^3=\theta Z^3.
\]
The curve is nonsingular because its partial derivatives are $X^2$, $Y^2$, and $\theta Z^2$, which have no common projective zero.  Hence $C_\theta$ has genus one.  Its three points at infinity are $[t:1:0]$ with $t\in E^*$.

\subsection{Pole divisors}
For an axis phase, say $u\ne0$ and $v=0$, the rational function $uX/Z$ has three simple poles at the points at infinity.  The Artin--Schreier conductor degree is therefore six.

Assume now $u,v\ne0$.  At the point $[v/u:1:0]$ the numerator $uX+vY$ vanishes.  The line $uX+vY=0$ is the flex tangent there: substituting $X=(v/u)Y$ into the cubic leaves $\theta Z^3$.  Thus the zero has intersection multiplicity three and cancels the pole completely at that point.  The remaining two points at infinity are simple poles, so the conductor degree is four.

A rational function with a simple pole cannot be of the form $g^2+g+c$, because every pole of $g^2+g$ has even order.  Hence all nonconstant phases satisfy the non-coboundary hypothesis for the Artin--Schreier bound.

\subsection{Estimate}
The Hasse error for the trivial phase is at most $2\sqrt q$.  There are six axis pairs, each bounded by $6\sqrt q$, and nine off-axis pairs, each bounded by $4\sqrt q$.  Accounting also for the affine/projective correction yields
\[
h_n(\theta)\ge \frac{q-2-(2+6\cdot6+9\cdot4)\sqrt q}{16}
=\frac{q-2-74\sqrt q}{16}.
\]
The right-hand side is positive for every even $n\ge14$.  Exact enumeration gives
\[
\min_{\theta\ne0}h_{10}(\theta)=46,
\qquad
\min_{\theta\ne0}h_{12}(\theta)=208.
\]
Consequently $h_n(\theta)>0$ for every even $n\ge10$ and every $\theta\ne0$.  By the section-product criterion, no nonzero member of the scalar relative-trace family is APN in these dimensions.
'''
(man/'SECTION_UNIFORM_NONEXISTENCE_DRAFT.tex').write_text(uniform)

skeleton=r'''\documentclass[11pt]{article}
\usepackage{amsmath,amssymb,amsthm,booktabs,longtable,hyperref}
\newtheorem{theorem}{Theorem}[section]
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{remark}[theorem]{Remark}
\newcommand{\F}{\mathbb F}
\newcommand{\Tr}{\operatorname{Tr}}
\newcommand{\tr}{\operatorname{tr}}

\title{Dimension Rigidity and Projective Geometry of Relative-Trace Switchings of the Gold Cube}
\author{Draft architecture -- authors omitted}
\date{}

\begin{document}
\maketitle

\begin{abstract}
We classify a fixed scalar relative-trace switching subfamily of the Gold APN function in every even dimension.  The family is nontrivial only in dimensions four, six, and eight; for every even dimension at least ten we prove uniform nonexistence.  In dimension eight we classify the exceptional rank-two extensions through the projective line over $\F_4$ and obtain a trace selector.  The resulting functions lie in known global EA/CCZ classes; the contribution is the coefficient classification, dimension rigidity, projective selector, and centre-relative organization.
\end{abstract}

\section{Introduction}
% Write only after the full Arshad-thesis check.  Explicitly cede priority for
% Edel--Pott rank-one switching, H-equivalence, and general codimension-two modifications.

\section{Preliminaries and low-rank derivative updates}
\section{The scalar relative-trace family}
\section{Dimension four}
\section{Dimension six}
\section{Dimension eight: rank one}
\input{SECTION_UNIFORM_NONEXISTENCE_DRAFT}
\section{Projective rank-two extensions in dimension eight}
\section{EA/CCZ partition and historical classes}
\section{Computation and reproducibility}
\section{Conclusion}

\appendix
\input{APPENDIX_ORBIT_CERTIFICATES_N6_N8}

\end{document}
'''
(man/'RELATIVE_TRACE_SWITCHING_MANUSCRIPT_SKELETON_2026-08-02.tex').write_text(skeleton)

nextp=r'''# Next-stage plan after Stage 5W/5X

Date: 2026-08-02

## Immediate priority: close the remaining bibliographic risk

1. Obtain the full PDF of Razi Arshad's 2018 thesis through an alternative mirror, library access, or direct request.
2. Search it for `Tr(x)L(x)`, codimension-two modifications, products of traces, order-nine coefficients, and dimension-eight examples.
3. Update only the affected priority-confidence fields; do not change mathematical results.

## Stage 5Y — manuscript-grade theorem sections

1. Write Sections 2 and 3 in full: low-rank update theorem, forbidden sets, and section-product criterion.
2. Write the complete $n=4$, $n=6$, and $n=8$ rank-one proofs around the compact certificates.
3. Integrate the completed uniform nonexistence section and independently line-check every character factor and affine/projective correction.
4. Write the projective rank-two theorem in coordinate-free form.
5. Add the exact Frobenius EA maps and historical IDs 15 and 18.

## Stage 5Z — internal referee pass

1. Produce a theorem-by-theorem checklist: hypotheses, normalization, field embeddings, equivalence type, and computer-assisted dependency.
2. Attempt adversarial counterexample tests for every statement in dimensions 4, 6, and 8.
3. Compile the LaTeX manuscript and verify appendix tables against canonical JSON.
4. Only after this pass, draft the introduction and select a journal.

## Deferred Stage 6

Build a centre-independent derivative-incidence atlas for non-Gold quadratic APN representatives, beginning with six-bit Banff/Kim classes and the cached quadratic Kasami-class centres.
'''
(chk/'NEXT_STAGE_PLAN_AFTER_STAGE5X.md').write_text(nextp)

# Validation report
files=[
 chk/'STAGE5W_DEEP_PRIORITY_AUDIT_2026-08-02.md',
 chk/'PAPER_CLAIMS_EVIDENCE_BIBLIOGRAPHY_MATRIX_V3_2026-08-02.md',
 chk/'STAGE5X_MANUSCRIPT_ARCHITECTURE_2026-08-02.md',
 chk/'MANUSCRIPT_THEOREM_DEPENDENCY_GRAPH_2026-08-02.md',
 man/'APPENDIX_ORBIT_CERTIFICATES_N6_N8.tex',
 man/'SECTION_UNIFORM_NONEXISTENCE_DRAFT.tex',
 man/'RELATIVE_TRACE_SWITCHING_MANUSCRIPT_SKELETON_2026-08-02.tex',
 res/'deep_priority_audit.json',
 chk/'NEXT_STAGE_PLAN_AFTER_STAGE5X.md'
]
val={"schema":"stage5wx-validation-v1","date":"2026-08-02","checks":{},"files":[]}
for p in files:
    b=p.read_bytes(); val['files'].append({"path":str(p.relative_to(root)),"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()})
val['checks']={
 "priority_audit_has_unresolved_thesis_flag": "HTTP 429" in audit_md,
 "matrix_forbids_new_class_claim": "We construct new eight-bit APN classes" in matrix,
 "appendix_n6_row_count": len(cert['n6']['rows']),
 "appendix_n8_row_count": len(cert['n8']['rows']),
 "expected_n6_rows":12,
 "expected_n8_rows":33,
 "uniform_constant_74_present":"74\\sqrt q" in uniform,
 "historical_ids_15_18_present":"15 and 18" in arch,
 "manuscript_skeleton_created":True
}
(chk/'STAGE5WX_VALIDATION_REPORT.json').write_text(json.dumps(val,indent=2)+"\n")
files.append(chk/'STAGE5WX_VALIDATION_REPORT.json')

# Delta zip preserving root-relative paths
zip_path=Path('/mnt/data/APN_STAGE5WX_PRIORITY_MANUSCRIPT_2026-08-02.zip')
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for p in files:
        z.write(p, arcname=str(p.relative_to(root)))
print(zip_path)
print('n6 rows',len(cert['n6']['rows']),'n8 rows',len(cert['n8']['rows']))
print('zip sha256',hashlib.sha256(zip_path.read_bytes()).hexdigest())
