# CY Landscape Explorer — état du projet

> **Ce document remplace intégralement les versions précédentes**, y compris
> l'addendum de session `CY_Landscape_Explorer_SESSION_EQUIVARIANCE.md`, qui est
> supprimé. La version antérieure décrivait l'équivariance de f comme bloquée par
> le cocycle de `#7669`, listait six candidats « réellement contraints », et
> annonçait 15 tests. Aucune de ces trois affirmations ne tient : le verrou est
> levé, la partition en « six contraints » mesurait la portée de l'ancien test et
> non une propriété des candidats, et la suite compte 39 tests.

---

## 1. Objet

Recherche de fibrés vectoriels stables à trois générations sur les variétés de
Calabi–Yau de la liste CICY d'Oxford (7 890 entrées), pour la compactification de
la corde hétérotique E₈×E₈.

Constructions employées : monades classiques et monades positives ; la branche
**extensions**, longtemps désactivée pour cause de défaut (§4.7), est rouverte
sur un chemin de calcul correct (§5.10) et constitue désormais la ligne
principale — la voie des monades positives paraissant fermée au rang 5 (§5.9).

---

## 2. Résultat principal à ce jour

Chaîne complète : monade → stabilité de Hoppe → symétrie librement agissante →
idéal Γ-covariant → f équivariante → stabilité **restreinte au sous-espace
équivariant** → surjectivité de f.

Point de départ : scan `scan_wilson2` sur les 194 CICYs à symétrie librement
agissante (classification de Braun, JHEP 1104 (2011) 005), 122 min sur 7 cœurs.
115 fibrés Hoppe-stables, 0 rejeté par l'audit, 108 passant le test nécessaire
sur les charges.

> **Reproductibilité (§5.23).** Ces deux candidats sortaient d'un générateur qui
> tirait **dix** configurations au hasard dans une famille de 2 201, sur un RNG
> partagé — ils ont d'ailleurs disparu du scan suivant, sans qu'aucun filtre ne
> les élimine. La famille est maintenant **énumérée** : leur présence est
> démontrée et non plus tirée au sort. Le scan qui fait foi est `scan_wilson4`,
> à relancer (§9) ; `scan_wilson2` reste la seule trace produite avant la
> correction.

Balayage d'équivariance sur ces 108 candidats (≈ 1 h, mono-cœur) : 71 dans le
domaine du modèle S/I, 3 878 couples (candidat, symétrie, λ).

| | |
|---|---|
| tués par h⁰(V) équivariant | 3 452 |
| indéterminés (certificat de surjectivité hors de portée) | 423 |
| **survivants** | **3 couples, soit 2 candidats** |

**Les deux candidats retenus :**

| CICY | jauge | rang | cohomologie | Γ | λ | n_gen amont | générations |
|---|---|---|---|---|---|---|---|
| **6890** | SO(10) | 4 | [0, 6, 0, 0] | ℤ₂ | **+1 seul** | 6 | 3 |
| **6947** | SO(10) | 4 | [0, 6, 0, 0] | ℤ₂ | **+1 seul** | 6 | 3 |

Trois remarques, dans l'ordre d'importance.

**`#6890` était invisible pour l'ancien test.** Γ n'y agit que par des phases,
donc `equivariance.py` le marquait `test_non_trivial = False` : le test sur les
charges ne disait **rien** à son sujet, ni pour ni contre. Le test au niveau des
polynômes, lui, dit quelque chose, et `#6890` passe. C'est toute la catégorie
laissée en suspens par l'ancienne partition qui se rouvre.

**Corollaire : la partition « 6 réellement contraints / 25 non concluants » de
l'ancienne version n'a plus lieu d'être.** Elle mesurait la puissance du test
disponible à l'époque, pas une propriété des candidats.

**La structure équivariante n'est pas libre.** Sur les deux candidats, une seule
des deux relèvements de ℤ₂ donne un fibré : λ = +1 est certifié surjectif,
λ = −1 présente un déficit de rang stable. Rien avant le §5.4 ne voyait cette
distinction.

**Le critère de Hoppe est vérifié en entier sur eux**, restreint au sous-espace
équivariant : h⁰(∧^p V) = 0 pour p = 1, 2, 3 — h³(V) inclus (§5.5). Plus la
surjectivité de f, certifiée à λ = +1 (§5.4).

**Et les trois générations sont établies par décomposition explicite** :
H¹(V) = 3 invariants + 3 anti-invariants sous ℤ₂ (§5.6), et non par division de
6 par 2.

**Le spectre est calculé** : 16 → H¹(V) = 3 + 3 (§5.6), 10 → H¹(∧²V) = 8, dont
2 + 6 sous ℤ₂ (§5.7). Avec une ligne de Wilson ℤ₂ en Pati–Salam, cela donne
**3 générations complètes (4,2,1) + (4̄,1,2)**, et 2 ou 6 bidoublets de Higgs
selon la corrélation choisie (§5.8).

**Ces deux-là ne sont pas des modèles du Modèle Standard, et ne peuvent pas
l'être.** Les lignes de Wilson préservent le rang ; SO(10) est de rang 5, le
groupe du MS de rang 4. Avec |Γ| = 2 on plafonne à Pati–Salam ou SU(5) flipped
(§5.8). Aller plus loin demande un Γ plus gros — la liste de Braun en offre — ou
un mécanisme de brisure supplémentaire.

---

## 3. Architecture

Le dépôt contient **huit** points d'entrée `main_*.py`. Un seul est maintenu.
Les sept autres sont antérieurs, n'ont reçu aucune des corrections du §4 ni du
§5 — **ni l'annulation d'anomalie**, qui est une condition physique — et
**refusent désormais de tourner** (§5.22). Deux d'entre eux, `main_full_scan` et
`main_monads`, étaient déjà morts : ils importent `cohomology_end_V_approx`,
supprimée par la réécriture de `monad_wedge` (§4.3).

```
cy_landscape/
├── main_optimized.py          scan principal (multiprocessing, checkpoint)
├── main*.py  (7 autres)       OBSOLÈTES — refusent de tourner (§5.22)
└── core/
    ├── intersection.py        nombres d'intersection d_ijk, c2(TY)·J
    ├── chi_exact.py           χ par Riemann-Roch — EXACT, sert de préfiltre
    ├── exact_cohomology.py    h^i(O(a)) par Koszul + d₁ + certification
    ├── koszul_exact.py        variante historique (idée de d₁, jamais branchée)
    ├── monads.py              MonadBundle, cohomologie de V, non-dégénérescence
    ├── positive_monads.py     générateur positif (élagage exact)
    ├── monad_wedge.py         ∧²V : χ exact + bornes rigoureuses
    ├── hoppe_fast.py          critère de Hoppe, phase 0 incluse
    ├── sections.py            anneaux de sections, rangs réels mod p
    ├── braun_symmetry.py      lecture COMPLÈTE de CicyQuotients.m           [+]
    ├── covariant_ring.py      idéal Γ-covariant, CovariantRing              [+]
    ├── equivariant_monad.py   f équivariante, h⁰(∧^p V(−H)) restreint     [+]
    │                          tout p (donc Hoppe complet, h³ inclus),
    │                          certificat de surjectivité,
    │                          décomposition de H¹(V) et H¹(∧²V) sous Γ
    ├── extensions.py          fibré d'extension : χ et rang corrects,       [+]
    │                          Hoppe par bornes sur les quotients gradués,
    │                          domaine ÉNUMÉRÉ (monotone en max_charge),
    │                          cohomologie par bornes rigoureuses,
    │                          pente : certificat d'instabilité exact
    ├── gamma_action.py        action de Γ sur les sections (obsolète, §5.2)
    └── cohomology.py          extraction du spectre (partiellement obsolète, §6)

racine/
├── tests_regression.py        39 tests — À LANCER AVANT CHAQUE SCAN
├── validate_cohomology.py     harnais de validation du socle
├── audit_results.py           triage 1 : cohérence interne
├── triage_clean.py            triage 2 : n_anti, familles, doublons
├── verify_hoppe.py            re-vérification a posteriori de la stabilité
├── wilson_match.py            croisement avec la liste de Braun
├── equivariance.py            test nécessaire sur les charges
└── equivariance_f.py          chaîne complète au niveau des polynômes       [+]
```

`[+]` = ajouté lors de la session « équivariance », ou réparé depuis (§5.10,
§5.12).

---

## 4. Défauts trouvés et corrigés

Chacun a été détecté par une **référence indépendante**, jamais par le code
lui-même. C'est la leçon centrale du projet.

### 4.1 `intersection.py` — c2 croisé compté deux fois

`compute_c2_tangent` posait `c2_amb[r,s] = (n_r+1)(n_s+1)` puis sommait sur les
couples **ordonnés**. Le coefficient hors diagonale doit être la moitié.

- Bicubique P²×P²[3,3] : donnait 63, valeur connue **36**.
- Quintique : 50 dans les deux cas — m = 1, aucun terme croisé.
- χ(L) entier : 34 % avant, **100 %** après (2 800 tirages).

**Le bug avait survécu à sa propre validation, faite sur la seule quintique.**

### 4.2 `exact_cohomology.py` — d₁ jamais calculé

La différentielle d₁ agit à q fixé : chaque ligne de la page E₁ est un complexe
dont il faut prendre l'homologie. Le code lisait directement les anti-diagonales.

- Accord avec Riemann-Roch : 40,0 % → **71,2 %**.
- Quintique : h⁰(O(5)) = **125** (l'écart de 1 venait de là).

Ajouts : χ **toujours exact** (somme alternée de E₁, indépendante de la
dégénérescence), et une **certification par degré**. Taux de certification de h¹
et h² sur `cicylist.txt` : **43,7 %**, avec accord HRR et Serre à **100 %** sur ce
domaine.

### 4.3 `monad_wedge.py` — hypothèse de rang maximal

La v1 supposait toutes les applications induites de rang maximal. Sur 456
monades, **79 %** des vecteurs de Betti contredisaient le χ calculé par la même
fonction. Sur les monades de rang 3, où ∧²V ≅ V* donne la réponse exacte,
**1 cas juste sur 34** (écarts jusqu'à 7 251 840).

v2 : χ(∧²V) exact et inconditionnel, bornes rigoureuses, drapeau de
détermination, raccourci exact au rang 3. **Ne retourne jamais un nombre
inventé.** Vérifications : χ(∧²V) = −χ(V) au rang 3 et = 0 au rang 4,
**218/218** chacune. Limite : les bornes ne déterminent h^i que dans ~0,1 % des
cas.

### 4.4 `monads.py` — même hypothèse sur H^i(V)

Le nombre de générations restait juste (il vaut |χ(V)|), mais **la répartition
entre h¹ et h² était fausse** — donc n_anti, donc tout le classement. Accord de
l'ancienne version avec la version rigoureuse : **0 / 11**.

`compute_monad_cohomology_ex` : χ exact, bornes rigoureuses, h⁰ = h³ = 0 imposé
par la stabilité, contrainte χ croisée avec les intervalles.

### 4.5 `hoppe_fast.py` — le critère n'était pas testé

Pour c₁(V) = 0, Hoppe s'énonce : V stable ⟺ h⁰(∧^p V) = 0 pour p = 1..rk−1, soit
**H = 0**. Le code ne testait que H = e_i, strictement **plus faible**.

Phase 0 ajoutée, avec le tableau des isomorphismes (det V = O, dualité de Serre) :

| rang | tests |
|---|---|
| 3 | h⁰(V), h³(V) |
| 4 | h⁰(V), h⁰(∧²V), h³(V) |
| 5 | h⁰(V), h⁰(∧²V), h³(∧²V), h³(V) |

**Bug associé** : `_wedge3_h0_twisted` implémente h⁰(∧³V) = h³(V), valable au
rang 4 seulement, alors que le code l'appelait dès rk ≥ 4.

### 4.6 `monads.py` — monades scindées

`check_map_exists` ne vérifiait que l'existence d'une entrée non nulle par ligne.
Trois situations passaient : un b_i égal à un c_j (f_i isomorphisme, monade
scindée) ; une colonne de f nulle (O(b_i) facteur direct) ; rang structurel de f
< rank_C. Cas réel, CICY 7669 : `C = [-2,0,-3]` et `b_1 = [-2,0,-3]`.

`check_monad_nondegenerate` teste les trois par comparaisons de charges. Effet :
**305 → 115** candidats sur le scan Wilson.

**Réserve, non levée** : ces contrôles portent uniquement sur les CHARGES. Ils ne
vérifient pas que f est surjective, donc pas que V est un fibré. Voir §5.4 et §6.

### 4.7 Branche `extension` — désactivée

Le pipeline construisait une pseudo-monade B = F₁⊕F₂, C = F₂. Le noyau est de
rang rank(F₁), alors que le fibré d'extension est de rang rank(F₁)+rank(F₂).
Mesure sur `test_v3` : **1 571 entrées sur 1 571** en incohérence de rang.
Désactivée par défaut à l'époque. Le chemin correct est décrit au §5.10 ;
l'option qui réactivait ce chemin cassé n'existe plus (§5.12).

### 4.8 Compteurs de spectre structurellement nuls

Dans `cohomology.py` : `sp.n_exotics = 0` en dur pour SO(10) ; pour SU(5),
`max(0, n_10 + n_10bar − n_gen − 2·n_anti)` vaut **identiquement zéro** puisque
|a−b| + 2·min(a,b) = a+b. Seul E₆ compte réellement ses anti-générations.

Le « zéro exotique » de tous les candidats SO(10) et SU(5) est donc une
constante, pas un résultat — et il valait 25 points gratuits dans le score.
**Corrigé (§5.19)** : ces quantités valent désormais `None`, et une quantité non
calculée ne rapporte plus rien.

### 4.9 `sections.py` — `reduce_vec` ne réduisait rien

Trouvé par la définition même d'une réduction : tout élément de I_a doit avoir un
reste nul.

`rref_mod` renvoie `(rang, pivots)` et travaille sur une **copie**
(`M = M % p` crée un nouveau tableau) : la matrice de l'appelant restait brute.
`Ring.quotient` faisait `Mred = M[:rank]`, gardant les `rank` premières lignes
**non réduites**, et `reduce_vec` soustrayait `v[c]·Mred[r]` en supposant un
pivot égal à 1.

- Avant correctif : **30/30** éléments de l'idéal ressortaient non nuls.
- `dimY` n'est pas touché (il ne dépend que du nombre de pivots) : toute la
  partie « dimensions » du pipeline était juste.
- `h0_V_explicit` et `h0_wedge2_V_explicit` l'étaient en revanche. Contre une
  référence construite sans aucun appel à `reduce_vec` — rang du composé
  ⊕S_{b_i} → S_c/I_c, moins Σ dim I_{b_i} — l'ancien code donne **7/9**, le
  nouveau **9/9**, dont 4 cas à h⁰ non nul. Les deux écarts vont dans le sens
  dangereux : l'ancien annonçait h⁰(V) = 0 là où il vaut 1, donc **déclarait
  stable un fibré qui ne l'est pas**.
- Effet sur le catalogue existant : **aucun**. Sur les 71 monades de
  `scan_wilson2` dans le domaine valide, h⁰(V) vaut 0 avant comme après.

Correctif : `sections.rref_mod_full` renvoie la forme échelonnée **réduite**, et
`Ring.quotient` l'utilise. `Ring` accepte par ailleurs un paramètre `p` (défaut :
P = 32003), nécessaire parce que l'anneau covariant travaille dans un GF(p)
contenant les racines de l'unité de Γ.

---

## 5. Équivariance — le verrou est levé

### 5.1 Polynômes Γ-covariants

`sections.Ring` tirait les coefficients des polynômes définissants **au hasard**.
Pour que l'action descende au quotient X/Γ, il faut que Γ préserve l'idéal, ce
qu'un tirage aléatoire ne fait pas : tout test bâti dessus portait sur une
variété n'admettant pas l'action annoncée.

Le second bloc de générateurs de Braun — l'action sur les K polynômes — donne la
donnée manquante :

```
p_α(g·x) = Σ_β N[α][β] p_β(x)
```

Système linéaire sur les coefficients, résolu par
`covariant_ring.resoudre_covariants`. `verifier_covariance` **recontrôle par
re-substitution**, sans réutiliser le système résolu.

**Mesure sur les 195 CICYs de Braun** : 1 689 symétries sur 1 695 résolues,
covariance revérifiée **1 689/1 689** avec écart exactement nul. Les 6 restantes :
4 sans convention non dégénérée, 1 illisible, 1 hors taille.

La convention `N` — et non sa transposée ni son inverse — est la bonne : elle est
la seule valide pour les groupes d'ordre 8 et 16 de `#6947`.

**Contrôle de non-dégénérescence** : la fonction de Hilbert de l'anneau covariant
est comparée à celle de l'anneau aléatoire. Accord **126/126**. Le choix covariant
est un point particulier de la famille, pas une dégénérescence.

`braun_symmetry.py` lit la structure entière par équilibrage d'accolades, sans
regex sur la structure — les noms de Braun contiennent `#` et `$`, les entrées
contiennent `rt[n]` et `Exp[...]`. 195/195 entrées, 1 695 symétries, blocs
coordonnées/polynômes appariés partout.

### 5.2 Existence d'un f équivariant

Pour V = ker(f : B → C), la substitution S_g envoie H⁰(O(a)) sur
H⁰(O(a∘σ⁻¹)), donc

```
S_g( f_{j,i} ) = λ_g · f_{ρ(j), π(i)}
```

Les degrés se recollent : (c_j − b_i)∘σ⁻¹ = c_{ρ(j)} − b_{π(i)}.

**L'ancienne `gamma_action.espace_f_equivariant` écrivait la relation à
l'envers** — elle appliquait S_g à f_{ρ(j),π(i)} pour la comparer à f_{j,i}, ce
qui donne un degré en σ⁻². Son propre garde-fou
`if deg_img != degres[j][i]: return None` s'en apercevait et abandonnait dès que
σ n'était pas une involution. C'est pourquoi elle ne concluait jamais. Elle reste
en place mais **ne doit plus être appelée**.

**Le relèvement projectif est calculé, pas supposé.** λ_g n'est pas supposé être
une racine de l'unité d'ordre |g| : on calcule le plus petit n tel que
T_g^n = c·Id, on lit c, et λ_g parcourt les n racines n-ièmes de c dans GF(p).
Cette énumération englobe le groupe de Heisenberg de `#7669` au lieu de supposer
qu'il n'existe pas. L'ancienne version disait qu'il fallait « travailler avec
l'extension centrale » : l'ordre projectif fait cela sans la construire.

Si aucune puissance jusqu'à n_max n'est scalaire, le module renvoie
`operateur sans puissance scalaire` et **ne conclut pas**.

`equivariant_monad.verifier_descente` contrôle explicitement que
S_g(I_a) ⊆ I_{a∘σ⁻¹} — faux avec des polynômes aléatoires, vrai avec l'idéal
covariant — et le test refuse de tourner sinon.

### 5.3 Le test qui mord : la stabilité restreinte

**Point de méthode, à ne pas perdre.** Quand Γ agit par phases, la contrainte sur
f est une condition de vecteur propre : le sous-espace équivariant représente
environ 1/|Γ| de l'espace total et **n'est jamais vide**. Conclure de sa non
vacuité que le fibré descend serait la même erreur que le « zéro exotique » du
§4.8 — un résultat identiquement vrai.

On recalcule donc h⁰(V) avec f tiré dans le sous-espace équivariant, et on le
compare à la valeur générique dans le **même** anneau covariant, ce qui isole
l'effet de la contrainte. h⁰ étant semi-continu supérieurement, on prend le
minimum sur plusieurs tirages. C'est ce test qui élimine 3 452 couples sur 3 878.

**h⁰(∧²V) sous contrainte** (seconde moitié de Hoppe au rang ≥ 4) est calculé de
même. Validation par deux références indépendantes : sans contrainte le chemin
redonne `sections.h0_wedge2_V_explicit` **71/71**, et ces valeurs sont **71/71** à
l'intérieur des bornes rigoureuses de `monad_wedge`.

Il **mord** — h⁰(∧²V) augmente sous contrainte dans 1 060 couples — mais il n'a
**éliminé personne à lui seul** :

| | h⁰(∧²V) = 0 | h⁰(∧²V) > 0 |
|---|---|---|
| h⁰(V) = 0 | 426 | **0** |
| h⁰(V) > 0 | 2 392 | 1 060 |

La case en haut à droite est vide. À prendre comme une **observation, pas comme
un théorème** : soit h⁰(V) = 0 entraîne h⁰(∧²V) = 0 dans ce régime, soit le
régime testé est trop étroit pour exhiber le contre-exemple.

### 5.4 Surjectivité de f

**Le trou.** Tout suppose que V = ker(f) est un *fibré*, donc f surjective en
tout point. Rien ne le garantit sur le sous-espace équivariant, qui pourrait être
contenu dans le lieu où f chute de rang.

**Le critère.** Pour rank_C = 1, f est surjective ⟺ les f_i n'ont aucun zéro
commun sur Y. Soit J = (f₁,…,fₙ) ⊂ R. S'il existe **un** multidegré d ≥ 0 avec
J_d = R_d, il n'y a pas de zéro commun : en un point y de P^{n₁}×…×P^{n_m}, au
moins un monôme de multidegré d est non nul, donc R_d ne peut pas s'annuler
entièrement en y, alors que J_d le devrait. **Suffisant, jamais faux positif** ;
un échec ne conclut pas.

**Résultat sur les deux candidats — la surjectivité dépend du caractère :**

| CICY | λ | espace entier | sous-espace équivariant |
|---|---|---|---|
| 6890 | **+1** | certifié | **certifié en (0,1,5,1,0)** |
| 6890 | −1 | certifié | **non certifié** — déficit stable |
| 6947 | **+1** | certifié | **certifié en (5,1,0,1,0)** |
| 6947 | −1 | certifié | **non certifié** — déficit stable |

Ces certificats sont obtenus à des degrés où la source est suffisante (24 = 24),
donc le contraste est réel et non un artefact.

**Une fausse piste, consignée pour ne pas la refaire.** J'avais d'abord lu un
déficit de rang constant sur `#7300` comme la signature d'un lieu de base de
dimension positive. **C'était faux.** La mesure qui manquait est la dimension de
la SOURCE : aucun degré testé ne vérifiait source ≥ cible, et le « déficit »
n'était que cible − source. `f_sans_point_base` calcule désormais la source et
écarte ces degrés avec la mention `source insuffisante` — le rang y est majoré
par la source, l'échec est arithmétique et ne dit rien.

**Portée actuelle du critère** : concluant au rang 4 et sur les petits N, **hors
de portée au rang 5**. Le long d'un seul axe, la pente de dim R_{a+t·e_k} dépend
du degré de base : la cible croît plus vite que chaque source et l'écart ne se
referme pas (mesuré jusqu'à t = 11). Seule la croissance uniforme
d = c + t·(1,…,1) marche asymptotiquement — le terme dominant (1/6)Σd_ijk·t³ ne
dépend pas du degré de base — mais le coût croît en t³ et sort du budget dès
t = 2. D'où les 423 indéterminés du §2, qui **ne disent rien** : ni que ces
monades sont des fibrés, ni le contraire.

### 5.5 h³(V), et le critère de Hoppe en entier

Le dernier tiers manquait, et il manquait pour une raison de chemin : `hoppe_fast`
traite h³(V) par les bornes rigoureuses de `compute_monad_cohomology_ex`, **sans
f explicite**, donc rien à restreindre à un sous-espace.

La sortie passe par l'énoncé même du critère. Pour c₁(V) = 0, Hoppe dit
h⁰(∧^p V) = 0 pour p = 1..rk−1, et **toutes** ces quantités s'obtiennent par le
même noyau : pour C de rang 1, la résolution de ∧^p V donne

```
0 → ∧^p V → ∧^p B → ∧^{p−1}B ⊗ C → …
h⁰(∧^p V) = dim ker( H⁰(∧^p B) → H⁰(∧^{p−1}B ⊗ C) )
```

l'application étant la contraction par f,
e_{i₁}∧…∧e_{i_p} ↦ Σ_k (−1)^{k−1} f_{i_k} · e_{I∖i_k}.

Et det V = O donne ∧^{rk−1}V ≅ V*, donc par Serre sur un CY3 :

```
h³(V)      = h⁰(∧^{rk−1} V)
h³(∧²V)    = h⁰(∧^{rk−2} V)        (rang 5)
```

`h0_wedgep_V_sur_espace` calcule le cas général, `hoppe_sur_espace` balaie
p = 1..rk−1. Le tableau des isomorphismes du §4.5 devient inutile : au lieu de
traduire h³ puis de chercher un chemin de calcul, on calcule directement tous les
h⁰(∧^p V). Le piège de `_wedge3_h0_twisted` — l'identité ∧³V = V*, vraie au rang 4
seulement — disparaît de lui-même puisque l'exposant est un paramètre.

**Validation.** Deux contrôles de cohérence — p = 1 et p = 2 redonnent
`h0_V_sur_espace` et `h0_wedge2_V_sur_espace`, écrites séparément — et surtout un
contrôle contre une **valeur connue d'avance** : à p = rk, ∧^{rk}V = det V = O,
donc h⁰ doit valoir exactement h⁰(O_Y) = **1**. Obtenu sur `#6890` et `#6947`.
C'est le seul des trois qui teste les **signes** de la contraction ; vérifié en
supprimant le (−1)^{k−1}, auquel cas h⁰(∧⁴V) tombe à 0 et le test échoue.

**Résultat** : h³(V) = 0 sur les deux candidats, pour λ = +1 comme pour λ = −1.
h³ ne les élimine donc pas, et le critère de Hoppe est désormais vérifié en
entier sous contrainte d'équivariance. C'est la surjectivité, et elle seule, qui
départage les deux caractères.

### 5.6 Décomposition de H¹(V) sous Γ — 3 générations, vérifiées

**H¹(V) est ici un conoyau explicite.** La suite longue de 0 → V → B → C → 0
donne H⁰(B) → H⁰(C) → H¹(V) → H¹(B). Si H¹(B) = 0, alors
H¹(V) = coker(H⁰(B) → H⁰(C)), c'est-à-dire R_c modulo l'image de f — calculable
dans l'anneau covariant. Γ agit sur R_c par S_g, et l'image est stable puisque f
est équivariante, donc l'action descend au conoyau.

**H¹(B) = 0 est vérifié, pas supposé.** `koszul_cohomology_ex` le certifie pour
chacun des cinq b_i des deux candidats : h¹ = 0, `certified_by_degree[1] = True`.
Sans cette vérification, H¹(V) ne serait pas ce conoyau et tout le calcul
porterait sur autre chose.

**Le [0, 6, 0, 0] est un vrai h¹.** Deux routes indépendantes concordent :
`compute_monad_cohomology_ex` donne des bornes **déterminantes** — (6,6) pour h¹
et (0,0) pour h² — et le comptage direct h⁰(C) − h⁰(B) = 16 − 10 = 6 par Koszul
certifié redonne la même valeur. Ce n'est donc pas un |χ| déguisé, ce qui était
la crainte légitime après le §4.4.

**Résultat.**

| | R_c | image de f | **H¹(V)** |
|---|---|---|---|
| invariants (λ = +1) | 8 | 5 | **3** |
| anti-invariants (λ = −1) | 8 | 5 | **3** |
| total | 16 | 10 | 6 |

Identique sur `#6890` et `#6947`.

**3 + 3 est une valeur connue d'avance**, pas une observation : Γ agissant
librement, n_gen(X/Γ) = n_gen(X)/|Γ|, donc la partie invariante doit valoir
exactement 3. Une décomposition en 4+2 ou 5+1 aurait signalé une erreur dans
l'action, dans le conoyau ou dans le fibré. **Les trois générations sont donc
obtenues par décomposition explicite, et non en divisant 6 par 2.**

Les multiplicités sont calculées en exploitant la semi-simplicité (p impair,
groupe fini) : elles sont additives le long de 0 → im f → R_c → H¹(V) → 0, ce
qui évite de choisir un supplémentaire. Que leurs sommes retombent sur 16 et 10
valide cette semi-simplicité au lieu de la supposer — et c'est ce contrôle qui
tombe en premier si l'on fausse le calcul.

### 5.7 H¹(∧²V) — les 10 de SO(10)

Pour V de rang 4, le groupe de structure est SU(4) et le commutant dans E₈ est
SO(10) : les générations sont les **16** (H¹(V), §5.6) et les **10** — d'où
sortent les doublets de Higgs après brisure — sont H¹(∧²V).

**Deux suites courtes, pas une résolution longue.** La filtration de ∧²B et la
suite de la monade tensorisée par C donnent

```
0 → ∧²V → ∧²B → V⊗C → 0
0 → V⊗C → B⊗C → C²  → 0
```

la seconde identifiant V⊗C au noyau de B⊗C → C². D'où

```
H⁰(V⊗C)   = ker( H⁰(B⊗C) → H⁰(C²) )
H¹(∧²V)   = coker( H⁰(∧²B) → H⁰(V⊗C) )
```

la seconde égalité exigeant H¹(∧²B) = 0 — **certifié** par
`koszul_cohomology_ex` sur les dix charges b_i + b_j de chaque candidat.

**Résultat, identique sur `#6890` et `#6947` :**

| | dim | invariants | anti-invariants |
|---|---|---|---|
| H⁰(V⊗C) | 45 | 22 | 23 |
| image de α | 37 | 20 | 17 |
| **H¹(∧²V)** | **8** | **2** | **6** |

Soit **8 copies du 10** en amont, dont **2 invariantes** sous ℤ₂.

**Trois contrôles, dont un croisé.** β∘α = 0 — c'est un complexe, et la
construction serait fausse sinon. Les multiplicités sont additives (45 = 22+23,
37 = 20+17). Et surtout h⁰(∧²V) = dim ker α **retombe sur 0**, la même valeur
que celle donnée par `h0_wedge2_V_sur_espace`, qui l'obtient par un tout autre
chemin (§5.3) : deux constructions indépendantes de la même quantité.

**L'étiquetage est bien défini — vérifié, pas argumenté.** La liberté résiduelle
du relèvement est M → μM avec μ² = 1, soit μ = −1 : sur R_a, S_g est alors
multiplié par μ^{|a|}, où |a| est le degré **total**. L'étiquetage ne bascule
donc que si |a| est impair. Ici |c| = 5 (impair) mais |b_i + c| = 6 (pair) :
H¹(V) bascule — invisible, sa décomposition 3 + 3 étant symétrique — et
H¹(∧²V) **ne bascule pas**. Recalculé avec −M : 2 + 6 à l'identique sur les deux
candidats. Le compte de Higgs repose donc sur une quantité bien définie.

Au passage, λ = −1 donnerait h¹(∧²V) = 10 réparti en 6 + 4 : les deux structures
équivariantes ont des spectres **différents**, ce qui est une raison de plus pour
que le certificat de surjectivité du §5.4, qui élimine λ = −1, ne soit pas
optionnel.

### 5.8 Ligne de Wilson — ce que ℤ₂ permet, et ce qu'il interdit

**Une limitation de principe, à énoncer d'abord.** Une ligne de Wilson est un
homomorphisme Γ → SO(10) ; le groupe non brisé est le commutant de son image.
Cette brisure **préserve le rang**. SO(10) est de rang 5, le groupe du Modèle
Standard de rang 4 : aucune ligne de Wilson, quel que soit Γ, ne peut mener
directement au Modèle Standard. Avec |Γ| = 2, les commutants d'un élément
d'ordre 2 sont, à conjugaison près :

| image de Γ | groupe non brisé |
|---|---|
| — | SU(5) × U(1) — SU(5) « flipped » |
| — | SO(6) × SO(4) ≅ SU(4) × SU(2)_L × SU(2)_R — Pati–Salam |
| — | SO(8) × SO(2), SO(2) × SO(8) |

**Ces deux candidats donnent donc au mieux Pati–Salam ou SU(5) flipped**, pas le
Modèle Standard. Descendre jusqu'à SU(3)×SU(2)×U(1) demanderait un mécanisme
supplémentaire, ou un Γ plus gros — ce que la liste de Braun offre sur d'autres
CICYs.

**Le spectre Pati–Salam.** Sous SU(4)×SU(2)_L×SU(2)_R :
16 → (4, 2, 1) + (4̄, 1, 2) et 10 → (6, 1, 1) + (1, 2, 2). L'élément d'ordre 2
vaut −1 sur la partie SO(6) et +1 sur la partie SO(4), donc il sépare les deux
morceaux de chaque représentation avec des signes opposés. Ne survivent au
quotient que les états dont la valeur propre sous Γ égale celle sous la ligne de
Wilson. Avec les décompositions calculées aux §5.6 et §5.7 :

| | amont | invariants + anti | survivants sur X/Γ |
|---|---|---|---|
| 16 = H¹(V) | 6 | 3 + 3 | **3 (4,2,1) + 3 (4̄,1,2)** |
| 10 = H¹(∧²V) | 8 | 2 + 6 | **2 (1,2,2) + 6 (6,1,1)** |

Les **trois générations Pati–Salam complètes** sortent indépendamment du signe
relatif choisi : la décomposition 3 + 3 étant symétrique, les deux corrélations
donnent le même compte. C'est le résultat robuste.

Le contenu en Higgs, lui, **dépend du choix** : (1,2,2) — qui contient la paire
H_u, H_d — hérite soit des 2, soit des 6, selon la corrélation retenue entre Γ et
la ligne de Wilson. **Ce choix est un intrant de construction de modèle, pas une
conséquence de la géométrie.** Deux bidoublets est la lecture
phénoménologiquement intéressante ; six est tout aussi permis par ce qui précède.

Ce paragraphe est de la théorie des groupes appliquée à des nombres calculés,
non un calcul supplémentaire : aucune ligne de Wilson n'a été construite
explicitement, et le code n'en manipule pas.

### 5.9 Monades positives au rang 5 — voie probablement fermée

Le Modèle Standard exige un fibré de **rang 5** : le groupe de structure est
alors SU(5), de rang 4 comme le groupe du MS, et une ligne de Wilson ℤ₂
— diag(1,1,1,−1,−1), disponible sur **166 des 194 CICYs** — a pour commutant
exactement SU(3)×SU(2)×U(1). Pas besoin d'un Γ d'ordre élevé. La cible est donc
un rang 5 avec |χ(V)| = 6.

**Ces monades existent.** Sur 95 430 monades de rang 5 engendrées, le |χ| minimum
non nul vaut exactement **6**, avec 13 cas à 12. L'absence observée dans les
scans venait du filtre : `--wilson` fixe la cible à 3·|Γ|, et
`wilson_gros_gamma.json` ne contenait que les ordres ≥ 4 — les |χ| = 6 étaient
rejetés avant tout test.

**Mais elles sont toutes instables.** Sur les 24 monades de rang 5, rank_C = 1,
|χ| ∈ {6, 12}, portées par une CICY à ℤ₂, réparties sur neuf variétés :

| motif d'élimination | cas |
|---|---|
| H³(V) ≥ 4 | 20 |
| H³(V) ≥ 2 | 2 |
| H³(V) ≥ 1 | 2 |
| **stables** | **0** |

Ce sont des **bornes inférieures rigoureuses** : une borne strictement positive
prouve h³(V) > 0, donc la non-stabilité. Ces 24 fibrés sont démontrés instables,
pas seulement non certifiés.

Vingt-quatre cas ne sont pas une démonstration générale. Mais le motif est
uniforme, sans exception ni cause alternative, et il a une lecture géométrique :
**au rang 5, un indice faible force une cohomologie qui casse la stabilité**. Ceci
explique aussi pourquoi les 832 candidats de rang 5 du scan `scan_su5` avaient
tous |χ| ∈ {24, 48} : les petits indices existent mais ne survivent pas à Hoppe.

**Conséquence pratique : ne pas relancer de scan de monades visant le MS.** Le
`scan_su5` complet représentait ~50 h et ne trouverait rien si ce constat tient.

### 5.10 Branche `extension` — rouverte sur des bases saines

Seule construction non explorée du dépôt, et désormais la ligne principale.

**Le défaut 4.7, chiffré.** Sur `#7669`, F1 = O(1,0,0)⊕O(0,1,0) et F2 = F1* :
la pseudo-monade annonce rang 2 et χ = 6, le fibré d'extension a rang 4 et
χ = **0**. Ce n'est pas un écart, c'est un autre objet.

**Le chemin correct** (`extensions.chi_extension`, `extensions.hoppe_extension`)
n'utilise plus la pseudo-monade :

- χ(V) = χ(F1) + χ(F2), exact par additivité sur les suites exactes ;
- Hoppe par **borne supérieure** : ∧^p V a une filtration de quotients gradués
  ∧^a F1 ⊗ ∧^b F2 (a+b = p), tous sommes de fibrés en droites, donc tous
  calculables exactement par Koszul. Si toutes les sommes h⁰ sont nulles pour
  p = 1..rk−1, alors h⁰(∧^p V) = 0 et V est stable. **Suffisant, jamais de faux
  positif** ; sinon `indéterminé`, jamais « stable ».

**Le taux de conclusion est de ~36 %**, contre **0,1 %** pour les bornes de
`monad_wedge` (§4.3). C'est ce qui rend la branche exploitable : F1 et F2 étant
des sommes de fibrés en droites, tous les gradués le sont aussi, là où ∧²V d'une
monade ne se laisse encadrer que grossièrement.

Le test `Ext¹(F2,F1) ≠ 0` élimine par ailleurs les cas où seule l'extension
scindée existe — et une somme directe n'est jamais stable.

**Ces quatre points sont faits** — énumération, test de régression, branchement
dans `process_cicy`, renommage de l'option. Voir §5.12.

### 5.11 Audits de générateurs — trois défauts, un motif commun

Aucun n'a été trouvé en lisant le code ; tous en lui posant trois questions :
énumère-t-il ou tire-t-il ? est-il monotone en ses propres paramètres ? le RNG
est-il partagé ?

**`_find_positive_B` échantillonne sans saturer.** 50 tirages parmi N ∈ [2·10⁴,
3·10⁵]. Cinquante graines successives ne découvrent que ~120 B distincts, soit
**0,04 %**, avec une croissance strictement linéaire — aucune saturation. Un
résultat d'ABSENCE obtenu ainsi ne porte sur rien. D'où `enumerer_positive_B` et
l'option `--exhaustif-max` (§9).

**La branche rank_C = 2 n'était pas monotone.** Elle puisait dans le RNG partagé,
donc ses tirages dépendaient de ce que rank_C = 1 avait consommé. Activer
`exhaustif_max`, qui ne touche pourtant que rank_C = 1, faisait **disparaître 20
candidats sur 42**. Corrigé par un RNG dérivé propre à la branche.

**`generate_extensions` ne l'était pas non plus, et bien plus gravement.**
Passer de 200 à 800 tirages perdait 40 extensions au rang 4 ; passer de
`max_charge` 2 à 3 en perdait **216 sur 222**, soit 97 %. Le premier point a été
corrigé par un RNG dérivé par couple (rk1, rk2) ; le second ne pouvait pas
l'être par un réglage de RNG — changer les bornes du tirage change la suite
tirée. Il fallait énumérer : c'est fait (§5.12).

**Deux pièges d'outillage rencontrés au passage.** `enumerer_positive_B`
matérialisait toute l'énumération en mémoire : `MemoryError` après 153 min de
scan, avec sept workers. Elle est devenue un générateur — le plafond bornait le
NOMBRE de B, pas la mémoire. Et `equivariance_f.py` n'écrivait que les lignes
`etat == 'ok'` : le JSONL montrait 0 couple sur un groupe d'ordre compatible
alors que 26 candidats en portaient un, écartés en amont pour une raison que le
fichier ne contenait pas. Toutes les lignes sont désormais persistées.


### 5.12 Branche `extension` — énumérée, testée, branchée

Les quatre points laissés ouverts au §5.10 sont traités.

**1. `generate_extensions` énumère.** Le domaine est le produit de boîtes
[−max_charge, max_charge]^m sur les charges de F1 et de F2, la dernière étant
fixée par c₁(V) = 0. `enumerer_extensions` le parcourt en entier ;
`compte_extensions` en donne le cardinal ordonné **par convolution, sans rien
construire**, ce qui permet de refuser d'énumérer avant d'essayer.

La monotonie en `max_charge` devient vraie **par construction** — des boîtes
emboîtées — et non par chance de RNG. Mesure : sur tous les couples testés,
E(q) ⊂ E(q+1) sans une seule perte, là où le tirage en perd 87 % en moyenne.

Trois points de mise en œuvre valent d'être notés :

- **Aucun ensemble de déduplication n'est conservé.** Les multiensembles sont
  produits à indices non décroissants, donc chacun une fois. C'est exactement
  le piège qui avait causé le `MemoryError` de `enumerer_positive_B` (§5.11) :
  là-bas le plafond bornait le nombre de B, pas la mémoire.
- **Élagage exact** dans la recherche de F2 à somme prescrite : les `t` vecteurs
  restants contribuent au plus `t·q` par composante, donc une cible hors de
  cette portée coupe la branche sans perte.
- **Le domaine du tirage a été aligné.** Il acceptait la dernière charge
  jusqu'à `max_charge + 2` alors que les autres étaient tirées dans
  [−max_charge, max_charge] — artefact de la construction « tirer n−1 vecteurs,
  déduire le dernier », sans justification géométrique, et qui faisait dépendre
  le domaine de l'**ordre** des facteurs de F2. Sans cet alignement, l'inclusion
  « tirage ⊂ énumération » n'aurait pas été testable.

Taille du domaine, en tuples ordonnés :

| | max_charge 1 | 2 | 3 |
|---|---|---|---|
| m = 2, rang 4 | 722 | 14 450 | 106 722 |
| m = 3, rang 4 | 13 718 | 1 228 250 | 24 652 782 |
| m = 3, rang 5 | 265 302 | 1,1·10⁸ | 6,1·10⁹ |

D'où le plafond `--ext-exhaustif-max` (défaut 200 000, comme
`enumerer_positive_B`) : au-delà, retour à l'échantillonnage, et le résultat
porte le champ `ext_mode = 'echantillonne'`. **Sans ce champ un résultat
d'absence serait ininterprétable.**

**2. Le test de régression (26ᵉ test).** Quatre références, dont une négative :

| volet | référence |
|---|---|
| comptage | `compte_extensions` (convolution, ne construit rien) contre un comptage par énumération explicite écrit dans le test |
| monotonie | E(q) ⊂ E(q+1), la propriété visée |
| inclusion | tout tirage se retrouve dans l'énumération — c'est ce point qui autorise « aucun survivant sur le domaine » |
| **contrôle négatif** | le tirage, lui, doit être **vu** non monotone : il perd 245 extensions sur 281 d'un cran au suivant |

Vérifié en cassant, dans les deux sens : remplacer l'énumération par
l'échantillonnage fait tomber le volet « monotonie » (40 pertes dès m = 2,
rang 3) ; rendre le tirage monotone fait tomber le contrôle négatif. Un test qui
n'exigerait que la monotonie passerait pour un générateur qui ne produit rien.

**3. Branchement dans `process_cicy`.** La branche a désormais sa propre boucle
— un fibré d'extension n'est pas un noyau de monade, et rien du chemin des
monades ne s'y applique. Filtres, du moins cher au plus cher :

1. χ(V) = χ(F1) + χ(F2), arithmétique pure, **exact** ;
2. Ext¹(F2, F1) ≠ 0, sinon seule l'extension scindée existe ;
3. Hoppe par borne supérieure sur les quotients gradués — **suffisant, jamais
   de faux positif** ; `indéterminé` n'est jamais inscrit ;
4. cohomologie par **bornes rigoureuses**.

Le point 4 méritait d'être écrit : `compute_extension_cohomology` pose les
morphismes de liaison « de rang maximal pour une extension générique », soit
exactement l'hypothèse qui avait faussé §4.3 et §4.4. `cohomology_extension_ex`
ne la fait pas. Avec h⁰(V) = h³(V) = 0 — conséquences de la stabilité, déjà
prouvée à l'étape 3 — la suite exacte longue donne

```
h¹(V) = h¹(F1) + h¹(F2) − r₁,     r₁ = rg( H¹(F2) → H²(F1) ) inconnu
h²(V) = h¹(V) + χ(V)
```

d'où un intervalle pour h¹, déterminé exactement quand h¹(F2) = 0 ou
h²(F1) = 0. `n_gen = |χ(V)|` reste **légitime** puisque la stabilité est
prouvée ; `cohomology` vaut `None` quand h¹ n'est pas déterminé, et non
[0,0,0,0] — un zéro de remplissage serait relu comme n_anti = 0, c'est-à-dire
comme un spectre propre, le piège du §4.8. Quand les quatre degrés sont
certifiés, χ par Riemann-Roch est confronté à la somme alternée des h^i.

**4. `--with-extensions` n'existe plus.** L'option est interceptée et le
programme s'arrête avec le motif : elle activait le chemin par pseudo-monade du
défaut 4.7. Un ancien script l'aurait rappelée en silence. Le nouveau nom est
`--extensions`.

**Trois effets de bord trouvés en branchant**, tous du même genre — un champ
absent relu comme un zéro :

- `deduplicate_results` indexait sur `(cicy, type, b_charges, c_charges)`.
  Les extensions n'ont pas de (B, C) : **2 647 candidats se repliaient sur 132**,
  silencieusement. La clé porte maintenant sur (F1, F2) quand (B, C) est vide.
  Même correction dans `audit_results.py` et `triage_clean.py`.
- `audit_results` ne vérifiait le rang que pour les monades. Le contrôle
  rk(V) = rk(F1) + rk(F2) est ajouté — c'est lui qui avait révélé le défaut 4.7.
- `triage_clean` lisait `higgs_fiable` de la seule absence d'avertissement
  `wedge2_heuristique`. Or H¹(∧²V) n'est pas calculé sur cette branche : le
  `higgs = 0` du résultat est un remplissage. Le champ `higgs_certifie` est
  désormais respecté.

**Premier scan de contrôle** (244 CICYs à m ≤ 3, `--max-charge 2`, 1,3 min) :
837 216 extensions énumérées, 30 021 à |χ(V)| = 3 (3,6 %), **2 647
Hoppe-stables** — contre 6 monades sur le même domaine. Toutes en mode
`exhaustif`, toutes de rang 3, donc E₆ ; h¹ déterminé dans 2 232 cas. Le taux de
passage est deux ordres de grandeur au-dessus de celui des monades, ce qui rend
la branche exploitable — et rend d'autant plus nécessaire un critère de tri
au-delà de Hoppe.

**Ce qui reste sur cette branche** : H¹(∧²V), donc le compte de Higgs et les
exotiques ; `verify_hoppe.py` ignore toujours les extensions
(`extension_ignoree`) ; et la chaîne d'équivariance (`equivariance.py`,
`equivariance_f.py`) ne lit que des monades — c'est le prolongement naturel,
puisque c'est elle qui a départagé les candidats du §2.


### 5.13 Pente — ce que le critère de Hoppe ne voit pas

**Une réserve de principe, qui n'avait pas été énoncée.** Le critère de Hoppe
s'écrit « c₁(V) = 0 ⟹ V stable ⟺ h⁰(∧^p V) = 0 pour p = 1..rk−1 ». Cette
**équivalence** suppose Pic(X) de rang 1. Sur une CICY à m > 1 facteurs, la
stabilité de pente dépend de la classe de Kähler J, et la condition reste
**nécessaire sans être suffisante** : elle est aveugle à la polarisation.

`stable: True` se lit donc « non éliminé par Hoppe », et non « stable ». Cela
vaut pour `hoppe_fast` et donc pour tout le catalogue de monades. La différence
avec les extensions est qu'une extension **exhibe** ses sous-faisceaux : dans
0 → F1 → V → F2 → 0, chaque sous-somme de F1 est un sous-faisceau de V, et la
préimage de toute sous-somme propre de F2 en est un autre. Comme μ(V) = 0, la
stabilité exige deg_J(W) < 0 pour chacun — de l'arithmétique pure sur les d_ijk.
Pour une monade, aucun sous-faisceau ne se présente à ce prix : **la réserve y
reste entière et non instruite.**

#### Le piège, et le chiffre qui n'existe pas

Première version : « aucun J de la grille [1,4]^m ne rend tous les degrés
négatifs ⟹ instable ». Elle annonçait **635 extensions déstabilisées sur
2 647, soit 24 %**. Ce chiffre est faux, et il ne mesurait que la grille :

| J_max | 3 | 6 | 12 | 24 |
|---|---|---|---|---|
| sans témoin | 1 748 | 1 299 | 1 042 | 925 |

Aucune saturation. C'est mot pour mot le faux lieu de base du §5.4, où la mesure
manquante était la dimension de la source ; ici c'est la saturation de la
recherche. **Un échec de recherche sur une grille finie ne démontre rien** et ne
doit jamais être inscrit comme une élimination.

#### Ce qui est démontré, et comment

L'élimination passe par un **certificat**, pas par une recherche. Avec
D_i(J) = Σ_jk d_ijk J_j J_k, on a deg_J(v) = Σ_i v_i·D_i(J). Les d_ijk d'une
CICY dans un produit d'espaces projectifs sont positifs ou nuls — **vérifié sur
les 7 890 entrées de `cicylist.txt`**, avec D_i(J) ≥ 0 et D_i(1,…,1) > 0
partout. Donc si des coefficients y_k ≥ 0 vérifient Σ_k y_k v_k ≥ 0 composante
par composante, alors

```
Σ_k y_k · deg_J(v_k) = (Σ_k y_k v_k) · D(J) ≥ 0
```

et l'un au moins des deg_J(v_k) est ≥ 0 — **pour toute classe de Kähler de
l'orthant**. Il suffit d'exhiber le y. C'est le sens facile du théorème de
transposition de Motzkin.

**Réserve** : la réciproque — l'absence d'un tel y prouverait l'existence d'un
p ≥ 0 avec v·p < 0 partout — demanderait un solveur de programmation linéaire,
et le dépôt ne dépend que de numpy. On n'explore donc que les y entiers à petits
coefficients. Mesure de saturation, elle : 105 certificats à taille ≤ 2, et
**pas un de plus** en montant à taille ≤ 3 et coefficients ≤ 3. La granularité
n'est pas le facteur limitant sur ce lot.

Trois issues, jamais confondues :

| verdict | sens |
|---|---|
| `False` | instable, **démontré** par certificat |
| `True` | un témoin J rend tous ces degrés négatifs — condition **nécessaire** satisfaite, pas une preuve de stabilité |
| `None` | ni certificat ni témoin. **N'élimine pas** |

#### Effet mesuré

Le test est branché **après le préfiltre χ et avant toute cohomologie** — il ne
coûte qu'un produit matriciel, les D_i(J) étant mis en cache par CICY (0,8 ms
par extension, contre plusieurs Koszul pour Hoppe).

Sur le scan de contrôle (244 CICYs à m ≤ 3, `--max-charge 2`) :

| | |
|---|---|
| extensions énumérées | 837 216 |
| passant \|χ(V)\| = 3 | 30 021 |
| **écartées par certificat d'instabilité** | **14 936 (50 %)** |
| retenues après Hoppe | 2 544 |

Et parmi celles que Hoppe seul laissait passer, **105 sont démontrées
instables** — 4,0 % du lot précédent, pas les 24 % annoncés par la version
fautive.

#### Le 27ᵉ test

Quatre volets, dont deux opposés et un qui fige précisément l'erreur ci-dessus :

| volet | référence |
|---|---|
| valeur connue d'avance | deg au point J = v, confronté à Riemann-Roch via `ChiCalculator` (24·χ = 4·cube + 2·c₂·v) |
| contrôle positif | un sous-faisceau à c₁ ≤ 0 doit trouver un témoin |
| contrôle négatif construit | un sous-faisceau à c₁ ≥ 0 doit être **certifié** instable, sur toute CICY |
| **non-élimination sur échec** | sur `#14`, un témoin existe (J trouvé hors d'une grille [1,2]³) : le verdict à budget insuffisant doit être `None`, jamais `False` |

Vérifié en cassant, de trois façons : convertir les `None` en `False` fait
tomber le quatrième volet ; renvoyer toujours `True` fait tomber le contrôle
négatif ; renvoyer toujours `False` fait tomber le contrôle positif.


### 5.14 Hoppe suffisant — la réserve du §5.13 est levée sur les deux candidats

Le §5.13 laissait le résultat principal en suspens : `#6890` et `#6947` vivent
dans P¹×P¹×P¹×P¹×P⁴, donc **h¹¹ = 5**, et l'équivalence de Hoppe suppose
Pic(X) de rang 1. Leur verdict ne valait que « non éliminé ».

**La forme suffisante.** V est μ_J-stable **dès que**

```
h⁰(∧^p V(−H)) = 0   pour p = 1..rk−1 et tout H avec deg_J(H) ≥ 0
```

Écrit ainsi l'ensemble des H paraît infini. Il ne l'est pas :

- ∧^pV(−H) ⊂ ∧^pB(−H), et dans le modèle S/I employé ici h⁰(O(a)) = 0 dès
  qu'une composante de a est négative. Donc **H ≤ une charge de ∧^pB** : borné
  au-dessus.
- deg_J(H) = Σ_k H_k·D_k(J) ≥ 0 avec tous les D_k(J) > 0 et H_k ≤ hi_k borne H
  **en dessous**.

Le polytope est compact, et petit : **110 twists pour `#6890`, 143 pour
`#6947`**, tous degrés p confondus. `h0_wedgep_V_sur_espace` gagne un paramètre
`twist` — la résolution est la même tordue par O(−H), et la contraction par f
est inchangée puisque f_i est de degré c − b_i, que décaler source et cible du
même −H ne modifie pas.

#### Résultat

| CICY | λ | twists testés | sources non vides | verdict |
|---|---|---|---|---|
| **6890** | +1 | 110 | **110 / 110** (dim max 78) | **stable à J = (1,1,1,1,1)** |
| 6890 | −1 | 110 | 110 / 110 | stable |
| **6947** | +1 | 143 | **143 / 143** (dim max 78) | **stable à J = (1,1,1,1,1)** |
| 6947 | −1 | 143 | 143 / 143 | stable |

Calcul fait sur le **sous-espace équivariant**, pas sur l'espace générique : h⁰
étant semi-continu supérieurement, une stabilité générique n'entraînerait pas
celle du f contraint. Zéro twist hors de portée.

**Le point qui rend ce vert défendable, c'est la colonne « sources non
vides ».** Un critère suffisant vérifié sur des sources vides serait vrai sans
rien démontrer — la faute des exotiques structurellement nuls du §4.8. Ici les
253 twists portent tous sur un calcul de rang réel, jusqu'à la dimension 78.

**Conclusion : `#6890` et `#6947` sont des fibrés stables**, au sens propre et
pour une classe de Kähler explicite, et non plus seulement « non éliminés ». Le
§2 tient. La distinction entre λ = +1 et λ = −1 reste portée par la seule
surjectivité (§5.4), la stabilité ne les départageant pas.

#### Ce que le twist n'apporte pas, mesuré

Sur **1 592 monades positives et 4 501 twists non nuls** échantillonnés sur la
liste, la marge `dim source − dim cible` vaut **au plus −2** — elle n'atteint
jamais 0. Un twist H ≠ 0 ne peut donc pas déstabiliser par comptage de
dimensions sur cette famille : seul H = 0 le peut. Autrement dit la réserve du
§5.13, réelle en principe, **ne mord pas sur le domaine des monades positives**.
C'est une mesure, pas un théorème — la déstabilisation par non-injectivité de
l'application reste possible, et c'est pourquoi les 253 twists ont bien été
calculés plutôt que supposés nuls.

#### Le 28ᵉ test

| volet | référence |
|---|---|
| valeur connue d'avance, négative | monade construite avec c − b₁ à composante négative : f₁ = 0, donc O(b₁) ⊂ V et h⁰(V) = dim H⁰(O(b₁)) = 3 exactement |
| **anti-vacuité** | les 110 twists de `#6890` doivent avoir une source **non vide** |
| **le twist doit agir** | `source_max` doit dépasser strictement la source à H = 0 (78 contre 64) |
| structure du polytope | H = 0 doit y figurer à tout p, sinon le critère ne contiendrait pas l'ancien ; un D_k nul doit rendre `None` |
| contenance | tout ce que Hoppe nu élimine, la forme suffisante doit l'éliminer |

Le troisième volet vient d'un cassage qui **passait** : neutraliser le paramètre
`twist` laissait le test vert, tous les h⁰ valant 0 de toute façon sur `#6890`.
La variation de la dimension de la source est la seule trace visible du fait que
le twist agit.


### 5.15 Le catalogue passé au crible — un faux positif, trois stabilités

Le §5.14 fermait la question sur les deux candidats et la laissait ouverte
partout ailleurs : `hoppe_fast` en restait à `max_H = 1`. Deux outils la
ferment, chacun dans un sens, et il faut les distinguer.

#### L'éliminateur : `hoppe_twists`, appliqué aux 115

Bon marché, il ne demande aucun f explicite. Pour chaque H du polytope
deg_J(H) ≥ 0, la borne `dim ker ≥ dim source − dim cible` est
inconditionnelle : une valeur strictement positive **prouve** h⁰(V(−H)) > 0,
donc un sous-faisceau de pente ≥ 0 dans un fibré de pente nulle.

**Chaque h⁰ doit être certifié.** Sans cela la borne porterait sur des nombres
faux dans ~30 % des cas (§4.2) et l'élimination ne serait pas démontrée.

Résultat sur `scan_wilson2` :

| | |
|---|---|
| **instable, démontré** | **1 — `#7484`** |
| aucun twist déstabilisant, tous certifiés | 53 |
| des twists non certifiés → aucun verdict | 61 |

**`#7484`** (rang 4, SO(10), cohomologie [0, 12, 0, 0]) tombe sur
H = (−2, 0, 1), de degré 4 :

```
h⁰(O(b_i − H)) = 6 + 0 + 3 + 4 + 0 = 13     tous certifiés
h⁰(O(c   − H)) = 12                          certifié
                        ⟹  h⁰(V(−H)) ≥ 1
```

Il figurait au catalogue comme Hoppe-stable. **Ni H = 0 ni H = e_i ne le
voient** : le témoin a des composantes de signes mélangés, hors de portée de
`max_H` quelle que soit sa valeur. La phase est branchée dans `hoppe_fast`
(paramètre `D`) et dans `process_cicy`.

#### Le critère exact, et où il s'arrête

Dans l'autre sens — prouver la stabilité — il faut f explicite et des rangs
exacts. `hoppe_suffisant_generique` monte l'anneau et la base pleine, puis
applique le §5.14. Couverture réelle sur les 115 :

| | |
|---|---|
| **prouvés stables à J = (1,…,1)** | **3 : `#6715`, `#6890`, `#6947`** (tous rang 4, SO(10)) |
| hors du domaine du modèle S/I | 44 |
| rang 5, hors budget | 68 |

`#6715` est un gain net : 221 twists, tous à source non vide. Les deux autres
confirment le §5.14, ici pour un f **générique** — l'énoncé équivariant du
§5.14 reste le plus fort des deux, h⁰ ne pouvant que monter en un point spécial.

**Le mur du rang 5 est mesuré, pas supposé.** Le polytope se calcule
instantanément (834 twists pour `#21`), mais chaque twist demande un rang sur
∧^p B pour p = 1..4, soit C(6,p) blocs : une seule entrée dépasse la dizaine de
minutes. Le coût est dans les rangs, pas dans l'énumération. Élargir demanderait
de mettre `dimY` en cache par degré, ou d'admettre un `maxdim` plus bas et donc
plus d'indéterminés.

#### Bilan de portée, à ce jour

| | statut |
|---|---|
| `#6890`, `#6947` | **stables**, sur le sous-espace équivariant, à J = (1,1,1,1,1) (§5.14) |
| `#6715` | **stable** pour un f générique, à J = (1,…,1) |
| `#7484` | **instable**, démontré — faux positif du catalogue |
| 53 entrées | non éliminées par les twists, tous certifiés — mais **non prouvées stables** |
| 61 entrées | twists non certifiés — aucun verdict, dans aucun sens |

#### Le 29ᵉ test

| volet | référence |
|---|---|
| valeur connue d'avance | les cinq h⁰ de `#7484` figés un par un (6, 0, 3, 4, 0 et 12), et la borne 13 − 12 = 1 |
| deux verdicts opposés | `#7484` doit tomber, `#6890` doit survivre |
| mesure du gain | `hoppe_fast` sans `D` ne doit **pas** éliminer `#7484`, sinon le cas ne mesure plus rien |
| garde de certification | `#21` a 5 twists non certifiés sur 45 : le verdict doit être `None`, jamais `False` |

Vérifié en cassant, trois fois. Le dernier volet vient d'un cassage qui
**passait** : ignorer la certification ne change rien sur `#7484`, dont tous les
h⁰ sont certifiés. Il fallait une entrée où la garde mord réellement.


### 5.16 Le balayage du §2, rejoué avec les filtres actuels

Le balayage du §2 datait d'avant h³ (§5.5), d'avant la phase des twists (§5.15)
et d'avant tout le travail sur la pente. `#7484` ayant prouvé que le catalogue
contenait au moins un faux positif, il fallait savoir si les survivants
tenaient. Le catalogue a d'abord été purgé par `hoppe_twists` — **115 → 114**,
`#7484` retiré — puis la chaîne entière relancée.

| | rejoué | §2 |
|---|---|---|
| passent le test nécessaire sur les charges | 108 | 108 |
| couples évalués | 4 076 | 3 878 |
| tués par h⁰(V) équivariant | 3 624 | 3 452 |
| sans f équivariant | 36 | — |
| indéterminés | 449 | 423 |
| **survivants** | **3 couples, 2 candidats** | **3 couples, 2 candidats** |

**Les trois couples sont les mêmes** : `#6890` (deux entrées du catalogue) et
`#6947`, tous à λ = +1, tous avec le certificat de surjectivité — en (0,1,5,1,0)
et (5,1,0,1,0) respectivement, aux degrés mêmes du §5.4.

**Le §2 tient.** Il repose maintenant sur un balayage qui inclut h³, la phase
des twists et un catalogue purgé, et sur deux fibrés dont la stabilité est
démontrée et non plus seulement non éliminée (§5.14).

**Les 449 indéterminés sont tous du même motif** — `surjectivité de f non
certifiée` — et à **420 sur 449 de rang 5**. C'est exactement le mur du §5.4 :
au rang 5 le certificat J_d = R_d n'est pas atteignable. Ces 449 lignes ne
disent rien, ni dans un sens ni dans l'autre, et c'est le plus gros bloc non
instruit du projet.


### 5.17 Le mur du rang 5 n'était pas géométrique — il était dans la liste des degrés

Les 449 indéterminés du §5.16 ont **un seul motif**, `surjectivité de f non
certifiée`, dont 420 au rang 5. Le §5.4 en concluait que le critère J_d = R_d
est « hors de portée au rang 5 ». **Ce constat mesurait la liste des multidegrés
essayés, pas la géométrie.**

`_degres_a_essayer` engendrait trois familles : marches **longues** sur un seul
axe (c + t·e_k), pas **courts** mixtes (|v|₁ ≤ 2), et croissance uniforme
c + t·(1,…,1). Elle ne combinait **jamais** longueur et mixité. Puis elle triait
par coût croissant et gardait les `n_degres` **premiers** — exactement le mauvais
sens, la marge s'améliorant avec la taille.

Balayage en directions mixtes longues sur `#21` (rang 5, m = 5) :

| d | source | cible | marge |
|---|---|---|---|
| [4, 4, 2, 5, 4] | 2 670 | 1 278 | **+1 392** |
| [2, 4, 4, 2, 7] | 2 446 | 1 074 | **+1 372** |
| [2, 2, 8, 5, 4] | 2 642 | 1 398 | **+1 244** |

Des cibles bien sous le plafond de 6 000. Sur ce même `#21`, les quatre degrés
retenus par l'ancienne sélection ont des marges de **−42, −43, −45, −46** : le
critère ne pouvait aboutir sur aucun, quel que soit f.

**Deux corrections.** Montée anisotrope (c + t·(1,…,1) + s·e_k, plus une
recherche locale gloutonne sur la marge prédite, qui trouve des optima à support
plein sans énumérer {0..q}^m). Et sélection **par viabilité** : les degrés à
marge prédite ≥ 0 d'abord, les moins chers parmi eux. La marge est estimée par
`dim_multi` — formule fermée dans l'ambiant — parce que `dimY` construirait le
quotient, soit l'opération qu'on cherche à éviter ; le test exact reste dans la
boucle.

#### Un bug trouvé en écrivant le test

Le filtre s'écrivait `(_marge_predite(...) or -1) >= 0`. **Une marge exactement
nulle est falsy en Python**, donc traitée comme −1 et écartée. Or les deux
certificats du §5.4 sont précisément à marge nulle — source = cible = 24 sur
`#6890` comme sur `#6947`. Le filtre écartait exactement les degrés qui
certifient les deux seuls candidats du projet. Corrigé, et figé par le volet (d)
du 30ᵉ test, qui exige qu'une marge nulle soit conservée.

#### Effet mesuré, et ce qu'il ne dit pas

Sur les douze premières minutes du rejeu : **12 lignes `SURVIT` contre 3 dans le
run complet précédent**, toutes sur `#21`, SU(5) de rang 5, Γ = ℤ₂, λ = +1,
certifiées surjectives en (1,3,1,5,3) et (1,3,3,1,5). Le régime que le §5.4
déclarait inatteignable rend désormais des verdicts.

**Ces survivants ne sont pas des modèles à trois générations, et il faut le dire
tout de suite.** `#21` a |χ(V)| ∈ {24, 48} ; avec |Γ| = 2 cela donne **12 ou 24
générations**, pas 3. Ce sont des fibrés stables, équivariants et surjectifs, au
mauvais indice pour ce groupe. `#21` porte aussi des groupes d'ordre 4, 8 et 16 —
et |χ| = 24 avec |Γ| = 8 donnerait 3 générations. C'est cette combinaison-là qu'il
faut regarder, et le rejeu ne l'a pas encore atteinte.

**Coût.** Les degrés viables sont bien plus gros que les degrés simplement bon
marché : le rejeu avance environ dix fois plus lentement (141 lignes en 32 min,
contre 4 130 en 95 min). Le balayage complet demandera plusieurs heures. C'est le
prix d'un critère qui peut aboutir, contre un qui ne le pouvait pas.


### 5.18 Le filtre d'indice retombait sur tous les groupes sans le dire

`equivariance_f.py` limitait les symétries testées à celles dont l'ordre est
compatible avec l'indice — sauf que :

```python
groupes = set(r.get('groupes_utiles') or [])
if groupes is not None and not groupes:
    groupes = set(r.get('equivariant_possible') or [])   # repli SILENCIEUX
```

Faute de groupe d'ordre compatible, il les essayait **tous**. Mesure sur le
balayage du §5.16 :

| couples évalués | 4 076 |
|---|---|
| indice compatible (\|χ\| = 3·\|Γ\|) | **184** |
| indice **incompatible** | **3 892 — 95,5 %** |
| candidats distincts avec au moins un groupe compatible | **12**, sur 8 CICYs |

Ces 3 892 couples ne peuvent donner trois générations avec ce Γ, quel que soit
le verdict. Le calcul les traitait quand même — et certains ressortaient
étiquetés `SURVIT`. C'est ainsi que le §5.17 a d'abord annoncé des survivants
SU(5) de rang 5 sur `#21` : ils sont stables, équivariants et surjectifs, et ils
donnent **12 générations** (|χ| = 24 avec |Γ| = 2). Un filtre qui devient vide
sans le dire, c'est le défaut du §4.8 sous une autre forme.

**Deux corrections.** Le repli est supprimé : le candidat est écarté et **la
raison persistée** dans le JSONL — un fichier de résultats doit dire pourquoi un
cas n'a pas été traité (§5.11). Et le verdict porte désormais
n_gen(X/Γ) = |χ(V)|/|Γ| : `SURVIT — 3 gen sur X/Gamma`. Sans ce nombre, un
`SURVIT` ne dit rien du contenu physique.

#### Balayage corrigé

| | §5.16 | corrigé |
|---|---|---|
| lignes | 4 130 | 298 |
| couples évalués | 4 076 | **184** |
| écartés faute de groupe compatible | — (repli silencieux) | **73, tracés** |
| indéterminés | 449 | **21** |
| **survivants** | 3 | **3** |
| durée | ~95 min | **~4 min** |

Les trois survivants sont les mêmes, et le verdict le dit maintenant en toutes
lettres : **`#6890` (deux entrées) et `#6947`, Γ = ℤ₂, λ = +1, 3 générations sur
X/Γ**, certifiés surjectifs en (0,1,5,1,0) et (5,1,0,1,0).

Les 449 indéterminés du §5.16 tombent à **21** — l'effet conjugué du §5.17 (les
multidegrés viables) et de la suppression du bruit. Ils portent sur `#5259` (16
lignes), `#2565`, et les λ non retenus de `#6890` et `#6947`. **160 couples sont
éliminés franchement**, dont les 144 de `#7884` (E₆, ℤ₃×ℤ₃).

Le §2 tient, pour la troisième fois de suite et sur une chaîne à chaque fois plus
exigeante.


### 5.19 Trois constantes déguisées en résultats

Trois nombres du classement n'en étaient pas. Ils ne faussaient pas la
sélection — les candidats retenus le sont sur d'autres critères — mais ils
faussaient le **tri**, et deux d'entre eux rapportaient des points.

**Les exotiques SU(5)** : `max(0, n_10 + n_10bar − n_gen − 2·n_anti)`, avec
n_gen = |a−b| et n_anti = min(a,b). Or |a−b| + 2·min(a,b) = a+b pour tous
a, b ≥ 0 : l'expression vaut **identiquement zéro**. Vérifié sur 144 couples
dans le test, pour que l'énoncé ne repose pas sur la seule algèbre.

**Les exotiques SO(10)** : codés en dur à 0.

**Les singlets** : lus dans `end_V`, qui était une valeur de remplissage —
rank_V² − 1 — puisque h¹(End V) n'est pas calculé.

Les deux premiers valaient **25 points** de score à tout SO(10) et tout SU(5),
le troisième jusqu'à 10. Tous trois valent maintenant `None`, et
`compute_sm_compatibility` n'accorde de points qu'à une quantité réellement
calculée. **E₆ conserve son compte d'exotiques** — n_anti y est effectivement
calculé, c'est le seul cas — et c'est le verdict opposé qui rend le test
discriminant.

**Un quatrième, sur les Higgs E₆** : `max(0, n_gen − 3) + n_anti`, avec un 3
**codé en dur**. En mode Wilson, n_gen est le compte **en amont** du quotient
— 6, 9, 27… — et le 3 est le compte **voulu en aval** : la soustraction mélange
deux étages. Avec n_gen = 6 et n_anti = 0, elle fabriquait 3 Higgs à partir de
rien. Les Higgs d'un E₆ viennent des paires 27 + 27̄, donc de n_anti seul ; avec
ligne de Wilson ils sortent de la décomposition des 27 sous Γ, qui n'est pas
calculée ici (§6).

Effet sur le scan de contrôle : les scores des E₆ passent de 92,5 à 87,5 — les
10 points de singlets inventés disparaissent, les 25 points d'exotiques restent
car ils y sont mérités.


### 5.20 Situation dans la littérature — les deux candidats sont hors de la classe balayée

Anderson, Gray, He et Lukas ont publié en novembre 2009 (`arXiv:0911.1569`,
JHEP 02(2010)054) « une analyse complète » de **7 118 monades positives** sur
les **4 515 CICYs favorables**, ranges 3/4/5, avec les mêmes groupes E₆ / SO(10)
/ SU(5), et concluent que **toute la classe est écartée sur des bases
phénoménologiques**. C'est le seul balayage systématique du même terrain. Il
fallait savoir où `#6890` et `#6947` s'y situent.

**Ils n'y sont pas, et pour une raison de définition.** Leur équation (2.11)
demande

```
b_i^r > 0   et   c_a^r > 0   pour TOUT r, i, a
```

— une positivité **stricte**. Nos deux candidats ont **20 zéros dans B et 2 dans
C** chacun. Ils sont *semi-positifs*, hors des 7 118. Le générateur du dépôt
n'impose que leur condition (2.6), c_a ≥ b_i composante par composante, qui est
beaucoup plus large.

Ce n'est pas une échappatoire technique : les auteurs désignent eux-mêmes cette
classe comme la frontière intéressante.

> « la condition de positivité, quoique favorable à la stabilité, n'est **pas
> nécessaire** pour X avec h¹¹(X) > 1. Cela signifie que des monades
> semi-positives ou même “légèrement négatives” peuvent être stables, et que les
> monades positives sont vraisemblablement **un petit sous-ensemble** de tous les
> fibrés monades stables. »

Leur propre nouveau modèle standard est d'ailleurs bâti sur une monade
**semi-positive**, sur le bicubique.

**Leur couverture en quotients est de cinq variétés.** Le critère χ(V) ∈ 3·S(X)
appliqué au petit ensemble positif ne laisse que 91 modèles, sur **cinq CICYs
seulement** : la quintique, [P⁵|3 3], le tétraquadrique [P⁷|2 2 2 2], le
bicubique, et le tétra-quadrique (P¹)⁴. Le dépôt travaille sur les **194 CICYs**
de la classification de Braun — publiée en mars 2010, soit **quatre mois après**
leur article.

**Leurs trois modèles SO(10) ne sont sur aucune de nos CICYs.** Un sur la
quintique, `0 → V → O(2)³ ⊕ O(1)⁴ → O(4) ⊕ O(3)² → 0`, écarté parce que
c₂(V) = −45 n'est pas divisible par 25 ; deux sur [P⁵|3 3], exigeant des groupes
d'ordre 18 et 12 que la variété ne réalise pas.

**Un point de méthode qu'ils laissent explicitement ouvert.** Leur note de bas
de page 2 dit que la condition (2.6) rend C*⊗B globalement engendré, donc V un
fibré par un théorème de Fulton–Lazarsfeld — puis :

> « Il est bien possible de relaxer la condition (2.6) et d'obtenir encore un
> fibré. Mais cela demande une **analyse au cas par cas que nous ne
> considérerons pas** dans le présent article. »

Or ce théorème porte sur un f **générique**. Le f équivariant est un point
spécial, et rien ne garantit qu'il évite le lieu où f chute de rang — c'est
exactement la réserve du §5.4. Le certificat J_d = R_d est cette analyse au cas
par cas, faite sur le f contraint.

**Ce que cela établit, et ce que cela n'établit pas.** Les deux candidats
tombent hors du seul balayage systématique de ce terrain, dans la classe que ses
auteurs désignent comme prometteuse et non traitée, sur des variétés que leur
couverture n'atteignait pas. Cela **n'établit pas** que personne ne l'a fait
depuis 2009 : une recherche n'a rien trouvé de systématique — le travail récent
le plus proche, avec apprentissage par renforcement (2021), ne porte que sur
**deux variétés** de nombre de Picard 2 et 3 — mais une absence de résultat de
recherche n'est pas une preuve d'absence.

**Une condition physique que le dépôt ne teste pas.** L'annulation d'anomalie
exige que c₂(TX) − c₂(V) soit une classe effective, leur équation (2.9). Le
pipeline ne la vérifie **nulle part**. Contrôle a posteriori sur les deux
candidats :

| | c₂(TX) − c₂(V) | effectif |
|---|---|---|
| **6890** | (10, 18, 22, 18, 28) | oui |
| **6947** | (22, 18, 10, 18, 36) | oui |

Ils passent — mais par chance, pas par construction. **Les autres candidats du
catalogue n'ont pas été contrôlés.**


### 5.21 Annulation d'anomalie — 61 % du catalogue n'était pas physique

Le §5.20 a mis au jour une condition que le pipeline ne testait **nulle part** :
pour préserver la supersymétrie, la classe duale à c₂(TX) − c₂(V) doit être
**effective**, ce qui sur une CICY favorable se lit composante par composante.
C'est l'équation (2.9) de `arXiv:0911.1569`, et c'est une contrainte de
cohérence de la théorie — pas un raffinement. **Un fibré qui la viole n'est pas
un modèle**, quelles que soient sa stabilité, sa cohomologie et son nombre de
générations.

**L'ampleur du trou :**

| catalogue `scan_wilson2` | 115 entrées |
|---|---|
| **violant l'annulation d'anomalie** | **70 — 60,9 %** |
| la satisfaisant | 45 |

Les cas les plus nets sont les SU(5) de rang 5 sur `#21` : déficits à
(−36, 14, −6, 16, 12) ou (−60, −22, 12, 8, 6). Ce sont précisément les entrées
que le §5.17 avait fait ressortir en `SURVIT` — elles étaient déjà à 12
générations (§5.18), elles ne sont même pas des théories cohérentes.

**Le calcul.** c(V) = c(B)/c(C) sur une monade, avec c₁(B) = c₁(C) puisque
c₁(V) = 0, donc c₂(V) = c₂(B) − c₂(C) ; en développant, on retombe exactement
sur la formule (2.9). Pour une extension, la multiplicativité donne
c(V) = c(F1)·c(F2), donc V a la classe de F1 ⊕ F2. C'est de l'arithmétique pure
sur les d_ijk, placée **avec le préfiltre χ**, avant toute cohomologie.

**Le §2 tient.** Les deux candidats passent, avec des déficits confortables :

| | c₂(TX) − c₂(V) |
|---|---|
| **6890** | (10, 18, 22, 18, 28) |
| **6947** | (22, 18, 10, 18, 36) |

Chaîne complète rejouée sur le catalogue assaini — purgé par la phase des twists
puis par l'anomalie, **114 → 44** :

| | avant | catalogue physique |
|---|---|---|
| candidats entrant dans la chaîne | 114 | **44** |
| passant le test nécessaire sur les charges | 108 | 38 |
| couples évalués | 184 | 184 |
| indéterminés | 21 | 21 |
| **survivants** | **3** | **3** |

Les mêmes : `#6890` (deux entrées) et `#6947`, Γ = ℤ₂, λ = +1, **3 générations**.
Le résultat principal repose désormais sur un catalogue dont **chaque entrée est
une théorie cohérente**, ce qui n'était pas le cas.

**Portée du chiffre de 61 %.** Il ne signale pas une erreur de calcul mais une
condition absente : le générateur produit des monades qui satisfont c₁(V) = 0 et
la cible d'indice, sans jamais regarder c₂. Toutes les statistiques du document
antérieures à cette section décrivent un catalogue dont trois entrées sur cinq
n'étaient pas des modèles.


### 5.22 Où la vérification d'anomalie est-elle réellement en place ?

Question posée après le §5.21, et la réponse était **« à un seul endroit »**.

| | anomalie testée |
|---|---|
| `main_optimized.py` (scan) | oui, deux appels — monades et extensions |
| `tests_regression.py` | oui, 33ᵉ test |
| **`audit_results.py`** | **non** |
| `triage_clean.py`, `verify_hoppe.py` | non |
| les sept autres `main_*.py` | non |

Le trou qui comptait est `audit_results.py`. C'est l'outil dont le rôle est
précisément de rattraper a posteriori ce qu'un scan a laissé passer — sa
docstring énumère les tests qu'il applique. Passé sur `scan_wilson2`, il
annonçait **« 115 retenus, 0 écarté »** alors que 70 de ces entrées ne sont pas
des modèles. La purge faite au §5.21 était un script jetable, hors du dépôt :
rien ne permettait de la reproduire.

**Corrigé.** Le drapeau `anomalie` est ajouté, avec le vecteur c₂(TX) − c₂(V)
dans le champ `anomalie_deficit` — une composante négative identifie la
direction fautive. Sur `scan_wilson2` :

```
  monad        anomalie=1
  pos_monad    anomalie=69
  Candidats retenus : 45 / 115  (39,1 %)
```

Exactement le chiffre du §5.21, désormais reproductible d'une commande. Et si
`cicylist.txt` est absent, le script **le dit** au lieu d'omettre le test en
silence — un catalogue pourrait sinon ressortir « propre » sans avoir été
contrôlé.

**Les sept points d'entrée hérités refusent de tourner.** Ils n'ont ni
l'annulation d'anomalie, ni la phase des twists, ni les correctifs du §4.
Les laisser exécutables sans le dire offrait un chemin produisant des résultats
d'apparence normale et faux. Ils sont conservés pour l'historique, avec un refus
explicite qui renvoie vers `main_optimized`.

---

### 5.23 Le résultat principal du projet avait été trouvé par une loterie à dix tirages

C'est le défaut le plus grave trouvé jusqu'ici, et il ne se voit dans aucun
résultat : il se voit dans une **absence**.

**Le symptôme.** `scan_wilson3` (194 CICYs, 152,8 min, 49 entrées) ne contient
**ni #6890, ni #6947, ni #6715** — les trois candidats sur lesquels reposent
le §2 et une bonne part du §5. Vérification faite un filtre après l'autre :
non-dégénérescence, |χ| = 6, annulation d'anomalie (déficits (10, 18, 22, 18, 28)
et (22, 18, 10, 18, 36), donc effectifs), `check_map_exists`, `hoppe_fast` avec
et sans D. **Ils passent tout.** Ils n'ont pas été éliminés : ils n'ont pas été
**engendrés**.

**La cause.** Le catalogue les étiquette `type: monad`, donc générateur
*classique*. Dans `generate_monads`, ils sortaient du bloc « monades
anti-symétriques » :

```python
for _ in range(min(m * 3, 10)):        # DIX tirages
    ...
    i1 = rng.randint(0, m)             # sur le RNG PARTAGE
```

Deux défauts en une ligne.

1. **Couverture.** La famille visée — des `b_i` de la forme e_a, ou e_a ± e_b —
   compte, pour m = 5 et r_B = 5, **101 multiensembles** dans sa strate pure et
   **2 100** en autorisant un vecteur perturbé. Dix tirages en voyaient dix.
2. **Reproductibilité.** Le RNG étant partagé avec `generate_positive_monads`,
   la valeur des dix tirages dépendait de **ce que le générateur positif avait
   consommé avant**. La correction du §5.11 — qui ne touchait que le générateur
   positif — a redistribué la loterie. Les candidats sont tombés du mauvais côté.

Mesure de la loterie, sur 5 000 scans simulés : le tirage à dix trouvait **au
moins un** des trois candidats dans **6 scans**, jamais deux, jamais trois.
Une chance sur mille.

> **Le résultat principal du projet n'était pas reproductible.** Il n'était pas
> faux — les trois candidats sont bien là, et passent tout — mais sa présence
> dans un scan tenait au tirage. C'est exactement la leçon du §5.11, jamais
> appliquée à ce générateur-ci.

**La correction.** La famille est **énumérée**, par strates selon le nombre de
vecteurs perturbés, parce que la famille complète ne l'est pas (pour m = 5,
r_B = 6, elle dépasse déjà 2 000 000 de multiensembles à |c₁(B)|∞ ≤ 3) :

| strate | contenu | plafond | état sur les 194 CICYs |
|---|---|---|---|
| **k = 0** | tous les `b_i` sont des e_a | `--unite-max` (200 000) | **exhaustive partout** (m ≤ 10) |
| **k = 1** | un `b_i` de la forme e_a ± e_b | `--unite-perturbe-max` (20 000) | exhaustive jusqu'à m = 6, soit **158 des 194** |
| **k ≥ 2** | — | — | jamais énumérée, **déclarée** comme telle |

#6890 = {e₁, e₂, e₂, e₂, e₃} et #6947 = {e₀, e₀, e₀, e₁, e₃} sont dans la strate
pure ; #6715 = {e₃, e₃, e₀, e₀+e₂, e₃} a un vecteur perturbé et vient de k = 1.
Les trois sont désormais **démontrés présents**, pas trouvés.

`generate_monads` reçoit en outre un **RNG dérivé** de (seed, m, max_charge,
rank_V, r_C), comme `generate_positive_monads` depuis le §5.11. Une modification
de l'un ne peut plus déplacer ce que produit l'autre.

**Coût.** Le générateur rend de 900 (m = 1) à ~47 000 (m = 6) monades par CICY,
contre ~110 avant. Vérifié de bout en bout sur les trois CICYs : 5,3 min
cumulées à trois cœurs, 367 fibres Hoppe-stables, **les trois candidats
retrouvés**. Sur 194 CICYs, compter environ une heure à huit cœurs.

**Ce que chaque résultat déclare maintenant.** Champ `unite_strates` :
`{'k0': [mode, produit, total], 'k1': [...], 'k2+': ['non_couvert', 0, None]}`.
Le générateur est un filtre comme un autre, et c'est celui qui a fait disparaître
les candidats sans qu'une ligne le signale. Sans ce champ, une absence dans un
scan reste ininterprétable — règle des filtres, §8.

**`max_charge = 3` est figé** pour ce générateur et ne suit pas `--max-charge`
(qui n'atteint que le générateur positif). C'est la borne sous laquelle la
famille est énumérée et sous laquelle les candidats ont été trouvés ; la rendre
variable changerait la famille énumérée sans que rien ne le déclare.

**Le test (34ᵉ) porte deux verdicts opposés**, et c'est le second qui lui donne
sa valeur : l'énumération contient les trois candidats, **et** l'ancien tirage à
dix ne les contient pas (2 scans sur 2 000 en trouvaient un, aucun les trois).
Sans le second, le test passerait aussi bien avec un générateur qui les aurait
trouvés par hasard. Il vérifie de plus l'indépendance au RNG amont, l'exactitude
du comptage préalable qui décide du plafond (36 cas), et qu'un plafond atteint
se **déclare**. Les trois sabotages essayés — strate k = 1 désactivée, RNG
partagé restauré, comptage faussé de +1 — le font échouer.

---

### 5.24 Cinquante heures de calcul sans un octet sur le disque

Conséquence directe du §5.23, et découverte en la subissant.

Le générateur énuméré fait passer le catalogue de **115 à 19 579 fibrés
Hoppe-stables**, et le lot d'entrée d'`equivariance_f.py` de 108 à **14 945
candidats**. L'étape qui prenait une heure en prend cinquante-cinq. Or :

```python
sortie = []
for r in rs:
    ...
    sortie.append({**ident, **L})

with open(dst, 'w', encoding='utf-8') as fh:      # <- une seule fois, à la fin
```

Tout était accumulé en mémoire et écrit après le **dernier** candidat. Ni
checkpoint, ni gestionnaire de Ctrl+C, ni reprise — ni dans `equivariance_f.py`,
ni dans `equivariance.py`. Mesure réelle sur `scan_wilson4` : **2 h 20 de calcul,
zéro octet récupérable**. Le fichier `results_equivariance_f.jsonl` n'existait
pas.

Ce n'était pas un défaut tant que le lot faisait 108 candidats. Il l'est devenu
parce qu'un autre correctif a multiplié le lot par 138 — et rien dans le code ne
reliait les deux.

**La correction** reprend le mécanisme de `main_optimized` : JSONL en écriture
*append-only*, fichier de progression écrit **après chaque candidat**.

Un compteur ne suffit pas. Si le processus meurt **au milieu** d'un candidat,
quelques-unes de ses lignes sont déjà sur le disque : reprendre au candidat
suivant les laisserait en double. On enregistre donc aussi l'**offset** du
fichier après le dernier candidat *complet*, et la reprise **tronque** à cet
offset. Ce qui est relu est alors exactement ce qui a été validé.

Trois gardes refusent une reprise incohérente, et **disent pourquoi** :

| situation | conséquence si acceptée |
|---|---|
| l'empreinte sha256 de l'entrée a changé | verdicts attribués aux **mauvais** candidats |
| `--cicy` diffère du checkpoint | idem — le filtre change la numérotation |
| le JSONL est plus court que l'offset | reprise dans le vide |

Repartir de zéro **en silence** ferait passer un recommencement complet pour une
reprise : le motif est affiché.

**Le test (35ᵉ)** vérifie les deux moitiés. Fidélité : un lot traité en morceaux,
avec coupures, produit un JSONL identique — ligne pour ligne, clé pour clé — à
celui du même lot traité d'un trait. Refus : un checkpoint qui ne correspond plus
est rejeté avec son motif.

Un point mérite d'être noté, parce qu'il a failli passer. La première version du
test attendait que le hasard du minutage produise une coupure *au milieu* d'un
candidat — et **retirer la troncature ne le faisait pas échouer** : la coupure
tombe presque toujours *entre* deux candidats, où il n'y a rien à tronquer. Le
test écrit maintenant lui-même des lignes au-delà de l'offset valide et exige
qu'elles disparaissent. C'est le cinquième cassage « évident » de ce dépôt à
passer au premier essai (§8).

Quatre sabotages le font échouer : troncature retirée, garde d'empreinte
neutralisée, garde qui refuse *tout* (elle est alors prise par « le lot n'avance
pas entre deux coupures »), et ouverture du JSONL en `'w'` au lieu de `'a'`.

**Un cinquième piège, découvert sur la machine de Franck et pas sur la mienne.**
Le test comparait les lignes **en ordre de fichier**. Or depuis la version
parallèle, `imap_unordered` rend les lots dans l'ordre où ils finissent — qui
dépend du nombre de cœurs. Sur 8 cœurs le test échouait avec
« reprise infidèle : **30 lignes contre 30** », deux runs ayant produit
exactement les mêmes lignes dans un ordre différent. La comparaison se fait
désormais sur les lignes **triées** : l'ordre n'est pas dans le contrat, le
multiensemble l'est. Le tri ne relâche rien — une ligne en trop, en moins ou
différente reste détectée, et les quatre sabotages le confirment.

**Et une troisième fois, dans la foulée.** Le volet « JSONL amputé » coupait le
fichier *à la moitié* après avoir exigé plus de quatre lignes — ce qui suppose
les lots gras. Sur 8 cœurs, `imap_unordered` rend d'abord les lots « aucun
groupe compatible », qui n'écrivent qu'**une** ligne : quatre lots, quatre
lignes, et le test échouait sur un message réduit à `4`. Le point de coupure est
maintenant lu **dans le fichier lui-même** — chaque ligne porte son `_lot`, donc
on sait quel lot a été écrit en dernier et combien de lignes il compte — et si ce
lot n'en a qu'une, le volet se déclare **non exercé** au lieu de couper au hasard
en croyant l'avoir fait.

Trois fois, donc, que ce test s'est révélé dépendant de l'environnement : le
chronomètre, l'ordre des lignes, la grosseur des lots. Chaque fois, la même
faute de ma part — écrire une vérification qui suppose ce qu'elle observe sur
*ma* machine. Un test dont le verdict dépend de la machine ne vaut rien, et
celui-ci a demandé trois corrections pour n'en plus dépendre.

---

### 5.25 Replier les orbites : facteur 4, et le contrôle qui empêche que ce soit un §5.23 de plus

Le générateur énuméré produit la famille **complète** des sommes de vecteurs
unité. Sur une CICY dont la matrice de configuration est symétrique, cette
famille contient les images les unes des autres. Les **12 monades survivantes de
#6947** sont les 12 arrangements d'un même motif de multiplicités (3, 1, 1, 0)
sur ses quatre facteurs P¹ — et son groupe d'automorphismes est d'ordre **24**,
de sorte que ces 12 forment **une seule orbite**. Douze fois le même calcul.

`Aut(config)` = les permutations des facteurs projectifs qui préservent les
dimensions et redonnent le même **multiensemble** de lignes — les équations ne
sont pas numérotées.

| CICY | \|Aut\| | ce que cela change |
|---|---|---|
| 7300 | 72 | 762 candidats → 24 orbites |
| 5302 | 48 | 564 → 20 |
| **6947** | **24** | 176 → 12 |
| 5 | 12 | 334 → ~28 |
| 6715 | 6 | 158 → 32 |
| **6890** | **1** | **aucun repli possible** |

Sur les 194 CICYs : **14 945 candidats → 3 688 orbites, facteur 4,05.** Soit
environ **14 h au lieu de 55**. Mesuré de bout en bout sur #6947 : 62,8 s sans
repli, 13,3 s avec (contrôle compris), **verdicts identiques sur les 176
candidats**, mêmes 12 survivants.

**Ce qui n'est PAS démontré.** Qu'une permutation préserve la matrice de
configuration n'implique **pas** qu'elle commute avec l'action de Γ lue chez
Braun, laquelle est attachée à des coordonnées précises. L'égalité des verdicts
sur une orbite est une **hypothèse**, vérifiée sur 919 lignes réelles (#6890,
#6947, #6715 : 0 discordance sur 42 orbites non triviales) et non prouvée.

C'est exactement la situation du §5.23 — un filtre qui peut faire disparaître des
candidats — sauf qu'ici on la choisit. D'où trois précautions :

1. **Le repli est optionnel.** Par défaut chaque candidat est évalué.
2. **Aucune ligne ne disparaît.** Le verdict du représentant est recopié sur
   chaque membre avec `verdict_replique: True` et l'identité du représentant. Le
   JSONL compte autant de lignes qu'un balayage complet ; l'aval ne voit aucune
   différence, sinon un champ qui dit d'où vient le verdict.
3. **Le repli se contrôle lui-même.** `--controle-orbites N` (défaut 20) évalue
   **pour de vrai** N membres non représentants et compare. Une discordance est
   affichée, comptée, et le run déclare le repli invalide.

**Le point qui a failli tout vider de sens.** Le contrôle tirait d'abord *N
orbites*, un membre chacune. J'ai saboté `canonique` pour qu'elle range **tous**
les candidats d'une CICY dans une seule orbite — un repli entièrement faux. Le
run a affiché : *176 candidats → 1 tâche, 1 contrôle, 0 discordance.* Le contrôle
validait le sabotage, parce qu'une orbite unique ne reçoit qu'un tirage. Il tire
maintenant **N couples (orbite, membre)** : une orbite géante reçoit une part des
contrôles proportionnelle à sa taille, et le sabotage produit six discordances.

**Le test (36ᵉ)** vérifie que `Aut(config)` est un groupe (identité, stabilité par
composition et inverse, sur 8 CICYs, avec des ordres allant de 1 à 72 — un ordre
constant trahirait une constante déguisée), que les verdicts sont invariants sur
les orbites des sorties réelles, que le JSONL replié a exactement les mêmes
candidats et verdicts que le complet, et que le repli abusif **déclenche** les
discordances.

---

### 5.26 Sept cœurs inutilisés, et une unité de travail mal choisie

Découvert en regardant le Gestionnaire des tâches pendant un balayage : **12 % de
CPU sur une machine à 8 cœurs**, c'est-à-dire *un seul cœur*. `main_optimized`
distribue son travail sur un `Pool` depuis toujours ; `equivariance_f.py` ne
l'avait jamais fait, parce qu'il tournait en une heure sur 108 candidats.

**Et l'unité de travail n'était pas la bonne.** Une tâche, ce n'est pas un
calcul : c'est un candidat confronté à **toutes** les réalisations que Braun
donne de son groupe.

```
#480  : 394 symetries -> {'Z2': 21, 'Z4': 5, 'Z2 x Z2': 368}
#6947 :   9 symetries -> {'Z2': 1, 'Z4': 1, 'Z2 x Z2': 1, ...}
```

**368 réalisations pour un seul candidat**, à ~24 s chacune : la tâche 823 dure
plus de deux heures, pendant lesquelles le script n'affiche rien et ne sauvegarde
rien — `analyser` ne contient aucun `print` et ne rend la main qu'à la fin. Sur
les 3 698 tâches, il y a **60 201 couples (tâche, réalisation)**, dont **19 461
pour #480** et 9 968 pour #2357 : deux CICYs portent la moitié du calcul.

L'unité est donc maintenant le **lot** = (tâche, tranche de réalisations),
`--taille-lot` (défaut 16). Trois bénéfices d'un seul changement : les grosses
tâches se répartissent entre workers, la granularité du checkpoint tombe de 2 h 28
à quelques minutes, et l'écran redevient vivant.

**Mesure** (conteneur à 2 cœurs, 400 candidats) : `-j 1` → 226,4 s, `-j 2` →
119,3 s, soit **1,90×**. À 7 workers, compter un facteur ~6.

**Le checkpoint a dû changer de nature.** Un compteur « les n premiers sont
faits » n'a plus de sens quand le lot 12 finit avant le lot 7. Il enregistre
désormais **l'ensemble des lots terminés, avec leur nombre de lignes**, et il n'y
a plus d'offset. À la reprise, **le fichier fait foi**, en deux passes : compter
les lignes par lot, puis ne garder que les lots dont le compte correspond
exactement. Un lot amputé est recalculé au lieu de rester marqué « fait » avec la
moitié de ses lignes.

Un checkpoint **séquentiel est migré** (`fait: 823` → les 823 tâches, tronquées à
leur offset). Sans cela, passer à la version parallèle aurait jeté six heures de
calcul déjà faites.

**Quatre scénarios donnent un résultat identique ligne pour ligne** : séquentiel,
parallèle, trois coupures successives, et migration depuis un checkpoint
séquentiel.

**Trois défauts trouvés en écrivant ces vérifications**, tous invisibles sans
elles :

| | ce qui se passait |
|---|---|
| **12 discordances d'orbite fictives sur 18** | le représentant est relu depuis le JSONL (`[1, 1]`), le membre de contrôle est encore en mémoire (`(1, 1)`) : `str()` les distingue. Une garde qui crie à tort est aussi inutile qu'une garde muette |
| **18 candidats étiquetés « hors domaine »** | au lieu de « aucun groupe d'ordre compatible » : le test de domaine passait avant celui du groupe |
| **une ligne en double après amputation** | on filtrait le fichier *avant* de restreindre le checkpoint, donc les lignes d'un lot à moitié écrit restaient **et** le lot était recalculé |

---

### 5.27 Les huit candidats à ℤ₂×ℤ₂ : ce sont les seuls qui survivent, et c'est le problème

Le balayage complet donne **33 099 lignes `SURVIT` → 2 857 candidats (B, C)
distincts → 691 orbites** sur **91 CICYs**, toutes SO(10) rang 4, toutes à
3 générations sur le quotient. Intégrité vérifiée : **zéro** ligne où `survit`
serait vrai sans que h⁰ générique, h⁰ équivariant, Hoppe complet et la
surjectivité soient tous passés. Les trois candidats connus retombent où ils
doivent : #6890 → 12 orbites (|Aut| = 1), #6947 → 1 (|Aut| = 24), #6715 → 3.

683 de ces orbites ont Γ = ℤ₂ et butent sur l'argument de rang du §5.8. Seules
**huit**, sur sept CICYs (#22, #480, #2357 ×2, #2534, #2568, #5421, #6829), ont
Γ = ℤ₂×ℤ₂ — donc deux lignes de Wilson, donc la seule route vers le rang 4.

**Elles sont toutes obstruées.** Les deux générateurs **anticommutent** sur
H¹(V) :

```
  #  22 : ordres projectifs (2,2) constantes (1,1)   AB = k.BA,  k = -1
  # 480 : ...                                        k = -1
  (les huit, sans exception)
```

Le relèvement de ℤ₂×ℤ₂ au fibré est une représentation **projective** de cocycle
non trivial. Rééchelonner les générateurs ne le supprime pas : le rapport k est
invariant. H¹(V) = 12 ne se décompose donc pas en 3+3+3+3 ; il porte des
irréductibles de dimension 2 de l'extension centrale. **Le comptage
n_gen = |χ|/|Γ| = 3 n'est pas établi pour ces huit.**

**La règle, mesurée** sur trois CICYs et onze degrés chacune : le cocycle vaut
**(−1)^|a|** — la parité du degré total.

| \|a\| | k |
|---|---|
| 1 | −1 |
| 2 | +1 |
| 5 (#6829) | −1 |
| 6 (#2534) | +1 |

Or les huit ont tous **\|c\| = 5**, et ce n'est pas un hasard : B est une somme
de **cinq** vecteurs unité, donc \|c₁(B)\| = 5, impair, structurellement.

**Et voici ce qui condamne l'interprétation.** Dans le catalogue, parmi les
candidats dont ℤ₂×ℤ₂ est un groupe utile et que le pipeline sait certifier
(rank_C = 1) :

| | nombre | survivants |
|---|---|---|
| \|c\| **pair** — non obstrués | **2 338** | **0** |
| \|c\| **impair** — obstrués | **14** | **8** |

Zéro survivant sur 2 338 candidats physiquement admissibles ; huit sur les
quatorze qui ne le sont pas. **La survie est anti-corrélée avec l'admissibilité.**
Ces huit ne survivent pas malgré le cocycle : tout indique qu'ils survivent
**à cause** de lui — l'espace « équivariant » calculé dans le cas projectif
n'est pas celui d'une structure équivariante véritable, et y annuler h⁰(V) ne
démontre pas ce qu'on croit.

C'est le §5.3 sous une autre forme : un sous-espace équivariant qui n'est pas
celui qu'on croit, et un test qui passe pour cette raison.

**Je retire donc l'annonce faite dans la foulée du scan** — « huit candidats qui
échappent à l'argument de rang ». Ce sont les huit que la machinerie traite mal.

**Ce que cela ouvre.** `decomposition_h1_V_abelien` (nouveau) traite un Γ abélien
à plusieurs générateurs et **refuse de rendre des multiplicités** quand les
générateurs ne commutent pas, au lieu de produire un tableau qui ne décrit rien.
Et la vraie question devient : pourquoi aucun des 2 338 candidats non obstrués
ne survit-il ? Là est le prochain chantier — pas dans les huit.

---

### 5.28 Une case nulle n'est pas une sortie de domaine — et elle fermait la seule route restante

Le §5.27 ferme la route ℤ₂×ℤ₂ : les 2 322 candidats admissibles sont tués par la
stabilité équivariante, et les 8 « survivants » sont des artefacts de relèvement
projectif. Restait à regarder les autres groupes d'ordre 4. Le catalogue en
contient **sept à ℤ₄** — et ℤ₄ est **cyclique**, donc structurellement à l'abri
du cocycle : un générateur unique ne peut pas ne pas commuter avec lui-même.

**Les sept étaient écartés « hors domaine », et jamais évalués.** Cause exacte,
mesurée : **une seule** charge négative chacun, sur les 37 que teste
`domaine_valide` — et toujours de la forme c_j − b_i.

Or ces degrés-là ne sont pas de même nature que les autres. Ce sont les **cases
de la matrice f**. Une case de degré négatif signifie H⁰(O(c_j − b_i)) = 0, donc
une case identiquement nulle — et toute la machinerie le traite déjà ainsi :

| | ce qui est déjà fait |
|---|---|
| `espace_f_equivariant` | calcule `actif = all(x >= 0 ...)` et saute la case |
| `h0_V_generique` | insère un bloc de zéros |
| `decomposition_h1_V` | saute les cases absentes de `offsets` |

La condition était donc **plus stricte que ce que le code consommateur demande**.
Même forme que la ligne `len(c_charges) != 1` du §6 : une condition unique qui
écartait en amont des candidats parfaitement traitables.

**Après relaxation** — les charges qui doivent *porter* des sections (b, c,
bᵢ+bⱼ, c+b) restent exigées positives et certifiées ; seules les cases c−b
peuvent être négatives — **quatre des sept passent** le domaine, et le résultat
est celui qu'on espérait sans oser l'annoncer :

| CICY | ambiant | anomalie | h⁰(V) générique | h⁰(V) équivariant |
|---|---|---|---|---|
| **#6826** | [1,1,4] | OK | 0 | **0** à λ = ±1 |
| **#6836** | [1,1,1,1,3] | OK | 0 | **0** à λ = ±1 |
| **#6836** (2ᵉ) | [1,1,1,1,3] | OK | 0 | **0** à λ = ±1 |
| **#7735** | [1,1,5] | OK | 0 | **0** à λ = ±1 |

Et le test **mord** : les deux autres relèvements λ = ±i donnent h⁰ = 4, 5 ou 6.
Ce ne sont pas des zéros obtenus parce que la contrainte serait vide — c'est
exactement le piège du §5.3, et il est écarté ici par mesure.

Ce sont les **premiers candidats d'ordre 4 sans cocycle** à passer l'étape qui
tue tout le reste.

> **⚠ CONCLUSION RETIRÉE — voir §5.29.** Trois de ces quatre calculent h⁰(V)
> dans un modèle qui **sous-compte les sections**. `domaine_valide` certifie le
> h⁰ de Koszul mais ne vérifie jamais que dim(S/I) lui est égal, et l'écart est
> réel. Ce qui reste acquis du §5.28 : la relaxation des cases nulles est juste
> et testée ; ce qui tombe : les quatre candidats.

**Ce qui manque encore, et c'est une seule chose.** Ils ont tous **rank_C = 2**.
Or `hoppe_sur_espace`, `f_sans_point_base` et `decomposition_h1_V` supposent
rank_C = 1 et se déclarent non calculables au-delà. Le verdict reste donc
`indéterminé`, et le nombre de générations sur le quotient — 12/4 = 3 attendu —
n'est pas encore établi.

> L'item « rank_C = 2 » du §6 change de statut. Ce n'était qu'une limite connue
> de l'outillage ; c'est désormais **la dernière porte** avant la seule route
> vers le rang 4 que ce projet ait jamais ouverte.

Réserve à ne pas perdre : #6826 a h¹(B₅) = 1 ≠ 0, donc H¹(V) n'y est pas
simplement le conoyau — ce candidat-là demandera un traitement à part.

**Le test (37ᵉ)** vérifie les deux sens : une case c−b négative ne fait plus
sortir du domaine, **et** un bᵢ ou un cⱼ négatif le fait toujours. Supprimer le
contrôle de signe restant le fait échouer.

---

### 5.29 Le domaine certifie une chose et en vérifie une autre

Trouvé en généralisant la décomposition de H¹(V) à rank_C = 2 pour trancher le
§5.28 — et c'est la généralisation qui a révélé la faute, pas une relecture.

Les quatre candidats ℤ₄ ont h¹(V) = 12 au catalogue. Le conoyau calculé donne
**11, 8, 8, 14**. Un conoyau qui vaut 8 quand h⁰(C) − h⁰(B) = 12 est impossible :
le rang de f ne peut pas dépasser h⁰(B). La cible était donc trop petite.

```
  #6836, charge (0,0,0,0,1) sur P1^4 x P3 :   dim(S/I) = 4   h0 Koszul = 8
  #6836, charge (0,0,1,1,1) :                 dim(S/I) = 16  h0 Koszul = 24
```

**Le modèle sous-compte d'un facteur deux** — H¹ d'un terme de Koszul contribue
à h⁰(Y) sans avoir d'antécédent polynomial dans l'anneau ambiant.

Et `domaine_valide` ne le voit pas : elle vérifie que le h⁰ de Koszul est
**certifié**, jamais que dim(S_a/I_a) lui est **égal** — ce qui est pourtant
tout ce que le modèle prétend. Elle certifie une chose et en vérifie une autre.

**Ce que cela retire.** Les quatre candidats ℤ₄ du §5.28 calculaient h⁰(V)
équivariant dans un espace qui n'est pas H⁰(Y, ·). Le résultat « h⁰ = 0 à
λ = ±1 » ne veut donc rien dire pour trois d'entre eux. **Le §5.28 est retiré**
sur ce point. Ce qui en reste, et qui est juste : la relaxation des cases c−b
négatives, défendue par le 37ᵉ test.

**Ce que cela ne retire pas.** Vérification faite sur **120 candidats survivants
tirés du scan : 120 cohérents, zéro écart**. Les 691 orbites ne sont pas
touchées. #6890 et #6947 sont exacts sur leurs 26 charges chacun ; **#6715 est
en écart sur une** — (3,0,2,3,0), 48 contre 52 — et mérite un réexamen.

**Ce que je n'ai PAS fait, et pourquoi.** Ajouter `dim(S/I) == h⁰` comme
condition de `domaine_valide` fait tomber **cinq** tests de non-régression, dont
le cas SU(5)/ℤ₂×ℤ₂ de #6947. C'est-à-dire que le corpus de référence du projet
contient des résultats obtenus hors du modèle. Changer le critère d'acceptation
de tout le pipeline dans le même geste que la découverte du problème serait
précisément la précipitation que le §8 proscrit.

La mesure est donc **exposée sans filtrer** : `charges_hors_modele(ambient,
config, charges)` rend la liste des charges en écart. La décision — resserrer le
domaine et réviser ce qui tombe, ou étendre le modèle pour qu'il représente
vraiment H⁰(Y, ·) — reste à prendre, les yeux sur les cinq tests concernés.

---

### 5.30 Où passent les sections manquantes — le critère exact, et ce qu'il reste à construire

Le §5.29 constate que le modèle S_a/I_a sous-compte. Voici **d'où** viennent les
sections qu'il ignore, et donc ce qu'il faudrait construire pour les représenter.

La résolution de Koszul de O_Y(a) sur l'ambiant A donne la suite spectrale
d'hypercohomologie

```
    E1^{-p,q} = H^q(A, ^p)  ==>  H^{q-p}(Y, O(a))
```

où ∧^p = ⊕_{|S|=p} O_A(a − Σ_{k∈S} d_k). **H⁰(Y) reçoit toute la diagonale
q = p**, pas seulement le coin (0, 0) :

```
    h0(Y, O(a)) = [ ligne q = 0 : ce que le modele represente ]
                + [ somme_{p >= 1} des termes q = p : ce qu'il ignore ]
```

**Le cas #6836, charge (0,0,0,0,1)** — quatre termes p = 1 de la forme
b = (0,0,0,−2,0) : un facteur P¹ en degré −2, donc H¹(P¹, O(−2)) = 1 chacun.
Quatre classes, qui portent h⁰ de 4 à 8.

**Ce que ces classes sont** : pas des monômes de l'anneau ambiant, mais des
classes de Čech, du type g_k /(x_j y_j) modulo la relation de Koszul.

`analyse_modele(ambient, config, a)` rend `{naif, manquant, exact, termes,
modele_exact}`.

**Validé, et dans la seule direction qui serve.** `modele_exact` est une
condition **suffisante** de fiabilité : sur l'échantillon, **54 charges déclarées
dans le modèle, 54 fois dim(S/I) = h⁰, zéro contre-exemple**. Au-delà du critère,
`naif + manquant` n'est qu'une **borne supérieure** — les différentielles
supérieures peuvent tuer des termes, et 14 écarts sur 144 le montrent, tous par
excès. Le dire est le point : un critère suffisant qui se présenterait comme une
égalité serait la même faute que celle du §5.29.

**Deux pièges rencontrés en l'écrivant**, tous deux des troncatures :

| | ce qui manquait |
|---|---|
| ligne q = 0 arrêtée à p = 1 | suppose ⊕ₖH⁰(a−d_k) → H⁰(a) **injective** — faux dès qu'il y a des syzygies. 17 désaccords sur 80 |
| diagonale plafonnée à p = 3 | #7293 reçoit sa contribution en **p = 5**, terme de tête, dual de Serre. Tronquer la diagonale, c'est reproduire le défaut qu'on mesure |

**Ce qui reste à construire, et c'est le vrai travail.** Corriger `dimY` ne sert
à rien seul : les consommateurs ont besoin de la **multiplication**
H⁰(O(a)) ⊗ H⁰(O(b)) → H⁰(O(a+b)) et de l'**action de Γ**, sur les classes de
Čech comme sur les monômes. Il faut donc

1. une base explicite des classes de Čech contribuant à chaque charge —
   accessible : ce sont des produits de Künneth de monômes et de classes
   H^{n_i}(P^{n_i}, O(d)) à base monomiale négative (Bott) ;
2. le produit de Čech H¹ × H⁰ → H¹, et la reconnaissance du cas où le produit
   redevient une section ordinaire (différentielle de Koszul) ;
3. l'action de Γ sur ces classes, qui permute les facteurs et les Koszul.

Tant que ce n'est pas fait, `modele_exact` doit **écarter** les charges
concernées plutôt que de les traiter dans un modèle qui les sous-compte — mais
poser ce filtre fait tomber cinq tests de non-régression (§5.29), et cette
révision-là reste à mener.

---

### 5.31 Les classes manquantes, construites — première des trois étapes

Le §5.30 dit **combien** de sections le modèle ignore et **d'où** elles viennent.
`cech.py` en construit maintenant la **base**.

Une classe manquante est un couple **(S, w)** : un sous-ensemble S de p équations,
et un élément w de la base de H^p(A, O(a − Σ_{k∈S} d_k)) — lui-même un produit de
Künneth où chaque facteur projectif porte soit un monôme ordinaire (q_i = 0), soit
un **monôme négatif**, exposants tous ≤ −1 (q_i = n_i, Bott).

Sur **#6836, charge (0,0,0,0,1)**, la base sort telle que la théorie l'annonce :

```
  S=(0,)  ((0,0), (0,0), (0,0), (-1,-1), (0,0,0,0))
  S=(1,)  ((0,0), (-1,-1), (0,0), (0,0),  (0,0,0,0))
  S=(2,)  ((0,0), (0,0), (-1,-1), (0,0),  (0,0,0,0))
  S=(3,)  ((-1,-1), (0,0), (0,0), (0,0),  (0,0,0,0))
```

Quatre classes, chacune **1/(x y) sur un P¹ différent** — les quatre qui portent
h⁰ de 4 à 8. Ce ne sont pas des monômes de l'anneau ambiant, et c'est exactement
pour cela que S/I ne les voyait pas.

**Validé contre deux calculs indépendants déjà présents dans le dépôt** :

| | résultat |
|---|---|
| \|base\| == `manquant` d'`analyse_modele` | **145 charges, 145 d'accord** |
| `cardinal_hq` (Künneth explicite) == `h_ambient` | **435 dimensions, 435 d'accord** |

Le 39ᵉ test vérifie en outre que les exposants sont bien **négatifs sur exactement
un facteur** et que les quatre classes de #6836 sont portées par quatre facteurs
**distincts** : une base entièrement positive, ou concentrée au même endroit,
signalerait qu'on a construit autre chose que ce qu'on croit.

**Restent les deux étapes qui mordent** : le produit de Čech H¹ × H⁰ → H¹, avec la
reconnaissance du cas où le produit redevient une section ordinaire ; puis l'action
de Γ, qui permute à la fois les facteurs projectifs et les générateurs de Koszul.
C'est là que les candidats ℤ₄ redeviendront calculables — pas avant.

---

## 6. Ce qui reste faux ou absent

| | état |
|---|---|
| voie « gros Γ » | **fermée, mesurée.** Le verrou n'était ni la certification Koszul ni les charges négatives, mais la ligne `len(c_charges) != 1` de `domaine_valide` : les 26 candidats à groupe d'ordre compatible sont tous E₆ à rank_C = 2, dont 24 satisfont tout le reste. Contrainte levée (`rank_c_max=None`, défaut 1 pour ne pas casser `hoppe_fast`) puis test relancé : **574 couples, 544 tués par h⁰(V) équivariant, 28 sans f équivariant, 0 survivant**. Aucun candidat à Γ d'ordre ≥ 4 ne passe la stabilité restreinte. Étendre ∧^p V et la surjectivité à rank_C = 2 est donc **inutile** : il n'y aurait rien à leur donner à manger |
| surjectivité au rang 5 | **bloquant ensuite** : les 40 candidats du scan « gros Γ » sont tous de rang 5, précisément le régime où le certificat J_d = R_d n'est pas atteignable (§5.4). Le critère de Hoppe les départagera, mais le verdict final restera `indéterminé` tant que ce point n'est pas traité |
| **branche extension — aucun discriminant en aval** | **mesuré, et structurel.** `hoppe_extension` n'est satisfait que si toutes les bornes sur les quotients gradués s'annulent, donc h⁰(F2) = 0 sur tout candidat retenu, donc le cup-produit ∪e part d'un espace nul : h⁰(V) = 0 **quelle que soit la classe d'extension**, vérifié sur **2 542 / 2 542**. Le test du §5.3 — celui qui tue 3 624 couples chez les monades — y est identiquement vide. Le test au niveau des charges l'est aussi : **2 792 / 2 792 passent, 0 via une permutation non triviale**. La branche produit un catalogue, pas une sélection. **Restent** aussi : H¹(∧²V) donc Higgs et exotiques ; `verify_hoppe.py` ignore les extensions |
| extensions au rang 5 | **non énumérables** : 2,9·10⁵ tuples ordonnés dès m = 2, max_charge 2, au-dessus du plafond de 200 000 ; 1,1·10⁸ à m = 3. Aucune extension de rang 5 n'a jamais été engendrée — la route SU(5) par les extensions est donc **non explorée, faute de pouvoir l'atteindre**, et non par un résultat négatif |
| **portée du verdict de Hoppe** | `stable: True` signifie « non éliminé » (§5.13). **Levé sur `#6890`, `#6947`** (équivariant, §5.14) et `#6715` (générique, §5.15). La phase des twists est branchée dans `hoppe_fast` et `process_cicy`, et a démontré un **faux positif du catalogue, `#7484`**. **Restent sans verdict** : 53 entrées non éliminées mais non prouvées stables, et 61 dont les twists ne sont pas certifiés |
| critère exact au rang 5 | **hors budget, mesuré** : le polytope est instantané, mais les rangs sur ∧^p B pour p = 1..4 dépassent la dizaine de minutes par entrée. 68 des 115 entrées du catalogue sont dans ce cas. Mettre `dimY` en cache par degré est le gain le plus évident |
| décision exacte de la pente | le certificat de Motzkin démontre l'instabilité quand il aboutit ; **la réciproque demanderait un solveur de programmation linéaire**, absent du dépôt (numpy seul). Les verdicts `None` ne sont donc ni des éliminations ni des validations |
| énumération de `generate_extensions` | **faite** (§5.12) : monotone en `max_charge` par construction, avec plafond `--ext-exhaustif-max` et champ `ext_mode` qui dit, résultat par résultat, si l'énoncé est démontré sur le domaine ou seulement sondé |
| ligne de Wilson explicite | **non construite** — le §5.8 est de la théorie des groupes appliquée aux nombres calculés ; le code ne manipule aucune ligne de Wilson. La corrélation entre Γ et la ligne de Wilson, qui décide de 2 ou 6 bidoublets de Higgs, reste un intrant |
| Modèle Standard hors de portée avec ℤ₂ | limitation de **principe** : les lignes de Wilson préservent le rang, SO(10) est de rang 5 et le MS de rang 4 (§5.8). Ces deux candidats plafonnent à Pati–Salam ou SU(5) flipped |

| balayage complet avec Hoppe complet | **fait** (§5.16) : catalogue purgé par la phase des twists (115 → 114), chaîne relancée en entier, mêmes 3 couples survivants. Le §2 tient |
| surjectivité au rang 5 | **le mur était dans la liste des multidegrés** (§5.17), et le balayage corrigé (§5.18) ramène les indéterminés de 449 à **21**. Reste ouverte, et indépendante de l'équivariance (§4.6), la question « le catalogue contient-il des monades non surjectives ? » |
| filtre d'indice | **corrigé** (§5.18) : plus de repli silencieux sur tous les groupes, raison persistée, et le verdict porte n_gen(X/Γ). 73 candidats sont écartés explicitement faute de groupe d'ordre compatible |
| domaine du modèle S/I | 37 candidats sur 108 hors domaine, faute de charges positives, dont `#5452`, `#6826`, `#7745`, `#7669`. Ni retenus ni éliminés. Élargir demande de passer par Koszul plutôt que par le quotient monomial |
| `end_V` (nombre de singlets) | **corrigé (§5.19)** : vaut `None`, et ne rapporte plus de points. h¹(End V) reste non calculé |
| exotiques SO(10) et SU(5) | **corrigé (§5.19)** : `None` au lieu d'un zéro structurel, plus de 25 points gratuits. E₆ conserve son compte réel |
| Higgs E₆ en mode Wilson | **corrigé (§5.19)** : le 3 codé en dur est supprimé. Le compte avec ligne de Wilson demanderait la décomposition des 27 sous Γ, non calculée |
| h^i hors certification | ~52 % des cas — d_r (r ≥ 2) ou ambiguïté de rang |
| annulation d'anomalie | **corrigée (§5.21)** : branchée avec le préfiltre χ, pour les monades comme pour les extensions. Elle écartait **70 des 115** entrées du catalogue. Les deux candidats du §2 passent. Également en place dans `audit_results.py`, donc reproductible sur un catalogue déjà produit (§5.22). **Toujours absente** de `triage_clean.py` et `verify_hoppe.py`, qui travaillent en aval d'un catalogue déjà audité |
| couplages de Yukawa | hors périmètre |

**Point de méthode sur les Higgs sans lignes de Wilson** : un E₆ avec n_anti = 0
n'est une impasse que *sans* quotient, où les Higgs viennent des paires 27 + 27̄.
Avec une ligne de Wilson, les 27 se décomposent selon les représentations de Γ et
les doublets sortent de cette décomposition. **n_anti = 0 est alors la propriété
recherchée, pas un défaut.**

---

## 7. Performance

| | v1 | v2 |
|---|---|---|
| Banc 7 CICYs, `--max-charge 4 --n-random 400` | 206 s | 15 s |
| Générateur positif seul (m=6, rang 4) | 13,2 s | 1,3 s |

Trois leviers :

**Préfiltre χ.** Pour V stable de pente nulle, h⁰ = h³ = 0, donc n_gen = |χ(V)|,
et χ(V) = χ(B) − χ(C) se calcule par arithmétique pure. Sélectivité mesurée :
**0,01 %** sur 797 027 monades — le travail cohomologique divisé par ~10⁴.

**Élagage exact du générateur positif.** Pour r_C = 1, un vecteur c ayant une
composante égale à 1 est infaisable. 15 624 vecteurs c → 4 095 retenus à m = 6,
max_charge = 4, espace de recherche **strictement identique**. Vérifié : 0 faux
rejet sur 10 160 vecteurs infaisables.

**Non-dégénérescence.** Rejette 29 % des monades par comparaisons de charges.

Extrapolation : le scan `--max-ps 6` qui prenait plusieurs **semaines** en v1 est
de l'ordre de la demi-journée en v2, sur un domaine plus large.

La chaîne d'équivariance (`equivariance_f.py`) coûte ≈ 1 h mono-cœur sur les 108
candidats de `scan_wilson2`. Elle n'est pas parallélisée et reconstruit un
`CovariantRing` par couple (candidat, symétrie) : un cache par (CICY, symétrie)
est le gain le plus évident si le besoin s'en fait sentir.

---

## 8. Discipline de validation

**`tests_regression.py` — 39 tests, ~4 min. À lancer après chaque modification,
avant chaque scan.**

Il rassemble toutes les références indépendantes utilisées : c2 sur la bicubique
et la quintique, intégralité de χ, χ du module contre Riemann-Roch, accord des
h^i certifiés, dualité de Serre sur les paires certifiées, χ(∧²V) aux rangs 3
et 4, les monades scindées réelles de `#7669`, l'élagage exact, l'action de Γ
d'ordre 3, la lecture des entrées symboliques `rt[n]`, l'appariement invariant par
permutation, les ordres de groupe, la cible d'indice en mode Wilson, et la
décomposition isotypique.

Huit ajouts de la session « équivariance » :

| test | référence indépendante | vérifié en cassant |
|---|---|---|
| `reduce_vec` annule l'idéal | définition d'une réduction, 90 tirages | oui : l'ancien `Mred = M[:rank]` fait tomber le test |
| idéal Γ-covariant | re-substitution sans le système résolu ; + Hilbert inchangée | — (écart exact 0) |
| équivariance de f mord | `#6947` SO(10)/ℤ₂ survit **et** SU(5)/ℤ₂×ℤ₂ tombe | oui : un test qui accepterait tout échoue sur le second cas |
| h⁰(∧²V) équivariant | accord avec `sections.h0_wedge2_V_explicit` sans contrainte **et** `#6715` 0 → 21 sous contrainte | oui : ignorer la base fait tomber le second volet |
| ∧^p V général | **valeur connue d'avance** : à p = rk, ∧^{rk}V = det V = O donc h⁰ = h⁰(O_Y) = 1 | oui : supprimer le signe (−1)^{k−1} de la contraction fait tomber h⁰(∧⁴V) à 0 |
| spectre sous Γ | **valeur connue d'avance** : Γ libre ⇒ partie invariante = 6/2 = 3, donc 3 + 3 ; + additivité des multiplicités ; + β∘α = 0 et h⁰(∧²V) concordant par deux chemins | oui : ignorer l'image de f casse d'abord le contrôle d'additivité |
| surjectivité de f | contrôle **négatif construit** : un f à facteur commun a un zéro commun, le critère doit refuser ; + λ=+1 certifie et λ=−1 non sur `#6890` | oui : une fonction renvoyant toujours `True` échoue |
| ordre projectif | racines n-ièmes construites comme puissances, comptage = pgcd(n, p−1) | oui : l'ancienne `racines_niemes` rendait la liste vide sur 7³ = 343 |
| phase des twists | les cinq h⁰ de `#7484` figés un par un, et la borne 13 − 12 = 1 ; **deux verdicts opposés** (`#7484` tombe, `#6890` survit) ; et `#21`, dont 5 twists sur 45 ne sont pas certifiés, doit rendre `None` et non `False` | oui, dont un cassage qui passait d'abord : ignorer la certification ne change rien sur `#7484`, tous ses h⁰ étant certifiés |
| Hoppe suffisant | valeur connue d'avance négative (f₁ = 0 ⟹ h⁰(V) = dim H⁰(O(b₁))) ; **anti-vacuité** : les 110 twists doivent avoir une source non vide ; et le twist doit **agir**, la source devant varier avec H | oui, dont un cassage qui passait d'abord : neutraliser `twist` laissait le test vert |
| multidegrés du certificat | marges figées sur `#21` (ancien : tout négatif ; nouveau : tout viable), vérification EXACTE par `dimY` et non par l'estimation, non-régression sur les degrés certifiants de `#6890` et `#6947`, et **le zéro falsy** : une marge nulle doit être conservée | oui, trois fois, dont le bug lui-même — réintroduire `marge or -1` fait tomber le volet des degrés certifiants |
| pente des sous-faisceaux | degré au point J = v contre Riemann-Roch ; **deux verdicts opposés construits** ; et un cas réel où un témoin existe hors grille, qui doit rendre `None` et non `False` | oui, de trois façons : convertir les `None` en `False`, renvoyer toujours `True`, renvoyer toujours `False` |
| extensions énumérées | comptage par convolution contre comptage par énumération ; **contrôle négatif construit** : le tirage doit être vu perdre 245/281 extensions d'un cran de `max_charge` au suivant | oui, dans les deux sens : revenir au tirage casse la monotonie ; rendre le tirage monotone casse le contrôle négatif |

Trois d'entre eux figent **deux verdicts opposés**, et trois confrontent le code à
une valeur connue d'avance : 125 pour la quintique, 1 pour det V = O, 3 + 3 pour
la décomposition du spectre. Un test qui n'exigerait que
la survie passerait pour un module qui accepte tout ; un test qui n'exigerait que
l'élimination passerait pour un module qui rejette tout. Les deux ensemble
n'admettent qu'un module qui discrimine.

### Principe

Un scan coûte deux heures ; un test coûte une seconde. Et surtout : **aucun des
neuf défauts n'a été trouvé par le code lui-même**. Tous l'ont été par une
référence externe — Riemann-Roch, dualité de Serre, une identité de rang 3, une
valeur de la littérature, la définition d'une réduction.

Corollaires pratiques :

- Ne jamais valider un module sur un seul exemple. Le bug du c2 avait été
  « validé » sur la quintique, où m = 1 le rendait invisible (§4.1).
- Ne jamais interpréter un chiffre sans mesurer ce qui le borne. Le « lieu de
  base » du §5.4 était un artefact de source insuffisante, et le contrôle qui l'a
  démenti coûtait trois lignes.

### La règle des filtres

**Un filtre doit déclarer combien il a laissé passer, et pourquoi.** Sans cela,
son silence se lit comme une sélection.

C'est le défaut le plus fréquent de ce dépôt — sept occurrences, toutes de la
même forme : une condition devient vide, ou universellement vraie, et rien ne
le signale. Le résultat continue de sortir, avec l'apparence d'un tri.

| | ce qui devenait vide | ce que le silence faisait croire |
|---|---|---|
| §4.8 | exotiques SO(10) et SU(5), identiquement nuls | « modèle sans exotiques », et 25 points de score |
| §5.3 | sous-espace équivariant, jamais vide quand Γ agit par phases | « le fibré descend au quotient » |
| §5.12 | `deduplicate_results` indexait sur (B, C), absent des extensions | 2 647 candidats repliés sur 132 |
| §5.13 | recherche de témoin J sur une grille trop petite | « 24 % des fibrés sont déstabilisés » |
| §5.17 | marge exactement nulle, falsy en Python | les degrés certifiants des deux candidats, écartés |
| §5.18 | `groupes_utiles` vide ⇒ repli sur **tous** les groupes | 95,5 % des couples calculés hors cible, certains marqués `SURVIT` |
| §5.25 | contrôle du repli tiré **par orbite** : une orbite géante n'en recevait qu'un | un repli entièrement faux validé par « 0 discordance » |
| **§5.23** | **le générateur lui-même : 10 tirages sur une famille de 2 201** | **« ces fibrés n'existent pas sur ces CICYs »** — alors qu'ils n'avaient pas été engendrés |
| **§5.27** | **l'espace « équivariant » d'un relèvement PROJECTIF** | **8 candidats à ℤ₂×ℤ₂** — les seuls survivants, et les seuls obstrués |
| **§5.29** | **`domaine_valide` certifie h⁰ sans vérifier dim(S/I) = h⁰** | un modèle qui sous-compte les sections d'un facteur 2, sans un mot |

La septième est la plus coûteuse : c'est un filtre qu'on n'avait pas identifié
comme tel. Un générateur incomplet ne se distingue d'un résultat d'absence par
aucune trace dans le fichier de sortie.

Deux d'entre eux ont produit un chiffre publiable qui n'existait pas. Aucun
n'aurait été trouvé en relisant le code : tous l'ont été en demandant au filtre
combien il avait laissé passer.

D'où trois exigences, tenues par les tests :

1. **Compter les deux côtés.** Un filtre rapporte ce qu'il retient *et* ce qu'il
   écarte, avec le motif. `hoppe_twists` rend `n_twists` et `non_certifies` ;
   `hoppe_suffisant_sur_espace` rend `sources_non_vides` — sans lui, un critère
   suffisant vérifié sur des sources vides serait vrai sans rien démontrer.
2. **Ne jamais remplacer par zéro ce qui n'est pas calculé.** `None` se voit,
   `0` se lit comme une qualité. Les exotiques, les singlets, la cohomologie
   non déterminée des extensions sont désormais `None`, et ne rapportent aucun
   point.
3. **Persister les lignes écartées, avec leur raison.** Un fichier de résultats
   doit dire pourquoi un cas n'a pas été traité (§5.11, §5.18), faute de quoi
   son silence se lit à tort comme une absence de candidats.

Et le contrôle qui met tout cela à l'épreuve : **casser le filtre dans les deux
sens.** Un module qui accepte tout et un module qui rejette tout doivent chacun
faire tomber un volet différent du test. Sur les quatorze tests ajoutés depuis, sept
cassages « évidents » se sont révélés **passants** au premier essai — c'est en
les voyant passer qu'on a trouvé le vrai angle mort.

---

## 9. Commandes usuelles

```bash
# avant tout
python tests_regression.py

# validation du socle sur les vraies CICYs
python validate_cohomology.py cicylist.txt --n-cicys 60 --max-charge 4

# vérification ciblée d'un candidat, sans attendre le balayage complet
#   equivariance_f.py --cicy N filtre à la LECTURE : quelques minutes au lieu
#   de plusieurs dizaines d'heures. Écrit dans le dossier passé en argument,
#   donc un dossier par candidat pour ne pas s'écraser l'un l'autre.
mkdir scan_w4_c6890 ; copy scan_wilson4\results_equivariant.jsonl scan_w4_c6890\
python -u equivariance_f.py cicyquotients.m cicylist.txt scan_w4_c6890 --cicy 6890
python resume_cible.py scan_w4_c6890 scan_w4_c6947 scan_w4_c6715

# balayage complet avec repli par orbite (§5.25) : ~14 h au lieu de ~55.
#   Aucune ligne ne disparait ; --controle-orbites verifie le repli en cours
#   de route et declare le run invalide s'il trouve une discordance.
python -u equivariance_f.py cicyquotients.m cicylist.txt scan_wilson4 `
       --replier-orbites -j 7 | Tee-Object -FilePath scan_wilson4_equiv_f.log

# reprise d'un balayage interrompu (§5.24) : relancer LA MEME commande.
#   Le checkpoint est dans <dossier>/progress_equivariance_f.json.
#   --reset pour repartir de zéro.
python -u equivariance_f.py cicyquotients.m cicylist.txt scan_wilson4

# scan ciblé Wilson (194 CICYs, |chi| = 3|Gamma|)
python wilson_match.py CicyQuotients.m cicylist.txt --results scan/results_clean.jsonl
python -m cy_landscape.main_optimized cicylist.txt -j 7 --output scan_wilson2 --reset \
       --wilson wilson_cicys.json --max-charge 5 --n-random 5000

# LE MÊME, générateur classique énuméré (§5.23) — c'est celui-ci qu'il faut
# relancer : scan_wilson3 a perdu #6890, #6947 et #6715 faute de les engendrer.
# Compter ~1 h à 8 cœurs. NE PAS écraser scan_wilson2, seule trace des trois.
# (PowerShell : accent grave en continuation de ligne)
python -u -m cy_landscape.main_optimized cicylist.txt -j 7 --output scan_wilson4 --reset `
       --wilson wilson_cicys.json --max-charge 5 --n-random 5000 |
       Tee-Object -FilePath scan_wilson4.log

# scan ordinaire
python -m cy_landscape.main_optimized cicylist.txt --max-ps 6 -j 7 --output scan --reset \
       --max-charge 5 --n-random 3000

# scan avec la branche extension (chemin propre, domaine enumere)
python -m cy_landscape.main_optimized cicylist.txt --max-ps 3 -j 7 \
       --output scan_ext --reset --extensions --max-charge 2 --n-random 200

# triage
python audit_results.py scan_wilson2
python triage_clean.py scan_wilson2
python equivariance.py CicyQuotients.m cicylist.txt scan_wilson2

# scan cible « gros Gamma », a charges positives donc DANS le domaine
#   wilson_gros_gamma.json = wilson_cicys.json filtre aux groupes d'ordre >= 4
#   (57 CICYs sur 194). --positive-only ecarte le generateur classique, qui
#   produit des charges negatives et donc des candidats non testables.
python -u -m cy_landscape.main_optimized cicylist.txt -j 7 --output scan_gros_gamma \
       --reset --wilson wilson_gros_gamma.json --positive-only \
       --max-charge 4 --n-random 3000
# equivariance.py D'ABORD : c'est lui qui produit results_equivariant.jsonl
# et le champ groupes_utiles dont equivariance_f.py se sert.
python -u equivariance.py   cicyquotients.m cicylist.txt scan_gros_gamma
python -u equivariance_f.py cicyquotients.m cicylist.txt scan_gros_gamma
#   Resultat : 3 732 couples, 3 376 tues par h0(V), 356 indetermines (rang 5),
#   0 survivant -- et AUCUN couple evalue sur un groupe d'ordre compatible,
#   ces 26 candidats-la etant hors domaine faute de certification Koszul.

# critere de Hoppe SUFFISANT sur un candidat (twists, §5.14)
#   hoppe_suffisant_sur_espace(anneau, b, c, base, offsets, dims, degres,
#                              p, rng, D)  avec D = hoppe_fast.vecteur_D(d, J)

# chaine complete au niveau des polynomes (~1 h)
python -u equivariance_f.py cicyquotients.m cicylist.txt scan_wilson2
python -u equivariance_f.py cicyquotients.m cicylist.txt scan_wilson2 --cicy 6890
```

Le scan extension écarte d'abord par certificat de pente (§5.13), puis par
Hoppe ; le champ `pente_verdict` de chaque résultat vaut `true` (témoin trouvé)
ou `null` (indéterminé) — jamais `false`, ces cas-là étant écartés.

Options utiles : `--sampling-threshold auto`, `--extensions` (branche extension,
chemin propre du §5.10), `--ext-exhaustif-max` (plafond d'énumération, défaut
200 000 ; 0 force l'échantillonnage), `--n-gen`, et pour `equivariance_f.py` :
`--cicy`, `--tous-groupes`, `--input`.

`--with-extensions` n'existe plus : il activait le chemin par pseudo-monade du
défaut 4.7. Le programme s'arrête avec ce motif plutôt que de l'ignorer.

**Sous PowerShell**, le séparateur est `;` et non `&&`, `tee` s'écrit
`Tee-Object`, et `python -u` est nécessaire pour voir la sortie progresser sur un
run long. `equivariance_f.py` bascule stdout en `errors='replace'` : les
étiquettes de jauge contiennent des indices Unicode (« E₆ ») que la console
cp1252 ne sait pas encoder, ce qui faisait planter le script à la première ligne
de résultat après plusieurs minutes de calcul. Le JSONL de sortie est écrit en
UTF-8 explicite et garde l'étiquette exacte.

---

## 10. Données externes

- **Liste CICY** : `cicylist.txt`, page CICY d'Oxford. Champs `Num`, `NumPs`,
  `NumPol`, `Eta`, `H11`, `H21`, `C2`, `Redun`, puis la matrice de configuration.
  **Ne contient pas les symétries.** `parse_oxford` transpose : `config` est
  K × m (polynômes × facteurs).
- **Quotients libres** : `cicyquotients.m`, page CicyQuotients d'Oxford.
  195 entrées, au format
  `{Num -> n, Conf -> matrice, H11, H21, Symmetries -> {...}}`. Chaque symétrie :
  `{gap-id, "nom", {générateurs sur les coordonnées}, {générateurs sur les
  polynômes}, H11, H21}`. Les entrées valent 0, 1, ou des puissances de `rt[n]`
  = exp(2iπ/n). Les noms contiennent eux-mêmes des `#` et des `$` pour distinguer
  des actions inéquivalentes — ne jamais les utiliser comme séparateur.
  **Le second bloc de générateurs est désormais exploité** (§5.1).

Attention : certaines entrées de `cicylist.txt` ne sont pas des 3-variétés
(dim(ambiant) − K ≠ 3). Les filtrer, sinon elles échouent à 100 % pour une raison
sans rapport avec le module testé.
