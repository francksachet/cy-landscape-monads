# CY Landscape Explorer — état du projet

> **Ce document remplace intégralement les versions précédentes**, y compris
> l'addendum de session `CY_Landscape_Explorer_SESSION_EQUIVARIANCE.md`, qui est
> supprimé. La version antérieure décrivait l'équivariance de f comme bloquée par
> le cocycle de `#7669`, listait six candidats « réellement contraints », et
> annonçait 15 tests. Aucune de ces trois affirmations ne tient : le verrou est
> levé, la partition en « six contraints » mesurait la portée de l'ancien test et
> non une propriété des candidats, et la suite compte 25 tests.

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

```
cy_landscape/
├── main_optimized.py          scan principal (multiprocessing, checkpoint)
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
    ├── equivariant_monad.py   f équivariante, h⁰(∧^p V) restreint pour     [+]
    │                          tout p (donc Hoppe complet, h³ inclus),
    │                          certificat de surjectivité,
    │                          décomposition de H¹(V) et H¹(∧²V) sous Γ
    ├── extensions.py          fibré d'extension : χ et rang corrects,       [~]
    │                          Hoppe par bornes sur les quotients gradués
    ├── gamma_action.py        action de Γ sur les sections (obsolète, §5.2)
    └── cohomology.py          extraction du spectre (partiellement obsolète, §6)

racine/
├── tests_regression.py        25 tests — À LANCER AVANT CHAQUE SCAN
├── validate_cohomology.py     harnais de validation du socle
├── audit_results.py           triage 1 : cohérence interne
├── triage_clean.py            triage 2 : n_anti, familles, doublons
├── verify_hoppe.py            re-vérification a posteriori de la stabilité
├── wilson_match.py            croisement avec la liste de Braun
├── equivariance.py            test nécessaire sur les charges
└── equivariance_f.py          chaîne complète au niveau des polynômes       [+]
```

`[+]` = ajouté lors de la session « équivariance ». `[~]` = réparé (§5.10).

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
Désactivée par défaut ; `--with-extensions` la réactive.

### 4.8 Compteurs de spectre structurellement nuls

Dans `cohomology.py` : `sp.n_exotics = 0` en dur pour SO(10) ; pour SU(5),
`max(0, n_10 + n_10bar − n_gen − 2·n_anti)` vaut **identiquement zéro** puisque
|a−b| + 2·min(a,b) = a+b. Seul E₆ compte réellement ses anti-générations.

Le « zéro exotique » de tous les candidats SO(10) et SU(5) est donc une
constante, pas un résultat — et il vaut 25 points gratuits dans le score.
**Non corrigé** (voir §6).

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

**Ce qui reste avant de scanner** : énumérer `generate_extensions` au lieu
d'échantillonner (§5.11), un test de régression sur sa monotonie, le branchement
dans `process_cicy`, et le renommage de `--with-extensions` pour que l'ancien
chemin cassé cesse d'être atteignable.

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
`max_charge` 2 à 3 en perdait **216 sur 222**, soit 97 %. Le premier point est
corrigé (RNG dérivé par couple (rk1, rk2)) ; **le second ne l'est pas** et ne
peut pas l'être par un réglage de RNG — changer les bornes du tirage change la
suite tirée. Seule une énumération le règlerait.

**Deux pièges d'outillage rencontrés au passage.** `enumerer_positive_B`
matérialisait toute l'énumération en mémoire : `MemoryError` après 153 min de
scan, avec sept workers. Elle est devenue un générateur — le plafond bornait le
NOMBRE de B, pas la mémoire. Et `equivariance_f.py` n'écrivait que les lignes
`etat == 'ok'` : le JSONL montrait 0 couple sur un groupe d'ordre compatible
alors que 26 candidats en portaient un, écartés en amont pour une raison que le
fichier ne contenait pas. Toutes les lignes sont désormais persistées.

---

## 6. Ce qui reste faux ou absent

| | état |
|---|---|
| voie « gros Γ » | **fermée, mesurée.** Le verrou n'était ni la certification Koszul ni les charges négatives, mais la ligne `len(c_charges) != 1` de `domaine_valide` : les 26 candidats à groupe d'ordre compatible sont tous E₆ à rank_C = 2, dont 24 satisfont tout le reste. Contrainte levée (`rank_c_max=None`, défaut 1 pour ne pas casser `hoppe_fast`) puis test relancé : **574 couples, 544 tués par h⁰(V) équivariant, 28 sans f équivariant, 0 survivant**. Aucun candidat à Γ d'ordre ≥ 4 ne passe la stabilité restreinte. Étendre ∧^p V et la surjectivité à rank_C = 2 est donc **inutile** : il n'y aurait rien à leur donner à manger |
| surjectivité au rang 5 | **bloquant ensuite** : les 40 candidats du scan « gros Γ » sont tous de rang 5, précisément le régime où le certificat J_d = R_d n'est pas atteignable (§5.4). Le critère de Hoppe les départagera, mais le verdict final restera `indéterminé` tant que ce point n'est pas traité |
| **branche extension** | chemin de calcul correct et testé (§5.10), taux de conclusion ~36 %. **Restent** : énumérer `generate_extensions` (non monotone en `max_charge`, §5.11), un test de régression sur cette monotonie, le branchement dans `process_cicy`, le renommage de `--with-extensions` |
| énumération de `generate_extensions` | **à faire** — c'est elle qui donnerait des énoncés démontrés plutôt que des sondages. L'espace est un produit de boîtes avec la dernière charge fixée par c₁(V) = 0, donc énumérable pour m et rk petits |
| ligne de Wilson explicite | **non construite** — le §5.8 est de la théorie des groupes appliquée aux nombres calculés ; le code ne manipule aucune ligne de Wilson. La corrélation entre Γ et la ligne de Wilson, qui décide de 2 ou 6 bidoublets de Higgs, reste un intrant |
| Modèle Standard hors de portée avec ℤ₂ | limitation de **principe** : les lignes de Wilson préservent le rang, SO(10) est de rang 5 et le MS de rang 4 (§5.8). Ces deux candidats plafonnent à Pati–Salam ou SU(5) flipped |

| balayage complet avec Hoppe complet | le balayage du §2 a été fait avant l'ajout de h³ (§5.5). Il ne peut que resserrer : `survit` exige maintenant strictement plus. Les deux retenus ont été revérifiés individuellement, le reste est à repasser |
| surjectivité au rang 5 | critère hors de portée (§5.4). La question « le catalogue contient-il des monades non surjectives ? » est **ouverte et non instruite**, et indépendante de l'équivariance (§4.6) |
| domaine du modèle S/I | 37 candidats sur 108 hors domaine, faute de charges positives, dont `#5452`, `#6826`, `#7745`, `#7669`. Ni retenus ni éliminés. Élargir demande de passer par Koszul plutôt que par le quotient monomial |
| `end_V` (nombre de singlets) | valeur de remplissage codée en dur — **sans aucune valeur** |
| exotiques SO(10) et SU(5) | structurellement nuls (§4.8) — fausse le classement, pas la sélection |
| colonne `H` en mode Wilson | calcule `max(0, n_gen − 3)` avec le 3 en dur ; en amont n_gen vaut 6, 9, 27… donc le chiffre affiché n'a pas de sens physique |
| h^i hors certification | ~52 % des cas — d_r (r ≥ 2) ou ambiguïté de rang |
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

**`tests_regression.py` — 25 tests, ~4 min. À lancer après chaque modification,
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

---

## 9. Commandes usuelles

```bash
# avant tout
python tests_regression.py

# validation du socle sur les vraies CICYs
python validate_cohomology.py cicylist.txt --n-cicys 60 --max-charge 4

# scan ciblé Wilson (194 CICYs, |chi| = 3|Gamma|)
python wilson_match.py CicyQuotients.m cicylist.txt --results scan/results_clean.jsonl
python -m cy_landscape.main_optimized cicylist.txt -j 7 --output scan_wilson2 --reset \
       --wilson wilson_cicys.json --max-charge 5 --n-random 5000

# scan ordinaire
python -m cy_landscape.main_optimized cicylist.txt --max-ps 6 -j 7 --output scan --reset \
       --max-charge 5 --n-random 3000

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

# chaine complete au niveau des polynomes (~1 h)
python -u equivariance_f.py cicyquotients.m cicylist.txt scan_wilson2
python -u equivariance_f.py cicyquotients.m cicylist.txt scan_wilson2 --cicy 6890
```

Options utiles : `--sampling-threshold auto`, `--with-extensions` (réactive la
branche cassée), `--n-gen`, et pour `equivariance_f.py` : `--cicy`,
`--tous-groupes`, `--input`.

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
