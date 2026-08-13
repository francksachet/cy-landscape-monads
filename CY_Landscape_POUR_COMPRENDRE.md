# CY Landscape Explorer — de quoi s'agit-il ?

*Document d'introduction. Aucune connaissance préalable en géométrie ou en
physique théorique n'est supposée. Pour l'état technique du projet, voir
`CY_Landscape_Explorer_PROJET.md`.*

---

## 1. La question de départ

Le Modèle Standard de la physique des particules décrit remarquablement bien ce
qu'on observe : trois « générations » de matière — l'électron et ses deux
cousins plus lourds, avec les quarks correspondants — et une poignée de forces.
Mais il ne dit pas *pourquoi* trois générations, ni d'où viennent les masses.

La théorie des cordes propose une réponse d'un genre inhabituel : ces nombres ne
seraient pas des constantes fondamentales, mais des **conséquences de la forme
d'un espace qu'on ne voit pas**.

L'idée est que l'univers a dix dimensions, dont six sont enroulées sur
elles-mêmes à une échelle si petite qu'elles nous échappent. L'image classique
est celle d'un tuyau d'arrosage vu de loin : il paraît être une ligne, à une
dimension, alors qu'il en a deux — la longueur, et le tour du tuyau. Les six
dimensions supplémentaires seraient ainsi « enroulées » en chaque point de notre
espace.

La forme précise de cet enroulement détermine la physique qu'on observe. Changez
la forme, vous changez le nombre de générations de particules.

**Le but du projet est de chercher, parmi des formes possibles, celles qui
donnent exactement trois générations.**

---

## 2. Les formes candidates : les variétés de Calabi–Yau

Toutes les formes ne conviennent pas. La cohérence de la théorie — en
particulier la supersymétrie, une propriété qu'on souhaite préserver — impose
une contrainte géométrique forte. Les espaces à six dimensions qui la
satisfont s'appellent des **variétés de Calabi–Yau**.

On ne sait pas les classer toutes. Mais on sait en construire des familles
entières par une recette simple : prendre un espace ambiant bien compris, et y
découper une forme par des équations polynomiales — comme une sphère est
l'ensemble des points vérifiant x² + y² + z² = 1.

La liste utilisée ici, dite **liste CICY d'Oxford**, contient **7 890** formes
obtenues ainsi. Chacune est décrite par une petite matrice d'entiers qui dit
combien d'équations, de quel degré, dans quel espace ambiant. Tout le calcul
part de cette matrice.

---

## 3. Le champ de jauge : les fibrés vectoriels

Une forme ne suffit pas. Il faut aussi préciser ce que fait le champ de force —
l'analogue du champ électromagnétique — dans ces dimensions cachées.

L'objet mathématique qui décrit cela s'appelle un **fibré vectoriel**. L'idée :
en chaque point de l'espace, on attache un petit espace de directions internes.
Une image imparfaite mais utile : imaginez qu'en chaque point du globe soit
planté un mât, orienté d'une certaine manière. L'ensemble de tous ces mâts,
avec la règle qui dit comment leur orientation varie d'un point à l'autre,
c'est le fibré. Le **rang** est le nombre de directions attachées à chaque
point : rang 3, 4 ou 5 dans ce projet.

Construire directement de tels objets est difficile. On les fabrique donc
indirectement, par une construction appelée **monade**. Le principe est celui
d'une soustraction : on prend deux objets simples, faciles à décrire, et on
définit le fibré cherché comme « ce qui, dans le premier, s'annule par une
application vers le second ». En notation du projet : deux collections de fibrés
en droites B et C, une application f entre elles, et V = noyau de f.

Cette application f est une matrice de polynômes. Le fibré est entièrement
déterminé par les degrés de B et C — appelés **charges** — et par le choix de f.

---

## 4. Les deux conditions à satisfaire

### La stabilité

Toutes les monades ne donnent pas un objet physiquement acceptable. La condition
requise s'appelle la **stabilité**, et elle garantit que la configuration résout
bien les équations du mouvement en préservant la supersymétrie. Sans elle, le
modèle est incohérent.

Le **critère de Hoppe** la traduit en quelque chose de calculable : une famille
de nombres, notés h⁰, doivent tous être **nuls**. Ces nombres comptent des
solutions d'équations ; qu'ils s'annulent signifie qu'aucune direction interne
ne « se détache » du reste, ce qui est exactement ce qu'on veut.

### Le nombre de générations

Ici intervient un fait remarquable. Le nombre de générations de particules se
lit sur un **invariant topologique** de la configuration : une quantité qui ne
dépend que de la forme globale, pas des détails, et qui se calcule par une
formule d'arithmétique pure.

C'est l'analogue du théorème selon lequel, sur n'importe quel polyèdre convexe,
sommets − arêtes + faces = 2, quelle que soit sa forme précise. Rien à mesurer :
on compte.

Ce calcul, appelé **caractéristique d'Euler**, coûte quelques multiplications,
là où les autres quantités demandent de l'algèbre linéaire lourde. Le projet
s'en sert comme **préfiltre** : sur 797 027 monades engendrées, 66 seulement
passent — soit 0,01 %. Le travail coûteux est ainsi divisé par un facteur
10 000. C'est ce qui rend l'exploration possible.

---

## 5. Replier l'espace : quotients et lignes de Wilson

Un problème subsiste. Les constructions naturelles donnent un groupe de forces
trop gros — appelé E₆, SO(10) ou SU(5) selon les cas — dont le Modèle Standard
n'est qu'un morceau. Il faut le « casser » en quelque chose de plus petit.

La méthode passe par un **repliement**. Si la forme possède une symétrie — comme
un anneau qu'on peut faire tourner d'un demi-tour sans le changer — on peut
identifier les points qui se correspondent, et travailler sur la forme repliée
plutôt que sur l'originale. C'est ce qu'on appelle prendre le **quotient** par un
groupe de symétrie Γ.

Ce repliement a deux effets.

D'abord il **divise le nombre de générations** par le nombre d'éléments du
groupe. Pour obtenir 3 générations après un repliement en deux, il faut donc en
viser 6 avant — d'où les cibles de 6, 12, 27 ou 75 générations qu'on voit passer
dans les calculs.

Ensuite il autorise une **ligne de Wilson** : un objet qui n'existe que sur
l'espace replié, et qui casse précisément le groupe de forces. C'est l'outil qui
doit mener au Modèle Standard.

Mais le fibré ne survit pas automatiquement au repliement. Il faut qu'il soit
compatible avec la symétrie — techniquement, **équivariant**. Vérifier cette
compatibilité, et surtout vérifier que le fibré reste stable une fois contraint
par elle, a été l'essentiel du travail récent. Ce n'est pas une formalité : sur
3 878 configurations testées, 3 452 sont éliminées à cette seule étape.

---

## 6. Ce que fait concrètement le programme

Une cascade de filtres, du moins cher au plus cher :

1. lire les 7 890 formes et calculer leurs invariants ;
2. engendrer des monades par millions ;
3. **préfiltre arithmétique** sur le nombre de générations (élimine 99,99 %) ;
4. tester la **stabilité** sur les survivantes ;
5. croiser avec la liste des formes possédant une symétrie exploitable ;
6. reconstruire les équations de la forme pour qu'elles soient **compatibles
   avec la symétrie** — sans quoi le repliement n'a pas de sens ;
7. chercher une application f compatible, puis **retester la stabilité** sous
   cette contrainte ;
8. vérifier que l'objet obtenu est bien un fibré ;
9. calculer le **spectre de particules** sur l'espace replié.

Un balayage complet représente quelques heures de calcul.

---

## 7. Où en est-on

Deux configurations passent aujourd'hui toute la chaîne, sur les formes
numérotées **6890** et **6947** de la liste d'Oxford. Chacune donne :

- un groupe de forces SO(10), repliement par une symétrie d'ordre 2 ;
- **3 générations**, obtenues non pas en divisant 6 par 2, mais en décomposant
  explicitement l'espace des solutions et en comptant ce qui survit au
  repliement ;
- un spectre calculé, incluant les particules dont sortiraient les bosons de
  Higgs.

**Ce ne sont pas des modèles du Modèle Standard, et elles ne peuvent pas
l'être.** Une symétrie d'ordre 2 est trop petite pour casser SO(10) jusqu'au
Modèle Standard ; on plafonne à un groupe intermédiaire, dit de Pati–Salam.
Aller plus loin demande des symétries plus grandes, et les candidats explorés
avec de telles symétries sont tous tombés à l'étape de stabilité.

C'est donc un résultat partiel et honnête : la chaîne fonctionne de bout en
bout, elle produit des configurations vérifiées, mais la cible finale n'est pas
atteinte.

---

## 8. Une remarque de méthode, qui vaut au-delà du sujet

Neuf défauts sérieux ont été trouvés dans ce code au fil du projet. **Aucun n'a
été détecté par le code lui-même.** Tous l'ont été en confrontant un résultat à
une **référence extérieure et indépendante** : une valeur connue de la
littérature, une identité mathématique que le résultat devait vérifier, un
second calcul de la même quantité par un chemin entièrement différent.

Un exemple parlant. Une erreur de facteur 2 avait été « validée » sur une forme
particulière — la quintique — où le terme fautif est absent. Elle a survécu à sa
propre vérification pendant des mois, jusqu'à ce qu'on teste une seconde forme
où il apparaît.

D'où deux règles que le projet s'applique :

- **ne jamais valider un calcul sur un seul exemple** ;
- **ne jamais interpréter un nombre sans mesurer ce qui le borne** — un résultat
  qui semble révéler un phénomène profond peut n'être que l'artefact d'une
  contrainte qu'on a oublié de regarder.

La suite de tests automatiques, aujourd'hui au nombre de 23, est construite sur
ce principe. Chacun s'appuie sur une référence indépendante, et chacun a été
validé en **réintroduisant volontairement le bug** qu'il est censé attraper. Un
test qui ne casse jamais ne prouve rien.
