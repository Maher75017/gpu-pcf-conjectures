# gpu-pcf-conjectures (en français)

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24485589-blue.svg)](https://doi.org/10.5281/zenodo.24485589)

Trois fractions continues polynomiales trouvées par une recherche sur GPU. Chacune coïncide avec une forme close sur environ 220 chiffres, et aucune ne figure dans les références et la base de données que nous avons consultées.

**Mise à jour (28 novembre 2026).** B et B' sont maintenant **démontrées** : ce sont les cas k = 2 et k = 3 d'une famille infinie qui vient de la fraction continue d'Euler d'une série simple (voir [plus bas](#mise-à-jour--b-et-b-sont-démontrées-une-famille-infinie)). **A, pour 1/(π − 3), reste une conjecture numérique** : elle n'a pas cette structure. Question sur MathOverflow : [Are these continued fractions known?](https://mathoverflow.net/questions/518279)

Auteur : Maher El Hmissi ([Maher75017](https://github.com/Maher75017)), novembre 2026. Code de recherche et analyse écrits avec l'aide de Claude (Anthropic). English: [README.md](README.md).

Convention (celle de la Ramanujan Machine et de LIReC) : V = a(0) + b(1)/(a(1) + b(2)/(a(2) + b(3)/(a(3) + ...))).

## Les trois formules

| | a(n) | b(n) | Valeur | Statut |
| --- | --- | --- | --- | --- |
| **A** | 7n² + 17n + 7 | 2n²(2n − 1)² | 1 / (π − 3) | conjecture |
| **B** | 10n² + 11n + 4 | −4n²(2n − 1)² | 2 / (3 ln 3 − 4 ln 2) = 2 / ln(27/16) | démontrée (famille, k = 2) |
| **B'** | 20n² + 26n + 9 | −9n²(2n − 1)² | 3 / (10 ln 2 − 6 ln 3) = 3 / ln(1024/729) | démontrée (famille, k = 3) |

`python3 verifier_formules.py` (besoin de `mpmath`) calcule chaque fraction en entiers exacts avec 1 500 puis 3 000 termes : les trois coïncident avec leur forme close sur 219 à 220 chiffres, et les deux profondeurs concordent.

**Correction.** Nous avions écrit qu'aucune n'est une série classique déguisée, parce que les dénominateurs q(n) des réduites ont de grands facteurs premiers. Ce test était incomplet : pour B et B', ce sont les numérateurs p(n) qui n'ont que de petits facteurs premiers, la signature de la fraction continue d'Euler d'une série, et c'est ainsi qu'elles sont démontrées plus bas. Pour A, p(n) et q(n) ont tous deux de grands facteurs premiers (p(3) = 3² × 5 × 7 × 6 311 ; q(3) = 3 × 93 827).

## Mise à jour : B et B' sont démontrées (une famille infinie)

Pour tout réel k > 1,

- a(n) = 2(k² + 1)n² + (3k² − 1)n + k², b(n) = −k²n²(2n − 1)², a(0) = k²,
- V = k / L_k, avec L_k = (k − 1) ln(k − 1) + (k + 1) ln(k + 1) − 2k ln k = Σ_{m≥1} 1 / (m(2m − 1) k^(2m−1)).

k = 2 donne B (L₂ = ln(27/16)), k = 3 donne B' (L₃ = ln(1024/729)), k = 4 donne 4 / ln(84375/65536) avec a(n) = 34n² + 47n + 16, et ainsi de suite.

**Démonstration.** Les termes t_m = 1 / (m(2m − 1)k^(2m−1)) de la série ont pour rapport r_m = t_(m+1)/t_m = m(2m − 1) / (k²(m + 1)(2m + 1)). La fraction continue d'Euler donne S = t₁ / (1 − r₁/(1 + r₁ − r₂/(1 + r₂ − ⋯))) avec t₁ = 1/k, donc V = k/S = k²(1 − r₁/(1 + r₁ − ⋯)). La transformation d'équivalence par c_n = k²(n + 1)(2n + 1) la met sous la forme a(n) = c_n(1 + r_n) et b(n) = −r_n c_(n−1) c_n, qui sont exactement les polynômes ci-dessus, avec a(0) = c₀ = k². La série converge pour k > 1 et la fraction continue gagne environ log₁₀(k²) chiffres par terme (pour k = 1, elle tend encore vers 1/ln 4, lentement).

`python3 verifier_famille.py` (moins d'une seconde) vérifie la valeur sur environ 128 chiffres pour k = 2 … 100 et pour k = 5/2, √2 et π, la structure d'Euler en arithmétique exacte pour k = 2 … 20 et n ≤ 200, et la série.

**Ce que cela veut dire.** B et B' sont des conséquences élémentaires d'une méthode classique (Euler, 1748) : elles ne sont pas nouvelles au sens profond, et nous n'avions pas repéré la famille avant cette mise à jour. A ne rentre pas dans ce cadre : pas de « queues » polynomiales de degré 2, numérateurs p(n) à grands facteurs premiers (les plus grands : 73, 347, 6 311, 27 953, … contre 3, 5, 7, 7, … pour B), et aucune combinaison u·p(n) + v·q(n) avec |u|, |v| ≤ 12 n'a que de petits facteurs premiers. A est la question ouverte. La famille et son identification avec la fraction continue d'Euler ont été trouvées avec l'aide de Claude.

## Vérifications de nouveauté (et leurs limites)

- **Henri Cohen, « Parametric continued fractions for π², ζ(3), and other constants » (arXiv 2304.11727)**, la liste de référence des familles de bas degré : A, B et B' n'y figurent pas. L'article ne donne aucune fraction explicite pour π (seulement la mention de familles paramétriques) ni pour log(a/b).
- **L'article de la Ramanujan Machine (Raayoni et al., Nature 2021 ; arXiv 1907.00205)** : elles ne sont pas parmi ses conjectures.
- **LIReC, la base publique de la Ramanujan Machine** (827 fractions continues, interrogée le 24 novembre 2026 avec `search/check_lirec.py`) :
  - même forme canonique : A, B, B' **absentes** ; contrôles positifs (leur conjecture 8/π² et la fraction d'Apéry pour ζ(3)) **présents** ;
  - fraction de la base dont la limite vaut (p + qV)/(r + sV), entiers jusqu'à 1 000 : **aucune** pour B et B'. Pour A, 24 fractions de la base sont dans ce cas, mais seulement parce que ce sont d'autres formules pour π, et 1/(π − 3) est une transformation de π.
- Formes canoniques LIReC (haut, bas) : A `[[8,-8,2,0,0],[49,140,79,-30,-21]]` ; B `[[-16,16,-4,0,0],[100,20,-29,-3,12]]` ; B' `[[-36,36,-9,0,0],[400,240,-124,-48,27]]`.

**Limites.** La base publique LIReC est une sélection, et la littérature sur les fractions continues est vaste. L'absence dans ces sources est un indice de nouveauté, pas une preuve. Corrections et références bienvenues. Mise à jour : B et B' ont ensuite été démontrées par la fraction continue d'Euler (voir plus haut) ; nous ne connaissons pas de démonstration pour A.

## Comment elles ont été trouvées

Voir le tableau de [README.md](README.md) : six recherches (environ 15 milliards de fractions sur GPU, puis des recherches ciblées sur processeur), chacune avec un contrôle positif (formules connues placées dans l'espace de recherche, par exemple celle d'Apéry pour ζ(3)) et un contrôle négatif (constantes aléatoires). A est apparue au 4ᵉ tir puis de nouveau au 5ᵉ, B au 5ᵉ, B' au 6ᵉ. Au 6ᵉ tir, nous avions écrit « aucun autre membre de cette forme » : **c'était faux**. Les cas k = 4, 5, … de la famille sont dans cet espace, mais l'étape de détection (un test en double précision que exp(m/V) est une fraction simple, m ∈ {1, 2, 3, 4, 6}) les a manqués : précision insuffisante pour k = 4, m = 5 non essayé pour k = 5, fractions trop grandes pour k ≥ 6. Une correspondance à 10 chiffres ne compte qu'après un recalcul en entiers exacts qui coïncide sur au moins 50 chiffres (80 au 6ᵉ tir, plus de 200 pour la vérification finale).

Les résultats négatifs font partie du dossier : aucune formule « à la Apéry » de degré 5 pour ζ(5) avec des coefficients jusqu'à (255, 511, 64) et |c| ≤ 32 (2ᵉ tir).

## Reproduire

Les commandes sont dans [README.md](README.md). Les messages des scripts sont en français.

## Voir aussi

[exact-exp32](https://github.com/Maher75017/exact-exp32) et [tt-erf-fp32](https://github.com/Maher75017/tt-erf-fp32) : la même approche de recherche exhaustive avec contrôles explicites.

## Licence

MIT (voir `LICENSE`).
