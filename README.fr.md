# gpu-pcf-conjectures (en français)

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24485589-blue.svg)](https://doi.org/10.5281/zenodo.24485589)

Trois fractions continues polynomiales trouvées par une recherche sur GPU. Chacune coïncide avec une forme close sur environ 220 chiffres, et aucune ne figure dans les références et la base de données que nous avons consultées. Ce sont des **conjectures numériques, pas des théorèmes**.

Auteur : Maher El Hmissi ([Maher75017](https://github.com/Maher75017)), novembre 2026. Code de recherche et analyse écrits avec l'aide de Claude (Anthropic). English: [README.md](README.md).

Convention (celle de la Ramanujan Machine et de LIReC) : V = a(0) + b(1)/(a(1) + b(2)/(a(2) + b(3)/(a(3) + ...))).

## Les trois conjectures

| | a(n) | b(n) | Valeur |
| --- | --- | --- | --- |
| **A** | 7n² + 17n + 7 | 2n²(2n − 1)² | 1 / (π − 3) |
| **B** | 10n² + 11n + 4 | −4n²(2n − 1)² | 2 / (3 ln 3 − 4 ln 2) = 2 / ln(27/16) |
| **B'** | 20n² + 26n + 9 | −9n²(2n − 1)² | 3 / (10 ln 2 − 6 ln 3) = 3 / ln(1024/729) |

`python3 verifier_formules.py` (besoin de `mpmath`) calcule chaque fraction en entiers exacts avec 1 500 puis 3 000 termes : les trois coïncident avec leur forme close sur 219 à 220 chiffres, et les deux profondeurs concordent.

Ce ne sont pas des séries classiques déguisées : les dénominateurs q(n) des réduites ont de grands facteurs premiers (pour A, q(3) = 3 × 93 827 ; pour B, q(3) = 2 × 3 × 11 × 2 557), contrairement à la fraction continue d'Euler d'une série hypergéométrique simple.

## Vérifications de nouveauté (et leurs limites)

- **Henri Cohen, « Parametric continued fractions for π², ζ(3), and other constants » (arXiv 2304.11727)**, la liste de référence des familles de bas degré : A, B et B' n'y figurent pas. L'article ne donne aucune fraction explicite pour π (seulement la mention de familles paramétriques) ni pour log(a/b).
- **L'article de la Ramanujan Machine (Raayoni et al., Nature 2021 ; arXiv 1907.00205)** : elles ne sont pas parmi ses conjectures.
- **LIReC, la base publique de la Ramanujan Machine** (827 fractions continues, interrogée le 24 novembre 2026 avec `search/check_lirec.py`) :
  - même forme canonique : A, B, B' **absentes** ; contrôles positifs (leur conjecture 8/π² et la fraction d'Apéry pour ζ(3)) **présents** ;
  - fraction de la base dont la limite vaut (p + qV)/(r + sV), entiers jusqu'à 1 000 : **aucune** pour B et B'. Pour A, 24 fractions de la base sont dans ce cas, mais seulement parce que ce sont d'autres formules pour π, et 1/(π − 3) est une transformation de π.
- Formes canoniques LIReC (haut, bas) : A `[[8,-8,2,0,0],[49,140,79,-30,-21]]` ; B `[[-16,16,-4,0,0],[100,20,-29,-3,12]]` ; B' `[[-36,36,-9,0,0],[400,240,-124,-48,27]]`.

**Limites.** La base publique LIReC est une sélection, et la littérature sur les fractions continues est vaste. L'absence dans ces sources est un indice de nouveauté, pas une preuve. Aucune démonstration n'est tentée ici (par exemple par transformations de Bauer–Muir ou contraction de fractions de Gauss). Corrections et références bienvenues.

## Comment elles ont été trouvées

Voir le tableau de [README.md](README.md) : six recherches (environ 15 milliards de fractions sur GPU, puis des recherches ciblées sur processeur), chacune avec un contrôle positif (formules connues placées dans l'espace de recherche, par exemple celle d'Apéry pour ζ(3)) et un contrôle négatif (constantes aléatoires). A est apparue au 4ᵉ tir puis de nouveau au 5ᵉ, B au 5ᵉ, B' au 6ᵉ. Une correspondance à 10 chiffres ne compte qu'après un recalcul en entiers exacts qui coïncide sur au moins 50 chiffres (80 au 6ᵉ tir, plus de 200 pour la vérification finale).

Les résultats négatifs font partie du dossier : aucune formule « à la Apéry » de degré 5 pour ζ(5) avec des coefficients jusqu'à (255, 511, 64) et |c| ≤ 32 (2ᵉ tir).

## Reproduire

Les commandes sont dans [README.md](README.md). Les messages des scripts sont en français.

## Voir aussi

[exact-exp32](https://github.com/Maher75017/exact-exp32) et [tt-erf-fp32](https://github.com/Maher75017/tt-erf-fp32) : la même approche de recherche exhaustive avec contrôles explicites.

## Licence

MIT (voir `LICENSE`).
