# APE007 - Analyse d'Enrichissement des Phénotypes

Ce projet combine des analyses en Julia et Python pour étudier l'enrichissement de différents phénotypes à partir de données de séquençage (fastq) d'APE007.

## Structure du Projet
```text
Serine_proteases_project/
├── data/
│ ├── raw/ # Fichiers fastq bruts
│ ├── processed/ # Données traitées
│ │ ├── count_matrix.csv
│ │ ├── normalized_counts.csv
│ │ └── sequence_quality_metrics.csv
│ └── phe_analysis/
│ ├── enrichissement.csv
│ ├── enrichment_with_controls.png
│ └── phenotype_statistics.csv
├── src/
│ ├── Julia/
│ │ ├── count_analysis.jl
│ │ ├── enrichment.jl
│ │ └── utils.jl
│ └── visualisation_enrichissement.ipynb
└── README.md
```
## Prérequis

### Julia
1. Julia (version 1.6 ou supérieure)

2. Packages Julia requis :
```julia
using Pkg
Pkg.add([
    "DataFrames",
    "CSV",
    "Statistics",
    "Plots",
    "StatsBase",
    "FastaIO",     # Pour la lecture des séquences
    "BioSequences" # Pour le traitement des séquences
])
```

### Python
Les dépendances Python sont listées dans requirements.txt. Installation :
```bash
pip install -r requirements.txt
```

Versions principales :
- pandas==1.4.2
- numpy==1.22.3
- matplotlib==3.5.1
- seaborn==0.11.2
- scipy==1.8.0
- jupyter==1.0.0

## Format des Données d'Entrée

### Données Brutes
- Format : fichiers FASTQ (.fastq ou .fq)
- Localisation : `data/raw/`
- Structure attendue :
  - Reads de séquençage avec scores de qualité
  - Headers contenant les identifiants uniques

## Analyse des Données

### Étape 1 : Traitement des Données avec Julia

1. Dans le dossier APE007_analysis :
```bash
julia
```

2. Exécution des analyses :
```julia
include("src/Julia/count_analysis.jl")
include("src/Julia/enrichment.jl")
```

### Étape 2 : Visualisation avec Python
```bash
jupyter notebook src/visualisation_enrichissement.ipynb
```

## Résultats

### Sorties Julia
1. Fichiers de comptage (`data/processed/`) :
   - `count_matrix.csv` : Matrice de comptage brut
   - `normalized_counts.csv` : Comptages normalisés
   - `sequence_quality_metrics.csv` : Métriques de qualité

2. Fichiers d'enrichissement (`data/phe_analysis/`) :
   - `enrichissement.csv` : Valeurs d'enrichissement calculées
   - Colonnes : mutant, phenotype, enrichissement

### Sorties Python
1. Visualisations (`data/phe_analysis/`) :
   - `enrichment_with_controls.png` : 
     * Boxplots des enrichissements par phénotype
     * Points individuels
     * Lignes de référence pour les contrôles
   
2. Statistiques :
   - `phenotype_statistics.csv` :
     * Moyennes par phénotype
     * Écarts-types
     * P-values vs contrôles
     * Nombre de séquences par phénotype

## Description des Analyses

### Analyse Julia
1. Traitement des données brutes :
   - Lecture des fichiers FASTQ
   - Contrôle qualité des séquences
   - Comptage des occurrences

2. Calcul des enrichissements :
   - Normalisation des comptages
   - Calcul des ratios d'enrichissement
   - Transformation log2
   - Comparaison avec les contrôles

### Analyse Python
Visualisation et statistiques :
1. Statistiques descriptives/ les phénotype attribué viennent de GenBank 
2. Visualisation des distributions
3. Tests statistiques vs contrôles

## Interprétation des Résultats

- Enrichissement > 0 : enrichissement positif
- Enrichissement < 0 : enrichissement négatif
- Références :
  - Vert (ref_tryp) : contrôle positif
  - Rouge (ref_chym) : contrôle négatif
  - Gris (0) : pas d'enrichissement

## Dépannage

Problèmes courants :
1. Erreur de lecture FASTQ :
   - Vérifier l'intégrité des fichiers
   - Vérifier les permissions

2. Erreur de mémoire :
   - Augmenter la mémoire disponible pour Julia
   - Traiter les fichiers par lots

3. Erreur de visualisation :
   - Vérifier la présence de tous les fichiers intermédiaires
   - Vérifier les versions des packages Python

## Contact

Paveyranne Amaury : amaury.paveyranne@espci.fr

## Licence

