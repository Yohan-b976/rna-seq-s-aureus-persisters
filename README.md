# RNA-seq Reproducibility Pipeline (Nextflow + Docker)

Ce dépôt propose un pipeline **Nextflow DSL2** entièrement reproductible pour l’analyse RNA-seq de *Staphylococcus aureus* en conditions de persistance, basé sur l’étude :

> [Nature Communications](https://www.nature.com/articles/s41467-020-15966-7) (2020) “Intracellular Staphylococcus aureus persisters upon antibiotic exposure.”  

Chaque étape s’exécute dans un conteneur Docker dédié pour garantir la reproductibilité et la portabilité.

**Résumé du workflow :**  
FASTQ → trimming → download du génome → indexation → alignement → comptage → analyse différentielle.
---

##  Structure du dépôt
```
rna-clean/
├── README.md
├── run.sh
├── main.nf
├── nextflow.config
├── reports
├── .gitignore
├── containers/
│   ├── bowtie/
│   │     └── Dockerfile
│   ├── cutadapt/
│   │     └── Dockerfile
│   ├── deseq2/
│   │     ├── Dockerfile
│   ├── sratoolkit/
│   │     └── Dockerfile
│   ├── subread/
│   │     └── Dockerfile
│   └── tidyverse/
│         └── Dockerfile
└── data/
│     ├── mapping_aureowiki.tsv
│     └── samples.tsv
├── results (créé automatiquement)
└── scripts/
      ├── run_deseq2.R
      ├── deseq2_plots.R
      └── gene_pathway_array.R
```

---

##  Conteneurs utilisés

Toutes les étapes utilisent des images Docker **construites localement** à partir des Dockerfiles du dossier `containers/`.

Versions exactes :

* **SRA Toolkit** : latest
* **Cutadapt** : 1.11
* **Bowtie (Bowtie1)** : 0.12.7
* **Samtools** : 0.1.19
* **Subread / featureCounts** : 1.4.6-p3
* **DESeq2** : 1.16 (R 3.4.1 via micromamba)
* **Tidyverse** : 4.3.2 (R, avec bioconductor-keggrest, tidyverse, venndiagram, ggrepel, scales)

Build des images :

```
cd containers/bowtie && docker build -t alantrbt/bowtie:latest .
cd containers/cutadapt && docker build -t alantrbt/cutadapt:1.11 .
cd containers/sratoolkit && docker build -t alantrbt/sratoolkit:latest .
cd containers/subread && docker build -t alantrbt/subread:latest .
cd containers/deseq2 && docker build -t alantrbt/deseq2:latest .
cd containers/tidyverse && docker build -t alantrbt/tidyverse:latest .
```

---

##  Données analysées

Les identifiants SRA proviennent de l'étude originale. Ils sont listés dans :

```
data/samples.tsv
```

Ce sont des **lectures single-end**.

---

##  Exécution du pipeline

Assurez-vous que Nextflow et Docker sont installés.

### 1. Lancer le pipeline

```
bash run.sh
```

Ou directement :

```
nextflow run main.nf -with-docker
```

### 2. Résultats produits

Tous les résultats seront écrits dans `results/` :

```
results/
├── raw_fastq/
├── trimmed/
├── index/
├── aligned/
├── counts/
└── deseq2/
```

Les productions finales principales sont :

```
results/deseq2/deseq2_results.csv   # Résultats tabulaires DESeq2
results/deseq2/plots/               # Graphiques et visualisations générés
```
> Remarque :
> Le nombre de CPU utilisé par chaque étape du pipeline est défini par le paramètre cpus dans le fichier nextflow.config (valeur par défaut : 2).
> Adaptez cette valeur selon les ressources disponibles sur votre machine ou serveur pour optimiser les performances du workflow.
---

##  Description du workflow

Le pipeline contient 7 étapes :

1. **Téléchargement FASTQ** (SRA Toolkit)
2. **Trimming** (Cutadapt 1.11)
3. **Téléchargement du génome & GFF3** (NCBI — CP000253.1)
4. **Indexation** (Bowtie1)
5. **Alignement single-end** (Bowtie1 + Samtools)
6. **Comptage** (featureCounts)
7. **Analyse différentielle** (DESeq2 1.16)

Le code du pipeline complet est dans `main.nf`.

---

##  Analyse différentielle (DESeq2)

Le script R utilisé est embarqué dans le conteneur :

```
scripts/
├── run_deseq2.R
├── deseq2_plots.R
└── gene_pathway_array.R
```

Trois scripts R sont utilisés pour l’analyse différentielle et la production des graphiques :

scripts/run_deseq2.R :

* Fusionne les fichiers `counts_*`
* Génère la matrice d'expression
* Applique DESeq2
* Écrit le fichier `deseq2_results.csv`

scripts/deseq2_plots.R :

* Génère les graphiques principaux (MA plots, visualisations des gènes liés à la traduction, etc.) à partir des résultats DESeq2
* Produit des fichiers PDF dans results/deseq2/plots/

scripts/gene_pathway_array.R :

* Crée la table d’association gènes / pathways à partir de KEGG
* Utilisée pour annoter les résultats et les graphiques
---
## Fichiers de métadonnées et d’annotation

* data/samples.tsv :
Contient la liste des échantillons utilisés dans l’analyse, avec leur nom, condition biologique (control/persister) et identifiant SRA. Ce fichier est utilisé pour associer chaque fichier de comptage à sa condition expérimentale lors de l’analyse différentielle.

* data/mapping_aureowiki.tsv :
Fichier d’annotation des gènes de Staphylococcus aureus, issu d’AureoWiki. Il associe chaque locus_tag à un symbole, une description du produit et un identifiant de gène. Ce fichier permet d’enrichir les résultats et les graphiques avec des informations fonctionnelles et des annotations biologiques.

> N’hésitez pas à modifier ou compléter ces fichiers selon vos besoins pour analyser d’autres conditions, échantillons ou jeux de gènes.

---

##  Reproductibilité

* Toutes les versions logicielles sont figées dans les Dockerfiles.
* Le pipeline est entièrement décrit par `main.nf`.
* Le workflow peut être relancé avec :

```
nextflow run main.nf -resume
```

---

##  Références

> [Nature Communications](https://www.nature.com/articles/s41467-020-15966-7) (2020) “Intracellular Staphylococcus aureus persisters upon antibiotic exposure.”  

---

## Equipe

- Turbot Alan
- Beaumatin Yohan
- Faramus François
- Bouteville Tristan
