### High throughput molecular dynamics analysis of protein units ###

Ce projet a pour objectif la création d'un outil permettant une analyse visuelle rapide et automatisée des unités protéiques (PU). A cet effet, le script analysis_PU.py a été développé, avec la version 3.7.9 de Python.

# Utilisation

Pour faire fonctionner le script analysis_PU.py, il faut tout d'abord que l'arborescence dans laquelle il se trouve soit correcte. Ainsi, le script analysis_PU.py soit être dans le même dossier qu'un dossier 'data', contenant les PUs. Dans ce dossier 'data', il faut que des dossiers de répliques soient présentes, sous le nom 'PUs_sele_pdb_Replicat_[numéro de le réplique]'. Chacun de ces dossiers doit contenir des sous-dossiers relatifs à chaque PU étudiée, et contenant à leur tour tous les fichiers de dynamiques moléculaires nécessaires. Cette arborescence peut être schématisée ainsi, dans le cas de l'utilisation de trois répliques :

Analysis_PU/
├── analysis_PU.py
└── data/
    ├── PUs_sele_pdb_Replicat_1/
    │   ├── PU_1/
    │   │   └── *contient tous les fichiers relatifs à la dynamique moléculaire de PU_1*
    │   ├── PU_2/
    │   │   └── *contient tous les fichiers relatifs à la dynamique moléculaire de PU_2*
    │   └── ...
    ├── PUs_sele_pdb_Replicat_2/
    │   ├── PU_1/
    │   │   └── *contient tous les fichiers relatifs à la dynamique moléculaire de PU_1*
    │   ├── PU_2/
    │   │   └── *contient tous les fichiers relatifs à la dynamique moléculaire de PU_2*
    │   └── ...
    └── PUs_sele_pdb_Replicat_3/
        ├── PU_1/
        │   └── *contient tous les fichiers relatifs à la dynamique moléculaire de PU_1*
        ├── PU_2/
        │   └── *contient tous les fichiers relatifs à la dynamique moléculaire de PU_2*
        └── ...
        
Pour lancer le script analysis_PU.py, la ligne de commande est la suivante :

    python analysis_PU.py -pu [nom de la PU à analyser] -nbr [nombre de répliques à traiter]

# Fichiers de sortie

Pour la PU traitée, et pour chaque réplique la concernant, neuf fichiers d'analyses sont générés dans le dossier 'outputs' :
- une image PyMOL de la structure tridimensionnelle de la PU,
- le graphe suivant l'évolution du RMSD,
- le graphe témoignant du RMSF par résidu,
- le graphe suivant l'évolution du rayon de giration,
- la carte de contact entre les résidus, avec la fréquence associée à ces contacts,
- la carte de fréquence des protein blocks le long de la séquence,
- le logo plot associé aux protein blocks,
- le graphe du Neq le long de la séquence,
- une image contenant toutes les sorties précédentes.

# Exemple de commande fonctionnelle

Par exemple, pour la PU 1F5QB_4_142_198 en trois répliques, la ligne de commande sera :

    python analysis_PU.py -pu 1F5QB_4_142_198 -nbr 3

