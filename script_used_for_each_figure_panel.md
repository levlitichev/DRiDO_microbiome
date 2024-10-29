| Figure | Panel | Description                                                                                                                     | Script                                        |
| ------ | ----- | ------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------- |
| 1      | a     | Experiment schematic.                                                                                                           | n/a                                           |
| 1      | b     | Lifespan.                                                                                                                       | lifespan.Rmd                                  |
| 1      | c     | Schematic of data generation.                                                                                                   | n/a                                           |
| 1      | d     | PCoA of all samples, genera.                                                                                                    | kraken.Rmd                                    |
| 2      | a     | Volcano plot of age coefficients.                                                                                               | kraken.Rmd                                    |
| 2      | b     | Bifido versus age.                                                                                                              | kraken.Rmd                                    |
| 2      | c     | Uniqueness versus age.                                                                                                          | kraken.Rmd                                    |
| 2      | d     | Age prediction, AL, genera.                                                                                                     | prediction.Rmd                                |
| 2      | e     | Most important genera for age prediction.                                                                                       | prediction.Rmd                                |
| 3      | a     | Schematic comparing DO AL mice, B6 mice, and human datasets.                                                                    | n/a                                           |
| 3      | b     | Percent genera associated with age in DO AL, B6, and human datasets.                                                            | compare_DO_to_B6_to_humans.Rmd                |
| 3      | c     | Correlation of age coefficients across DO AL, B6, and human datasets.                                                           | compare_DO_to_B6_to_humans.Rmd                |
| 3      | d     | Uniqueness across DO AL, B6, and human datasets.                                                                                | compare_DO_to_B6_to_humans.Rmd                |
| 3      | e     | Schematic of cohousing experiment.                                                                                              | n/a                                           |
| 3      | f     | Cohousing PCoA.                                                                                                                 | cohousing.Rmd                                 |
| 3      | g     | Distance to uncohoused control.                                                                                                 | cohousing.Rmd                                 |
| 3      | h     | Age prediction, cohousing experiment.                                                                                           | cohousing.Rmd                                 |
| 3      | i     | Uniqueness in cohousing experiment.                                                                                             | cohousing.Rmd                                 |
| 4      | a     | Heritability volcano plot.                                                                                                      | kraken.Rmd                                    |
| 4      | b     | Percent heritable taxa across studies.                                                                                          | perc_heritable_taxa_by_study.Rmd              |
| 4      | c     | PVE for all experimental variables.                                                                                             | kraken.Rmd                                    |
| 4      | d     | Genome-wide QTL mapping.                                                                                                        | QTL.Rmd                                       |
| 5      | a     | Barplots of DR coefficients.                                                                                                    | kraken.Rmd                                    |
| 5      | b     | UMGS1815 increased by diets.                                                                                                    | kraken.Rmd                                    |
| 5      | c     | Ligilactobacillus increased by diets.                                                                                           | kraken.Rmd                                    |
| 5      | d     | Magnitude of DR coefficients.                                                                                                   | kraken.Rmd                                    |
| 5      | e     | Similarity of coefficients across diets.                                                                                        | kraken.Rmd                                    |
| 5      | f     | Genera affected by CR or fasting exclusively.                                                                                   | kraken.Rmd                                    |
| 5      | g     | Emergencia versus diet.                                                                                                         | kraken.Rmd                                    |
| 5      | h     | Roseburia versus diet.                                                                                                          | kraken.Rmd                                    |
| 5      | i     | Predicting binary DR status.                                                                                                    | prediction.Rmd                                |
| 5      | j     | Predicting dietary group.                                                                                                       | prediction.Rmd                                |
| 5      | k     | Predicting dietary group, stratified by diet.                                                                                   | prediction.Rmd                                |
| 5      | l     | Predicting age, trained on AL.                                                                                                  | prediction.Rmd                                |
| 5      | m     | UBA11957 versus age.                                                                                                            | kraken.Rmd                                    |
| 5      | n     | Ligilactobacillus versus age.                                                                                                   | kraken.Rmd                                    |
| 6      | a     | Schematic of association and mediation analysis.                                                                                | n/a                                           |
| 6      | b     | Microbiome-phenotype associations by domain.                                                                                    | pheno_assoc_and_mediation.Rmd                 |
| 6      | c     | Examples of genera associated with body weight.                                                                                 | pheno_assoc_and_mediation.Rmd                 |
| 6      | d     | Paramuribaculum v. percent fat.                                                                                                 | pheno_assoc_and_mediation.Rmd                 |
| 6      | e     | Overlap of microbiome-phenotype association and mediation hits.                                                                 | pheno_assoc_and_mediation.Rmd                 |
| 6      | f     | Heatmap of microbiome-phenotype associations.                                                                                   | pheno_assoc_and_mediation.Rmd                 |
| 6      | g     | Akkermansia v. EE.                                                                                                              | pheno_assoc_and_mediation.Rmd                 |
| 6      | h     | Pathway associated with CO2 production.                                                                                         | pheno_assoc_and_mediation.Rmd                 |
| 6      | i     | Microbiome-phenotype associations by domain, cross-sectional analysis.                                                          | pheno_assoc_and_mediation.Rmd                 |
| ED1    | a     | Number of read-pairs per sample, stratified by sample type.                                                                     | qc.Rmd                                        |
| ED1    | b     | PCoA including control samples.                                                                                                 | kraken_unaggregated.Rmd                       |
| ED1    | c     | Stacked barplots of positive controls.                                                                                          | metaphlan_poscons.Rmd                         |
| ED1    | d     | PCoA without control samples (before aggregating by stool_ID).                                                                  | kraken_unaggregated.Rmd                       |
| ED1    | e     | Same as d, but highlighting samples sequenced from the same DNA, different libraries.                                           | kraken_unaggregated.Rmd                       |
| ED1    | f     | Same as d, but highlighting samples sequenced from the same library.                                                            | kraken_unaggregated.Rmd                       |
| ED2    | n/a   | Sample mix-ups: code not included (requires a file with 3M data points).                                                        | n/a                                           |
| ED3    | a     | Histogram of pairwise sample distances, highlighting 13 outliers.                                                               | qc.Rmd                                        |
| ED3    | b     | Kraken v. MPA, % unclassified.                                                                                                  | kraken_v_metaphlan.Rmd                        |
| ED3    | c     | Kraken v. MPA, mean relab per genus.                                                                                            | kraken_v_metaphlan.Rmd                        |
| ED3    | d     | Example of community-wide and specialized pathways.                                                                             | humann.Rmd                                    |
| ED3    | e     | PCA of all samples, pathways.                                                                                                   | humann.Rmd                                    |
| ED4    | a     | Day of experiment versus age, colored by cohort.                                                                                | create_cross_sectional_slices.Rmd             |
| ED4    | b     | Number of age-associated features with time as a random effect, fixed effect, or omitted.                                       | effect_of_chronological_time.Rmd              |
| ED4    | c     | Correlation of age coefficients between longitudinal (time as fixef, time as ranef, no time) models and cross-sectional models. | effect_of_chronological_time.Rmd              |
| ED5    | a     | Uniqueness versus age, for different numbers of mice per age.                                                                   | kraken.Rmd                                    |
| ED5    | b     | Alpha diversity versus age.                                                                                                     | kraken.Rmd                                    |
| ED5    | c     | Proportion age-associated features, genera versus species.                                                                      | kraken.Rmd                                    |
| ED5    | d     | Uniqueness versus age, species-level data.                                                                                      | kraken.Rmd                                    |
| ED5    | e     | Kraken v. MPA, age coefficients.                                                                                                | kraken_v_metaphlan.Rmd                        |
| ED5    | f     | Volcano plot of age coefficients, pathways.                                                                                     | humann.Rmd                                    |
| ED5    | g     | Example pathway that decreases with age.                                                                                        | humann.Rmd                                    |
| ED5    | h     | Example of community-wide age-associated pathway.                                                                               | humann.Rmd                                    |
| ED5    | i     | Example of specialized age-associated pathway.                                                                                  | humann.Rmd                                    |
| ED5    | j     | Example pathway that increases with age.                                                                                        | humann.Rmd                                    |
| ED5    | k     | Functional uniqueness versus age.                                                                                               | humann.Rmd                                    |
| ED5    | l     | Age prediction, all mice, genera.                                                                                               | prediction.Rmd                                |
| ED5    | m     | Age prediction, AL mice, species.                                                                                               | prediction.Rmd                                |
| ED5    | n     | Age prediction, AL mice, pathways.                                                                                              | prediction.Rmd                                |
| ED5    | o     | Most important pathways for age prediction.                                                                                     | prediction.Rmd                                |
| ED6    | a     | Percent pathways associated with age in DO AL, B6, and human datasets.                                                          | compare_DO_to_B6_to_humans.Rmd                |
| ED6    | b     | Uniqueness across human studies.                                                                                                | compare_DO_to_B6_to_humans.Rmd                |
| ED6    | c     | Blautia across human studies.                                                                                                   | compare_DO_to_B6_to_humans.Rmd                |
| ED6    | d     | ɑ-diversity across human studies.                                                                                               | compare_DO_to_B6_to_humans.Rmd                |
| ED6    | e     | Correlation of pathway age coefficients across DO AL, B6, and human datasets.                                                   | compare_DO_to_B6_to_humans.Rmd                |
| ED6    | f     | Flavin biosynthesis pathway across DO AL, B6, and human datasets.                                                               | compare_DO_to_B6_to_humans.Rmd                |
| ED6    | g     | Correlations with genera for universal age-associated pathways.                                                                 | humann.Rmd                                    |
| ED6    | h     | Schematic of germ-free cohousing experiment.                                                                                    | n/a                                           |
| ED6    | i     | PCoA of germ-free cohousing experiment.                                                                                         | germfree_cohousing.Rmd                        |
| ED7    | a     | Volcano plot of pathway heritability.                                                                                           | humann.Rmd                                    |
| ED7    | b     | Correlations with genera for example heritable pathways.                                                                        | humann.Rmd                                    |
| ED7    | c     | Proportion heritable features, genera versus species.                                                                           | kraken.Rmd                                    |
| ED7    | d     | Kraken v. MPA, heritability.                                                                                                    | kraken_v_metaphlan.Rmd                        |
| ED7    | e     | Heritability, ASReml versus lme4qtl.                                                                                            | asreml_v_lme4qtl.Rmd                          |
| ED7    | f     | Heritability, our study versus Schlamp 2021.                                                                                    | compare_to_schlamp.Rmd                        |
| ED7    | g     | Cross-sectional versus longitudinal heritability.                                                                               | asreml_longitudinal_v_cross_sectional_age.Rmd |
| ED7    | h     | PVE for all experimental variables, pathways.                                                                                   | humann.Rmd                                    |
| ED7    | i     | Allele effects for top 6 QTL.                                                                                                   | QTL.Rmd                                       |
| ED8    | a     | Barplots of DR coefficients, pathways.                                                                                          | humann.Rmd                                    |
| ED8    | b     | L-lysine biosynthesis II increased by DR.                                                                                       | humann.Rmd                                    |
| ED8    | c     | L-lysine biosynthesis II is a specialized pathway.                                                                              | humann.Rmd                                    |
| ED8    | d     | Urea cycle decreased by DR.                                                                                                     | humann.Rmd                                    |
| ED8    | e     | Urea cycle is a community-wide pathway.                                                                                         | humann.Rmd                                    |
| ED8    | f     | Magnitude of DR coefficients, pathways.                                                                                         | humann.Rmd                                    |
| ED8    | g     | Similarity of DR coefficiencts, pathways.                                                                                       | humann.Rmd                                    |
| ED8    | h     | Pathways affected by CR or fasting exclusively.                                                                                 | humann.Rmd                                    |
| ED8    | i     | Predicting DR status, pathways.                                                                                                 | prediction.Rmd                                |
| ED8    | j     | Most important genera for binary DR prediction.                                                                                 | prediction.Rmd                                |
| ED8    | k     | Most important pathways for binary DR prediction.                                                                               | prediction.Rmd                                |
| ED8    | l     | Predicting dietary group, pathways.                                                                                             | prediction.Rmd                                |
| ED8    | m     | Predicting dietary group, pathways, stratified by diet.                                                                         | prediction.Rmd                                |
| ED8    | n     | Proportion DR-associated features, genera versus species.                                                                       | kraken.Rmd                                    |
| ED8    | o     | Kraken v. MPA, DR coefficients.                                                                                                 | kraken_v_metaphlan.Rmd                        |
| ED9    | a     | Predicting age, trained on AL, pathways.                                                                                        | prediction.Rmd                                |
| ED9    | b     | Predicting age, trained on 40% CR.                                                                                              | prediction.Rmd                                |
| ED9    | c     | PCoA of just AL and 40% CR 10 and 28-month-old samples.                                                                         | kraken.Rmd                                    |
| ED10   | n/a   | Histogram of microbiome-phenotype assocation p-values.                                                                          | pheno_assoc_and_mediation.Rmd                 |
| ED10   | n/a   | Histograms of microbiome-phenotype p-values.                                                                                    | pheno_assoc_and_mediation.Rmd                 |
