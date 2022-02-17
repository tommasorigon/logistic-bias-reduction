# Bias reduction methods for logistic regression via Diaconis-Ylvisaker conjugate prior
This repository is associated with the article "Rigon, Aliverti (2022). [*Conjugate priors and bias reduction for logistic regression models*"](https://arxiv.org/abs/AXID) and contains code implementing bias-reduction via
Diaconis-Ylvisaker conjugate prior and tutorials reproducing the results of the article.

- [BIRTHWEIGHT](./BIRTHWEIGHT) focuses on the `birthweight` dataset available in the `MASS` package, and contains code to reproduce results from Section 4.1 of the paper. 
- [HIGH-DIMENSIONAL-SYNTHETIC](./HIGH-DIMENSIONAL-SYNTHETIC) focuses on a high-dimensional setting
inspired by Sur and Candes (2019, PNAS) and reproduces results from Section 4.2 of the paper and some additional considerations on computational time in settings with large *n*, large *p* and correlated design matrix.

- [ENDOMETRIAL](./ENDOMETRIAL) focuses on the `endometrial` dataset studied in Heinze and Schemper (2002, Stat. in Med.) and contains code to reproduce results from the Supplementary Materials of the paper.
