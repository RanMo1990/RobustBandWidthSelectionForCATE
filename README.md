📖 Robust Plug-in Bandwidth Selection for CATE Estimation
This repository contains the theoretical development of a robust plug-in bandwidth selection method for nonparametric estimation of Conditional Average Treatment Effects (CATE) under inverse probability weighting (IPW) and robust loss functions.

📌 Project Overview
Derived new asymptotic expressions for bias and variance of weighted and robust kernel estimators.
Developed a plug-in bandwidth formula accounting for:
1.Treatment assignment mechanisms via IPW terms (Ti / π(Xi)).
2.Heavy-tailed distributions and outlier-robust loss functions.
3.Extended classical kernel regression plug-in methods to causal inference settings with weighted samples.

🔧 Requirements
This project uses the following R packages:
  1.stats	Core R package for statistical modeling (glm)
  2.MASS	For robust regression models via rlm() and psi.bisquare
  3.robustbase	For functions like mad() used in robust scale estimation
  4.KernSmooth	For kernel smoothing and plug-in bandwidth selection (dpill, dpik, cpblock, blkest)

🧠 Key Technical Challenges
Handling IPW weights inside kernel smoothing.
Managing bias terms induced by robust losses and nonsymmetric error structures.
Deriving explicit, tractable forms for optimal bandwidth selection under these complexities.

📈 Outcomes
A complete asymptotic bias-variance expansion customized for CATE estimators with robustness.
A plug-in bandwidth selection rule balancing robustness and statistical efficiency.
Mathematical proofs verified through detailed theoretical derivations.

🔧 Repository Structure
/theory/ — Full derivations, formulas, and LaTeX source.
/simulation/ — Simulation scripts to validate theoretical results.
/application/ — Application examples on real-world datasets such as NHANES.


🚀 Future Work
Implement simulation studies to verify theoretical predictions.
Develop Python/R packages for automated robust bandwidth selection in causal inference workflows.

🧑‍💻 About
Ph.D. in Statistics | Specialized in robust inference, nonparametric methods, and causal inference.
