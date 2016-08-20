Algorithms  
-------
sparse_gmcp: MCP-based bi-level selection methods for integrative analysis of multiple datasets under the heterogeneity model.

Maintainer
-------
Jin Liu   <jin.liu@duke-nus.edu.sg>


Publication
-------
Liu, J., Huang, J., & Ma, S. (2013). Integrative analysis of multiple cancer genomic datasets under the heterogeneity model. Statistics in medicine, 32(20), 3509-3521.


Description
-------
Marker selection under the heterogeneity model calls for bi-level selection to determine whether a covariate is associated with response in any study at all as well as in which studies it is associated with responses. In this study, we consider two minimax concave penalty (MCP) based penalization approaches, sparse group MCP and  composite MCP, for marker selection under the heterogeneity model.

Usage
-------
1. sparse group MCP for linear model.
    - source("CallC_int.r")
2. sparse group MCP for Logistic model.
    - source("Logistic_SGM_MM.r")
3. Composite MCP for linear model.
    - source("bi_level.r")
4. Composite MCP for logistic model.
    - source("logistic_bi_level.r")