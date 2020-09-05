# OSEM
Title: Learning Bayesian networks from ordinal data

Date: 2020.09.02

Author: Xiang Ge Luo

Supervisor: Dr. Jack Kuipers

Abstract: Bayesian networks are powerful frameworks for studying the dependency structure of variables in a complex system. The problem of learning Bayesian networks is tightly associated with the given data types. Ordinal data, such as stages of cancer, Likert scale survey questions, and letter grades for exams, is ubiquitous in applied research. However, existing solutions are mainly for continuous and categorical data. In this thesis, we propose an iterative score-and-search method - called the Ordinal Structural EM (OSEM) algorithm - for learning Bayesian networks from ordinal data. Unlike traditional approaches with the multinomial distribution, we explicitly respect the ordering amongst the categories. More precisely, we assume that the ordinal variables originate from marginally discretizing a set of Gaussian variables, which follow in the latent space a directed acyclic graph. Then, we adopt the structural EM algorithm and derive closed-form scoring functions for efficient graph searching. Through simulation studies, we demonstrate the superior performance of our method compared to the alternatives and analyze various factors that may influence the learning results.

Description: The OSEM algorithm is contained in the R file `ordinalScore.R`. Other source code is directly taken from the BiDAG package (https://github.com/cran/BiDAG), with a few modifications made in order to utilize the new algorithm. 

References: 

            J. Kuipers, P. Suter and G. Moffa (2018) <arXiv:1803.07859v2>

            M. Kalisch et al.(2012) <doi:10.18637/jss.v047.i11>
