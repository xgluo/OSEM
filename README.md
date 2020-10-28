# OSEM
Title: Learning Bayesian networks from ordinal data

Date: 2020.10.27

Authors: Xiang Ge Luo, Giusi Moffa, Jack Kuipers

Abstract: Bayesian networks are a powerful framework for studying the dependency structure of variables in a complex system. The problem of learning Bayesian networks is tightly associated with the given data type. Ordinal data, such as stages of cancer, rating scale survey questions, and letter grades for exams, are ubiquitous in applied research. However, existing solutions are mainly for continuous and nominal data. In this work, we propose an iterative score-and-search method - called the Ordinal Structural EM (OSEM) algorithm - for learning Bayesian networks from ordinal data. Unlike traditional approaches designed for nominal data, we explicitly respect the ordering amongst the categories. More precisely, we assume that the ordinal variables originate from marginally discretizing a set of Gaussian variables, whose structural dependence in the latent space follows a directed acyclic graph. Then, we adopt the Structural EM algorithm and derive closed-form scoring functions for efficient graph searching. Through simulation studies, we illustrate the superior performance of the OSEM algorithm compared to the alternatives and analyze various factors that may influence the learning accuracy. Finally, we demonstrate the applicability and practicality of our method with a real-world application on psychological survey data from 408 patients with co-morbid symptoms of obsessive-compulsive disorder and depression.

Description: The OSEM algorithm is contained in the R file `ordinalScore.R`. Other source code is directly taken from the BiDAG package (https://github.com/cran/BiDAG), with a few modifications made in order to utilize the new algorithm. 

References: 

            J. Kuipers, P. Suter and G. Moffa (2018) <arXiv:1803.07859v2>

            M. Kalisch et al.(2012) <doi:10.18637/jss.v047.i11>
