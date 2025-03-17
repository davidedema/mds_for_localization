<p align='center'>
    <h1 align="center">MDS for localization</h1>
    <p align="center">
    Project for Distributed Robot Perception at the University of Trento A.Y.2024/2025
    </p>
    <p align='center'>
    Developed by:<br>
    De Martini Davide <br>
    </p>   
</p>

----------

- [Project Description](#project-description)
- [Installation](#installation)
- [Running the project](#running-the-project)


## Project Description
The aim of this project is to implement Multidimensional Scaling (MDS) techinques for the localization problem. Centralized and distributed solution are implemented accounting also for ambiguities.


This codebase implements the algorithms explained in the following papers:

- [Extending the Classical Multidimensional Scaling Algorithm Given Partial Pairwise Distance Measurements](https://ieeexplore.ieee.org/document/5419066) by Amar et al.
- [Solving Ambiguities in MDS Relative Localization](https://ieeexplore.ieee.org/document/7251461) by Di Franco et al.
- [Robust Distributed Network Localization with Noisy Range Measurements](https://dl.acm.org/doi/10.1145/1031495.1031502) by Moore et al.

## Installation

Clone the repository
```
git clone https://github.com/davidedema/mds_for_localization.git
```

## Running the project
Open the project in matlab and run the `main.m` file, a menu will appear giving the possibility to test:

1) Partial connectivity scenario (Amar et al.)
2) Fully connectivity scenario (Di Franco et al.) 
3) Distributed scenario (Moore et al.)
4) Classical MDS
