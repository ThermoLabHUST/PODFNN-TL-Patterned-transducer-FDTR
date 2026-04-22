# PODFNN-TL-Patterned-transducer-FDTR

MATLAB codes for COMSOL-based dataset generation, POD-FNN surrogate training, transfer learning, and inverse fitting in patterned transducer FDTR.

## Overview

This repository contains the MATLAB source codes developed for the proposed POD-FNN-TL framework for patterned FDTR analysis. The supplied codes cover the main computational steps of the present study, including:

1. COMSOL-based high-fidelity forward signal generation for dataset construction  
2. Training-dataset generation through Latin hypercube sampling  
3. POD-based dimensionality reduction and POD-FNN surrogate training  
4. Transfer-learning-based fine-tuning for expanded parameter domains  
5. Main analysis routines for signal generation, sensitivity evaluation, and inverse fitting  

## Code organization

The codes are organized as separate scripts according to their functions.

- One script is used for high-fidelity dataset generation, including COMSOL model construction, parameter sampling, and FDTR signal extraction.
- A second script is used as the main analysis program, which supports full-order or surrogate-based signal generation, sensitivity analysis, and automatic fitting.
- Additional scripts are provided for POD-FNN surrogate training and transfer-learning-based model updating.

## Software requirements

- MATLAB
- COMSOL Multiphysics 6.3
- LiveLink for MATLAB

The COMSOL-controlled parts of the workflow were implemented using COMSOL Multiphysics 6.3 with LiveLink for MATLAB. To run the COMSOL-based scripts properly, MATLAB should be launched from the COMSOL with MATLAB interface, or otherwise in a LiveLink-enabled COMSOL environment so that the COMSOL-MATLAB connection is available.

For other compatible COMSOL versions, minor modifications to the scripts may be required.

## Scope of the released code

These scripts correspond to the computational workflow used in this study and are intended to support reproduction of the reduced-order modeling, surrogate prediction, transfer-learning adaptation, sensitivity evaluation, and inverse-analysis procedures reported in the associated manuscript.
