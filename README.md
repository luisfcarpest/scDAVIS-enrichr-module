# EnrichR Functional Enrichment Module for scDAVIS

This repository contains the implementation of a functional enrichment analysis module developed for integration into the **scDAVIS** single-cell transcriptomics analysis platform.

## Project Description

The module extends the analytical capabilities of scDAVIS by incorporating native functional enrichment analysis based on the **EnrichR** framework, enabling downstream biological interpretation of gene signatures derived from single-cell RNA sequencing workflows.

The implementation supports multiple gene input strategies, interactive enrichment exploration, and integrated graphical visualization within the existing reactive Shiny architecture of scDAVIS.

## Implemented Features

- Functional enrichment analysis using EnrichR
- Support for multiple gene input modes:
  - User-selected gene sets from scDAVIS
  - Uploaded custom gene lists (`.txt`)
  - Differential expression-derived gene signatures
- Dynamic annotation database selection
- Interactive enrichment result tables
- Reactive graphical visualization of enrichment signatures
- Export functionalities for plots and tabular outputs

## Technical Stack

- **R**
- **Shiny**
- **EnrichR**
- **ggplot2**
- **DT**

## Repository Structure
- module_enrichR.R   # EnrichR module implementation
- server.R           # Main server integration (changes made from the original server.R for the new module integration)
- ui.R               # User interface integration (changes made from the original ui.R for the new module integration)
- EnrichRShinyApp.R  # Application prototype developed prior to module implementation to perform standalone functional enrichment analysis.
