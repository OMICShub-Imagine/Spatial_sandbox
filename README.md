# Spatial Sandbox

Internal collection of scripts, workflows and utilities developed by the Bioinfo Core at Institut Imagine for the analysis of spatial omics datasets.

This repository serves as our team's sandbox: a place to develop, test and share analysis code before it is incorporated into production pipelines or dedicated repositories.

## Current contents

### CODEX
- Integration and batch correction workflows
- Neighborhood analysis
- Morphology analysis
- Visualization utilities

### Xenium
- Quality control workflows
- Filtering and preprocessing scripts
- Exploratory analyses

### QuPath
- Groovy scripts for exporting image measurements

## Notes

- Scripts are under active development and may change without notice.
- Code is provided as reusable examples rather than a standardized pipeline.
- Individual scripts may assume specific directory structures or metadata conventions used within our projects.
- Documentation will be expanded as workflows mature.

## Repository organization

```
CODEX/
Xenium/
QuPath/
utils/
```

*(Folders may evolve as new analyses are added.)*

## Requirements

Most analyses are written in **R** and make extensive use of:

- Seurat
- tidyverse
- ggplot2
- patchwork

Additional package requirements are documented within individual scripts.

## Contributing

This is an internal collaborative repository. Team members are encouraged to:

- keep scripts modular and well commented;
- document assumptions and input requirements;
- avoid committing large datasets;
- add a short description when introducing new workflows.

## Disclaimer

This repository is intended for research and method development. Scripts have not necessarily been optimized, standardized, or extensively validated for general use.
