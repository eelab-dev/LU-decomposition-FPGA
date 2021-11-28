# EE4 Project Report Template

## Prerequisites

This template requires, as a minimum, the following packages:
- babel
- isodate
- ifdraft
- setspace
- hyperref
- biblatex
- glossaries
- fancyhdr
- caption
- subcaption

The following is a list of optional packages. The relevant example usage is
commented out.
- microtype
- pgfplots
- amsmath
- listingsutf8
- graphicx

Compiling also requires the following executables to be present on your system:
- pdflatex
- biber
- makeglossaries

## Compiling

To successfully compile this template, you need to run this series of commands:
    pdflatex ee4thesis
    makeglossaries ee4thesis
    biber ee4thesis
    pdflatex ee4thesis
    pdflatex ee4thesis

This will most likely differ slightly on a Windows host.
