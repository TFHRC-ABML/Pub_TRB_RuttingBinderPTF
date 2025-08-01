# Data and Analysis Code for TRB Manuscript

## Laboratory Evaluation of the Rutting Performance of Different Lanes at Pavement Testing Facility (PTF) at Asphalt Binder Level

[![GitHub License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/TFHRC-ABML/Pub_TRB_RuttingBinderPTF/blob/main/LICENSE)

📄 **[Paper Link (under review, will be out soon)](XXXX)**

### Authors and Contributors
- *S. Farhad Abdollahi (farhad.abdollahi.ctr@dot.gov)*
- *Behnam Jahangiri (behnam.jahangiri.ctr@dot.gov)*
- *Aaron Leavitt (aaron.leavitt@dot.gov)*
- *Adrian Anderiescu (adrian.anderiescu.ctr@dot.gov)*

---

This repository contains the data and analysis code accompanying the TRB manuscript (currently under review) titled:  
**"Laboratory Evaluation of the Rutting Performance of Different Lanes at Pavement Testing Facility (PTF) at Asphalt Binder Level."**

---

## Table of Contents

1. [Requirements](#requirements)
2. [Dataset](#dataset)
3. [Analysis Code](#analysis-code)
4. [Acknowledgments](#acknowledgments)
5. [Citation](#citation)

---

## Requirements

The analysis is implemented in Python and utilizes several external libraries. For full compatibility, it is recommended to install the specific 
configuration listed below. However, stable upgraded versions of these libraries are also expected to work properly:

- Python 3.8.20  
- `numpy` 1.24.3  
- `scipy` 1.10.1  
- `pandas` 2.0.3  
- `matplotlib` 3.7.2  

---

## Dataset

The dataset is provided as an Excel file, `Data.xlsx`, which contains multiple sheets corresponding to different laboratory tests and material properties:

- **HTPG**: Measured HTPG of tank and recovered binders.
- **MSCR**: Measured recovery (e.g., $R_{3.2}$) and non-recoverable creep compliance (e.g., $J_{nr,3.2}$) parameters for tank and recovered binder at 64°C. 
- **Binder Frequency Sweep**: Detailed binder frequency sweep results at different temperatures ranging from 10°C to 82°C and loading frequencies ranging from 0.1 rad/s to 100 rad/s for all the tank and recovered binders at different aging levels.  

**NOTE**: **Tank binders** refer to the binders collected from the asphalt binder tank in the asphalt plant during the production of the mixtures for construction. However, **recovered binders** refer to the binders extracted and recovered from the loose mixtures sampled right behind the asphalt paver. 

---

## Analysis Code

The analysis is organized into four Jupyter Notebooks:

- `Task01_Binder_HTPG.ipynb`: This notebook includes the codes to read, analyze, and plot the High-Temperature Performance Grade (HTPG) asphalt binders.
- `Task02_Binder_MSCR.ipynb`: This notebook includes the codes to read, analyze, and plot the Multiple Stress Creep Recovery (MSCR) asphalt binders.
- `Task03_Binder_FrequencySweep.ipynb`: This notebook includes the codes to read, analyze, and plot the Zero-Shear Viscosity (ZSV) and Shenoy parameter for asphalt binders, based on the frequency sweep test results.
- `Task04_Binder_CorrelationStudy.ipynb`: This notebook includes the codes to evaluate the correlation between different parameters and plot a correlation graph.

---

## Acknowledgments

We gratefully acknowledge Bethel La Plana, Steve Portillo, Scott Parobeck, and Frank Davis for their efforts in preparing and testing the asphalt binders used in this study.

---

## Citation

If you use this code or dataset in your research, please cite the following:

```bibtex
@misc{AbdollahiRutBinder2025,
  title     = {Laboratory Evaluation of the Rutting Performance of Different Lanes at Pavement Testing Facility (PTF) at Asphalt Binder Level},
  author    = {Abdollahi, Seyed Farhad and Jahangiri, Behnam and Leavitt, Aaron and Anderiescu, Adrian},
  year      = {TBD},
  month     = {TBD},
  publisher = {Transportation Research Board},
  doi       = {XXXXX},
  urldate   = {2025-07-27},
  archiveprefix = {XXX},
  langid    = {english},
  keywords  = {asphalt binder, rutting, high-temperature performance grade (HTPG), MSCR, ZSV, reclaimed asphalt pavement (RAP), recycling agent (RA)}
}
```

---

📬 For questions or feedback, please contact **Farhad Abdollahi** at *farhad.abdollahi.ctr@dot.gov*.