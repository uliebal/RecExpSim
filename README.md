# BioLabSim: Data Analysis for Biotechnology

## Introduction
This Jupyter book is BioLabSim: a collection of workflows to simulate different steps of strain engineering and fermentation in industrial biotechnology. The data is generated from a virtual organism simulation with various models of microbial metabolism, genetics and physiology {cite:p}`Liebal2023`. 

### How Complex Data Permeates Biotechnology
<!-- ```{margin}
<img src='Figures/Jupyter/FermProSim_Objective.png' alt='Variables in a fermentation'  width='250'/><br>
``` -->
Data analysis plays a crucial role in biotechnology by providing insights essential for optimizing fermentation processes. Effective fermentation necessitates precise control of variables such as nutrient concentrations, pH levels, and temperature. Monitoring and measuring these variables using advanced sensors and analytical tools generate data that informs decision-making and process adjustments. For instance, in microbial fermentation, understanding the dynamics of metabolite production and consumption is pivotal for achieving high yields and product quality {cite:p}`nielsen2016engineering`. Accurate data analysis enables biotechnologists to uncover patterns and correlations in the data, guiding the modification of fermentation conditions to enhance productivity and efficiency {cite:p}`Narayanan2020`. Without robust data analysis, harnessing the potential of biotechnological processes would be severely limited.
<!-- 
```{margin}
<img src='Figures/Jupyter/multi-omics.png' alt='Data diversity'  width='250'/><br>
``` -->
The integration of high-throughput measurement techniques {cite:p}`Wehrs2020` and omics technologies, such as genomics, proteomics, and metabolomics, has significantly elevated the complexity of data in biotechnology {cite:p}`Pinu2019`. These advanced methodologies capture an immense volume of molecular information, offering a comprehensive view of intricate biological processes. However, this expanded scope comes with a challenge – the data becomes more intricate and multifaceted. Traditional analytical methods struggle to cope with the sheer volume and intricacy of omics data, which encompasses intricate molecular interactions across various biological networks. This complexity demands the development and application of novel computational tools and algorithms to extract meaningful insights {cite:p}`Volk2020`. The advent of high-throughput and omics technologies has thus transformed biotechnological research, requiring a paradigm shift in data analysis approaches to fully exploit the potential of these data-rich techniques.

### Improving Biotech with New Computational Methods
<!-- ```{margin}
<img src='Figures/Jupyter/MetaEngSim_Segment_DNAMetabol.png' alt='New methods and analysis types'  width='250'/><br>
``` -->
As biotechnological data becomes increasingly complex and intricate, the demand for sophisticated analysis methods and advanced statistical techniques grows. The intricate nature of modern data, stemming from diverse sources such as omics technologies and high-throughput experiments, requires analytical approaches that can unravel intricate patterns and correlations. Techniques like machine learning, network analysis, and multivariate statistics are now crucial for extracting meaningful insights from such complex datasets  {cite:p}`Oliveira2019`. Moreover, the adoption of FAIR (Findable, Accessible, Interoperable, and Reusable) standards is essential to ensure that intricate biotechnological data remains comprehensible and shareable {cite:p}`Rehnert2022`. These standards facilitate data integration and collaboration, enabling researchers to collectively address intricate challenges and unlock new frontiers in biotechnology.

Paragraph about design-build-test-learn cycles and the automated strain design workflow.

:::{dropdown} Summary
* Data analysis is crucial in biotechnology, guiding effective fermentation through precise measurement of variables, such as nutrient levels and metabolite concentrations, to optimize microbial production systems.
* Integration of high-throughput measurement and omics technologies elevates data complexity, challenging traditional methods and demanding novel computational approaches for meaningful insights.
* Industrial biotechnology's intricate data necessitates advanced analysis methods like machine learning and network analysis, while FAIR standards ensure accessible and collaborative data management.
:::

---

## User Notes

To run BioLabSim directly choose among the following JupyterHubs:

- RWTHjupyter: RWTH Aachen University JupyterHub: [![](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/uliebal/RecExpSim/HEAD)

- RWTHjupyter: RWTH Aachen University JupyterHub: [![](https://jupyter.pages.rwth-aachen.de/documentation/images/badge-launch-rwth-jupyter.svg)](https://jupyter.rwth-aachen.de/hub/spawn?profile=biolabsim)

It is possible to download the examples in BioLabSim and run them locally. This requires the installation of packages to do the simulations, see [Developer Notes](DevopNotes).


(DevopNotes)=
## Developer Notes

### Setup project for local development

```bash

# Setup the python virtual environment next to it. (use Python 3.9)
python3.9 -m venv py39-env

# Activate your environment. (Broad topic that depends on what software and OS is used)
source py39-env/bin/activate

# Clone the repository to a nearby folder.
git clone https://git.rwth-aachen.de/ulf.liebal/biolabsim.git repo-biolabsim

# Enter the newly cloned repository.
cd repo-biolabsim

# Install all required python libraries.
pip install -r requirements.txt

# See the Notebook for examples on how to use the library.
```



## Contacts

*Ulf Liebal*

Institute of Applied Microbiology-iAMB, Aachen Biology and Biotechnology-ABBT, RWTH Aachen University, Worringerweg 1, 52074 Aachen Germany



<!-- Last update: 8 June, 2022

Contact: ulf.liebal@iamb.rwth-aachen.de -->

Licence: See LICENCE file @https://git.rwth-aachen.de/ulf.liebal/biolabsim, or @https://github.com/uliebal/BioLabSim

## References
:::{bibliography}
:::