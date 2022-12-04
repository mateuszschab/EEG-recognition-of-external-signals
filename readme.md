# EEG recognition of external signals

![MATLAB](https://img.shields.io/badge/MATLAB-R2021b-blue.svg)
[![GitHub Issues](https://img.shields.io/github/issues/mateuszschab/EEG-recognition-of-external-signals.svg)](https://github.com/mateuszschab/EEG-recognition-of-external-signals/issues)
![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)

## Basic Overview
While collecting the experiment, the user stared at the LED, which blinked at different frequencies. The brain in certain areas synchromizes its electrical activity with the frequency of the diode. The task is to process the signal in different ways and determine the frequency of the blinking diode, and ultimately determine the effectiveness of each method. 

Methods used:
1) the spectrum of the EEG signal
2) Features determined from the Butterworth filter
3) Features determined from single frequency bars
4) Features determined from the CCA.

A detailed description of the task in Polish can be found in the file. The efficiency, as well as the conclusions can be found in the MATLAB program code. 

**Acknowledgements**
---

+ [West Pomeranian University of Technology in Szczecin](https://www.wi.zut.edu.pl/en/) Faculty of Computer Science and Information Technology - place of data collection and processing