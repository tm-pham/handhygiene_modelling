# The potential impact of intensified community hand hygiene interventions on respiratory tract infections: <br> A Modelling Study
Thi Mui Pham<sup>1,*</sup>, Yin Mo<sup>2,3,4,5</sup>, Ben S. Cooper<sup>2,3</sup>

<sup>1</sup> Julius Center for Health Sciences and Primary Care of the UMC Utrecht, Utrecht University, Utrecht, The Netherlands <br>
<sup>2</sup> Centre for Tropical Medicine and Global Health, Nuffield Department of Medicine, University of Oxford, Oxford, United Kingdom <br>
<sup>3</sup> Mahidol-Oxford Tropical Medicine Research Unit, Mahidol University, Bangkok, Thailand <br>
<sup>4</sup> Division of Infectious Diseases, Department of Medicine, National University Hospital, Singapore <br>
<sup>5</sup> Department of Medicine, National University of Singapore, Singapore <br>
<sup>*</sup> <thi.mui.pham@posteo.de>

## Abstract
Increased hand hygiene amongst the general public has been widely promoted as one of the most important non-pharmaceutical interventions for reducing transmission during the ongoing COVID-19 pandemic and is likely to continue to play a key role in long-term efforts to suppress transmission before a vaccine can be deployed. For other respiratory tract infections community hand hygiene interventions are supported by evidence from randomised trials, but information on how effectiveness in reducing transmission scales with achieved changes in hand hygiene behaviour is lacking. This information is of critical importance when considering the potential value of substantially enhancing community hand hygiene frequency to help suppress COVID-19. Here, we developed a simple model-based framework for understanding the key determinants of the effectiveness of changes in hand hygiene behaviour in reducing transmission and use it to explore the potential impact of interventions aimed at achieving large-scale population-wide changes in hand hygiene behaviour. Our analyses show that the effect of hand hygiene is highly dependent on the duration of viral persistence on hands and that hand washing needs to be performed very frequently or immediately after hand contamination events in order to substantially reduce the probability of infection. Hand washing at a lower frequency, such as every 30 minutes or with a delay of 15 minutes after contamination events, may be adequate to reduce the probability of infection when viral survival on hands is longer, such as when hands are contaminated with mucus. Immediate hand washing after contamination is more effective than hand washing at fixed-time intervals even when the total number of hand washing events is similar. This event-prompted hand washing strategy is consistently more effective than fixed-time strategy regardless of hand contamination rates and should be highlighted in hand hygiene campaigns.

## Model code
The simulation code was built by Thi Mui Pham using R (version 4.0.1).

### Most important functions

handhygiene_functions.R: Functions to compute the probability of infection with respect to hand hygiene frequency/delay as a function of half-life of prob. of persistence 

handhygiene_main.R: Main code to produce Figure 3 in main text of the manuscript

handhygiene_fig4.R: Code to produce Figure 4 in main text of the manuscript

handhygiene_figS1_halflife_CMI.R: Code to produce Figure S1

handhygiene_num_time_hw.R: Code to produce Figure S2

handhygiene_heatmap.R and handhygiene_plot_combine_heatmap.R: Code to produce Figure S20 and Figure S21

handhygiene_plotting_functions.R: Functions for plotting


