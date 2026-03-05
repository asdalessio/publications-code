# Novel g-computation algorithms for time-varying actions with recurrent and semi-competing events

## Alena Sorensen D’Alessio, Lucas M. Neuroth, Jessie K Edwards, Chantel L. Martin, Paul N Zivich

**Abstract**

Background: A core aspect of epidemiology is determining the impacts of potential public health interventions over time. With long follow-up periods, epidemiologists may need to consider semi-competing events, in which a terminal event, like death, precludes a non-terminal event, like hypertension. Time-varying confounding poses an additional challenge when studying time-varying interventions or actions. Existing methods do not simultaneously address semi-competing events and time-varying confounding.

Methods: We propose two novel g-computation algorithms for causal effects with semi-competing events and time-varying actions. To explore performance of our novel g-computation estimators, we conducted a Monte Carlo simulation study. We then applied our estimator to investigate how cigarette smoking prevention throughout young and middle adulthood might impact prevalent hypertension using data from Waves III (aged 18-26 years) - VI (aged 39-51 years) of the National Longitudinal Study of Adolescent to Adult Health.

Results: Our simulations show that the novel g-computation estimators had little bias and appropriate confidence interval coverage. They outperformed existing alternative estimators across sample sizes. In the illustrative application, the novel estimator identified a small reduction in prevalence of hypertension and risk of death in midlife had all cigarette smoking been prevented across follow-up compared to the observed smoking patterns. 

Conclusion: As long-running cohorts progress in age, death within the study sample will become an increasing concern for studies of aging-related outcomes, life course analyses, and investigations into chronic disease development. Our novel g-computation estimators provide a simultaneous solution.

**File Manifesto**

`data/`

- `single_sim_observed_data.csv`: Single simulated data set using the data generating mechanism described in the paper

`Rcode/`

- Estimators described in the paper
  - `semicomp_gcomp_estimator_functions.R`: functions for the proposed g-computation estimators (ICE and Standard g-computation)
  - `alternative_estimator_functions.R`: functions for the alternative estimators
- Code to replicate the simulation results in the paper
  - `multistate_DGM.R`: data generating mechanism described in the paper
  - `sim.R`: full simulation experiement comparing proposed estimators to alternative estimators
  - `run_sim_job.R`: batch submission to run `sim.R` in parallel
  - `sim_functions.R`: functions to assist with processing simulation results
  - `evaluate_sim.R`: calculation of simulation evaluation metrics as described in the paper
- Code to evaluate `single_sim_observed_data.csv`
  - `single_sim.R`: Application of proposed estimators using single simulated data set
    - Requires `semicomp_gcomp_estimator_functions.R`,`alternative_estimator_functions.R`,`multistate_DGM.R`, and `sim_functions.R`
