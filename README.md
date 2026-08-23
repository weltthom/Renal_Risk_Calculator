Renal Risk Calculator

An interactive Shiny app that estimates the 5- and 10-year risk of eGFR falling below 60 mL/min/1.73 m² in patients with eGFR ≥60 mL/min/1.73 m² at baseline, based on a Cox proportional hazards model developed on pooled Freiburg and Zurich cohort data.

Live app: https://thomas-welte.shinyapps.io/Risk_Calcultor_ShinyApp/

What it does

Given a patient's sex, age, HbA1c, eGFR, triglycerides, and total cholesterol, the calculator returns the model's estimated probability of eGFR dropping below 60 mL/min/1.73 m² within 5 and 10 years.

The calculator is restricted to patients with eGFR ≥60 mL/min/1.73 m² at baseline, matching the population in which the model was developed and validated. Predictions outside the supported input ranges (shown in the app) are not available.

This tool is intended to illustrate the published model and to support clinical judgment, not to replace it. Individual predictions should be interpreted alongside a patient's full clinical picture.

Repository contents
File	Description
app.R	Shiny application (UI + server + risk calculation)
fit_dev_final_pooled.rds	Fitted Cox model (pooled Freiburg + Zurich development cohort)
Running locally
r
install.packages(c("shiny", "survival", "dplyr"))
shiny::runApp()
Deploying as a static site (shinylive)

This app can be exported to a fully static, serverless site using shinylive, and hosted on GitHub Pages with no backend R process required:

r
install.packages("shinylive")
shinylive::export(".", "site")

Preview locally before publishing:

r
httpuv::runStaticServer("site/")

A GitHub Actions workflow can automate the export-and-deploy step on every push to main — see the r-shinylive deployment example.

Model

Cox proportional hazards model, pooled Freiburg + Zurich development cohort (n = 27,059; 7,818 events). Predictors: sex, age, HbA1c, eGFR, triglycerides, total cholesterol at baseline. Endpoint: time to eGFR <60 mL/min/1.73 m².

Publications
Development (Freiburg cohort): Arnold F et al., "HbA1c-dependent projection of long-term renal outcomes," J Intern Med 2024. doi:10.1111/joim.13736
External validation (Zurich cohort): add citation / DOI once published
Citing this tool

Add a suggested citation (e.g. app DOI via Zenodo, or a reference to the validation publication) once available.
