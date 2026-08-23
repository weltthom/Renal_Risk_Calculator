# ============================================================
# RENAL RISK CALCULATOR
# Freiburg + Zurich pooled development model
# ============================================================

library(shiny)
library(survival)
library(dplyr)
library(ggplot2)


# ============================================================
# 1. LOAD MODEL
# ============================================================

fit_dev_final <- readRDS("fit_dev_final_pooled.rds")


# ============================================================
# 1b. LIPID UNIT CONVERSION CONSTANTS
# ============================================================
# Model was fit using mg/dL for triglycerides and cholesterol.
# Standard SI conversion factors, used to convert mmol/L inputs
# to mg/dL before prediction.

TG_MGDL_PER_MMOL <- 88.57
CHOL_MGDL_PER_MMOL <- 38.67


# ============================================================
# 2. PREDICTION FUNCTION
# ============================================================
# Adapted from the pooled-model script. Two changes from that
# version, both to keep the app consistent with its original
# clinical framing (calculator is for patients with eGFR >=60,
# predicting risk of dropping below that threshold):
#   - eGFR validation floor raised from 15 to 60 (was matched to
#     the app's UI restriction, per your choice to keep >=60)
#   - eGFR validation ceiling left at 150, matching the function
#     as given (UI numericInput max is separately capped at 200 -
#     see note below if you want these aligned too)

renal_risk <- function(
    sex,
    age,
    HbA1c,
    eGFR,
    triglycerides,
    cholesterol,
    times = c(5, 10)
) {

  # ----------------------------------------------------------
  # Input validation
  # ----------------------------------------------------------

  if (!sex %in% c("Male", "Female", "M", "F", "male", "female", "m", "f")) {
    stop("sex must be Male or Female.")
  }

  if (age < 18 || age > 100) {
    stop("Age is outside the supported range.")
  }

  if (HbA1c < 4 || HbA1c > 12) {
    stop("HbA1c is outside the supported range.")
  }

  if (eGFR < 60 || eGFR > 150) {
    stop("This calculator is validated for eGFR >=60 mL/min/1.73 m2 (up to 150).")
  }

  if (triglycerides <= 0) {
    stop("Triglycerides must be > 0.")
  }

  if (cholesterol <= 0) {
    stop("Cholesterol must be > 0.")
  }

  if (any(times <= 0)) {
    stop("Prediction times must be positive.")
  }

  # ----------------------------------------------------------
  # Sex coding
  # ----------------------------------------------------------

  sex_model <- case_when(
    sex %in% c("Male", "M", "male", "m") ~ "M",
    sex %in% c("Female", "F", "female", "f") ~ "F",
    TRUE ~ NA_character_
  )

  # ----------------------------------------------------------
  # New patient
  # ----------------------------------------------------------

  newdata <- data.frame(
    sex = factor(
      sex_model,
      levels = fit_dev_final$xlevels$sex
    ),
    age_y1 = age,
    HBA1_y1 = HbA1c,
    eGFR_y1 = eGFR,
    TG_y1_imputed = triglycerides,
    CHOL_y1_imputed = cholesterol
  )

  if (anyNA(newdata)) {
    stop("One or more input variables are invalid.")
  }

  # ----------------------------------------------------------
  # Prediction (log-log CI: stays within [0,1] on the survival
  # scale, standard choice for survival probabilities)
  # ----------------------------------------------------------

  sf <- survfit(
    fit_dev_final,
    newdata = newdata,
    conf.type = "log-log"
  )

  times_days <- times * 365.25

  surv <- summary(
    sf,
    times = times_days,
    extend = TRUE
  )

  predicted_risk <- 1 - surv$surv
  # survival CI flips to become the risk CI: the upper bound on
  # survival is the lower bound on risk, and vice versa
  risk_ci_lower <- 1 - surv$upper
  risk_ci_upper <- 1 - surv$lower

  data.frame(
    time_years = times,
    predicted_risk = predicted_risk,
    predicted_risk_percent =
      round(100 * predicted_risk, 1),
    ci_lower_percent =
      round(100 * risk_ci_lower, 1),
    ci_upper_percent =
      round(100 * risk_ci_upper, 1)
  )
}


# ============================================================
# 2b. SURVIVAL CURVE FUNCTION (for plotting)
# ============================================================
# Shares the same validation and newdata construction as renal_risk()
# above, but returns a fine time grid of predicted survival
# probabilities instead of just two timepoints, for plotting the
# patient's individualized predicted survival curve. Note this is
# the *model's predicted* survival curve for this covariate profile,
# not an empirical Kaplan-Meier curve (KM curves are nonparametric
# cohort-level estimates, not something a Cox model produces for a
# single hypothetical patient) - labelled accordingly in the plot.

renal_survival_curve <- function(
    sex,
    age,
    HbA1c,
    eGFR,
    triglycerides,
    cholesterol,
    max_years = 10,
    n_points = 100
) {

  if (!sex %in% c("Male", "Female", "M", "F", "male", "female", "m", "f")) {
    stop("sex must be Male or Female.")
  }

  if (age < 18 || age > 100) {
    stop("Age is outside the supported range.")
  }

  if (HbA1c < 4 || HbA1c > 12) {
    stop("HbA1c is outside the supported range.")
  }

  if (eGFR < 60 || eGFR > 150) {
    stop("This calculator is validated for eGFR >=60 mL/min/1.73 m2 (up to 150).")
  }

  if (triglycerides <= 0) {
    stop("Triglycerides must be > 0.")
  }

  if (cholesterol <= 0) {
    stop("Cholesterol must be > 0.")
  }

  sex_model <- case_when(
    sex %in% c("Male", "M", "male", "m") ~ "M",
    sex %in% c("Female", "F", "female", "f") ~ "F",
    TRUE ~ NA_character_
  )

  newdata <- data.frame(
    sex = factor(
      sex_model,
      levels = fit_dev_final$xlevels$sex
    ),
    age_y1 = age,
    HBA1_y1 = HbA1c,
    eGFR_y1 = eGFR,
    TG_y1_imputed = triglycerides,
    CHOL_y1_imputed = cholesterol
  )

  if (anyNA(newdata)) {
    stop("One or more input variables are invalid.")
  }

  sf <- survfit(
    fit_dev_final,
    newdata = newdata,
    conf.type = "log-log"
  )

  grid_years <- seq(0, max_years, length.out = n_points)
  grid_days <- grid_years * 365.25

  surv <- summary(
    sf,
    times = grid_days,
    extend = TRUE
  )

  data.frame(
    time_years = grid_years,
    risk = 1 - surv$surv,
    ci_lower = 1 - surv$upper,
    ci_upper = 1 - surv$lower
  )
}


# ============================================================
# 3. UI
# ============================================================

ui <- fluidPage(

  tags$head(

    tags$title(
      "Renal Risk Calculator"
    ),

    tags$style(
      HTML("

      body {
        background: #f4f7fa;
        font-family: Arial, sans-serif;
        color: #263238;
      }

      .header {
        background: #173f5f;
        color: white;
        padding: 30px;
        margin-bottom: 30px;
      }

      .header h1 {
        margin: 0;
        font-size: 30px;
      }

      .header p {
        margin-top: 8px;
        margin-bottom: 0;
        opacity: 0.9;
      }

      .card {
        background: white;
        padding: 25px;
        margin-bottom: 20px;
        border-radius: 12px;
        border: 1px solid #e0e6ea;
        box-shadow: 0 2px 10px rgba(0,0,0,0.06);
      }

      .card-title {
        color: #173f5f;
        font-size: 20px;
        font-weight: bold;
        margin-top: 0;
        margin-bottom: 22px;
      }

      .risk-card {
        text-align: center;
        background: white;
        padding: 25px;
        border-radius: 12px;
        border: 1px solid #e0e6ea;
        box-shadow: 0 2px 10px rgba(0,0,0,0.06);
        margin-bottom: 20px;
      }

      .risk-label {
        color: #607d8b;
        font-size: 16px;
        font-weight: bold;
      }

      .risk-value {
        color: darkred;
        font-size: 48px;
        font-weight: bold;
        margin-top: 12px;
      }

      .risk-subtitle {
        color: #90a4ae;
        font-size: 13px;
      }

      .interpretation {
        background: #eef5fb;
        border-left: 4px solid #20639b;
        padding: 16px;
        border-radius: 5px;
        line-height: 1.6;
        margin-top: 20px;
      }

      .note {
        color: #78909c;
        font-size: 12px;
        line-height: 1.6;
      }

      .btn-calc {
        width: 100%;
        background: #20639b;
        color: white;
        border: none;
        border-radius: 8px;
        font-size: 17px;
        font-weight: bold;
        padding: 13px;
      }

      .btn-calc:hover {
        background: #173f5f;
        color: white;
      }

      .footer {
        text-align: center;
        color: #90a4ae;
        font-size: 12px;
        padding: 25px;
      }

      ")
    )
  ),


  # ==========================================================
  # HEADER
  # ==========================================================

  div(
    class = "header",

    div(
      class = "container",

      h1(
        "Renal Risk Calculator"
      ),

      p(
        "Individualized prediction of the risk for eGFR <60 ml/min/1.73 m² at 5 and 10 years",
        "in patients with eGFR ≥60 ml/min/1.73 m²"
      )

    )
  ),


  # ==========================================================
  # MAIN
  # ==========================================================

  div(
    class = "container",

    fluidRow(

      # ========================================================
      # INPUT PANEL
      # ========================================================

      column(
        width = 5,

        div(
          class = "card",

          h2(
            class = "card-title",
            "Patient characteristics"
          ),


          radioButtons(
            inputId = "sex",

            label = "Sex",

            choices = c(
              Male = "Male",
              Female = "Female"
            ),

            selected = "Male",

            inline = TRUE
          ),


          numericInput(
            "age",
            "Age (years)",
            value = 55,
            min = 18,
            max = 100,
            step = 1
          ),


          numericInput(
            "HbA1c",
            "HbA1c (%)",
            value = 6.5,
            min = 4,
            max = 12,
            step = 0.1
          ),


          numericInput(
            "eGFR",
            "eGFR (mL/min/1.73 m²)",
            value = 80,
            min = 60,
            max = 150,
            step = 1
          ),

          div(
            class = "note",
            style = "margin-top:-10px; margin-bottom:15px;",
            "This calculator was developed in patients with eGFR >=60 ",
            "mL/min/1.73 m². Predictions below this range are not validated."
          ),


          radioButtons(
            inputId = "lipid_unit",

            label = "Triglyceride / cholesterol units",

            choices = c(
              "mg/dL" = "mgdl",
              "mmol/L" = "mmol"
            ),

            selected = "mgdl",

            inline = TRUE
          ),


          numericInput(
            "triglycerides",
            "Triglycerides (mg/dL)",
            value = 150,
            min = 0.1,
            max = 5000,
            step = 1
          ),


          numericInput(
            "cholesterol",
            "Total cholesterol (mg/dL)",
            value = 200,
            min = 0.1,
            max = 1000,
            step = 1
          ),


          br(),


          actionButton(
            "calculate",
            "Calculate risk",
            class = "btn-calc"
          )

        )

      ),


      # ========================================================
      # RESULTS
      # ========================================================

      column(
        width = 7,

        div(
          class = "card",

          h2(
            class = "card-title",
            "Predicted risk of eGFR < 60 ml/min/1.73 m²"
          ),


          fluidRow(

            column(
              width = 6,

              div(
                class = "risk-card",

                div(
                  class = "risk-label",
                  "5-year risk"
                ),

                div(
                  class = "risk-value",

                  textOutput(
                    "risk5",
                    inline = TRUE
                  )

                ),

                div(
                  class = "risk-subtitle",
                  textOutput(
                    "risk5_ci",
                    inline = TRUE
                  )
                ),

                div(
                  class = "risk-subtitle",
                  "Probability of eGFR < 60 ml/min/1.73m²"
                )

              )

            ),


            column(
              width = 6,

              div(
                class = "risk-card",

                div(
                  class = "risk-label",
                  "10-year risk"
                ),

                div(
                  class = "risk-value",

                  textOutput(
                    "risk10",
                    inline = TRUE
                  )

                ),

                div(
                  class = "risk-subtitle",
                  textOutput(
                    "risk10_ci",
                    inline = TRUE
                  )
                ),

                div(
                  class = "risk-subtitle",
                  "Probability of eGFR < 60 ml/min/1.73m²"
                )

              )

            )

          )

        ),


        # ======================================================
        # PREDICTED SURVIVAL CURVE
        # ======================================================

        div(
          class = "card",

          h2(
            class = "card-title",
            "Risk curve"
          ),

          plotOutput(
            "survival_curve",
            height = "320px"
          ),

          div(
            class = "note",
            style = "margin-top:12px;",
            "This curve shows the model's predicted cumulative probability ",
            "that this patient's eGFR falls below 60 mL/min/1.73 m\u00b2 over time, ",
            "with the shaded band showing the 95% confidence interval. It is ",
            "the model's individualized prediction, not an empirical ",
            "Kaplan-Meier curve estimated directly from cohort data."
          )

        ),


        # ======================================================
        # MODEL INFORMATION
        # ======================================================

        div(
          class = "card",

          h2(
            class = "card-title",
            "Model information"
          ),

          p(
            "The calculator is based on the final Cox proportional ",
            "hazards model developed jointly on pooled Freiburg and ",
            "Zurich cohort data. ",
            a(
              "Original publication",
              href = "https://onlinelibrary.wiley.com/doi/10.1111/joim.13736",
              target = "_blank"
            ),
            " | ",
            a(
              "Multicenter validation publication",
              href = "SECOND_PUBLICATION_URL",
              target = "_blank"
            ), 
            " | ",
            a(
              "Source code",
              href = "https://github.com/weltthom/Renal_Risk_Calculator",
              target = "_blank"
            ),
          ),

          p(
            "The model was developed on pooled Freiburg and Zurich data ",
            "(n = 27,059; 7,818 events); ",
            "see the linked publications for development and validation details."
          ),

          h3(
            style = "color:#173f5f; font-size:16px; margin-top:22px; margin-bottom:10px;",
            "How this calculator works"
          ),

          p(
            class = "note",
            style = "font-size:13px; line-height:1.7; color:#455a64;",
            "This tool implements a Cox proportional hazards model relating ",
            "sex, age, HbA1c, eGFR, triglycerides, and total cholesterol at ",
            "baseline to the hazard of eGFR falling below 60 mL/min/1.73 m². ",
            "Entering a patient's baseline values returns the model's estimated ",
            "probability of crossing that threshold within 5 and 10 years, ",
            "based on the fitted hazard ratios and baseline survival function ",
            "from the pooled development cohort."
          ),

          p(
            class = "note",
            style = "font-size:13px; line-height:1.7; color:#455a64;",
            "The calculator is restricted to patients with eGFR \u226560 ",
            "mL/min/1.73 m\u00b2 at baseline, matching the population in which ",
            "the model was developed and validated; predictions outside the ",
            "input ranges shown are not supported. This tool is intended to ",
            "illustrate the published model and to support clinical judgment, ",
            "not to replace it \u2014 individual predictions should be interpreted ",
            "alongside a patient's full clinical picture."
          )

        )

      )

    ),


    div(
      class = "footer",

      "Renal Risk Calculator · Freiburg + Zurich pooled model"

    )

  )

)


# ============================================================
# 4. SERVER
# ============================================================

server <- function(
    input,
    output,
    session
) {

  # ==========================================================
  # LIPID UNIT TOGGLE
  # ==========================================================
  # Converts displayed values/bounds when the user switches units.
  # Actual mg/dL conversion for prediction happens separately below,
  # at calculate time, using whatever unit is currently selected.

  previous_lipid_unit <- reactiveVal("mgdl")

  observeEvent(input$lipid_unit, {

    old_unit <- isolate(previous_lipid_unit())
    new_unit <- input$lipid_unit

    if (identical(old_unit, new_unit)) {
      return()
    }

    if (new_unit == "mmol" && old_unit == "mgdl") {

      updateNumericInput(
        session, "triglycerides",
        label = "Triglycerides (mmol/L)",
        value = round(input$triglycerides / TG_MGDL_PER_MMOL, 2),
        min = 0.05, max = 30, step = 0.1
      )

      updateNumericInput(
        session, "cholesterol",
        label = "Total cholesterol (mmol/L)",
        value = round(input$cholesterol / CHOL_MGDL_PER_MMOL, 2),
        min = 0.05, max = 30, step = 0.1
      )

    } else if (new_unit == "mgdl" && old_unit == "mmol") {

      updateNumericInput(
        session, "triglycerides",
        label = "Triglycerides (mg/dL)",
        value = round(input$triglycerides * TG_MGDL_PER_MMOL, 1),
        min = 0.1, max = 5000, step = 1
      )

      updateNumericInput(
        session, "cholesterol",
        label = "Total cholesterol (mg/dL)",
        value = round(input$cholesterol * CHOL_MGDL_PER_MMOL, 1),
        min = 0.1, max = 1000, step = 1
      )

    }

    previous_lipid_unit(new_unit)

  })


  # ==========================================================
  # CALCULATION
  # ==========================================================
  # renal_risk() now does its own eGFR >=60 validation internally
  # (via stop()), so this just needs to catch that error rather
  # than pre-checking eGFR itself, avoiding the previous duplicate
  # ">=60" logic living in two places at once.

  result <- eventReactive(

    input$calculate,

    {

      tryCatch(
        {
          # Convert lipid inputs to mg/dL (what the model was fit on)
          # regardless of which unit the user selected
          tg_mgdl <- if (input$lipid_unit == "mmol") {
            input$triglycerides * TG_MGDL_PER_MMOL
          } else {
            input$triglycerides
          }

          chol_mgdl <- if (input$lipid_unit == "mmol") {
            input$cholesterol * CHOL_MGDL_PER_MMOL
          } else {
            input$cholesterol
          }

          out <- renal_risk(
            sex = input$sex,
            age = input$age,
            HbA1c = input$HbA1c,
            eGFR = input$eGFR,
            triglycerides = tg_mgdl,
            cholesterol = chol_mgdl,
            times = c(5, 10)
          )

          curve <- renal_survival_curve(
            sex = input$sex,
            age = input$age,
            HbA1c = input$HbA1c,
            eGFR = input$eGFR,
            triglycerides = tg_mgdl,
            cholesterol = chol_mgdl
          )

          list(
            risk5 = out$predicted_risk[out$time_years == 5],
            risk10 = out$predicted_risk[out$time_years == 10],
            risk5_ci_lower = out$ci_lower_percent[out$time_years == 5],
            risk5_ci_upper = out$ci_upper_percent[out$time_years == 5],
            risk10_ci_lower = out$ci_lower_percent[out$time_years == 10],
            risk10_ci_upper = out$ci_upper_percent[out$time_years == 10],
            curve = curve
          )
        },

        error = function(e) {

          showNotification(
            conditionMessage(e),
            type = "error",
            duration = 8
          )

          NULL
        }

      )

    },

    ignoreInit = TRUE
  )


  # ==========================================================
  # 5 YEAR
  # ==========================================================

  output$risk5 <- renderText({

    r <- result()

    if (is.null(r)) {
      return("—")
    }

    sprintf(
      "%.1f%%",
      r$risk5 * 100
    )

  })


  output$risk5_ci <- renderText({

    r <- result()

    if (is.null(r)) {
      return("")
    }

    sprintf(
      "95%% CI: %.1f%%\u2013%.1f%%",
      r$risk5_ci_lower,
      r$risk5_ci_upper
    )

  })


  # ==========================================================
  # 10 YEAR
  # ==========================================================

  output$risk10 <- renderText({

    r <- result()

    if (is.null(r)) {
      return("—")
    }

    sprintf(
      "%.1f%%",
      r$risk10 * 100
    )

  })


  output$risk10_ci <- renderText({

    r <- result()

    if (is.null(r)) {
      return("")
    }

    sprintf(
      "95%% CI: %.1f%%\u2013%.1f%%",
      r$risk10_ci_lower,
      r$risk10_ci_upper
    )

  })


  # ==========================================================
  # SURVIVAL CURVE PLOT
  # ==========================================================

  output$survival_curve <- renderPlot({

    r <- result()

    if (is.null(r)) {
      return(NULL)
    }

    curve <- r$curve

    marker_points <- data.frame(
      time_years = c(5, 10),
      risk = c(r$risk5, r$risk10)
    )

    ggplot(curve, aes(x = time_years, y = risk)) +
      geom_ribbon(
        aes(ymin = ci_lower, ymax = ci_upper),
        fill = "#F08080",
        alpha = 0.15
      ) +
      geom_line(color = "darkred", linewidth = 1.1) +
      geom_point(
        data = marker_points,
        aes(x = time_years, y = risk),
        color = "darkred",
        size = 3
      ) +
      geom_text(
        data = marker_points,
        aes(
          x = time_years,
          y = risk,
          label = sprintf("%.1f%%", risk * 100)
        ),
        vjust = -1.1,
        color = "darkred",
        fontface = "bold",
        size = 4.2
      ) +
      geom_vline(
        xintercept = c(5, 10),
        linetype = "dashed",
        color = "#90a4ae"
      ) +
      scale_y_continuous(
        labels = function(x) paste0(x * 100, "%"),
        limits = c(0, 1)
      ) +
      scale_x_continuous(
        breaks = seq(0, 10, by = 2)
      ) +
      labs(
        x = "Years since follow up",
        y = "Probability of eGFR < 60 ml/min/1.72 m²"
      ) +
      theme_classic(base_size = 13) + 
      theme(
        panel.grid.minor = element_blank(),
        axis.title = element_text(color = "#455a64"),
        axis.text = element_text(color = "#607d8b")
      )

  })


  # ==========================================================
  # INTERPRETATION
  # ==========================================================

  output$interpretation <- renderText({

    r <- result()

    if (is.null(r)) {

      return(
        "Enter the patient characteristics and click 'Calculate risk'."
      )

    }

    paste0(

      "For this patient, the model estimates a ",

      sprintf(
        "%.1f%%",
        r$risk5 * 100
      ),

      " probability of eGFR <60 ml/min/1.73 m² within 5 years and a ",

      sprintf(
        "%.1f%%",
        r$risk10 * 100
      ),

      " probability within 10 years."

    )

  })

}


# ============================================================
# 5. RUN
# ============================================================

shinyApp(
  ui = ui,
  server = server
)

