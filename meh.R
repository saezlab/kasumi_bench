library(readr)
library(pROC)
library(ggplot2)

# MIBI  ----
ct <- read_rds("rocs/dcis.ct.rds")
expr <- read_rds("rocs/dcis.expr.rds")
ws <- read_rds("rocs/dcis.ws.rds")

ggroc(ct, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("MIBI DCIS") + theme_classic()

ggsave("roc.mibi.pdf")

ggroc(list(sm.ct = ct$sm, ws.ct = ws$ct, sm.expr = expr, ws.expr = ws$expr), legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("MIBI DCIS") + theme_classic()

ggsave("roc.misty.mibi.pdf")

# CODEX ----

ct <- read_rds("rocs/ctcl.ct.rds")
expr <- read_rds("rocs/ctcl.expr.rds")
ws <- read_rds("rocs/ctcl.ws.rds")

ggroc(ct, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("CODEX CTCL") + theme_classic()

ggsave("roc.codex.pdf")

ggroc(list(sm.ct = ct$sm, ws.ct = ws$ct, sm.expr = expr, ws.expr = ws$expr), legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("CODEX CTCL") + theme_classic()

ggsave("roc.misty.codex.pdf")

# IMC ----

ct <- read_rds("rocs/bc.ct.rds")
expr <- read_rds("rocs/bc.expr.rds")
ws <- read_rds("rocs/bc.ws.rds")

ggroc(ct, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("IMC BC") + theme_classic()

ggsave("roc.imc.pdf")

ggroc(list(sm.ct = ct$sm, ws.ct = ws$ct, sm.expr = expr, ws.expr = ws$expr), legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  labs(color = NULL) +
  ggtitle("IMC BC") + theme_classic()

ggsave("roc.misty.imc.pdf")
