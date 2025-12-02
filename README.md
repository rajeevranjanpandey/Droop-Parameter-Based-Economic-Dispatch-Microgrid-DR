# Droop Parameter-Based Economic Dispatch for Islanded Microgrid Considering Demand Response

This repository contains the implementation, dataset, and reproducible results for the work:

> **Droop Parameter-Based Economic Dispatch for Islanded Microgrid Considering Demand Response**  
> Authors: *Ganesh Kumar Pandey*, *Debapriya Das*  
> International Journal of Ambient Energy, 2024  
> DOI: 10.1080/01430750.2024.2305326

---

## 🔍 Overview

This project optimises day-ahead economic dispatch in an **islanded AC microgrid** using:

| Module | Purpose |
|--------|---------|
| Droop-based DG control | Decentralised power sharing |
| Modified Newton–Raphson | MG load flow without slack bus |
| Particle Swarm Optimisation | Optimal droop + DR scheduling |
| Incentive-based Demand Response | Peak load curtailment |

Key benefits:
- Reduced operating cost
- Better fairness in DG sharing
- Peak demand trimming via DR
- Voltage & frequency within IEEE 1547.7 limits

---

## 🧠 Optimisation Problem

### 🎯 Objective

$$
\min
\sum_{t=1}^{24}
\left[
\sum_{i=1}^{N_c}
\left(
y_i^t - \gamma_i^t x_i^t
\right)
+
\frac{
\sum_{j=1}^{N_{\mathrm{DG}}}
\left(
\lambda_j^t - \bar{\lambda}^t
\right)^2
}{
N_{\mathrm{DG}} - 1
}
\right]
$$

Where:
- $\lambda_j^t$ = incremental cost of DG $j$  
- $x_i^t$ = curtailed power by customer $i$  
- $y_i^t$ = incentive paid to customer $i$  
- $\gamma_i^t$ = interruptibility factor  

---

### 🔒 Constraints

#### DG output limits
$$
P_{j}^{\min}
\le
P_{j}^t
\le
P_{j}^{\max}
$$

#### Droop coefficient limits
$$
m_{p,j}^{\min}
\le
m_{p,j}^t
\le
m_{p,j}^{\max}
\quad,\quad
n_{q,j}^{\min}
\le
n_{q,j}^t
\le
n_{q,j}^{\max}
$$

#### Power balance
$$
\sum_{j=1}^{N_{\DG}}
P_{DG,j}^t
=
P_{\load}^t
-
P_{\DR}^t
+
P_{\loss}^t
$$

#### Demand Response Feasibility

**Individual rationality**
$$
y_i^t - C_i(\theta_i, x_i^t) \ge 0
$$

**Incentive compatibility**
$$
y_i^t - C_i(\theta_i, x_i^t)
\ge
y_{i-1}^t - C_{i-1}(\theta_i, x_{i-1}^t)
$$

#### Utility daily budget
$$
\sum_{t=1}^{24}
\sum_{i=1}^{N_c}
y_i^t
\le
U_B
$$

---

## 📂 Repository Structure

