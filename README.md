# SSRT Calculation

SSRTs can either be calculated:

- Using a Mean estimation (M<sub>RT</sub> - M<sub>SSD</sub>)
- Utilising the RT distribution of correct Go Trials to estimate the Stop-Signal Reaction Time.

---

Since the 2nd approach uses an integrative procedure that estimates the SSRT
from the RT distribution using p(response | signal); in cases where p(response
| signal) == 50%, both estimates should produce similar results.

(see Verbruggen et al. 2019)

---

## Installation

devtools::install_github("sonntasa/SSRT")
