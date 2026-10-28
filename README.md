# SSRT Calculation

SSRTs can either be calculated:

- Using a Mean estimation (M<sub>RT</sub> - M<sub>SSD</sub>)
- Utilising the RT distribution of correct Go Trials to estimate the Stop-Signal Reaction Time.

---

Since the 2nd approach uses an integrative procedure that estimates the SSRT
from the RT distribution using p(response | signal); in cases where p(response
| signal) == 50%, both estimates should produce identical results.

(see Verbruggen et al. 2019)

---

# SSD Calculation

The Stop Signal Delay can be calculated as either the mean Value of all SSDs, or, more robustly,
using the mid-run estimates method which is more suited to Staircase paradigms. This method
computes the mean value by averaging the midpoints of every second run (see Levitt, 1970).

## Installation

pak::pak("sonntasa/sstFun")
