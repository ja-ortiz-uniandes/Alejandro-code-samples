# Village sampling simulation  <!-- omit in toc --> 

- [Contact](#contact)
- [How to view the report](#How-to-view-the-fial-report)
- [Known issues](#known-issues)
- [License](#license)

 Thank you for looking at my simulation code sample in `R`! This is a report for a project I thought of while doing a questionnaire for the World Bank's Development Impact Evaluation (DIME). The question was: "*There is a program that is implemented at the village level. Households within the same village are very similar but households between villages are not. To maximize the likelihood of detecting the programs effect is it better to sample more households within each village or to sample more villages?*"

This problem is closely related to the standard *frequentist* sample size formula for detecting a treatment effect in a clustered design:

$$
n = \frac{\sigma^2 \cdot (z_{1-\alpha/2} + z_{1-\beta})^2}{D^2} \cdot \frac{1 + \rho (m - 1)}{p(1 - p)}
$$

Where:
- $n$ — required total sample size
- $\sigma^2$ — variance of the outcome variable
- $D$ — minimum detectable effect size
- $z_{1-\alpha/2}$, $z_{1-\beta}$ — critical values from the standard normal distribution for the chosen significance level $\alpha$ and power $1 - \beta$
- $\rho$ — intraclass correlation coefficient (ICC), measuring similarity of observations within the same cluster
- $m$ — number of observations per cluster
- $p$ — proportion of the sample assigned to the treatment group (so $1-p$ is the proportion in the control group)

My simulation visualices this logic to compare trade-offs between increasing the number of clusters (villages) versus the number of observations per cluster (households).

In this directory you will find a fully replicable and portable project in R that uses simulation techniques to answer this question. I decided to do this project as I believe is a good demonstration of my skills in `R`. It covers function writing, data management, visualization, parallelization, memory usage, error handling, local vs. global variables, planning for edge cases, use of many different data types from vectors to data.tables through matrices and lists; and brings it all together using a real-world question about sampling techniques. It demonstrates the use of many different and popular packages like `tidyverse`, `data.table`, `furrr`, and `plotly`. Lastly, this report itself shows the use of Rmarkdown with HTML, CSS, and YAML elements mixed in.

However, it's worth saying that this is only a sample of my abilities in R. If there is any particular demonstration you would like to see, feel free to contact me.

## How to view the final report
Please download the `R-simualtion-report.html` file, then open it with your browser.

## Contact

Please don't hesitate to write me an email if you find any issues, problems, or if you have any comments, suggestions, or questions. You can find my email on my GitHub profile.

## Known issues

None at the time of writing. If you know of any please create an issue or email me directly.

## License

[![CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND-lightgrey)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

This work is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](https://creativecommons.org/licenses/by-nc-nd/4.0/).

[![CC BY-NC-ND 4.0](https://licensebuttons.net/l/by-nc-nd/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc-nd/4.0/)
