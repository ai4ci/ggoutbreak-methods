CC: tkl@lanl.gov, "Leon Danon" leon.danon@bristol.ac.uk

PCOMPBIOL-D-26-00077
Instantaneous Reproduction Number Estimation From Modelled Incidence
PLOS Computational Biology

Dear Dr. Challen,

Thank you for submitting your manuscript to PLOS Computational Biology. As with all papers, your manuscript was reviewed by members of the editorial board. Based on our assessment, we have decided that the work does not meet our criteria for publication and will therefore be rejected. If external reviews were secured, reviewers' comments will be included at the bottom of this email.

We are sorry that we cannot be more positive on this occasion. We very much appreciate your wish to present your work in one of PLOS's Open Access publications. Thank you for your support, and we hope that you will consider PLOS Computational Biology for other submissions in the future.

Yours sincerely,

Sasikiran Kandula
Academic Editor
PLOS Computational Biology

Thomas Leitner
Section Editor
PLOS Computational Biology


Additional Editor Comments (if provided):

The manuscript by Drs. Challen and Danon has been reviewed by three reviewers. I believe the concerns raised by reviewers 1 and 3 are substantial, and if meaningfully addressed, would result in a distinctly different manuscript that might need to be evaluated separately. In addition to the more theoretical issues, a stronger justification for a comparison with EpiEstim (Reviewer 1, 1e) as well as support for the proposed method's advance relative to EpiNow2 (Reviewer 3) should be considered.

All reviewers note that the authors provided incorrect link(s) (for example: https://ai4ci.github.io/ggoutbreak-paper/) in the data availability section, which likely prevented replicability checks.


[Note: HTML markup is below. Please do not edit.]

Reviewers' Comments (if peer reviewed):

Reviewer's Responses to Questions

Comments to the Authors:
Please note here if the review is uploaded as an attachment.

Reviewer #1: See PDF

Reviewer #2: This is a nice paper, but needs some minor refinements throughout.

1. Introduction, lines 12-15: this is a rather confusing definition of the instantaneous reproduction number, and I am not convinced it is a correct one. Consider providing a more standard definition.
2. Methods, lines 60-65: the jump between equations (1) and (2) needs more explanation. In particular, equation (2) is described as holding “in expectation”, but the expression in equation (2) is a random variable. In addition, is some uncertainty lost in any averaging that does take place, and how do we know approximating the expectation of a quotient by the quotient of expectations is reasonable here?
3. Methods, lines 97-102: I am unsure about the validity of the method in the scenario with negative serial intervals, given that the renewal equation itself (working directly with symptomatic case counts) doesn’t make sense in that case. I would recommend removing these claims (similarly with the related text in the Discussion, lines 298-303).
4. Methods, line 120: specify that the infectivity profile estimates are equally weighted.
5. Methods, “Validation”, and Supplementary 1: is this branching process model not entirely equivalent to the renewal equation (which would be easier to simulate directly)?
6. Methods, lines 142-144 (“In each simulation run, the degree of outward edges in the network of realised infections at any given time is an instantaneous Rt”): I again found this a confusing and possibly incorrect characterisation of R_t.
7. Results, Fig 1: does panel A show case counts for the GAM or Locfit model (not currently stated in the caption)? It may actually be nice to show both (if that doesn’t break journal figure size requirements).
8. Results, caption to Fig 2B: I would specify that a bias closer to zero is better (since more negative could be considered “lower”).
9. Discussion, lines 322-324: this information would be useful earlier in the manuscript – when I first saw the figures, I was concerned about a possible impact of bad initial estimates for EpiEstim and Locfit on Fig 2. As an aside, I wondered if the early Locfit estimates would be substantially improved by only fitting to the data from day 2 onwards (since the model seemingly can’t capture the somewhat artificial drop-off in incidence between days 1 and 2 caused by initiating the model with many cases on a single day).
10. Supplementary 1, section 2 (“Simulation”): notwithstanding point 5 above, I didn’t understand the notation X_j,i,t (specifically, what the subscript j refers to).
11. Supplementary 1, section 3 (“Validation”): “The minimum root mean squared error between prediction and known R_t is taken as the estimate lag”. Presumably this should be the delay leading to the minimum error?
12. Supplementary 1, section 3 (“Validation”): was “the replicate resulting in the highest and lowest R_t estimates for each scenario” selected based on the highest/lowest estimates across time points for each replicate, or based on the average of R_t across time points, or was a different replicate selected at each time point, or something else?
13. Supplementary 1, Fig 4: the x-axis label is given as “Dispersion”, which I think is an error. Also, the caption does not refer to panel F.
14. Supplementary 2, equation (2): should this be g^{-1}() rather than g’()?
15. Supplementary 2, section 3 (“Conclusion”): the first line here (“In this paper…”) reads as though this is a standalone piece rather than a supplement.
16. Supplementary 3, section 2 (“Bias from assuming independence”): where does the lower bound in the second line of equation (9) come from? I found this section rather hard to follow in general, and more explanation would be useful.
17. I noticed several typos: line 128 (“this this”), caption to Fig 1 (“esimtates”), caption to Fig 2 (“cirtical”), Supplementary 1 section 5 (“‘R_t Locfit’m”), Supplementary 2 section 1 (“properly of the estimator”), Supplementary 2 section 2 (“PIT values (also sometimes described as universal residuals [7]) u_i”). Please check carefully throughout.
18. As far as I can tell, the linked GitHub repository doesn’t exist, and the Zenodo doi links to a different project.

Reviewer #3: In this paper, Challen and colleagues propose a statistical method for estimation of time-varying instantaneous reproduction numbers from a time series of infection incidence estimates. This method is meant as a downstream analysis tool that can be coupled with a statistical method of choice for estimating incidence curves to obtain reproduction number estimates. The key rationale of this approach is to improve reproduction number estimation when case observations do not represent the infection incidence curve well (e.g. due to biological/reporting delays). In these cases, classical methods such as by Cori et al. are misspecified. Instead, the authors propose to 1) estimate the underlying incidence curve with a statistical method that accounts for the observation process and 2) subsequently plug in these estimates into the method proposed here to estimate Rt.

I find the general motivation of the problem addressed in this paper well justified and relevant for epidemiological practice. However, the authors do not discuss an obvious and attractive alternative to their framework that is already used in practice: joint estimation of reproduction numbers and incidence. For example, the authors have mentioned Abbott et al., whose model (EpiNow2) indeed provides infection incidence estimates that could be used as an input to the method presented here. However, EpiNow2 also readily provides fully Bayesian Rt estimates itself. I think a discussion of how the presented method compares to this alternative type of approach would be needed. I expect that lower complexity/computational burden of estimation could be a potential argument (of course this will mostly depend on the method used for incidence estimation).

In the following, I have further remarks on the methodology and discussion.

I find that the core assumption of log-normality of the estimated Poisson rates needs to be justified in more detail. First, it should be clarified that this assumption holds only asymptotically for sufficiently large I_t, and the underlying theory on asymptotic normality of MLE and Bayesian posterior distributions (under smooth, non-degenerate priors) for model parameters could be shortly described. Second, without wanting to be overly pedantic, this paper makes a bit of a philosophical leap from frequentist to Bayesian interpretations when modeling rate estimates from MLE as random variables to estimate R_t. From a frequentist perspective, \lambda_t in Eq (2) is a fixed ground truth value. Using a confidence distribution for \lambda_t based on MLE assumes that the estimation error is consistent with the distribution of the true rates in the data generating process. Also, the standard error of \hat{\lambda_t} itself is estimated but this uncertainty is not propagated here. Overall, my impression is that the approach presented here conceptually fits better to a Bayesian inference setting where the posterior for \lambda_t is approximated by a (multivariate) normal. I am not arguing that the use of MLE estimates is not pragmatically acceptable with the present method, but the theoretical assumptions implied by doing so could be developed a bit more.

The authors highlight decoupling of incidence modelling and Rt estimation as a main feature of their approach, arguing that this allows for correction of biases. I agree, however this decoupling can also be a source of bias. Specifically, if incidence estimation methods use non-parametric smoothing (e.g. via P-splines, ARIMA etc.), this can introduce additional assumptions about the shape of the incidence curve and therefore about transmission dynamics. A nice feature of EpiEstim, with all its limitations, is that it makes a very explicit assumption about the smoothness of Rt in terms of a smoothing window. The method presented here offers more flexibility w.r.t. to the smoothing approach used, but there is also a risk that the modeled smoothness of the incidence curve has "hidden" implications for the smoothness of Rt (and, in my experience, these effects can be quite strong). I therefore think it is important to highlight this risk to potential users of this method.

Related to the above point, one practical question is: how does the smoothing hyperparameter used in EpiEstim (window size) relate to smoothing hyperparameters in the methods used to estimate the incidence curve? For example, the authors have here used 14 days both for the EpiEstim window size as well as for the hyperparameters in the incidence method, i.e. local polynomial bandwith (locfit) or the knot distance (GAM). But are these parameters comparable, do they imply similar smoothness? I believe this could also have consequences for the performance of the different methods in simulation-based validation study.

Also related to the above point, but more conceptually: in a step-wise procedure, there will probably be a contradiction between i) the first step (estimation of a smooth incidence curve), which typically assumes that incidence evolves according to some general state space model (e.g. ARIMA) or follows a higher-order polynomial and ii) the second step (Rt estimation), which assumes that incidence evolves according to a renewal process. This is of course mostly a conceptual disagreement, but could lead to practical problems if the method used in step 1 significantly contradicts assumptions of the renewal process.

I appreciate that different degrees of ascertainment noise were tested in the simulation analysis. To account for overdispersion in case reporting, many methods also model case counts as negative binomial distributed instead of Poisson distributed. I am curious whether the method proposed here can also be used with such incidence models, and if not, it would be good to highlight this as a potential limitation. Another concern is overdispersion in the offspring distribution - not simply in the observation process - which can bias Rt estimates by EpiEstim (see Němcová et al. https://doi.org/10.1101/2025.07.31.25332479). It would be good to comment on this issue as well.

One important application of Rt estimation is real-time disease monitoring, but this scenario (real-time Rt estimation) is not studied here. However, this is exactly the setting in which uncertainty quantification and smoothing assumptions are most important. I think it would be very helpful to include a small simulation study to assess the performance of the method in a real-time setting with some kind of epidemiological delay (e.g. incubation period for symptom onset data or reporting delay for case data). On the one hand, this is a setting in which the strength of this method compared to the default EpiEstim could really be demonstrated. If the model used for infection incidence curve estimation can quantify uncertainty near the present due to reporting delays, then I expect the Rt estimation method proposed here to also reflect this uncertainty in Rt. It would be very valuable to demonstrate this strength. On the other hand, a real-time setting would also nicely highlight potential pitfalls of the approach presented here, e.g. false extrapolation of the incidence curve by polynomial smoothing.

The authors have mentioned that the method also works with infectivity profiles with negative support, but I could not find a demonstration of this. Could an example of this be included in the appendix?

It is reassuring that the optimal lags for the Rt estimates under the author's method is 0. For EpiEstim, the lag of half the smoothing window size is already a well known characteristic (see also Gostic et al.). However, I guess for the approach proposed here there should be no lag as long as the incidence estimation method is temporally accurate? If not, how could the optimal lag even be identified without a ground truth?

It is great that the authors have tested different metrics as part of their method validation. I don't have the expertise to assess the faithfulness of the "novel" PIT-based calibration metric proposed. PIT histograms have been around for a long time and I wonder to what extend this or similar metrics derived from PIT histograms might already exist (or if they don't exist, maybe for good reasons). The simulation results in the appendix appear promising, but typically a formal analysis of the characteristics of a metric is needed before it can be generally trusted. I am not demanding the authors go down that route, but if this metric is here introduced as novel, a more detailed discussion on previous work around PIT histograms and the need for future, more extensive study of this metric would be important. Also, it might be helpful to still show examples of actual PIT histograms from the simulation-based comparison of methods, as readers will be used to interpreting these instead of some derived metric.
## Minor points

Eq (1): The specification here seems to imply that the \lambda_t are independently log-normal distributed, whereas the remainder of the paper rather seems to assume that they are multivariate log-normally distributed. This could be made more clear here.

Line 25: Note that EpiEstim has also been extended to account for some of these issues, such as data aggregation (using EM-style methods). See Nash, Rebecca K., Samir Bhatt, Anne Cori, und Pierre Nouvellet. 2023. „Estimating the Epidemic Reproduction Number from Temporally Aggregated Incidence Data: A Statistical Modelling Approach and Software Tool“. _PLOS Computational Biology_ 19 (8): e1011439. [https://doi.org/10.1371/journal.pcbi.1011439](https://doi.org/10.1371/journal.pcbi.1011439).

Line 97: symptomatic case counts. Does this mean case counts by date of symptom onset?

Validation/simulation study, line 142: "...the degree of outward edges in the network of realised infections at any given time is an instantaneous Rt." I would disagree, since Rt described only the expected outward degree, not the realized one. In my understanding, this means that the "simulation reproduction number" is not only close to the parameterised Rt but identical to it.

Have the authors made all data and (if applicable) computational code underlying the findings in their manuscript fully available?
The PLOS Data policy requires authors to make all data and code underlying the findings described in their manuscript fully available without restriction, with rare exception (please refer to the Data Availability Statement in the manuscript PDF file). The data and code should be provided as part of the manuscript or its supporting information, or deposited to a public repository. For example, in addition to summary statistics, the data points behind means, medians and variance measures should be available. If there are restrictions on publicly sharing data or code —e.g. participant privacy or use of data from a third party—those must be specified.

Reviewer #1: No: The link https://ai4ci.github.io/ggoutbreak-paper/ is dead.

Reviewer #2: No: It seems the authors have provided the wrong links in error

Reviewer #3: No: I believe the authors intend to make all data and code available, however the link to the study github repository provided in the mansuscript does not resolve (the repository might be private?). Moreover, the provided zenodo DOI seems to resolve to a different project not related to this paper. These issues need to be fixed.

PLOS authors have the option to publish the peer review history of their article (what does this mean?). If published, this will include your full peer review and any attached files.


If you choose “no”, your identity will remain anonymous but your review may still be made public.


Do you want your identity to be public for this peer review? For information about this choice, including consent withdrawal, please see our Privacy Policy.

Reviewer #1: No

Reviewer #2: No

Reviewer #3: No




[NOTE: If reviewer comments were submitted as an attachment file, they will be attached to this email and accessible via the submission site. Please log into your account, locate the manuscript record, and check for the action link "View Attachments". If this link does not appear, there are no attachment files.]

