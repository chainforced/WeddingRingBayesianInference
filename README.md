[![ru](https://img.shields.io/badge/lang-ru-red.svg)](https://github.com/chainforced/WeddingRingBayesianInference/blob/main/README.RU.md)


I'm using here bayesian inference to find expected probabilities of 2 inductive arguments:<br />
 **P(not married | no ring & my data) = 62%** <br />
 **P(married | ring & my data) = 92%** <br />
In other words, given the data I have, if one sees that a person has NO wedding ring, one should rationally expect that this person is NOT married with only 62% confidence level. <br />
And if one sees that a person has wedding ring, one should rationally expect that this person is married with 92% confidence level.<br />
Prior distributions are such, that prior expectations of these probabilities are: <br />
**P(not married | no ring) = 96%** <br />
**P(married | ring) = 99%** <br />
Hence, my analysis falsifies hypothesis, that the first inductive argument is accurate, but didn't change estimation of accuracy of the second argument much. <br />

Here's a little article of mine, describing what data I use and some basic logic of how it's connected to the probabilities I estimate
https://medium.com/@chainforced/the-absence-of-a-wedding-ring-is-an-uninformative-sign-of-singleness-of-random-people-fbf1f924088f
