---
layout: post
title: "Hierarchical Data Management and Visualization"
output: html_document
tags: [r, mixed models, random effects, nested data, data science, dplyr, ggplot2]
---

This blog has a few goals. First, you'll see how to simulate a nested data set using the assumptions of a linear mixed-effects model. Then, you'll learn about *R* packages that can help to summarize and visualize similar hierarchical data in fast, reproducible, and easily generlizable ways. Finally, the user can change the simulation parameters to visualize the emergent effects of variability at different hierarchical scales. 

## Data Generation

Working with hierarchical data can be a pain, but there is a suite of *R* packages in the `tidy` family that are unbelievably helpful. What I mean is that these packages will change your scientific life forever. For real. 

Hyperbole aside, let's start by generating some realistic, nested data. We'll be testing Bergmann's rule that body size increases with elevation (as a proxy for temperature) by sampling along different mountains. We'll sample multiple individuals of many species, representing many genera, sampled across multiple mountains. Lots of nestedness here. 

There's a lot of non-independence going on with this type of nested sampling. You'd expect body size to vary more among genera, than among species within a given genera, and there is probably variation among mountains. This all effects the intercept of body size. There might be similar non-independence in the relationship between body size and elevation (i.e. the slope). We'll code this type of non-independence as random effects on the interecept and slope, below. 

This blog is focused on managing and visualizing hierarchical data. In a later post, I'll use this same generative example to illustrate proper statstical models to handle all of this non-independence.  


{% highlight %}
## Sample sizes:
mount_N <- 10 # Number of mountains sampled
obs_N_perMount <- sample(c(100:300), size=mount_N, replace=T)
obs_N <- sum(obs_N_perMount)
# Note: the same species could be found on different mountains
genera_N <- 15 # Number of genera
species_N_perGenus <- sample(c(2:8), size=genera_N, replace=T) # Number of species per genus
species_N <- sum(species_N_perGenus)

## Create placeholders
# Make sure these are factors, which will be used later. 
Mount <- factor(rep(c(1:mount_N), times=obs_N_perMount))
Genera <- factor(sample(c(1:genera_N), size=obs_N, replace=T))
Species <- factor(sample(c(1:species_N), size=obs_N, replace=T))
  
## Now assign parameters that determine the relationship between elevation and body size
beta_mean <- 2.0 # Slope: Body size declines with elevation
alpha_mean <- 0.0 # Intercept: Mean body size (centered and scaled)

# Random effects on slope
species_sd_beta <- 0.5 
genera_sd_beta <- 2.5
mount_sd_beta <- 1.0

species_rand_beta <- rnorm(species_N, 0, species_sd_beta)
genera_rand_beta <- rnorm(genera_N, 0, genera_sd_beta)
mount_rand_beta <- rnorm(mount_N, 0, mount_sd_beta)

# Random effects on intercept
species_sd_alpha <- 2.5
genera_sd_alpha <- 10.0
mount_sd_alpha <- 2.5

species_rand_alpha <- rnorm(species_N, 0, species_sd_alpha)
genera_rand_alpha <- rnorm(genera_N, 0, genera_sd_alpha)
mount_rand_alpha <- rnorm(mount_N, 0, mount_sd_alpha)

## Elevation (centered and scaled)
Elevation <- rnorm(obs_N, 0, 1)

## Generate body weights
# Note: Because we generated the data above in a smart way, we can vectorize this process!
Weight <- vector(mode="numeric", length=obs_N)
# Add the intercept with random effects:
Weight <- alpha_mean + species_rand_alpha[Species] + genera_rand_alpha[Genera] + mount_rand_alpha[Mount]
# Add the slope with random effects:
Weight <- Weight + (beta_mean + species_rand_beta[Species] + genera_rand_beta[Genera] + mount_rand_beta[Mount]) * Elevation

# Put this into a big data frame:
berg <- data.frame(Mount, Genera, Species, Elevation, Weight)
head(berg)
{% endhighlight %}



{% highlight text %}
##   Mount Genera Species   Elevation    Weight
## 1     1      1      54  0.24974103 -1.104280
## 2     1      1      56 -1.41061910 -5.014356
## 3     1     15      29  0.49633511 11.954253
## 4     1     15      65  0.08548281  9.427705
## 5     1     10      27  0.85287980  6.672971
## 6     1     12       1 -0.64909766 -4.222328
{% endhighlight %}

Very crudely, let's look at the overall pattern in the data, across all species.

{% highlight %}
plot(berg$Weight ~ berg$Elevation, pch=20, xlab="Elevation", ylab="Weight")
{% endhighlight %}

<img src="/figs/2016-9-10-hierarchical-data-management-and-visualization/all_data_plot-1.png" title="center" alt="center" style="display: block; margin: auto;" />

Well, that's a lot of variation, which makes it hard to see a clear pattern. Does the inherent hierarchy in the data obscure the body size - elevation relationship? What are the patterns among mountains? Among genera? Among species? 

## Data Management and Synthesis

I'm going to focus on the `tidy` family of packages that can help us look at the raw data in more manageable chunks, and the `ggplot2` package

{% highlight r %}
library(tidyverse) # This launches all the best packages
{% endhighlight %}

Brad Boehmke does a great job highlighting the functions of the `dplyr` and `tidyr` packages in his [R publication](https://rpubs.com/bradleyboehmke/data_wrangling), and you should definitely read it. I'll just highlight a few useful functions relevant to this dataset. 

Here's a simple question: What's the average body weight on each mountain? We could write a `for`-loop that partitions the data frame and then calculates an average, but yuck. Instead, use `dplyr` and the `summarize()` function.

{% highlight r%}
berg %>% # This symbol pipes the result to the next function
  group_by(Mount) %>% # Essentially cluster all the data for each mountain
  summarize(Weight_avg = mean(Weight), Weight_sd = sd(Weight)) # Create new columns that summarize Weight
{% endhighlight %}



{% highlight text %}
## # A tibble: 10 × 3
##     Mount Weight_avg Weight_sd
##    <fctr>      <dbl>     <dbl>
## 1       1  4.9636721  10.56421
## 2       2  4.0290075  10.91626
## 3       3  2.1202840  10.74022
## 4       4  0.1002801  11.49470
## 5       5  2.8156702  11.54933
## 6       6  1.9584231  12.57281
## 7       7  4.3874073  12.08083
## 8       8  9.3074673  10.30616
## 9       9  5.8189615  11.38677
## 10     10  0.5453171  12.10775
{% endhighlight %}

What about the mean weights for genera on different mountains? 

{% highlight r%}
berg %>% 
  group_by(Mount, Genera) %>% # Cluster all the data for each mountain AND genus
  summarize(Weight_avg = mean(Weight), Weight_sd = sd(Weight)) %>%
  print(n=20) # Print 20 rows, instead of the default 10
{% endhighlight %}



{% highlight text %}
## Source: local data frame [150 x 4]
## Groups: Mount [?]
## 
##     Mount Genera  Weight_avg Weight_sd
##    <fctr> <fctr>       <dbl>     <dbl>
## 1       1      1  -2.7573401  3.307295
## 2       1      2   8.5972132  4.181637
## 3       1      3   2.9854196  2.335625
## 4       1      4  21.0834808 11.479890
## 5       1      5  15.7930094  4.096622
## 6       1      6   7.4557598  4.285728
## 7       1      7 -19.7661844  3.677072
## 8       1      8  -4.3076198  4.861927
## 9       1      9  15.5273034  4.934868
## 10      1     10   0.7746332  4.438193
## 11      1     11  -5.6647755  2.666068
## 12      1     12   0.4239082  2.928701
## 13      1     13   0.6039859  5.401860
## 14      1     14  13.7291090  3.167712
## 15      1     15  11.2012089  2.628091
## 16      2      1  -3.3664910  1.907666
## 17      2      2   8.9181201  3.174841
## 18      2      3   1.2582477  2.300102
## 19      2      4  20.8658468 10.419646
## 20      2      5  15.1743084  3.173604
## # ... with 130 more rows
{% endhighlight %}

Perhaps simpler, we want to know how many species we found on each mountain.

{% highlight r %}
berg %>% 
  group_by(Mount) %>% 
  summarize(species_N_perMount = length(unique(Species))) 
{% endhighlight %}



{% highlight text %}
## # A tibble: 10 × 2
##     Mount species_N_perMount
##    <fctr>              <int>
## 1       1                 71
## 2       2                 73
## 3       3                 65
## 4       4                 74
## 5       5                 65
## 6       6                 73
## 7       7                 61
## 8       8                 67
## 9       9                 63
## 10     10                 74
{% endhighlight %}

## Data Visualization

We can even pipe these data frames right into `ggplot2`! Let's take a look at the variability in average weights across the different mountains. Here we'll average to the species level, because we have multiple individuals sampled from each species.


{% highlight r %}
berg %>% 
  group_by(Mount, Species) %>% # Cluster all the data for each mountain AND SPECIES
  summarize(Weight_avg = mean(Weight)) %>%
  ggplot(aes(y=Weight_avg, x=Mount)) + 
  geom_boxplot()
{% endhighlight %}

![center](/figs/2016-9-10-hierarchical-data-management-and-visualization/boxplot-1.png)

Now, more relevant to the question at hand, what is the relationship between body size and elevation? First, is it consistent across mountains? Let's look at this visually. We'll use the `facet_wrap()` function from the `ggplot2` package.

{% highlight r %}
ggplot(berg, aes(x=Elevation, y=Weight))+
  geom_point(shape=20)+
  facet_wrap(~Mount, ncol=5) # We'll see all 10 mountains separately
{% endhighlight %}

![center](/figs/2016-9-10-hierarchical-data-management-and-visualization/by_mount-1.png)

Upon each mountain it seems like there isn't a clear picture. What's going on? 

Let's use `facet_wrap()` to look at the each genus separately, and let's highlight the species with different colors. 

{% highlight r %}
ggplot(berg, aes(x=Elevation, y=Weight, color=Species))+
  geom_point(shape=20)+
  facet_wrap(~Genera, ncol=5)+
  guides(color=F) # There are many species, so suppress the color legend
{% endhighlight %}

![center](/figs/2016-9-10-hierarchical-data-management-and-visualization/by_genus-1.png)

Now it becomes clear why Bergmann's rule was being obscured when we looked across the whole data set and when we looked at mountains individually. There is a lot of variation among genera in the body size - elevation relationship. This makes sense, because we coded a large random variation in slope among genera. 

Go back and change the random variation in slopes and intercepts at the different hierarchical levels. How does this effect the overall pattern among mountains? This type of simulation can help you understand how much data you might need to collect in order to find statistically meaningful patterns. 
