---
output: github_document
---

# Argo Canada Development Blog

<!-- badges: start -->
[![Render & Deploy Site](https://github.com/ArgoCanada/blog/actions/workflows/render-distill.yaml/badge.svg?branch=master)](https://github.com/ArgoCanada/blog/actions/workflows/render-distill.yaml)
<!-- badges: end -->

The goal of the [Argo Canada Development Blog](https://argocanada.github.io/blog/) is to provide a venue for communicating chunks of code related to Argo development to facilitate synergy between various Argo development projects. To contribute to the blog:

1. [Fork this repository](https://github.com/ArgoCanada/blog/fork)
2. Clone to create a local copy.
3. Create a post from R using `distill::create_post("The Title Of My Post")`
4. Edit the .Rmd file. When you're done, click "Knit" in RStudio (`R -e 'rmarkdown::render("_posts/my-post-dir/the-title-of-my-post.Rmd")` for command-line R holdouts)
4. Commit (only commit the post folder...maybe Dewey can find a better workflow than keeping the rendered blog in the docs/ folder...) 
5. push, and create a [pull request](https://github.com/ArgoCanada/blog/pulls)

The blog is built by [distill for RMarkdown](https://rstudio.github.io/distill/), which is like blogdown but optimized for scientific publishing.
