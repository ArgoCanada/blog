
# Argo Canada Development Blog

<!-- badges: start -->
[![Render & Deploy Site](https://github.com/ArgoCanada/blog/actions/workflows/render-distill.yaml/badge.svg?branch=master)](https://github.com/ArgoCanada/blog/actions/workflows/render-distill.yaml)
<!-- badges: end -->

The goal of the [Argo Canada Development Blog](https://argocanada.github.io/blog/) is to provide a venue for communicating chunks of code related to Argo development to facilitate synergy between various Argo development projects. To contribute to the blog

1. [Fork this repository](https://github.com/ArgoCanada/blog/fork)
2. Clone to create a local copy.
3. Create a post from R using `distill::create_post("The Title Of My Post")`
4. Edit the .Rmd file. When you're done, click "Knit" in RStudio (`R -e 'rmarkdown::render("_posts/my-post-dir/the-title-of-my-post.Rmd")` for command-line R holdouts). You can also run `rmarkdown::render_site()` to generate a preview of the entire site.
5. Commit, push, and create a [pull request](https://github.com/ArgoCanada/blog/pulls)

For Python-based posts, you can use the Python engine for RMarkdown (which uses [reticulate](https://rstudio.github.io/reticulate/) under the hood and might have to be configured in the first R chunk). Alternatively, you can create your post as an IPython Notebook (e.g., in VSCode using 'Create New Blank notebook' from the command palette) and convert it to Markdown using `python[3] -m nbconvert --to markdown path/to/post.ipynb`. You can then copy the markdown content to the .Rmd file you just created and click 'Knit' (or `rmarkdown::render("path/to/post.Rmd")` from the console).

The blog is built by [distill for RMarkdown](https://rstudio.github.io/distill/), which is like blogdown but optimized for scientific publishing.
