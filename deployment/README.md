
# Canadian Argo Deployments

The page serves for Canadian Argo groups to submit deployment plans is as much or as little detail as possible. This may be done using any of the following methods: 

- Pull request (preferred)
- Form submission
- Email

## Submit Floats via Pull Request (PR)

The best way to submit a float is by forking this repository, making the required changes, and then creating a pull request. The PR will then be reviewed and merged by Chris at DFO, and the site will automatically update. This method is beneficial as it encourage good version control as well as tracking of who has submitted which floats. 

The changes can be:

1. Directly editing `canada_deployments.csv`
2. Submitting a file to the directory `file_submission`

The file `canada_deployments.csv` is the workhorse of the [Canadian deployments page](https://argocanada.github.io/blog/deployment.html). Editing this directly will update the website upon merging to the main branch of the repository. 

Submitting a file via the `file_submission` folder may be more convenient as groups may be able to simply submit files they are using for their own internal data management. In this case, the submitter just needs to ensure that the headings in their file match the headings used in `canada_deployments.csv`. It is ok for extra columns to be present, but columns that do not match will not be synthesized into the deployment table.

Column headings will follow OceanOps syntax. Headings in `canada_deployments.csv` are as follows:

PROGRAM, PRINCIPAL INVESTIGATOR, INSTITUTE, STATUS, MODEL, DEPLOYMENT DATE, DEPLOYMENT LAT, DEPLOYMENT LON, DEPLOYMENT SHIP, IMEI, REF, SERIAL NUMBER

Remember to submit a PR when you are done!!

## Form Submission

Floats may be submitted using the following [google form](https://docs.google.com/forms/d/e/1FAIpQLScsgQ2kJjSdpOTmVxt89-o6vQZ-7mUyqXlRuPBjgJOXiavfgQ/viewform?usp=sf_link). The form will only work for 1 float at a time, so may be convenient for single or few submissions. 

## Email Submission

Finally, floats may be submitted by simply emailing the file and/or float information to [chris.gordon@dfo-mpo.gc.ca](mailto:chris.gordon@dfo-mpo.gc.ca). If the float is already submitted to OceanOps, you can forward the notification email from submission.
