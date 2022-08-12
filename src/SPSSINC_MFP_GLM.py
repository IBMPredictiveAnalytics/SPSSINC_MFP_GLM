# module to run an R program as an Extension command
from __future__ import with_statement
#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 1989, 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

__author__ = "SPSS, JKP"
__version__ = "1.0.2"

# history
# 10-apr-2009 original version
#  12-jun-2009 add ask=FALSE parameter to plot for interactive mode
# 12-aug-2022 update for Python3
###import wingdbstub    ### debug

import spss, spssaux
from extension import Template, Syntax, checkrequiredparams, processcmd
import sys, inspect, tempfile, os, csv, codecs, textwrap

helptext="""SPSSINC MFP GLM  DEPENDENT=dependent variable 
COVAR = scale independent variables
FACTOR = categorical independent variables
VARIABLEALPHA=sig or list of sig, one for each variable
FPALPHA=sig or list of sig
DF= df or list of dfs, one for each scale variable
FAMILY = {gaussian* | binomial | poisson | Gamma | inverse.gaussian | quasi}
LINK = link-name
[/OPTIONS MISSING={LISTWISE*, FAIL}] [PLOTS]
[/SAVE [PROGRAMFILE=filespec] [RESIDUALSDATASET=datasetname]
[COEFSDATASET=datasetname].

DEPENDENT specifies the dependent variable.  Continuous predictors are
specified as COVAR and factors as FACTOR.

VARIABLEALPHA specifies the significance level for inclusion of variables,
and FPALPHA specifies the significance level for polynomial terms.  Both can be
a single number applied to all relevant variables or a list of numbers.  If a list,
there must be as many values as the total number of predictors (VARIABLEALPHA)
or the number of covariates.  Use 1 for non-polynomial covariates.

DF, also a single number or a list of numbers, specifies the parameters for
the polynomials.  It may be 1: linear model, 2: FP model with maximum degree 1,
4: FP model with maximum degree 2.

FAMILY specifies the model, defaulting to gaussian.

LINK specifies the link function.  Valid values depend on the family.
mu stands for 1/mu^2
gaussian:  identity* | log | inverse
binomial:   logit* | probit | cauchit | cloglog
gamma:    inverse* | identity | log
poisson:    log* | identity | sqrt
inversegaussian: mu | inverse | identity | log
quasi:        logit* | probit | cloglog | identity | inverse |log | mu | sqrt


Split files and weight are not honored by this command.

SPSSINC MFP GLM /HELP prints this information and does nothing else.

Example:
SPSSINC MFPGLM DEPENDENT=mpg COVAR=engine weight FACTOR=origin
DF=3 3.

Estimate a regression with fractional polynomials for the scale variables

DEPENDENT  specifies the dependent variable name.  

Categorical independent variables are automatically converted
appropriately to factors.  A constant term is automatically included.



MISSING = LISTWISE causes listwise deletion of missing values.  FAIL 
stops the procedure if missing values are encountered.

/SAVE PROGRAMFILE causes the R code for this procedure to be written to 
the specified file.  The generated program can be a useful starting point 
for additional specifications.

RESIDUALSDATASET causes a dataset containing the residuals to be created.
The dataset name must not already be in use.
The case number is included as cases will only be written for input cases with no
missing data.  Data filtered by SPSS Statistics are not passed to R and will not
generate cases in the residuals dataset.

COEFSDATASET causes a new dataset containing the coefficients to be created.
The dataset name must not already be in use.

This extension command requires the Python and R programmability plug-ins and 
the R car package.
"""
familylinks = {
    "gaussian" : ["identity", "log", "inverse"],
        "binomial" : ["logit","probit","cauchit","cloglog"],
        "gamma" : ["inverse", "identity","log"],
        "poisson" : ["log","identity","sqrt"],
        "inversegaussian" : ["mu","inverse","identity","log"],
        "quasi" : ["logit","probit","cloglog","identity","inverse","log","mu","sqrt"]}

def fpreg(dep, covar,factor=None, variablealpha=None, fpalpha=None, 
          df=[4], family="gaussian", link=None, missing="listwise", 
       programfile=None, residualsdataset=None, coefsdataset=None,
       plots=False):
    """run fractional polynomial regression via mfp package

    dep is the dependent variable name.


    missing is "listwise" or "fail" according to the desired missing value treatment.
    programfile, if specified, causes the generated R program to be written to that file.
    residualsdataset causes a dataset of residuals to be created.
    coefdataset causes a dataset of terms and coefficients to be created.
"""

    # This function generates an R program and submits it via INSERT.
    if factor is None:
        factor = []
    if variablealpha is None:
        variablealpha = [.05]
    if len(variablealpha) == 1:
        variablealpha = len(covar)*variablealpha
    alpha0 = variablealpha[0]
    if fpalpha is None:
        fpalpha = [.05]
    if len(fpalpha) == 1:
        fpalpha = len(covar) * [fpalpha[0]]
    if len(df) == 1:
        df = len(covar) * [df[0]]
    if 3 in df:
        raise ValueError("df must be 1, 3, or 4")
    if len(fpalpha) != len(covar):
        raise ValueError("The number of fp alpha values is different from the number of covariates")
    if len(variablealpha) != len(covar):
        raise ValueError("The number of variable alpha values is different from the number of variables")
    if len(df) != len(covar):
        raise ValueError("The number of df specifications is different from the number of covariates")
    if link is not None and not link in familylinks[family]:
        raise ValueError("Link function not valid for family.  Valid choices are " + ", ".join(familylinks[family]))
    if family == "inversegaussian":
        family = "inverse.gaussian"
    if link == "mu":
        link = "1/mu^2"
    if not link is None:

        family = family + "(" + link + ")"
    missing = missing == "listwise" and "na.exclude" or "na.fail"
    plots = plots and "TRUE" or "FALSE"
    # check the independent variables measurement levels and code accordingly
    dfex = []
    catlistnames = ['"' + item + '"' for item in factor]
    catlistrindex = [str(i+len(covar) + 2) for i in range(len(factor))]
    allvars = "c(" + ", ".join(['"' + v + '"' for v in [dep] + covar + factor]) + ")"
    allvars = "\n".join(textwrap.wrap(allvars, width=100))
    covarfuncs = []
    for i, v in enumerate(covar):
        covarfuncs.append("fp(%s, df=%s, select=%s, alpha=%s)" % \
                                  (v, df[i], variablealpha[i], fpalpha[i]))
    model = dep + "~" + \
            "\n".join(textwrap.wrap(" +".join(covarfuncs + factor)))

    catlistnames = "c(" + ", ".join(catlistnames) + ")"
    catlistnames = "\n".join(textwrap.wrap(catlistnames, width=100))
    catlistrindex =  "c(" + ", ".join(catlistrindex) + ")"
    catlistrindex = "\n".join(textwrap.wrap(catlistrindex, width=100))

    if residualsdataset:
        residuals = """dict<- spssdictionary.CreateSPSSDictionary(c("caseNumber", "Case Number", 0, "F8.0", "nominal"),
c("mfpResiduals", as.character(res$call[2]), 0, "F8.2", "scale"))
tryCatch({
spssdictionary.SetDictionaryToSPSS("%(residualsdataset)s", dict)
df = data.frame(reslm$residuals)
spssdata.SetDataToSPSS("%(residualsdataset)s", data.frame(row.names(df), reslm$residuals))},
error=function(e) {print(e)
print("Failed to create residuals dataset.  Dataset name must not already exist: %(residualsdataset)s")}
)
""" % locals()
    else:
        residuals= ""
    if coefsdataset:
        coefs = """dict<- spssdictionary.CreateSPSSDictionary(c("term", "Variable or Factor Value", 100, "A100", "nominal"),
c("coefficient", "Estimated Coefficient", 0, "F10.3", "scale"))
tryCatch({
spssdictionary.SetDictionaryToSPSS("%(coefsdataset)s", dict)
spssdata.SetDataToSPSS("%(coefsdataset)s", data.frame(row.names(res$coef), res$coef[,1]))},
error=function(e) {print(e)
print("Failed to create coefficients dataset.  Dataset name must not already exist: %(coefsdataset)s")
})
""" % locals()
    else:
        coefs = ""

    pgm = r"""BEGIN PROGRAM R.
library(mfp)
makefactor = function(var, vallabels, ordfac) {
# return a labeled factor for variable var and label set vallabels
# if vallabels is NULL, the raw values will label the factor.
# if ordfac is TRUE, returns an ordered factor
values = sort(unique(var))
lbls = values
lmatch = match(vallabels$values, values, nomatch=0)
lmatch2 = match(values[lmatch], vallabels$values)
lbls[lmatch] = vallabels$labels[lmatch2]
if (ordfac) {
            return (ordered(var, levels=values, labels=lbls))}
else {
            return (factor(var, levels=values, labels=lbls))
}
}
loopfactor = function(vars, indexes) {
if (length(indexes) == 0) {return(dta)}
for (v in 1:length(indexes)) {
            vl = spssdictionary.GetValueLabels(vars[v])
            dta[indexes[v]] = makefactor(dta[[indexes[v]]], vl, FALSE)
}
return(dta)
}

catlistnames = %(catlistnames)s
catlistrindex = %(catlistrindex)s
plots = %(plots)s

dta<-spssdata.GetDataFromSPSS(%(allvars)s)
# Convert NaN values to NA
is.na(dta)<- is.na(dta)

dta=loopfactor(catlistnames, catlistrindex)

res <- tryCatch(
            summary(
            reslm <- mfp(%(model)s, 
            data=dta, na.action=%(missing)s,
            family=%(family)s, alpha=%(alpha0)s)
            ),
error=function(e) 
            {return(c("ERROR:",e))}
)

if (!is.null(res$message)) {print(res$message);break} else {
            caption = paste("Dependent variable:", strsplit(as.character(reslm$call[2]),"~")[[1]][1])
            f1=res$family$family  # does not work directly as sprintf argument
            f2=res$family$link
            f3=res$dispersion  
            caption = paste(caption, sprintf("Family: %%s.  Link function: %%s.  Dispersion parameter: %%.3f",
            f1, f2, f3), sep="\n")
            colnames(res$coefficients)[4] = "Sig."
            spsspivottable.Display(coefficients(res), 
            title="Multiple Fractional Polynomial Linear Model: Coefficients", "SPSSINCMFPMODEL",
            caption=caption,
            isSplit=FALSE,
            format = formatSpec.Coefficient)

            colnames(reslm$fptable)[3] = "Fpalpha"
            colnames(reslm$fptable)[2] = "Variable alpha"
            spsspivottable.Display(reslm$fptable, title="Multiple Fractional Polynomial Linear Model: Fp Transformations", "SPSSINCMFPFPTRANS",
            caption=reslm$formula[3], isSplit=FALSE)

            spsspivottable.Display(reslm$scale, title="Multiple Fractional Polynomial Linear Model: Predictor Scaling",
            "SPSSINCMFPSCALE", isSplit=FALSE)

            rowl=c("Null deviance", "Residual deviance", "AIC")
            c1 = c(round(res$null.deviance, 3), round(res$deviance, 3), round(res$aic,3))
            c2 = c(res$df.null, res$df.residual, "")
            dfr = data.frame(cbind(c1, c2), row.names=rowl)
            colnames(dfr)[1:2] = c("Statistic", "D. f.")
            spsspivottable.Display(dfr, 
            title="Multiple Fractional Polynomial Regression: Deviance",    
            "SPSSINCMFPSTATS", isSplit=FALSE)

            if (plots) {
            plot(reslm, main="Multiple Fractional Polynomial Linear Model", col="blue", lwd=2, ask=FALSE)
            }

            %(residuals)s
            %(coefs)s
            spssdictionary.EndDataStep()

}

rm(dta)
rm(res)
rm(reslm)
END PROGRAM   .
""" % locals()
    # write program to a temporary or user-specified file
    if programfile:
        cmdfile = programfile.replace("\\", "/")
    else:
        cmdfile = (tempfile.gettempdir() + os.sep + "pgm.R").replace("\\", "/")
    f = codecs.open(cmdfile, "wb", encoding="utf_8_sig")
    f.write(pgm)
    f.close()


    spss.Submit("INSERT FILE='%s'" % cmdfile)
    if not programfile:
        os.remove(cmdfile)



def Run(args):
    """Execute the SPSSINC MFP GLM command"""

    args = args[list(args.keys())[0]]

    oobj = Syntax([
            Template("DEPENDENT", subc="",  ktype="existingvarlist", var="dep", islist=False),
                Template("COVAR", subc="",  ktype="existingvarlist", var="covar", islist=True),
                Template("FACTOR", subc="",  ktype="existingvarlist", var="factor", islist=True),
                Template("VARIABLEALPHA", subc="", ktype="float", var="variablealpha", vallist=[0.,1.], islist=True),
                Template("FPALPHA", subc="", ktype="float", var="fpalpha", vallist=[0.,1.], islist=True),
                Template("DF", subc="", ktype="int", var="df", vallist=[1,4], islist=True),
                Template("FAMILY", subc="", ktype="str", var="family", 
                         vallist=["gaussian", "binomial", "poisson", "gamma", "inversegaussian", "quasi"], islist=False),
                Template("LINK", subc="", ktype="str", var="link"),
                Template("MISSING", subc="OPTIONS",ktype="str", var="missing"),
                Template("PLOTS", subc="OPTIONS", ktype="bool", var="plots"),
                Template("PROGRAMFILE", subc="SAVE", ktype="literal", var="programfile"),
                Template("RESIDUALSDATASET", subc="SAVE", ktype="literal", var="residualsdataset"),
                Template("COEFSDATASET", subc="SAVE", ktype="literal", var="coefsdataset"),
                Template("HELP", subc="", ktype="bool")])

    # A HELP subcommand overrides all else
    if "HELP" in args:
        #print helptext
        helper()
    else:
        processcmd(oobj, args, fpreg, vardict=spssaux.VariableDict())

def helper():
    """open html help in default browser window

    The location is computed from the current module name"""

    import webbrowser, os.path

    path = os.path.splitext(__file__)[0]
    helpspec = "file://" + path + os.path.sep + \
        "markdown.html"

    # webbrowser.open seems not to work well
    browser = webbrowser.get()
    if not browser.open_new(helpspec):
        print("Help file not found:" + helpspec)
try:    #override
    from extension import helper
except:
    pass